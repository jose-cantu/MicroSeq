from __future__ import annotations  # will help me with using typehints in my code w/o throwing out immediate errors 

import os 
import stat 
import subprocess
from pathlib import Path 

import pytest

pytest.importorskip("pandas")
pytest.importorskip("Bio") 

# --- Test Suite Summary ---
# This pytest suite validates functions within the 'microseq_tests.assembly'
# package that handle paired-end DNA sequence assembly logic. 
# 
# Key functions under test include filename parsing (extract_sid_orientation), 
# file grouping logic (group_pairs), the main assembly orchestrator (assemble_pairs),
# and the command-line interface (CLI) wrapper for the paired assembly mode.
# 
# Tests employ mocking (monkeypatching) to simulate the behavior of external tools like 
# CAP3 without requiring them to be installed or run in the test environment.
# --------------------------

from microseq_tests.assembly import paired_assembly 
from microseq_tests.assembly.pairing import (
        DupPolicy,
        extract_sid_orientation,
        group_pairs,
        _extract_well
) 

CLI = ["python", "-m", "microseq_tests.microseq"] 

def test_extract_sid_orientation_basic(): 
    """
    Filename detectors should recover sample ID and orientation tokens 
    from standardized input filenames. 
    """ 
    # Assert specific input patterns map to expected Sample IDs (SID) and orient:
    assert extract_sid_orientation("A3_27F_x.fasta") == ("A3", "F")
    assert extract_sid_orientation("A3-1492R.fa") == ("A3", "R") 

def test_extract_sid_orientation_honors_custom_detectors():
    """Custom detectors passed in should take precedence over the defaults."""

    def custom(name: str) -> tuple[str, str | None]:
        if "custom-token" in name:
            return "sid-custom", "F"
        return Path(name).stem, None

    assert extract_sid_orientation(
        "sample_custom-token_read.fasta", detectors=[custom]
    ) == ("sid-custom", "F")

def test_group_pairs_error_on_dups(tmp_path: Path):
    """
    Verifies that duplicate forward reads raise an error under the strict 
    (DupPolicy.ERROR) policy when grouping files. 
    """ 
    # Create input files with duplicate forward reads for the same sample 'S1'. 
    (tmp_path / "S1_27F.fasta").write_text(">x\nA\n", encoding="utf-8")
    (tmp_path / "S1_8F.fasta").write_text(">y\nA\n", encoding="utf-8")
    (tmp_path / "S1_1492R.fasta").write_text(">z\nA\n", encoding="utf-8")

    # Assert that calling group_pairs with DupPolicy.Error raises the expected Valuue error 
    with pytest.raises(ValueError):
        group_pairs(tmp_path, dup_policy=DupPolicy.ERROR)


def test_group_pairs_enforce_well_requires_match(tmp_path: Path):
    (tmp_path / "sample_27F_A01.fasta").write_text(">x\nA\n", encoding="utf-8")
    (tmp_path / "sample_27F_A01.fasta").write_text(">x\nA\n", encoding="utf-8")

    pairs = group_pairs(tmp_path, enforce_same_well=True)

    assert pairs == {}

def test_group_pairs_report_missing_wells(tmp_path: Path):
    (tmp_path / "sample_27F.fasta").write_text(">x\nA\n", encoding="utf-8")
    (tmp_path / "sample_1492R.fasta").write_text(">x\nA\n", encoding="utf-8")
    
    pairs, meta = group_pairs(tmp_path, enforce_same_well=True, return_metadata=True)

    assert pairs == {}
    assert meta.get("_missing_well", {}).get("files")

def test_extract_well_validates_plate_range():
    assert _extract_well("sample_A1.fastq") == "A01"
    assert _extract_well("sample_h12.fastq") == "H12"
    assert _extract_well("sample_I01.fastq") is None
    assert _extract_well("sample_A13.fastq") is None


def test_assemble_pairs_creates_contigs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Running assemble_pairs should surface CAP3 artefacts for each sample."""

    in_dir = tmp_path / "input"
    out_dir = tmp_path / "output"
    in_dir.mkdir()

    (in_dir / "S1_27F.fa").write_text(">a\nA\n", encoding="utf-8")
    (in_dir / "S1_1492R.fa").write_text(">b\nA\n", encoding="utf-8")

    monkeypatch.setattr(
        paired_assembly,
        "load_config",
        lambda: {"tools": {"cap3": "/bin/true"}},
    )

    def fake_run(cmd, check, cwd, stderr, text):
        """Stub CAP3 invocation by writing the expected contig output."""

        contig = Path(cwd, f"{cmd[1]}.cap.contigs")
        contig.write_text(">contig\nA\n", encoding="utf-8")

        class Result:
            returncode = 0

        return Result()

    monkeypatch.setattr(paired_assembly.subprocess, "run", fake_run)

    paths = paired_assembly.assemble_pairs(in_dir, out_dir)

    expected = out_dir / "S1" / "S1_paired.fasta.cap.contigs"
    assert expected in paths


def test_assemble_pairs_keep_separate(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Keep-separate should fan out duplicates into individual assemblies."""

    in_dir = tmp_path / "input"
    out_dir = tmp_path / "output"
    in_dir.mkdir()

    forward1 = in_dir / "S1_27F.fa"
    forward2 = in_dir / "S1_8F.fa"
    reverse = in_dir / "S1_1492R.fa"

    forward1.write_text(">a\nA\n", encoding="utf-8")
    forward2.write_text(">b\nA\n", encoding="utf-8")
    reverse.write_text(">c\nA\n", encoding="utf-8")

    monkeypatch.setattr(
        paired_assembly,
        "load_config",
        lambda: {"tools": {"cap3": "/bin/true"}},
    )

    calls: list[tuple[tuple[str, ...], Path]] = []

    def fake_run(cmd, check, cwd, stderr, text):  # noqa: D401
        """Record CAP3 invocations and synthesize contig outputs."""

        calls.append((tuple(cmd), Path(cwd)))
        contig = Path(cwd, f"{cmd[1]}.cap.contigs")
        contig.write_text(">contig\nA\n", encoding="utf-8")

        class Result:
            returncode = 0

        return Result()

    monkeypatch.setattr(paired_assembly.subprocess, "run", fake_run)

    paths = paired_assembly.assemble_pairs(
        in_dir,
        out_dir,
        dup_policy=DupPolicy.KEEP_SEPARATE,
    )

    expected_dirs = [out_dir / "S1_1", out_dir / "S1_2"]
    expected_contigs = [
        directory / f"{directory.name}_paired.fasta.cap.contigs" for directory in expected_dirs
    ]

    assert paths == expected_contigs
    assert all(path.exists() for path in expected_contigs)

    # Two distinct CAP3 runs should have been executed in their respective directories.
    assert {cwd for _, cwd in calls} == set(expected_dirs)



def test_cli_paired_mode(tmp_path: Path):
    """
    An end to end integration test verifying the main command line interface for 
    paired assembly works correctly. 

    It assserts that the CLI wrapper executates the necessary steps and produces the final output file 
    in the correct directory structure. 
    """

    in_dir = tmp_path / "input"
    out_dir = tmp_path / "output" 
    in_dir.mkdir() 

    # Create dummy input files 
    (in_dir / "S1_27F.fa").write_text(">a\nA\n", encoding="utf-8")
    (in_dir / "S1_1492R.fa").write_text(">b\nA\n", encoding="utf-8")

    # Create a real, exectuble dummy 'cap3' script using standard shell command. 
    fake_cap3 = tmp_path / "cap3" 
    fake_cap3.write_text(
        "#!/bin/sh\n" "printf '>contig\\nA\\n' > \"$1.cap.contigs\"\n",
        encoding="utf-8",
    ) 
    # make the script execuatble usign file permissiions (chmod +x) 
    fake_cap3.chmod(fake_cap3.stat().st_mode | stat.S_IEXEC) 

    # Modify environment variables so test subprocess can fine the fake 'cap3'.
    env = os.environ.copy() 
    env["PATH"] = f"{tmp_path}{os.pathsep}{env['PATH']}"

    # Run the actual CLI command via subprocess. 
    res = subprocess.run( 
       CLI + [ "assembly", "--mode", "paired", "-i", str(in_dir), "-o", str(out_dir), ],
       check=False,# Check return code manually instead below 
       capture_output=True, 
       text=True,
       env=env,
    )

    # Assert the command ran succesfully (return code 0) and expected output file exists 
    assert res.returncode == 0, res.stderr 
    assert (out_dir / "S1" / "S1_paired.fasta.cap.contigs").exists() 












