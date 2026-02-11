from __future__ import annotations  # will help me with using typehints in my code w/o throwing out immediate errors 

import os 
import stat 
import subprocess
import types
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
    assert _extract_well("sample_-B07.fastq") == "B07" 
    assert _extract_well("sample__A10.fastq") == "A10" 
    assert _extract_well("sample_I01.fastq") is None
    assert _extract_well("sample_A13.fastq") is None


def test_assemble_pairs_creates_contigs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Running assemble_pairs should surface CAP3 artefacts for each sample."""

    in_dir = tmp_path / "input"
    out_dir = tmp_path / "output"
    in_dir.mkdir()

    (in_dir / "S1_27F.fasta").write_text(">a\nA\n", encoding="utf-8")
    (in_dir / "S1_1492R.fasta").write_text(">b\nA\n", encoding="utf-8")

    monkeypatch.setattr(
        paired_assembly,
        "load_config",
        lambda: {"tools": {"cap3": "/bin/true"}},
    )

    def fake_run(cmd, check, cwd=None, **kwargs):
        """Stub CAP3 invocation by writing the expected contig output."""

        if "--version" in cmd:
            class Result:
                returncode = 0
                stdout = "CAP3 test"
                stderr = ""

            return Result()

        contig = Path(cwd, f"{cmd[1]}.cap.contigs")
        contig.write_text(">contig\nA\n", encoding="utf-8")

        class Result:
            returncode = 0
            stdout = ""
            stderr = ""

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

    forward1 = in_dir / "S1_27F.fasta"
    forward2 = in_dir / "S1_8F.fasta"
    reverse = in_dir / "S1_1492R.fasta"

    forward1.write_text(">a\nA\n", encoding="utf-8")
    forward2.write_text(">b\nA\n", encoding="utf-8")
    reverse.write_text(">c\nA\n", encoding="utf-8")

    monkeypatch.setattr(
        paired_assembly,
        "load_config",
        lambda: {"tools": {"cap3": "/bin/true"}},
    )

    calls: list[tuple[tuple[str, ...], Path]] = []

    def fake_run(cmd, check, cwd=None, **kwargs):  # noqa: D401
        """Record CAP3 invocations and synthesize contig outputs."""

        if "--version" in cmd:
            class Result:
                returncode = 0
                stdout = "CAP3 test"
                stderr = ""

            return Result()

        calls.append((tuple(cmd), Path(cwd)))
        contig = Path(cwd, f"{cmd[1]}.cap.contigs")
        contig.write_text(">contig\nA\n", encoding="utf-8")

        class Result:
            returncode = 0
            stdout = ""
            stderr = ""

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


def test_assemble_pairs_non_singleton_merge_fallback_logs_counts(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    in_dir = tmp_path / "input"
    out_dir = tmp_path / "output"
    in_dir.mkdir()

    (in_dir / "S1_27F.fasta").write_text(">a\nA\n>b\nA\n", encoding="utf-8")
    (in_dir / "S1_1492R.fasta").write_text(">c\nA\n>d\nA\n", encoding="utf-8")

    monkeypatch.setattr(
        paired_assembly,
        "load_config",
        lambda: {
            "tools": {"cap3": "/bin/true"},
            "overlap_eval": {
                "min_overlap": 100,
                "min_identity": 0.8,
                "min_quality": 20.0,
                "quality_mode": "warning",
                "cap3_validate_pair_support": True,
            },
        },
    )

    def fake_run(cmd, check, cwd=None, **kwargs):
        if "--version" in cmd:
            class Result:
                returncode = 0
                stdout = "CAP3 test"
                stderr = ""

            return Result()

        contig = Path(cwd, f"{cmd[1]}.cap.contigs")
        contig.write_text(">contig\nA\n", encoding="utf-8")

        class Result:
            returncode = 0
            stdout = ""
            stderr = ""

        return Result()

    monkeypatch.setattr(paired_assembly.subprocess, "run", fake_run)

    paired_assembly.assemble_pairs(in_dir, out_dir)
    meta = (out_dir / "cap3_run_metadata.txt").read_text(encoding="utf-8")
    assert "merge_two_reads_skipped_non_singleton" in meta
    assert "f_n=2" in meta
    assert "r_n=2" in meta



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
    (in_dir / "S1_27F.fasta").write_text(">a\nA\n", encoding="utf-8")
    (in_dir / "S1_1492R.fasta").write_text(">b\nA\n", encoding="utf-8")

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
    repo_src = Path(__file__).resolve().parents[1] / "src"
    env["PYTHONPATH"] = f"{repo_src}{os.pathsep}{env.get('PYTHONPATH', '')}".rstrip(os.pathsep)

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












def test_assemble_pairs_runs_cap3_fallback_on_nonmerged_fast_merge(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    in_dir = tmp_path / "input"
    out_dir = tmp_path / "output"
    in_dir.mkdir()

    (in_dir / "S1_27F.fasta").write_text(">a\nA\n", encoding="utf-8")
    (in_dir / "S1_1492R.fasta").write_text(">b\nA\n", encoding="utf-8")

    monkeypatch.setattr(
        paired_assembly,
        "load_config",
        lambda: {
            "tools": {"cap3": "/bin/true"},
            "overlap_eval": {
                "min_overlap": 100,
                "min_identity": 0.8,
                "min_quality": 20.0,
                "quality_mode": "warning",
            },
        },
    )

    calls: list[tuple[str, ...]] = []

    def fake_run(cmd, check, cwd=None, **kwargs):
        calls.append(tuple(cmd))
        if "--version" in cmd:
            class Result:
                returncode = 0
                stdout = "CAP3 test"
                stderr = ""

            return Result()

        contig = Path(cwd, f"{cmd[1]}.cap.contigs")
        contig.write_text(">contig\nA\n", encoding="utf-8")

        class Result:
            returncode = 0
            stdout = ""
            stderr = ""

        return Result()

    monkeypatch.setattr(paired_assembly.subprocess, "run", fake_run)

    paths = paired_assembly.assemble_pairs(in_dir, out_dir)
    assert paths
    assert any(len(cmd) > 1 and cmd[1].endswith("_paired.fasta") for cmd in calls)


def test_assemble_pairs_blocking_quality_low_skips_cap3_fallback(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    in_dir = tmp_path / "input"
    out_dir = tmp_path / "output"
    in_dir.mkdir()

    (in_dir / "S1_27F.fasta").write_text(">a\n" + "A" * 120 + "\n", encoding="utf-8")
    (in_dir / "S1_1492R.fasta").write_text(">b\n" + "A" * 120 + "\n", encoding="utf-8")

    monkeypatch.setattr(
        paired_assembly,
        "load_config",
        lambda: {
            "tools": {"cap3": "/bin/true"},
            "overlap_eval": {
                "min_overlap": 100,
                "min_identity": 0.8,
                "min_quality": 20.0,
                "quality_mode": "blocking",
            },
        },
    )

    calls: list[tuple[str, ...]] = []

    def fake_run(cmd, check, cwd=None, **kwargs):
        calls.append(tuple(cmd))
        if "--version" in cmd:
            class Result:
                returncode = 0
                stdout = "CAP3 test"
                stderr = ""

            return Result()

        raise AssertionError("CAP3 fallback should not run for blocking quality_low")

    monkeypatch.setattr(paired_assembly.subprocess, "run", fake_run)

    paths = paired_assembly.assemble_pairs(in_dir, out_dir)
    assert paths == []
    assert all("--version" in cmd for cmd in calls)


def test_assemble_pairs_runs_cap3_fallback_on_ambiguous_overlap(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    in_dir = tmp_path / "input"
    out_dir = tmp_path / "output"
    in_dir.mkdir()

    (in_dir / "S1_27F.fasta").write_text(">a\n" + "A" * 120 + "\n", encoding="utf-8")
    (in_dir / "S1_1492R.fasta").write_text(">b\n" + "A" * 120 + "\n", encoding="utf-8")

    monkeypatch.setattr(
        paired_assembly,
        "load_config",
        lambda: {
            "tools": {"cap3": "/bin/true"},
            "overlap_eval": {
                "min_overlap": 100,
                "min_identity": 0.8,
                "min_quality": 20.0,
                "quality_mode": "warning",
                "ambiguity_identity_delta": 0.0,
                "cap3_validate_pair_support": False,
            },
        },
    )

    calls: list[tuple[str, ...]] = []

    def fake_merge_two_reads(**_kwargs):
        report = types.SimpleNamespace(
            orientation="forward",
            overlap_len=120,
            identity=0.99,
            merge_status="ambiguous_overlap",
            high_conflict_mismatches=0,
        )
        return None, report

    def fake_run(cmd, check, cwd=None, **kwargs):
        calls.append(tuple(cmd))
        if "--version" in cmd:
            class Result:
                returncode = 0
                stdout = "CAP3 test"
                stderr = ""

            return Result()
        contig = Path(cwd, f"{cmd[1]}.cap.contigs")
        contig.write_text(">contig\n" + "A" * 120 + "\n", encoding="utf-8")

        class Result:
            returncode = 0
            stdout = ""
            stderr = ""

        return Result()

    monkeypatch.setattr(paired_assembly, "merge_two_reads", fake_merge_two_reads)
    monkeypatch.setattr(paired_assembly.subprocess, "run", fake_run)

    paths = paired_assembly.assemble_pairs(in_dir, out_dir)
    assert paths
    assert any(len(cmd) > 1 and cmd[1].endswith("_paired.fasta") for cmd in calls)


def test_assemble_pairs_high_conflict_route_runs_cap3(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    in_dir = tmp_path / "input"
    out_dir = tmp_path / "output"
    in_dir.mkdir()

    fwd_path = in_dir / "S1_27F.fasta"
    rev_path = in_dir / "S1_1492R.fasta"
    fwd_seq = "A" * 120
    rev_seq = "A" * 119 + "T"
    fwd_path.write_text(f">a\n{fwd_seq}\n", encoding="utf-8")
    rev_path.write_text(f">b\n{rev_seq}\n", encoding="utf-8")

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO
    fq = SeqRecord(Seq(fwd_seq), id="a", description="")
    rq = SeqRecord(Seq(rev_seq), id="b", description="")
    fq.letter_annotations["phred_quality"] = [40] * len(fwd_seq)
    rq.letter_annotations["phred_quality"] = [40] * len(rev_seq)
    SeqIO.write([fq], Path(f"{fwd_path}.qual"), "qual")
    SeqIO.write([rq], Path(f"{rev_path}.qual"), "qual")

    monkeypatch.setattr(
        paired_assembly,
        "load_config",
        lambda: {
            "tools": {"cap3": "/bin/true"},
            "overlap_eval": {
                "min_overlap": 100,
                "min_identity": 0.8,
                "min_quality": 20.0,
                "quality_mode": "warning",
                "high_conflict_q_threshold": 30,
                "high_conflict_action": "route_cap3",
                "cap3_validate_pair_support": False,
            },
        },
    )

    calls: list[tuple[str, ...]] = []

    def fake_run(cmd, check, cwd=None, **kwargs):
        calls.append(tuple(cmd))
        if "--version" in cmd:
            class Result:
                returncode = 0
                stdout = "CAP3 test"
                stderr = ""

            return Result()

        contig = Path(cwd, f"{cmd[1]}.cap.contigs")
        contig.write_text(">contig\n" + "A" * 121 + "\n", encoding="utf-8")

        class Result:
            returncode = 0
            stdout = ""
            stderr = ""

        return Result()

    monkeypatch.setattr(paired_assembly.subprocess, "run", fake_run)

    paths = paired_assembly.assemble_pairs(in_dir, out_dir)
    assert paths
    assert any(len(cmd) > 1 and cmd[1].endswith("_paired.fasta") for cmd in calls)


def test_assemble_pairs_cap3_validation_failed_marks_unverified(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    in_dir = tmp_path / "input"
    out_dir = tmp_path / "output"
    in_dir.mkdir()

    (in_dir / "S1_27F.fasta").write_text(">a\n" + "A" * 120 + "\n", encoding="utf-8")
    (in_dir / "S1_1492R.fasta").write_text(">b\n" + "C" * 120 + "\n", encoding="utf-8")

    monkeypatch.setattr(
        paired_assembly,
        "load_config",
        lambda: {
            "tools": {"cap3": "/bin/true"},
            "overlap_eval": {
                "min_overlap": 200,
                "min_identity": 1.0,
                "min_quality": 20.0,
                "quality_mode": "warning",
                "cap3_validate_pair_support": True,
            },
        },
    )

    def fake_run(cmd, check, cwd=None, **kwargs):
        if "--version" in cmd:
            class Result:
                returncode = 0
                stdout = "CAP3 test"
                stderr = ""

            return Result()

        contig = Path(cwd, f"{cmd[1]}.cap.contigs")
        contig.write_text(">contig\n" + "A" * 120 + "\n", encoding="utf-8")

        class Result:
            returncode = 0
            stdout = ""
            stderr = ""

        return Result()

    monkeypatch.setattr(paired_assembly.subprocess, "run", fake_run)

    paths = paired_assembly.assemble_pairs(in_dir, out_dir)
    assert paths == []
    marker = out_dir / "S1" / "S1_paired.cap3_validation.txt"
    assert marker.exists()
    assert marker.read_text(encoding="utf-8").strip() == "failed"
    assert (out_dir / "S1" / "S1_paired.fasta.cap.singlets").exists()


def test_assemble_pairs_cap3_validation_verified(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    in_dir = tmp_path / "input"
    out_dir = tmp_path / "output"
    in_dir.mkdir()

    fwd = "A" * 120
    rev = "T" * 120
    (in_dir / "S1_27F.fasta").write_text(f">a\n{fwd}\n", encoding="utf-8")
    (in_dir / "S1_1492R.fasta").write_text(f">b\n{rev}\n", encoding="utf-8")

    monkeypatch.setattr(
        paired_assembly,
        "load_config",
        lambda: {
            "tools": {"cap3": "/bin/true"},
            "overlap_eval": {
                "min_overlap": 200,
                "min_identity": 1.0,
                "min_quality": 20.0,
                "quality_mode": "warning",
                "cap3_validate_pair_support": True,
            },
        },
    )

    def fake_run(cmd, check, cwd=None, **kwargs):
        if "--version" in cmd:
            class Result:
                returncode = 0
                stdout = "CAP3 test"
                stderr = ""

            return Result()

        contig = Path(cwd, f"{cmd[1]}.cap.contigs")
        # include fwd and revcomp(rev)=A*120 so both reads are represented
        contig.write_text(">contig\n" + "A" * 240 + "\n", encoding="utf-8")

        class Result:
            returncode = 0
            stdout = ""
            stderr = ""

        return Result()

    monkeypatch.setattr(paired_assembly.subprocess, "run", fake_run)

    paths = paired_assembly.assemble_pairs(in_dir, out_dir)
    assert len(paths) == 1
    marker = out_dir / "S1" / "S1_paired.cap3_validation.txt"
    assert marker.exists()
    assert marker.read_text(encoding="utf-8").strip() == "verified"
