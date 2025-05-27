"""
tests/test_trim_pipeline.py

❖ Goal
   ─────
   • Start at raw ABI traces (tests/fixtures/*.ab1)
   • microseq trim (‑‑sanger) should:
       1. convert AB1 → FASTQ               (ab1_to_fastq.py)
       2. run BioPython sliding‑window QC   (biopy_trim.py)
       3. write a single combined FASTA     (fastq_folder_to_fasta.py)
   • We assert that the final FASTA exists and is non‑empty.

This is a smoke test: if any of the three steps break, the CLI
will exit non‑zero and/or the FASTA won’t appear, so the test fails.
"""

# ───────────────────────── imports ──────────────────────────
from pathlib import Path          # pathlib makes path maths ( / , .resolve() … ) easy & OS‑agnostic
import subprocess                 # run the CLI exactly as a user would
import sys                        # gives us sys.executable → the Python in the current env
import uuid                       # create a unique temp workdir per test run
import pytest
pytest.importorskip("pandas")
pytest.importorskip("Bio")

# pytest automatically supplies a `tmp_path` fixture: a unique temporary
# directory that is cleaned up after the test finishes.
def test_ab1_trim_to_fasta(tmp_path: Path) -> None:
    """
    AB1 → FASTQ → BioPython trim → FASTA
    """
    # create a unique sub‑folder inside pytest’s tmp dir
    work = tmp_path / f"run_{uuid.uuid4().hex}"

    # tests/fixtures/ holds four *.ab1 chromatograms shipped with the repo
    fixtures = Path(__file__).parent / "fixtures"

    # ───── invoke the CLI (subprocess) ────────────────────
    result = subprocess.run(
        [
            sys.executable,             # same Python that’s running pytest
            "-m", "microseq_tests.microseq",
            "--workdir", str(work),
            "trim",                     # sub‑command
            "-i", str(fixtures),        # input folder of AB1
            "--sanger",                 # switch on Sanger mode
            # no -o: workdir layout writes …/qc/trimmed.fastq(a)
        ],
        capture_output=True,            # collect stdout / stderr for debugging
        text=True,                      # decode bytes → str
        check=True,                     # raise CalledProcessError if exit≠0
    )

    # ───── assertion section ──────────────────────────────
    fasta = work / "qc" / "trimmed.fasta"

    # If CLI crashed, subprocess.run above would already have raised.
    # Here we make sure the expected artifact exists **and** isn't empty.
    assert fasta.exists(), (
        f"FASTA {fasta} not created.\n"
        f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"
    )

    # 100 bytes is arbitrary but guarantees we didn’t get an empty stub file
    assert fasta.stat().st_size > 100, (
        f"FASTA unexpectedly small ({fasta.stat().st_size} bytes)\n"
        f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}"
    )

