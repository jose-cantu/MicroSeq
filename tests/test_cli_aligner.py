# tests/tests_cli_aligner.py 

import subprocess, sys, pathlib, tempfile

def test_cli_blastn_runs():
    out = pathlib.Path(tempfile.gettempdir()) / "hits.tsv"
    res = subprocess.run([
        sys.executable, "-m", "microseq_tests.microseq",
        "blast", "-i", "tests/data/demo.fasta", "-d", "gg2",
        "-o", out, "--aligner", "blastn"
    ])
    assert res.returncode == 0
    assert out.exists()

