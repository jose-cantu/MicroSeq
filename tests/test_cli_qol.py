# tests/test_cli_qol.py

"""
CLI quality of life tests runned here 

What this test covers....
default_workdir fallback via config patch I made 
--link-raw symlink creation 
-v/-vv logging levels 
"""

import subprocess, logging, os
from pathlib import Path
import pytest
pytest.importorskip("pandas")
pytest.importorskip("Bio")
from microseq_tests.utility import utils
from microseq_tests.trimming import ab1_to_fastq

CLI = ["python", "-m", "microseq_tests.microseq"] # invoke module in-process 

def test_default_workdir(tmp_path, monkeypatch):
    monkeypatch.setattr(utils, "load_config",
                        lambda: {"default_workdir": str(tmp_path)}) 

    fasta = tmp_path / "demo.fasta"
    fasta.write_text(">a\nACGT\n")

    out = tmp_path / "hits.tsv"
    res = subprocess.run(
            CLI + ["blast", "-i", str(fasta),
                   "-d", "gg2", 
                   "-o", str(out),
                   "--identity", "50", "--qcov", "10", "--max_target_seqs", "1"],
            capture_output=True)
    assert res.returncode == 0 
    assert out.exists() 

# test function for --link-raw 
def test_link_raw_symlink(tmp_path, monkeypatch):
    # Create a fake AB1 file
    (tmp_path / "trace1.ab1").write_bytes(b"ABIF") 
    workdir = tmp_path / "run" 

    # skip actual conversion ot avoid BioPython parsing here 
    monkeypatch.setattr(ab1_to_fastq, "ab1_folder_to_fastq",
                      lambda *a, **k: None) 

    res = subprocess.run(
            CLI + ["--workdir", str(workdir),
                   "trim", "-i", str(tmp_path),
                   "--link-raw"],
            capture_output=True)
    assert res.returncode == 0 
    assert (workdir / "raw_ab1").is_symlink() 


# verbose level test function 
def test_verbose_level(tmp_path):
    fasta = tmp_path / "demo.fasta" 
    fasta.write_text("@id\nACGT\n+\n!!!!\n")

    res = subprocess.run(
    CLI + ["-vv", "--workdir", str(tmp_path),
           "blast", "-i", str(fasta),
           "-d", "gg2", "-o", str(tmp_path / "h.tsv"),
           "--identity", "50", "--qcov", "10", "--max_target_seqs", "1"],
    capture_output=True, text=True)

    assert res.returncode == 0
    assert "INFO" in res.stderr              









