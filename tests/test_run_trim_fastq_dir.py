from __future__ import annotations

import pytest
from pathlib import Path
import microseq_tests.pipeline as mp
import microseq_tests.trimming.quality_trim as qt

# --- Test Suite Summary ---
# This pytest suite validates the 'trim_fastq_inputs' and 'run_trim' functions from the
# 'microseq-tests' package. It uses monkeypatching to mock external dependencies
# (like the actual Trimmomatic execution and file conversions) to test the Python logic
# in isolation. Key behaviors tested include handling directories of files, creating
# combined output, generating summary statistics, cleaning up temporary files,
# and ensuring functions call their appropriate helpers.
# --------------------------


# Use pytest markers to ensure required optional dependencies (pandas, Biopython) are installed
pytest.importorskip("pandas") 
pytest.importorskip("Bio")

def test_trim_fastq_inputs_directory_fastq(monkeypatch, tmp_path):
    """
    Tests the 'trim_fastq_inputs' function when given an input directory 
    containing multiple FASTQ files. 

    Mocks the 'quality_trim' function to simulate trimming behavior w/o running 
    external tools. 

    Key things:
        1) All input files are processed and combined into a single output file.
        2) Temporary trimmed files are correctly unlinked/deleted after merging. 
        3) A summary TSV file is correctly generated, containing stats for each indivisual file,
        a header row, and a final '_combined' row with the aggreagted 
        weighted averages. 
        4) Running the function multiple times appends to the summary file w/o 
        duplicating the header. 
        5) The final output file is located in the specified 'trim_dir' 
    """
    input_dir = tmp_path / "input"  # isolated fixture root for inputs 
    input_dir.mkdir() # create a test directory 

    fastq_content = "@r\nACGT\n+\n!!!!\n" # 1-read FASTQ -> Q=0 per base 
    for name in ("a.fastq", "b.fastq"): # two sources to exercise concat + per-file 
        (input_dir / name).write_text(fastq_content)

    trimmed_calls: list[tuple[Path, Path]] = [] # remember temp outputs to assert cleanup 

    def fake_quality_trim(src: Path, dest: Path, **_: object) -> Path: 
        dest = Path(dest) 
        dest.write_text(fastq_content)
        trimmed_calls.append((Path(src), dest))  
        return dest

    monkeypatch.setattr(qt, "quality_trim", fake_quality_trim) # isolate from Trimmomatic subprocess 

    summary = tmp_path / "summary.tsv" 
    out = qt.trim_fastq_inputs(input_dir, tmp_path / "qc", summary_tsv=summary) # invoke unit test 

    assert out == tmp_path / "qc" / "trimmed.fastq" 
    combined = out.read_text() 
    assert combined.count("@r") == 2 # both records concatenate  

    # temporary trimmed files are cleaned up once merged 
    assert all(not dest.exists() for _, dest in trimmed_calls) 

    # Summary rows per source file plus combined aggregate with header 
    summary_lines = summary.read_text().strip().splitlines() 
    assert summary_lines[0] == "file\treads\tavg_len\tavg_q\tavg_mee\tavg_mee_per_kb\tavg_qeff\tmee_qc_label" 

    # Because inputs are lexicographically sorted ("a", "b") etc. then just assert exact row order 
    assert summary_lines[1] == "a.fastq\t1\t4.0\t0.00\t4.000\t1000.000\t0.00\treview"
    assert summary_lines[2] == "b.fastq\t1\t4.0\t0.00\t4.000\t1000.000\t0.00\treview" # per file row 
   
    # Global weighted row: R=2, B=1*4 + 1*4 = 8 ⇒ L̄=8/2=4.0; Q̄ = (0*8)/8 = 0.00 
    assert summary_lines[3] == "_combined\t2\t4.0\t0.00\t4.000\t1000.000\t0.00\treview" 

    # Creating combined file does not recreate source fastqs 
    assert sorted(p.name for p in input_dir.iterdir()) == ["a.fastq", "b.fastq"] 

    # Combine outputs in qc directory 
    assert out.parent == tmp_path / "qc" 

    # Running again appends rows w/o duplicating header 
    qt.trim_fastq_inputs(input_dir, tmp_path / "qc", summary_tsv=summary) 
    summary_lines = summary.read_text().strip().splitlines() 
    assert summary_lines.count("file\treads\tavg_len\tavg_q\tavg_mee\tavg_mee_per_kb\tavg_qeff\tmee_qc_label") == 1 
    assert summary_lines.count("_combined\t2\t4.0\t0.00\t4.000\t1000.000\t0.00\treview") == 2

    # verify concatenation order is stable 
    assert combined.startswith("@r") and combined.count("@r") == 2 

def test_run_trim_uses_helper(monkeypatch, tmp_path):
    input_dir = tmp_path / "input" 
    input_dir.mkdir() 
    (input_dir / "a.fastq").write_text("@r\nACGT\n+\n!!!!\n")

    helper_calls: list[tuple[Path, Path, Path | None]] = [] 

    def fake_trim_fastq_inputs(input_path: Path, trim_dir: Path, summary_tsv=None):
        summary_arg = Path(summary_tsv) if summary_tsv is not None else None 
        helper_calls.append((Path(input_path), Path(trim_dir), summary_arg)) 
        trim_dir = Path(trim_dir) 
        trim_dir.mkdir(parents=True, exist_ok=True) 
        out = trim_dir / "trimmed.fastq" 
        out.write_text("data") 
        return out 

    monkeypatch.setattr(mp, "trim_fastq_inputs", fake_trim_fastq_inputs)

    fasta_calls: list[tuple[Path, Path]] = [] 

    def fake_fastq_to_fasta(src_dir: Path, out_fasta: Path) -> Path:
        out_fasta = Path(out_fasta)
        out_fasta.write_text(">seq\nACGT\n")
        fasta_calls.append((Path(src_dir), out_fasta))
        return out_fasta 

    # FASTA conversion triggered on qc directory 
    monkeypatch.setattr(mp, "fastq_to_fasta", fake_fastq_to_fasta)

    summary = tmp_path / "summary.tsv" 

    rc = mp.run_trim(input_dir, tmp_path, summary_tsv=summary)

    assert rc == 0 
    assert helper_calls == [(input_dir, tmp_path / "qc", summary)]
    assert fasta_calls == [(tmp_path / "qc", tmp_path / "qc" / "trimmed.fasta")] 

def test_trim_fastq_inputs_with_mee_telemetry(monkeypatch, tmp_path):
    input_dir = tmp_path / "input"
    input_dir.mkdir()

    raw_content = (
        "@good\nACGT\n+\nIIII\n"
        "@bad\nACGT\n+\n!!!!\n"
    )
    (input_dir / "a.fastq").write_text(raw_content)

    def fake_quality_trim(src: Path, dest: Path, **_: object) -> Path:
        dest = Path(dest)
        dest.write_text(raw_content)
        return dest

    monkeypatch.setattr(qt, "quality_trim", fake_quality_trim)

    summary = tmp_path / "summary.tsv"
    out = qt.trim_fastq_inputs(
        input_dir,
        tmp_path / "qc",
        summary_tsv=summary,
    )

    combined = out.read_text()
    assert combined.count("@good") == 1
    assert combined.count("@bad") == 1

    summary_lines = summary.read_text().strip().splitlines()
    assert summary_lines[0] == "file\treads\tavg_len\tavg_q\tavg_mee\tavg_mee_per_kb\tavg_qeff\tmee_qc_label"
    assert summary_lines[1] == "a.fastq\t2\t4.0\t20.00\t2.000\t500.050\t3.01\treview"
    assert summary_lines[2] == "_combined\t2\t4.0\t20.00\t2.000\t500.050\t3.01\treview"

