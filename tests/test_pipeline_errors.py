from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("pandas")
pytest.importorskip("Bio")

from microseq_tests import pipeline


def test_paired_error_includes_suggestions(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Paired pipeline errors should surface the suggestion text."""

    infile = tmp_path / "input"
    infile.mkdir()

    fastq_dir = tmp_path / "out" / "passed_qc_fastq"
    fastq_dir.mkdir(parents=True)

    fastq_data = "@r1\nA\n+\n!\n"
    (fastq_dir / "S1_27F.fastq").write_text(fastq_data, encoding="utf-8")
    (fastq_dir / "S1_1492R.fastq").write_text(fastq_data, encoding="utf-8")

    def fake_run_trim(*_args, **_kwargs):
        return 0

    def fake_run_fastq_to_fasta(*_args, **_kwargs):
        return None

    def fake_fastq_to_paired_fastas(_fastq_dir: Path, fasta_dir: Path):
        fasta_dir.mkdir(parents=True, exist_ok=True)
        fwd = fasta_dir / "S1_27F.fasta"
        rev = fasta_dir / "S1_1492R.fasta"
        fwd.write_text(">fwd\nA\n", encoding="utf-8")
        rev.write_text(">rev\nA\n", encoding="utf-8")
        return [fwd, rev]

    monkeypatch.setattr(pipeline, "run_trim", fake_run_trim)
    monkeypatch.setattr(pipeline, "run_fastq_to_fasta", fake_run_fastq_to_fasta)
    monkeypatch.setattr(pipeline, "_fastq_to_paired_fastas", fake_fastq_to_paired_fastas)
    monkeypatch.setattr(pipeline, "assemble_pairs", lambda *_, **__: [])

    with pytest.raises(ValueError) as excinfo:
        pipeline.run_full_pipeline(
            infile,
            "nt",
            tmp_path / "out",
            mode="paired",
        )

    message = str(excinfo.value)
    assert "detected 1 forward and 1 reverse reads" in message
    assert "Example flag" in message
