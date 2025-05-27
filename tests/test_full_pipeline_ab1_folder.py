from pathlib import Path
import pytest
pytest.importorskip("pandas")
pytest.importorskip("Bio")
import microseq_tests.pipeline as mp


def test_full_pipeline_ab1_folder(tmp_path: Path, monkeypatch) -> None:
    fixtures = Path(__file__).parent / "fixtures"
    monkeypatch.setattr(mp, "run_blast_stage", lambda *a, **k: 0)
    monkeypatch.setattr(mp, "run_add_tax", lambda *a, **k: 0)
    paths = mp.run_full_pipeline(fixtures, "gg2", out_dir=tmp_path)
    fasta = tmp_path / "reads.fasta"
    assert paths["fasta"] == fasta
    assert fasta.exists()
    assert fasta.stat().st_size > 0
