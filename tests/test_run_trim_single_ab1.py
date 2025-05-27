from pathlib import Path
import pytest
pytest.importorskip("pandas")
pytest.importorskip("Bio")
import microseq_tests.pipeline as mp


def test_run_trim_single_ab1(tmp_path: Path) -> None:
    fixtures = Path(__file__).parent / "fixtures"
    ab1 = fixtures / "39764260.ab1"
    ret = mp.run_trim(ab1, tmp_path, sanger=True)
    assert ret == 0
    fasta = tmp_path / "qc" / "trimmed.fasta"
    assert fasta.exists()
    assert fasta.stat().st_size > 0
