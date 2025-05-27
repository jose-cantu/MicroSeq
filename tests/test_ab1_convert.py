import pathlib
import pytest
pytest.importorskip("Bio")
from microseq_tests.trimming.ab1_to_fastq import ab1_folder_to_fastq

def test_ab1_to_fastq(tmp_path: pathlib.Path):
    fixtures = pathlib.Path(__file__).parent / "fixtures" 
    fastqs = ab1_folder_to_fastq(fixtures, tmp_path)

    # should procedure one FASTQ per AB1 file 
    assert len(fastqs) == len(list(fixtures.glob("*ab1"))) 
    for fq in fastqs:
        assert fq.exists()
        assert fq.stat().st_size > 100 # tiny but non-empty 
