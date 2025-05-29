# tests/test_vsearch_aligner.py
from microseq_tests.aligners import load
from microseq_tests.aligners._parse import COLS
import pytest
from microseq_tests.aligners import load, _parse




def test_vsearch_cols():
    df = load("vsearch").run("tests/data/demo.fasta", "gg2", threads=1)
    assert list(df.columns) == COLS


@pytest.mark.parametrize("eng", ["blastn", "vsearch"])
def test_engines_shape(eng):
    df = load(eng).run("tests/data/demo.fasta","gg2",threads=1)
    assert list(df.columns) == _parse.COLS

