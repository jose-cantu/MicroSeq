# tests/test_vsearch_aligner.py
from microseq_tests.aligners import load, _parse
from microseq_tests.aligners._parse import COLS
import pytest
import shutil
from pathlib import Path
from microseq_tests.utility.utils import load_config, expand_db_path

if shutil.which("vsearch") is None:
    pytest.skip("vsearch not installed", allow_module_level=True)

cfg = load_config()
tmpl = cfg["databases"]["gg2"]["blastdb"]
DB_FASTA = Path(expand_db_path(tmpl) + ".fasta")
if not DB_FASTA.is_file():
    pytest.skip(f"Database FASTA not found: {DB_FASTA}", allow_module_level=True)




def test_vsearch_cols():
    df = load("vsearch").run("tests/data/demo.fasta", "gg2", threads=1)
    assert list(df.columns) == COLS


@pytest.mark.parametrize("eng", ["blastn", "vsearch"])
def test_engines_shape(eng):
    df = load(eng).run("tests/data/demo.fasta","gg2",threads=1)
    assert list(df.columns) == _parse.COLS

