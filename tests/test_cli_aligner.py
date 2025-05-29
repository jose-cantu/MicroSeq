# tests/tests_cli_aligner.py 

from click.testing import CliRunner 
from microseq_tests.microseq import cli 

def test_cli_blastn_runs():
    res = CliRunner().invoke(cli, [
        "blast",
        "--aligner", "blastn",
        "-i", "tests/data/demo.fasta",
        "-d", "gg2",
        "-o", "/tmp/out.tsv",
    ])
    assert res.exit_code == 0
