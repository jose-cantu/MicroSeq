import pytest
pytest.importorskip("pandas")
import pandas as pd
from pathlib import Path

biom = pytest.importorskip("biom")
from microseq_tests import pipeline as mp

DEMO = Path(__file__).parent / "data" / "demo.fasta"
META = Path(__file__).parent / "data" / "meta.tsv"
TAX  = Path.home() / ".microseq_dbs" / "gg2" / "taxonomy.tsv"

def test_e2e_wrappers(tmp_path: Path):
    # 1 trim skipped – demo is already fasta
    blast_tsv = tmp_path / "hits.tsv"
    tax_tsv   = tmp_path / "hits_tax.tsv"
    biom_out  = tmp_path / "profile.biom"

    # synthetic BLAST hits (avoids running the real BLAST binary)
    blast_tsv.write_text(
        "qseqid\tsseqid\tpident\tevalue\tbitscore\n"
        "S1\tGG2_0001\t99.0\t1e-10\t200\n"
        "S2\tGG2_0002\t98.0\t1e-9\t180\n"
    )

    tax = tmp_path / "taxonomy.tsv"
    tax.write_text(
        "sseqid\ttaxonomy\n"
        "GG2_0001\tk__Bacteria; p__Firmicutes\n"
        "GG2_0002\tk__Bacteria; p__Proteobacteria\n"
    )

    mp.run_add_tax(blast_tsv, tax, tax_tsv)
    mp.run_postblast(tax_tsv, META, biom_out, taxonomy_col="taxonomy")

    table = biom.load_table(str(biom_out))
    assert table.shape[0] >= 1
    assert table.shape[1] >= 2
    assert (biom_out.with_suffix(".csv")).is_file()
