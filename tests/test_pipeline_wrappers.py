import biom, pandas as pd
from pathlib import Path
from microseq_tests import pipeline as mp

DEMO = Path(__file__).parent / "data" / "demo.fasta"
META = Path(__file__).parent / "data" / "meta.tsv"
TAX  = Path.home() / ".microseq_dbs" / "gg2" / "taxonomy.tsv"

def test_e2e_wrappers(tmp_path: Path):
    # 1 trim skipped – demo is already fasta
    blast_tsv = tmp_path / "hits.tsv"
    tax_tsv   = tmp_path / "hits_tax.tsv"
    biom_out  = tmp_path / "profile.biom"

    mp.run_blast_stage(DEMO, "gg2", blast_tsv,
                       identity=50, qcov=10, max_target_seqs=50)
    mp.run_add_tax(blast_tsv, TAX, tax_tsv)
    mp.run_postblast(tax_tsv, META, biom_out)

    table = biom.load_table(str(biom_out))
    assert table.shape[0] >= 1
    assert table.shape[1] >= 2
    assert (biom_out.with_suffix(".csv")).is_file()
