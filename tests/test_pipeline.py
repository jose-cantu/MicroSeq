# tests/test_pipeline.py 

from pathlib import Path 
import subprocess, pandas as pd, biom 

DEMO = Path(__file__).parent / "data" / "demo.fasta" 
META = Path(__file__).parent / "data" / "meta.tsv" 
TAX_PATH = Path.home() / ".microseq_dbs" / "gg2" / "taxonomy.tsv" 
# adding as small function here to not repeat subprocess.check_call 
def run(cmd: list[str]) -> None:
    """ Exectute cmd (a list of strings) and abort if it fails. """ 
    subprocess.check_call([str(c) for c in cmd]) 

# writing actual test here 
def test_postblast_pipeline(tmp_path):
    blast = tmp_path / "hits.tsv"
    tax = tmp_path / "hits_tax.tsv"
    base = tmp_path / "profile" 

    run(["microseq", "blast",
         "-i", DEMO,
         "-d", "gg2",
         "-o", blast,
         "--identity", "50",
         "--qcov", "10",
         "--max_target_seqs", "50", 
         ]) 

    # add taxnonomy step 
    run(["microseq", "add_taxonomy",
         "--hits", str(blast),
         "--taxonomy", str(TAX_PATH),  
         "--out", str(tax),
         ])

    run(["microseq", "postblast",
         "-b", tax,
         "-m", META,
         "--sample-col", "sample_id",
         "-o", str(base.with_suffix(".biom")),
         ]) 


    biom_f = base.with_suffix(".biom")
    csv_f = base.with_suffix(".csv")

    assert biom_f.is_file() and csv_f.is_file()

    # ─── BIOM table ───────────────────────────────────────────────────────
    table = biom.load_table(str(biom_f))
    assert table.shape[0] >= 1          # ≥ 1 taxon  (rows)
    assert table.shape[1] >= 2          # ≥ 2 samples (columns)

    # ─── CSV mirror ───────────────────────────────────────────────────────
    df = pd.read_csv(csv_f, index_col=0)
    assert df.shape[0] >= 1             # ≥ 1 taxon  (rows)
    assert df.shape[1] >= 2             # ≥ 2 samples (columns) 
