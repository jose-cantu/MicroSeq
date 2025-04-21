# microseq_tests/tests/test_postblast.py 

from pathlib import Path 
from microseq_tests.post_blast_analysis import run
import pandas as pd 
from biom import load_table 


def test_postblast(tmp_path: Path):
    work = tmp_path 
    blast = Path(__file__).with_suffix("").parent / "fixtures" / "blast.tsv"
    meta = Path(__file__).with_suffix("").parent / "fixtures" / "meta.tsv" 
    out_biom = work / "sample.biom" 

    run(blast, meta, out_biom, write_csv=True) 

    csv = out_biom.with_suffix(".csv")
    assert out_biom.exists() and csv.exists() 

    df = pd.read_csv(csv, index_col=0)
    table = load_table(str(out_biom))
    assert table.shape[1] == 1 # one column 
    assert df.shape[0] == table.shape[0] # same taxa rows 



