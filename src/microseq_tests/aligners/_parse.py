# ------ src/microseq_tests/aligners/_parse.py -------------

from __future__ import annotations 
import pandas as pd 
from pathlib import Path 

COLS ="qseqid sseqid pident qlen qcovhsp length evalue bitscore stitle".split()

def parse_blast_tab(tsv: Path) -> pd.DataFrame:
    df = pd.read_csv(tsv, sep="\t", names = COLS, comment="#", skiprows=1)

    return df[COLS] 
 



