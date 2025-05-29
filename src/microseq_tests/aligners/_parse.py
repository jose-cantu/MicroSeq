# ------ src/microseq_tests/aligners/_parse.py -------------

from __future__ import annotations 
import pandas as pd 
from pathlib import Path 

COLS ="qseqid sseqid pident qlen qcovhsp length evalue bitscore stitle".split()

def parse_blast_tab(tsv: Path) -> pd.DataFrame:
    """Read BLAST/VSEARCH tabular output with typed columns."""
    df = pd.read_csv(tsv, sep="\t", names=COLS, comment="#", dtype=str)

    if not df.empty and df.iloc[0, 0] == "qseqid":
        df = df.iloc[1:]

    num_cols = ["pident", "qlen", "qcovhsp", "length", "evalue", "bitscore"]
    for col in num_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df[COLS].reset_index(drop=True)
 



