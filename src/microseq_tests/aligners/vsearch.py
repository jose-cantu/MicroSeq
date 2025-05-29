# src/microseq_tests/aligners/vsearch.py ---------------

from __future__ import annotations 
from pathlib import Path 
import subprocess, os, pandas as pd 
from .base import BaseAligner 
from ._parse import parse_blast_tab 
from .schema import validate 

class VsearchAligner(BaseAligner):
    name = "vsearch"

    def run(self, query_fa: str, db_key: str, 
            threads: int = 1, **kw) -> pd.DataFrame:

        out_tsv = Path(query_fa).with_suffix(".vsearch.tsv")
        db_fa = f"{os.environ['MICROSEQ_DB_HOME']}/{db_key}/greengenes2_db.fasta"

        subprocess.run([
            "vsearch", "--usearch_global", query_fa,
            "--db", db_fa,
            "--blast6out", out_tsv,
            "--id", "0.8",
            "--threads", str(threads),
        ], check=True)

        df = parse_blast_tab(out_tsv)
        if not df.empty:
            df["qcovhsp"] = 100 * df["length"] / df["qlen"] 
            df = validate(df) 

        return df 




