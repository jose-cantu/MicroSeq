# -------- src/microseq_tests/aligners/blast.py ----------

from __future__ import annotations 
from pathlib import Path
import pandas as pd 

from .base import BaseAligner 
from ._parse import parse_blast_tab 
from .schema import validate 
from microseq_tests.blast.run_blast import run_blast, BlastOptions 

class BlastAligner(BaseAligner):
    name = "blastn"
    
    def run(self,
            query_fa: str,
            db_key: str, 
            threads: int = 1,
            fast: bool = True,
            **kw,
            ) -> pd.DataFrame:
        ident_cut = kw.pop("identity_cutoff", 90.0)
        qcov_cut = kw.pop("qcov_cutoff", 90.0) 

        out_tsv = Path(query_fa).with_suffix(".blast.tsv")
        task = "megablast" if fast else "blastn" 

        # first pass 
        run_blast(
            query_fa,
            db_key,
            out_tsv,
            options=BlastOptions(task=task),
            threads=threads,
            **kw,
        )
        df = parse_blast_tab(out_tsv) 
        if not df.empty:
            df = validate(df) 
        

        if fast:
           needs_retry = df.empty or (
            (not df.empty) and (
               df.iloc[0]["pident"]   < ident_cut or
               df.iloc[0]["qcovhsp"] < qcov_cut
            )
           ) 
       
           if needs_retry:
                   run_blast(
                       query_fa,
                       db_key,
                       out_tsv,
                       options=BlastOptions(task="blastn"),
                       threads=threads,
                       **kw,
                   )
                   df = parse_blast_tab(out_tsv)
                   if not df.empty:
                       df = validate(df) 

        return df


