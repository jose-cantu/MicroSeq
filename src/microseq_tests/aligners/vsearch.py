# src/microseq_tests/aligners/vsearch.py ---------------

from __future__ import annotations 
from pathlib import Path
import subprocess
import os
import pandas as pd
from .base import BaseAligner
from ._parse import parse_blast_tab
from .schema import validate
from microseq_tests.utility.utils import load_config, expand_db_path

class VsearchAligner(BaseAligner):
    name = "vsearch"

    def run(self, query_fa: str, db_key: str, 
            threads: int = 1, **kw) -> pd.DataFrame:

        out_tsv = Path(query_fa).with_suffix(".vsearch.tsv")

        cfg = load_config()
        try:
            tmpl = cfg["databases"][db_key]["blastdb"]
        except KeyError as e:
            valid = ", ".join(cfg.get("databases", {}).keys())
            raise KeyError(
                f"[VsearchAligner] unknown DB key '{db_key}'. Valid keys: {valid}"
            ) from e

        db_fa = expand_db_path(tmpl) + ".fasta"
        if not Path(db_fa).is_file():
            raise FileNotFoundError(
                f"Database FASTA not found: {db_fa}"
            )

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




