#!/usr/bin/env python3 

"""
add_taxonomy.py - append GG2 taxon names to a MicroSeq BLAST table. =) 

Usage case here...
python add_taxonomy.py \
        --hits data/biom/ocular_isolates.csv \
        --taxonomy ~/.microseq_dbs/gg2/taxonomy.tsv \
        --out data/biom/ocular_isolates_with_tax.csv 
"""
from __future__ import annotations 
import re 
from pathlib import Path 
import pandas as pd

def run_taxonomy_join(hits_fp: Path, taxonomy_fp: Path, out_fp: Path) -> None:
    # ------------------------------------------------------------------ #
    hits = pd.read_csv(hits_fp, sep=None, engine="python")
    
    tax = (
        pd.read_csv(taxonomy_fp, sep="\t")                # original columns
            .rename(                                        # normalise names
                columns={
                    "Feature ID": "sseqid",
                    "Taxon": "taxonomy"
                }
            )
    )
    # ------------------------------------------------------------------ #
    # sanity-check
    
    merged = hits.merge(tax, on="sseqid", how="left")
    n_unmatched = merged["taxonomy"].isna().sum()
    print(f"[add_taxonomy] {n_unmatched}/{len(merged)} rows unmatched") 

    # --- reoder for post-BLAST 
    cols   = ["sample_id", "taxonomy", "evalue", "bitscore", "pident", "qcov"]
    merged = merged[cols]
    
    sep = "\t" if out_fp.suffix.lower() in {".tsv", ".txt"} else ","
    merged.to_csv(out_fp, sep=sep, index=False)
    print(f"[add_taxonomy] wrote {out_fp}") 
