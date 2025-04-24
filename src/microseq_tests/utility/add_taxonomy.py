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

_ACC_RE = re.compile(r"(GC[AF])[-_](\d+)")           # GCF_012345678


def strip_accession(s: str) -> str | None:
    """
    RS-GCF-000771715.1-NZ-ADEY01000022.1-2  ->  GCF_000771715
    GCF-014326775.1-NZ-CP053957.1--5        ->  GCF_014326775
    """
    m = _ACC_RE.search(s)
    return f"{m.group(1)}_{m.group(2)}" if m else None 


def run_taxonomy_join(hits_fp: Path, taxonomy_fp: Path, out_fp: Path) -> None:
    # ------------------------------------------------------------------ #
    hits = pd.read_csv(hits_fp, sep=None, engine="python")
    
    tax = (
        pd.read_csv(taxonomy_fp, sep="\t")                # original columns
            .rename(                                        # normalise names
                columns={
                    "Feature ID":  "accession",
                    "seqid":       "accession",   # future-proof
                    "accession_id":"accession",
                    "Taxon":       "taxon_name",
                    "lineage":     "taxon_name",
                }
            )
    )
    # ------------------------------------------------------------------ #
    # sanity-check
    required = {"accession", "taxon_name"}
    missing  = required - set(tax.columns)
    if missing:
        raise ValueError(f"[add_taxonomy] taxonomy.tsv missing {missing}")
        
    if "sseqid" not in hits:
        raise ValueError("hits file lacks 'sseqid' column")
        
    hits["accession"] = hits["sseqid"].map(strip_accession)
    
    merged = hits.merge(
        tax[["accession", "taxon_name"]],
        on="accession",
        how="left",
    )
    
    # report unmatched rows (optional)
    n_missing = merged["taxon_name"].isna().sum()
    if n_missing:
        print(f"[add_taxonomy] {n_missing}/{len(merged)} rows unmatched")
        
    # rename + reorder for post-BLAST
    merged = merged.rename(columns={"taxon_name": "taxonomy"})
    cols   = ["sample_id", "taxonomy", "evalue", "bitscore", "pident", "qcov"]
    merged = merged[cols]
    
    sep = "\t" if out_fp.suffix.lower() in {".tsv", ".txt"} else ","
    merged.to_csv(out_fp, sep=sep, index=False)
