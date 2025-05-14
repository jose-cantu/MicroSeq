#!/usr/bin/env python3

# src/microseq_tests/utility/add_taxonomy.py 

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
from microseq_tests.utility.io_utils import normalise_tsv 

# -------- helpers -----------------------------------------------------------

# Extract plain genome accession (GCF_012345678) from RS-GCF-... strings
_ACC_RE = re.compile(r"(GC[AF])[-_](\d+)")
def _strip_accession(s: str) -> str | None:
    m = _ACC_RE.search(str(s))
    return f"{m.group(1)}_{m.group(2)}" if m else None




def run_taxonomy_join(hits_fp: Path, taxonomy_fp: Path, out_fp: Path, fill_species: bool = False) -> None: 
    # ------------------------------------------------------------------ #
    # added a new revision in here auto-convert space or comma-delimited files to real TABs moving forward
    sep = "\t" if hits_fp.suffix.lower() == ".tsv" else None
    hits = pd.read_csv(
            normalise_tsv(hits_fp),
            sep=sep,
            engine="python", # auto-deatect header row 
            dtype=str, 
            ).rename(columns={"qseqid": "sample_id"})
   # canonicalize the subject ID so it can match the taxonomy table for each of the databases 
    hits["sseqid"] = ( 
       hits["sseqid"].astype(str)
           .str.rstrip("|") 
           .str.split("|").str[-1] # drop gi|.....|ref| or silva|..... 
           .str.split(" ").str[0] # drop everything after the first space
           .str.replace(r"\.\d+\.\d+$", "", regex=True)   # Silva special case           
                      ) 
    tax = (
        pd.read_csv(normalise_tsv(taxonomy_fp), sep="\t")                # original columns
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
    
    if fill_species:
        needs_fill = (merged["taxonomy"].str.count(";") < 6) & \
                     (merged["pident"].astype(float) >= 97.0) 

        species = (
            merged.loc[needs_fill, "stitle"]
                  .str.split(";").str[-1]
                  .str.strip()
                  .str.replace(" ", "_") 
                  )
        
        merged.loc[needs_fill, "taxonomy"] = (
            merged.loc[needs_fill, "taxonomy"].str.rstrip(";") + ";" + species 
            ) 

        # sanity-check 
        assert (merged.loc[needs_fill, "taxonomy"].str.count(";") >= 6).all() 

    # keep every original blast column, then tack taxonomy on the end 
    blast_cols = hits.columns.tolist() # order from input file 
    merged = merged[blast_cols + ["taxonomy"]] 

    out_sep = "\t" if out_fp.suffix.lower() in {".tsv", ".txt"} else ","
    merged.to_csv(out_fp, sep=out_sep, index=False)
    print(f"[add_taxonomy] wrote {out_fp}") 
