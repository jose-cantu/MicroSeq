#!/usr/bin/env python3 

"""
add_taxonomy.py - append GG@ taxon names to a MicroSeq BLAST table. =) 

Usage case here...
python add_taxonomy.py \
        --hits data/biom/ocular_isolates.csv \
        --taxonomy ~/.microseq_dbs/greengenes2/taxonomy.tsv \
        --out data/biom/ocular_isolates_with_tax.csv 
"""
import re 
import sys 
import argparse 
import pathlib 
import pandas as pd 

def _strip_accession(s: str) -> str } None: 
    """
    Convert an sseqid like ....
    RS-GCF-000771715.1-NZ-ADEY01000022.1-2
    -> GCF_000771715 
    """
    m = re.search(r"(GC[AF])[-_](\d+)", s)
    return f"{m.group(1)}_{m.group(2)}" if m else None

def _parse_args() -> argparse.Namespace:
    ap = argparse
