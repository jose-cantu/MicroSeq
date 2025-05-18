"""
microseq_tests.utility.taxonomy_utils
Embed full-length taxonomy strings from a metadata DataFrame intoa BIOM Table. 
"""

from __future__ import annotations 
from typing import Mapping, Dict 
from biom import Table 
import pandas as pd
import re   

def parse_lineage(lineage: str) -> Dict[str, str]:
    """
    Parse a semicolon-delimited lineage like
        d__Bacteria; p__Firmicutes; c__Bacilli; …
    into the seven canonical ranks.  Missing ranks → ''.
    """
    lineage = lineage.strip()
    rank_specs = [
        ("Domain",  r"d__"),
        ("Phylum",  r"p__"),
        ("Class",   r"c__"),
        ("Order",   r"o__"),
        ("Family",  r"f__"),
        ("Genus",   r"g__"),
        ("Species", r"s__"),
    ]
    out: Dict[str, str] = {}
    for col, prefix in rank_specs:
        m = re.search(re.escape(prefix) + r"([^;]+)", lineage)
        out[col] = m.group(1).strip() if m and m.group(1).strip() else ""
    return out


def embed_taxonomy_from_metadata(
    tbl,
    df,
    *,
    col: str = "taxonomy",   # column in df that already holds the lineage strings
):
    """
    Embed a seven-rank taxonomy dict *and* the BIOM-spec 'taxonomy' list
    onto every observation in `tbl`.

    Assumes each observation-ID is the full lineage string (as created by
    pivot_table(index="taxonomy", …) earlier in the pipeline).
    """
    # Build lookup   lineage-string → {'Domain': 'Bacteria', …}
    lut = {row[col]: parse_lineage(row[col]) for _, row in df.iterrows()}

    def _obs_meta(obs_id):
        parsed = lut.get(obs_id, {})
        return {"taxonomy": list(parsed.values()), **parsed}

    tbl.add_metadata(_obs_meta, axis="observation", inplace=True)
    return tbl


# small helper to unit-test independently 
def split_tax(tax_string: str) -> list[str]:
    """
    Robust splitter: accepts semicolon, pipe, tab, or space delimiters and strips extra 
    whitespace in the process. 
    """
    if ";" in tax_string:
        parts = tax_string.split(";")
    elif "|" in tax_string:
        parts = tax_string.split("|")
    else:
        parts = tax_string.split() 

    return [p.strip() for p in parts if p.strip()] 



