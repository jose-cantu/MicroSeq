"""
microseq_tests.utility.taxonomy_utils
Embed full-length taxonomy strings from a metadata DataFrame intoa BIOM Table. 
"""

from __future__ import annotations 
from typing import Mapping, Dict
from biom import Table 
import pandas as pd
import re   

def parse_lineage(lineage: str, fmt: str = "auto") -> Dict[str, str]:
    """Return a mapping of the canonical seven ranks."""
    if not isinstance(lineage, str):
        lineage = ""

    if fmt in {"auto", "gg2", "silva"}:
        parts = [p.split("__", 1)[-1] for p in lineage.rstrip(";").split(";")]
    else:  # ncbi or unknown -> assume no prefixes
        parts = [p.strip() for p in lineage.rstrip(";").split(";")]

    parts = [p.split(" strain", 1)[0] if "strain" in p else p for p in parts]
    parts = [p.strip() for p in parts]
    parts = (parts + [""] * 7)[:7]

    cols = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    return {c: v for c, v in zip(cols, parts)}


def embed_taxonomy_from_metadata(
    tbl: Table,
    df: pd.DataFrame,
    *,
    col: str = "taxonomy",   # column in df that already holds the lineage strings
    fmt: str = "auto",
) -> Table:
    """
    Embed a seven-rank taxonomy dict *and* the BIOM-spec 'taxonomy' list
    onto every observation in `tbl`.

    Assumes each observation-ID is the full lineage string (as created by
    pivot_table(index="taxonomy", …) earlier in the pipeline).
    """
    # Build lookup   lineage-string → {'Domain': 'Bacteria', …}
    lut = {row[col]: parse_lineage(row[col], fmt) for _, row in df.iterrows()}

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



