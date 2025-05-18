"""
microseq_tests.utility.taxonomy_utils
Embed full-length taxonomy strings from a metadata DataFrame intoa BIOM Table. 
"""

from __future__ import annotations 
from typing import Mapping, Dict, List, Union
from biom import Table 
import pandas as pd
import re   

def parse_lineage(
    line: str,
    fmt: str = "auto",
    *,
    as_dict: bool = False,
) -> Union[List[str], Dict[str, str]]:
    """Return the canonical seven taxonomic ranks from ``line``.

    Parameters
    ----------
    line
        The raw lineage string to parse.
    fmt
        Expected format of the lineage. ``"auto"``/``"gg2"``/``"silva"`` will
        strip ``d__`` style prefixes; anything else is treated as already
        prefix-free.
    as_dict
        If ``True`` return a mapping of rank name to value, otherwise return a
        list of rank values in order.
    """

    if not isinstance(line, str):
        line = ""

    if fmt in {"auto", "gg2", "silva"}:
        parts = [p.split("__", 1)[-1] for p in line.rstrip(";").split(";")]
    else:  # ncbi or unknown -> assume no prefixes
        parts = [p.strip() for p in line.rstrip(";").split(";")]

    parts = [
        p.split(" strain", 1)[0] if "strain" in p else p or "Unclassified"
        for p in parts
    ]
    parts = (parts + ["Unclassified"] * 7)[:7]

    ranks = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    if as_dict:
        return dict(zip(ranks, parts))
    return parts


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
    lut = {row[col]: parse_lineage(row[col], as_dict=True) for _, row in df.iterrows()}

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



