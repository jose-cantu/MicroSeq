"""
microseq_tests.utility.taxonomy_utils
Embed full-length taxonomy strings from a metadata DataFrame intoa BIOM Table. 
"""

from __future__ import annotations 
from typing import Mapping 
from biom import Table 
import pandas as pd 

def embed_taxonomy_from_metadata(
    tb1: Table,
    meta: pd.DataFrame,
    *,
    col: str = "gg_taxon",
    ) -> Table:
    """
    For every OTU row in tb1 find the first sample where it is present (>0),
    lookup meta.loc[sample, col], split on delimiters, and attach as the 
    'taxonomy' metadata field expected by ATIMA/QIIME.

    The function mutates the Table in place and also reutnrs it so you write 
    tb1 = embed_taxonomy_from_metadata(tb1, meta).
    """
    # --------------------- building a Lookup SampleID -> taxonomy string ------ 
    lookup: Mapping[str, str] = meta.set_index("SampleID")[col].to_dict() 

    tax_meta = {} 
    for otu in tb1.ids(axis="observation"):
        sample_idx = tb1.data(otu, axis="observation").nonzero()[0][0] # first hit 
        sample_id = tb1.ids(axis="sample")[sample_idx]
        tax_str = lookup.get(sample_id) 
        if tax_str is None: 
            raise KeyError(f"Sample {sample_id} missing {col} column") 
        tax_meta[otu] = {"taxonomy": split_tax(tax_str)} 

    tb1.add_metadata(tax_meta, axis="observation")
    return tb1 

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



