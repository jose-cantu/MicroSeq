# src/microseq_tests/utility/id_normalisers.py 

import re 

# specific regex for cleaning up sample_id names to match BIOM for ATIMA 

_SUFFIX_RE = re.compile(r"_(?:\d{4}-\d{2}-\d{2}_)?[A-H]\d{2}_trimmed$")

def none(x: str) -> str: 
    return x 

def strip_suffix(x: str) -> str: 
    """Remove data + well + '_trimmed' tail from a sample name.""" 
    return _SUFFIX_RE.sub("", x) 

NORMALISERS = {"none": none, "strip_suffix": strip_suffix} 
