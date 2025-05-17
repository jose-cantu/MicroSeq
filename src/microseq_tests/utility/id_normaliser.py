# src/microseq_tests/utility/id_normaliser.py 

# import re 
#
# # specific regex for cleaning up sample_id names to match BIOM for ATIMA 
#
# _SUFFIX_RE = re.compile(r"_(?:(?:\d{4}-\d{2}-\d{2})|\d{8})?_[A-H]\d{2}_trimmed$")
#
# def none(x: str) -> str: 
#     return x 
#
# def strip_suffix(x: str) -> str: 
#     """Remove data + well + '_trimmed' tail from a sample name.""" 
#     return _SUFFIX_RE.sub("", x) 
#
# NORMALISERS = {"none": none, "strip_suffix": strip_suffix}
#
# # r"_(?:\d{4}-\d{2}-\d{2}_)?[A-H]\d{2}_trimmed$"
#
# def strip_suffix_simple(x: str) -> str:
#     """Always drop the trailing _B07_trimmed (date optional)."""
#      return _SUFFIX_RE.sub("", x)

# NORMALISERS.update({"strip_suffix_simple": strip_suffix_simple})
# src/microseq_tests/utility/id_normaliser.py
# src/microseq_tests/utility/id_normalisers.py
import re

# ── regexes ────────────────────────────────────────────────
_WELL_RE   = re.compile(r'_(?:(?:\d{4}-\d{2}-\d{2})|\d{8})?_[A-H]\d{2}_trimmed$')
_PRIMER_RE = re.compile(r'-\d{4}R.*$')          # catches 1492R, 1429R, …

# ── helpers ────────────────────────────────────────────────
def none(x: str) -> str: return x

def strip_suffix(x: str) -> str:
    """Drop _well_trimmed and trailing -1492R… primer chunk."""
    x = _WELL_RE.sub('', str(x))
    return _PRIMER_RE.sub('', x)

# keep the simple one around because other pipelines already use it
def strip_suffix_simple(x: str) -> str:
    return _WELL_RE.sub('', str(x))

NORMALISERS = {
    'none': none,
    'strip_suffix': strip_suffix,          # ← NOTE: now the *two-step* version
    'strip_suffix_simple': strip_suffix_simple,
}

