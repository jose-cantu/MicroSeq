# src/microseq_tests/utility/id_normaliser.py 

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

_ORIG_WELL_RE = re.compile(r'_(?:\d{4}-\d{2}-\d{2}_)?[A-H]\d{2}_trimmed$', re.I)

def strip_suffix_legacy(x: str) -> str:
    """Original 2024 pattern: drops _date_B07_trimmed (must have dashes)."""
    return _ORIG_WELL_RE.sub('', str(x))

NORMALISERS.update({
    'strip_suffix_legacy': strip_suffix_legacy,
})
