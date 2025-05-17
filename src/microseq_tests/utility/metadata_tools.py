# src/microseq_tests/utility/metadata_tools.py 

from __future__ import annotations
import pandas as pd, logging

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------
def _rows_identical(df: pd.DataFrame) -> bool:
    """Return True if every row in *df* is exactly identical."""
    return df.nunique(dropna=False, axis=0).le(1).all()

# ---------------------------------------------------------------------
def resolve_duplicates(
    meta: pd.DataFrame,
    policy: str = "error",
    col: str = "SampleID",
) -> pd.DataFrame:
    """
    Deduplicate *meta* on the given *col*.

    policy
    ------
    error       → raise if non-identical duplicates are found  
    keep-first  → keep the first row, drop the rest (warn)  
    merge       → keep one row if all duplicates are identical or differ
                  only in Media/replicate-style columns; drop conflicts

    The caller (postblast) decides *col* by the --sample-col flag so the
    same logic works whether you match on 'SampleID', 'sample_id', etc.
    """
    assert policy in {"error", "keep-first", "merge"}

    dup_mask = meta.duplicated(subset=col, keep=False)      # fixed call
    if not dup_mask.any():
        return meta

    dup_df = meta.loc[dup_mask]
    good, bad = [], []

    for sid, grp in dup_df.groupby(col):
        if _rows_identical(grp):
            good.append(sid)
        elif grp["Media"].nunique() == len(grp):
            good.append(sid)                   # technical replicates
        else:
            bad.append(sid)                    # conflicting rows

    # ── policy handlers ─────────────────────────────────────────────
    if policy == "error" and bad:
        raise ValueError(
            f"Conflicting rows for {len(bad)} {col}(s): {bad[:5]} ..."
        )

    if policy == "keep-first":
        log.warning("Dropped %d duplicate rows (keep-first).", dup_mask.sum())
        return meta.drop_duplicates(subset=col, keep="first")

    # merge
    log.warning("Merged %d duplicate %s rows.", len(good), col)
    merged = (
        meta.loc[~dup_mask]
        .append(dup_df.groupby(col).first())
        .reset_index(drop=True)
    )
    if bad:
        log.warning("Conflict rows removed: %s", bad[:5])
    return merged


