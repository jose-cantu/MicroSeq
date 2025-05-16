# src/microseq_tests/utility/metadata_tools.py 

from __future__ import annotations 
import pandas as pd, logging 

log = logging.getLogger(__name__) 

def _rows_identical(df: pd.DataFrame) -> bool:
    return df.nunique(dropna=False, axis=0).le(1).all() 

def resolve_duplicates(meta: pd.DataFrame, policy: str = "error") -> pd.DataFrame:
    """
    De-duplicate rows after ID normalisation. 

    policy = error -> raise on conflicts. 
             keep first -> keep first row, warn 
             merge -> keep one row if rows identical or site/media etc. differs 

    """ 
    assert policy in {"error", "keep-first", "merge"}
    dupes = meta.duplicates("SampleID", keep=False) 

    if not dupes.any(): 
        return meta 

    dup_df = meta[dupes] 
    good, bad = [], [] 
    for sid, grp in dup_df.groupby("SampleID"):
        if _rows_identical(grp):
            good.append(sid) 
        elif grp["Media"].nunique() == len(grp):
            good.append(sid) # replicates 
        else:
            bad.append(sid) # conflict 

    if policy == "error" and bad: 
        raise ValueError(f"Conflicting rows for {len(bad)} SampleID(s): " 
                         f"{bad[:5]} ...") 

    if policy == "keep-first":
        log.warning("Dropped %d duplicate rows (keep-first).", dupes.sum()) 
        return meta.drop_duplicates("SampleID", keep="first")

    # merge policy 
    log.warning("Merged %d duplicate SampleID rows.", len(good))
    merged = (meta[~dupes]
              .append(dup_df.groupby("SampleID").first())
              .reset_index(drop=True))

    if bad:
        log.warning("Conflict rows removed: %s", bad[:5]) 
    return merged 

