# src/microseq_tests/utility/cutoff_sweeper.py 

"""
cutoff_sweeper.py
Explore identity / coverage thresholds that yield ~N PASS rows
from the hits_full / all_hits TSV written by microseq blast --relaxed.
Also created a variant for postblast predicting what survies after threshold cutoff. 
"""

from __future__ import annotations
from pathlib import Path
from itertools import product
import pandas as pd


# ---------- helper function -----------

def _load_hits(path: Path | str, 
               meta_cols: list[str] | None,
               need_bitscore: bool):
    """Return a DataFrame with columns:
        pident, qcovhsp, [bitscore], + meta_cols (renamed & uniform) """
    cols = pd.read_csv(path, sep="\t", nrows=0).columns 
    pid_col = "best_pident" if "best_pident" in cols else "pident"
    qcov_col = "best_qcov" if "best_qcov" in cols else "qcovhsp" 

    use = [pid_col, qcov_col]
    if need_bitscore:
        use.append("bitscore")
    if meta_cols:
        use += meta_cols 

    df = pd.read_csv(path, sep="\t", usecols=use)
    df = df.rename(columns={pid_col: "pident", qcov_col: "qcovhsp"})
    return df 


def suggest(path: str | Path,
            target: int,
            *,
            id_min: int = 80,
            id_max: int = 100,
            cov_min: int = 0,
            cov_max: int = 100,
            step: int = 1,
            top: int = 10,
) -> list[tuple[int, int, int]]:
    """
    Return the `top` (identity, qcov, pass_count) triples whose pass_count
    is closest to `target`.

    The search space:
        identity ∈ [id_min, id_max]
        qcov     ∈ [cov_min, cov_max]
        step     : both stepped by `step` %
    """
    path = Path(path)
    df = _load_hits(path, meta_cols=None, need_bitscore=False)
    best = df[["pident", "qcovhsp"]].to_numpy()

    results: list[tuple[int, int, int]] = []
    for ident, qcov in product(range(id_min, id_max + 1, step),
                               range(cov_min, cov_max + 1, step)):
        mask = (best[:, 0] >= ident) & (best[:, 1] >= qcov)
        count = int(mask.sum())
        results.append((ident, qcov, count))

    # sort by |count - target|  then by higher identity (tie-break)
    results.sort(key=lambda t: (abs(t[2] - target), -t[0], -t[1]))
    return results[:top]



# Collapse-aware variant that predicts rows/samples that survive MicroSeq postblast
def suggest_after_collapse(path, meta_cols, target, *, step=1,
                           id_min=80, id_max=100,
                           cov_min=0, cov_max=100, top=10):
    """
    Predict identity / qcov pairs that yield ≈ target PASS rows
    after collapsing to one hit per sample (postblast logic). Returns the top (identity, qcov, pass_count) tripples whose PASS samples/rows after per-sample collapse are closest to desired target value. 

    Parameters
    ----------
    path : str or Path
        hits_full.tsv written by `microseq blast --relaxed`
    meta_cols : list[str]
        Column(s) that define a unique sample, usually ["sample_id"].
    target : int
        Desired number of PASS rows after collapse.
    step : int or float
        Grid step size in percent (default 1).
    Returns
    -------
    list[tuple[int, int, int]]
        (identity, qcov, pass_count) sorted by |pass_count-target|.
    """
    # Load only nmeeded columns here 
    df = _load_hits(path, meta_cols=meta_cols, need_bitscore=True)  

    # keep highest identity (then the highest bitscore == best hit) per sample 
    best = (df.sort_values(["pident", "bitscore"], ascending=False).drop_duplicates(subset=meta_cols, keep="first")) 

    results = []
    for ident in range(id_min, id_max + 1, step):
        for qcov in range(cov_min, cov_max + 1, step):
            mask = (best.pident >= ident) & (best.qcovhsp >= qcov) 
            results.append((ident, qcov, int(mask.sum()))) 

    results.sort(key=lambda t: (abs(t[2] - target), -t[0], -t[1]))
    return results[:top] 


if __name__ == "__main__":       # quick CLI use:
    import argparse, textwrap
    pa = argparse.ArgumentParser(
        description="Suggest identity / qcov pairs to hit a target PASS count",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""
            Example:
            python cutoff_sweeper.py blast/all_hits.tsv 530
        """))
    pa.add_argument("table", help="Path to all_hits / hits_full TSV")
    pa.add_argument("target", type=int, help="Desired PASS count")
    pa.add_argument("--step", type=int, default=1, help="Percent step size")
    args = pa.parse_args()

    for ident, qcov, n in suggest(args.table, args.target, step=args.step):
        print(f"{ident:>3}%  {qcov:>3}%   {n} PASS")


