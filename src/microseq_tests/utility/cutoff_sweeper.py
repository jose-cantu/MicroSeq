# src/microseq_tests/utility/cutoff_sweeper.py 

"""
cutoff_sweeper.py
Explore identity / coverage thresholds that yield ~N PASS rows
from the hits_full / all_hits TSV written by microseq blast --relaxed.
"""

from __future__ import annotations
from pathlib import Path
from itertools import product
import pandas as pd

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
    df = pd.read_csv(path, sep="\t", usecols=["best_pident", "best_qcov"])
    best = df[["best_pident", "best_qcov"]].to_numpy()

    results: list[tuple[int, int, int]] = []
    for ident, qcov in product(range(id_min, id_max + 1, step),
                               range(cov_min, cov_max + 1, step)):
        mask = (best[:, 0] >= ident) & (best[:, 1] >= qcov)
        count = int(mask.sum())
        results.append((ident, qcov, count))

    # sort by |count - target|  then by higher identity (tie-break)
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


