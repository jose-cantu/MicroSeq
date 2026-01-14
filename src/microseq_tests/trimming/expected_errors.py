# microseq_tests/src/microseq_tests/trimming/expected_errors.py 

from __future__ import annotations 

from typing import Iterable

def expected_errors(phred: Iterable[int]) -> float:
    """
    Compute per-read expected errors (EE; often called MEE / maxee) from Phred scores.

    EE = Σ_i 10^(-Q_i / 10), where Q_i is the Phred quality score for base i.

    Interpretation:
        Under the Phred error model (p_i = 10^(-Q_i/10)), EE is the expected number of
        incorrect bases in the read. This is the quantity filtered by tools like VSEARCH
        via `--fastq_maxee`.

    Parameters
    ----------
    phred : Iterable[int]
        Per-base Phred quality scores as integers (already decoded from FASTQ; not ASCII chars).

    Returns
    -------
    float
        Expected errors for the read. Lower is better (e.g., EE <= 1.0 ≈ “≤1 expected error/read”).

    Notes
    -----
    - Time: O(L) for read length L; Memory: O(1).
    - Empty iterable returns 0.0 (sum over empty set).
    """
    return sum(10 ** (-q / 10) for q in phred) 
