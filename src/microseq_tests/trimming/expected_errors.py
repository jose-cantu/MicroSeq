# microseq_tests/src/microseq_tests/trimming/expected_errors.py 

from __future__ import annotations 

from typing import Iterable
import math 

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

def qeff_from_mee_per_kb(mee_per_kb: float) -> float:
    """
    Convert length-normalized expected errors (EE/kb) into a Phred equivalent quality score.

    Definitions: 
        - Per base error probability: p_i = 10^(-Q_i/10) 
        - Expected errors for a read: EE = Σ_i p_i 
        - Length normalized: mee_per_kb = 1000 * EE / L 
    
    NOtes: 
        - mee_per_kb must be > 0 (log10 domain). For mee_per_kb <= 0, return a high cap
        (e.g. 60) to represent “effectively perfect” and avoid math domain errors. 

    Examples:
      mee_per_kb = 1   -> qeff = 30
      mee_per_kb = 2   -> qeff ≈ 27
      mee_per_kb = 5   -> qeff ≈ 23
      mee_per_kb = 10  -> qeff = 20
    """ 
    if mee_per_kb <= 0: 
        return 60.0 
    return 30 - 10 * math.log10(mee_per_kb) 
