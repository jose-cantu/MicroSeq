from __future__ import annotations

BLAST_IDENTITY_MIN = 50.0
BLAST_IDENTITY_MAX = 100.0
BLAST_IDENTITY_DEFAULT = 97.0

BLAST_QCOV_MIN = 10.0
BLAST_QCOV_MAX = 100.0
BLAST_QCOV_DEFAULT = 80.0

BLAST_MAX_HITS_MIN = 1
BLAST_MAX_HITS_MAX = 500
BLAST_MAX_HITS_DEFAULT = 5

BLAST_THREADS_MIN = 1
BLAST_THREADS_MAX = 32
BLAST_THREADS_DEFAULT = 4


def validate_identity(value: float | str) -> float:
    x = float(value)
    if not (BLAST_IDENTITY_MIN <= x <= BLAST_IDENTITY_MAX):
        raise ValueError(
            f"identity must be between {BLAST_IDENTITY_MIN} and {BLAST_IDENTITY_MAX}; got {x}"
        )
    return x


def validate_qcov(value: float | str) -> float:
    x = float(value)
    if not (BLAST_QCOV_MIN <= x <= BLAST_QCOV_MAX):
        raise ValueError(
            f"qcov must be between {BLAST_QCOV_MIN} and {BLAST_QCOV_MAX}; got {x}"
        )
    return x


def validate_max_target_seqs(value: int | str) -> int:
    x = int(value)
    if not (BLAST_MAX_HITS_MIN <= x <= BLAST_MAX_HITS_MAX):
        raise ValueError(
            f"max_target_seqs must be between {BLAST_MAX_HITS_MIN} and {BLAST_MAX_HITS_MAX}; got {x}"
        )
    return x


def validate_threads(value: int | str) -> int:
    x = int(value)
    if not (BLAST_THREADS_MIN <= x <= BLAST_THREADS_MAX):
        raise ValueError(
            f"threads must be between {BLAST_THREADS_MIN} and {BLAST_THREADS_MAX}; got {x}"
        )
    return x


def validate_blast_params(
    *,
    identity: float | str,
    qcov: float | str,
    max_target_seqs: int | str,
    threads: int | str,
) -> tuple[float, float, int, int]:
    return (
        validate_identity(identity),
        validate_qcov(qcov),
        validate_max_target_seqs(max_target_seqs),
        validate_threads(threads),
    )
