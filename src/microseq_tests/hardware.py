from __future__ import annotations
import os
from math import floor

try:
    import psutil  # type: ignore
except Exception:  # pragma: no cover - missing psutil
    psutil = None


def _total_mem_gb() -> float:
    """Return total system memory in gigabytes.

    Falls back to 8 GB when psutil is unavailable.
    """
    if psutil is not None:
        try:
            total_bytes = psutil.virtual_memory().total
        except Exception:
            total_bytes = 8 * 1024 ** 3
    else:
        total_bytes = 8 * 1024 ** 3
    return total_bytes / (1024 ** 3)


def recommend_threads() -> int:
    """Suggest a safe default for ``--threads`` on this machine."""
    logical = os.cpu_count() or 1
    total_mem_gb = _total_mem_gb()
    max_by_mem = floor(total_mem_gb / 0.25)
    return max(1, min(logical, max_by_mem, 16))


def recommend_threads_cli() -> None:
    """CLI entry point: print ``recommend_threads()``."""
    print(recommend_threads())
