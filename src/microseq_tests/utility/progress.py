# microseq_tests/src/microseq_tests/utility/progress.py 

from __future__ import annotations 
from contextlib import contextmanager 
from tqdm import tqdm 
from typing import Iterator
import threading, logging 

_tls = threading.local()  # module-level, one per thread 
L = logging.getLogger(__name__) 

@contextmanager 
def stage_bar(total: int, *, desc: str = "", unit: str = "") -> Iterator[tqdm]:
    """
    Context manager that yields a tqdm and remembers the last one in a thread-local so nested calls 
    can find their parent automatically. This is used in a nested helper function such as run_blast that calls parent.update(1)
    without the caller having to pass the bar explicitly. 
    """ 
    outer = getattr(_tls, "current", None) 
    bar = tqdm(total=total, desc=desc, unit=unit, leave=False, ncols=80, bar_format=("{l_bar}{bar}| " "{n_fmt}/{total_fmt} " "[elapsed: {elapsed} < remaining: {remaining}]"),) # leave = False here so finished bar leaves and only see latest one 
    _tls.current = bar # make this bar the new "current one" 

    try:
        yield bar # hand the bar to the "with" block so caller does bar.update() 
    finally:
        bar.close() # tidy up 
        _tls.current = outer # restore previous parent (or None here) 


# Legacy shim - old tests import merge_hits through this module. 
# keeping legacy import at very bottom to avoid circular dependency. 
from microseq_tests.utility.merge_hits import merge_hits # noqa: E402, F401  
