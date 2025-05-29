# -----src/microseq_tests/aligners/base.py -------------------

from __future__ import annotations  
from typing import Protocol 
import pandas as pd 

class BaseAligner(Protocol):
    name: str

    def run(self,
            query_fa: str,
            db_key: str,
            threads: int = 1,
            **kw 
    ) -> pd.DataFrame: ...


