# --------- src/microseq_tests/aligners/schema.py ----------------

from __future__ import annotations 
import pandera as pa # dataframe validation library 
import pandas as pd # write pd.DataFrame in type hints 
from microseq_tests.aligners._parse import COLS 

schema = pa.DataFrameSchema(
    { # might do dictionary comprehesion instead later on....
     "qseqid":   pa.Column(str),
     "sseqid":   pa.Column(str),
     "pident":   pa.Column(float, coerce=True),
     "qlen":     pa.Column(int),
     "qcovhsp":  pa.Column(float, coerce=True),
     "length":   pa.Column(int),
     "evalue":   pa.Column(float),
     "bitscore": pa.Column(float, coerce=True),
     "stitle":   pa.Column(str),
    }
)

def validate(df: pd.DataFrame) -> pd.DataFrame:
    """Raise SchemaError if columns or dtypes deviate; otherwise return df unchanged.""" 
    return schema.validate(df, lazy=True) 
