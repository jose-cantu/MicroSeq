# ------- tests/test_schema.py ------------- 

from microseq_tests.aligners._parse import COLS, parse_blast_tab 
import pandas as pd, pathlib, tempfile 
tmp = pathlib.Path(tempfile.gettempdir()) / "dummy.tsv"
tmp.write_text("a\tb\t99\t100\t99\t100\t1e-5\t200\ttitle\n")

from microseq_tests.aligners.schema import validate 
df = parse_blast_tab(tmp) 
validate(df) 
