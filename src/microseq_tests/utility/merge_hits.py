# microseq_tests/src/microseq_tests/utility/merge_hits.py 

from __future__ import annotations 
from pathlib import Path
import logging, glob
from typing import Sequence  
from .progress import stage_bar 


L = logging.getLogger(__name__)
__all__ = ["merge_hits"] 

def merge_hits(input_specs: Sequence[str], out_tsv: str | Path) -> Path: 
    """ 
    Concatenate one or more BLAST TSVs into output single file. 
    input_specs can be literal paths or shell-type patterns. 
    Only the header lines (starting with '#') fro the first file are kept. 

    Returns the absolute :class: `pathlib.Path` of merged TSV. 
    """ 
    files: list[str] = []
    for spec in input_specs:
        matches = glob.glob(spec)
        files.extend(matches if matches else [spec])

    if not files:
        raise FileNotFoundError("merged_hits: no TSVs matched")

    out = Path(out_tsv).expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)

    L.info("Merging %d TSV -> %s", len(files), out)

    # ------------- concatenate while updating progress bar --------- 
    with stage_bar(len(files), desc="merge", unit="file") as bar, \
        out.open("w") as w:
        for idx, fp in enumerate(files):
            with open(fp) as r:     
                for line in r: 
                   if idx and line.startswith("#"):
                        continue
                   w.write(line)
            bar.update(1) 

    return out
