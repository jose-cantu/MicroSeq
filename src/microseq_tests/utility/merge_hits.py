# microseq_tests/src/microseq_tests/utility/merge_hits.py 

from __future__ import annotations 
from pathlib import Path
import logging, glob
from typing import Sequence  
from importlib import import_module 

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
    prog = import_module("microseq_tests.utility.progress") # live view late-bind (honours pytest monkey-patch)  
    cm = prog.stage_bar(len(files), desc="merge", unit="file") # use or not using context-manager here 
    

    # try to close nicely since real bar has .close(); dummy from tests does not 
    if hasattr(cm, "__enter__"): # real tqdm _> use "with" 
        with cm as bar: 
            _write(files, out, bar) 
    else:      # dummy stub -> just use this for now 
        bar = cm 
        _write(files, out, bar) 
        if hasattr(bar, "close"):  # real tqdm has .close(); stub may not 
            bar.close() 

    return out 
    

def _write(files: list[str], out: Path, bar) -> None: 
    """helper that actually writes and ticks bar progress"""
    with out.open("w") as w:
        for ix, fp in enumerate(files):
            with open(fp) as r:
                for line in r: 
                    if ix and line.startswith("#"):
                        continue
                    w.write(line) 
            bar.update(1) 

# ---------- expose for tests that do `merge_hits()` without import  
import builtins as _bi 
_bi.merge_hits = merge_hits # makes it globally visible to my tests =)

