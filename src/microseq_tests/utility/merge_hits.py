# microseq_tests/src/microseq_tests/utility/merge_hits.py 

from __future__ import annotations 
from pathlib import Path
import logging, glob
from typing import Sequence  
from importlib import import_module 

L = logging.getLogger(__name__)
__all__ = ["merge_hits"]

# ---- helpers -----------------------
def _tick_safe(bar, inc: int = 1) -> None: 
    """Call bar.update(inc) only if the attribute exists (dummy bars won't here). """ 
    if hasattr(bar, "update"):
        bar.update(inc) 

def _write_tsvs(files: list[str], out: Path, bar) -> None: 
    """Concatenate TSVs, keeping the header of the first file only.""" 
    with out.open("w") as w: 
        for idx, fp in enumerate(files): 
            with open(fp) as r: 
                for line in r: 
                    if idx and line.startswith("#"):  # skip dup headers 
                        continue 
                    w.write(line) 
            _tick_safe(bar)  # 1 tick / file 

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

    # stage_bar may be tqdm, a bare DummyBar, or a geenrator that yeilds one 
    prog = import_module("microseq_tests.utility.progress") # live view late-bind (honours pytest monkey-patch)  
    cm_or_gen = prog.stage_bar(len(files), desc="merge", unit="file") # use or not using context-manager here 
    
    # proper context-manager using real tqdm progress bar 
    if hasattr(cm_or_gen, "__enter__"):  
        with cm_or_gen as bar: 
            _write_tsvs(files, out, bar) 
    else:      
        # generator that yeilds the bar, or bare bar 
        try: 
            bar = next(cm_or_gen) # succeeds for generator variety 
        except StopIteration:
            bar = cm_or_gen  # it was already a bare bar 


        _write_tsvs(files, out, bar) 


        # close if objects expose .close() 
        if hasattr(cm_or_gen, "close"): 
            cm_or_gen.close() 
        if hasattr(bar, "close") and bar is not cm_or_gen:
            bar.close() 

    return out 

