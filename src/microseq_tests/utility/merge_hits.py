# microseq_tests/src/microseq_tests/utility/merge_hits.py 

from __future__ import annotations 
import glob, logging, sys 
from pathlib import Path 
from microseq_tests.utility.progress import stage_bar, _tls 

log = logging.getLogger(__name__) 

def merge_hits(in_specs: list[str | Path], out_path: str | Path) -> Path: 
    """ 
    Concatenate many MicroSeq Blast TSV files into one. 

    Parameters:
    -------
    in_specs
             Either explicit TSV paths *or* shell-style globs / directories.
    out_path
        Destination file (will be overwritten).
        
    Returns
    -------
    pathlib.Path
        The absolute path to the merged TSV.
    """ 
    # --- resolve globs/dir ------------------------- 
    paths: list[str] = []
    for spec in in_specs:
        p = Path(spec).expanduser()
        if p.is_dir():
            paths.extend(sorted(str(f) # deterministic order 
                for f in p.glob("*.tsv"))) 
        else: 
            g = glob.glob(str(p)) # expand *.tsv globs 
            paths.extend(g if g else [str(p)]) 

    if not paths: 
        raise FileNotFoundError("No TSV files matched the given parameters")

    out = Path(out_path).expanduser().resolve() 
    out.parent.mkdir(parents=True, exist_ok=True) 

    log.info("Merging %d TSV -> %s", len(paths), out) 

    # --- concatenate with the progress bar ----------- 
    with stage_bar(len(paths), desc="merge", unit="file") as bar, \
        out.open("w") as w:
        for i, src in enumerate(paths):
            with open(src) as r: 
                for line in r:
                    if i and line.startswith("#"):
                         continue # keep header only from first file 
                    w.write(line) 
            bar.update() 

            # pardent (outer) bar tick, if any 
            parent = getattr(_tls, "current", None)
            if parent:
                parent.update(1)

    return out 
