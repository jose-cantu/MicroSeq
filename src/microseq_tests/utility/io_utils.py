from __future__ import annotations 
import logging
L = logging.getLogger(__name__)

from pathlib import Path 
import re, shutil, logging 

__all__ = ["normalise_tsv"] 

_TAB_RX = re.compile(r"( {2,}|,)") # 2 + spaces or comma 
_NEEDS_RX = re.compile(r"\t") # have at least one TAB char 

def normalise_tsv(path: str | Path) -> Path:
    """
    Ensure path is a true tab-separated file. 
    Detect if lines contain TABS already then do nothing please. 
    If not, then collapse 2 space or single commas into one TAB. 
    Write to a tmp file, then atomically replace the original with the TSV/TAB version
    """
    p = Path(path)
    text = p.read_text()

    if "\t" in text: # already fine here 
        return p

    # ------- header repair -----------
    lines = text.splitlines() 
    if lines and "\t" not in lines[0]:
        # if the header still has no TAB, but has >1 word, 
        # split on any run of space/comma 
        header_parts = re.split(r"[ ,]+", lines[0].strip())  
        if len(header_parts) >= 2:
            lines[0] = "\t".join(header_parts) 
            text = "\n".join(lines) 

    L.info("[fix-tsv] converting spaces / commas to TABS in %s", p.name)
    fixed = re.sub(r"( {2,}|,)", "\t", text) # 2+ spaces or comma gets turned into TAB format.

    tmp = p.with_suffix(".tsv.tmp")
    tmp.write_text(fixed)
    shutil.move(tmp, p) # atomic replacement here
    return p 

def cli():
    import sys 
    for f in sys.argv[1:]:
        normalise_tsv(f)
    print("files are normalised now and ready to run in MicroSeq", *sys.argv[1:])
