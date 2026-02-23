"""
ab1_to_fastq.py

convert Sanger chomatograms (.ab1) to single-end FASTQ Phred-33. 

Usage

>>> from microseq_tests.trimming.ab1_to_fastq import ab1_folder_to_fastq 
>>> fastqs = ab1_folder_to_fastq("data/raw", "data/raw_fastq") 

Returns 
list[pathlib.Path]
one path per created FASTQ file. 
"""

from __future__ import annotations 
from pathlib import Path
from typing import List
import hashlib
from collections import Counter
from Bio import SeqIO 


def ab1_rel_key(ab1_path: Path, input_root: Path) -> str:
    """Return a collision-safe key based on AB1 path relative to the staging/input root.

    Important: do not resolve symlinks here, because `--link-raw` may point outside
    the staging tree and `.resolve()` would break `relative_to(...)`.
    """
    rel = ab1_path.relative_to(input_root)
    return "__".join(rel.with_suffix("").parts)



def _short_rel_hash(ab1_path: Path, input_root: Path, size: int = 8) -> str:
    rel = ab1_path.relative_to(input_root)
    return hashlib.sha1(rel.as_posix().encode("utf-8")).hexdigest()[:size]


def build_ab1_output_key_map(input_dir: Path) -> dict[Path, str]:
    """Map each AB1 path to a stable output key.

    Policy: keep legacy stem-based IDs when unique; only disambiguate collisions
    by falling back to staging-relative path keys.
    """
    ab1_paths = sorted(input_dir.rglob("*.ab1"))
    stem_counts = Counter(p.stem for p in ab1_paths)

    keys: dict[Path, str] = {}
    used: dict[str, Path] = {}
    for ab1 in ab1_paths:
        key = ab1.stem if stem_counts[ab1.stem] == 1 else ab1_rel_key(ab1, input_dir)
        if key in used and used[key] != ab1:
            key = f"{key}__{_short_rel_hash(ab1, input_dir)}"
        used[key] = ab1
        keys[ab1] = key
    return keys

def ab1_folder_to_fastq(
        input_dir: str | Path, 
        output_dir: str | Path, 
        *, 
        overwrite: bool = False,
        ) -> List[Path]:
    """
    Convert every ``*.ab1`` file in ``input_dir`` to FASTQ and write them to
    ``output_dir``.  ``input_dir`` is scanned recursively, so any
    sub‑directories containing traces are also processed.

    Paramteres 
    
    input_dir : str | Path
        Folder containing ``.ab1`` chromatograms.  Sub‑directories are
        searched recursively.
    output_dir : str | Path 
        Where FASTQ files will be written (created if missing). 
    overwrite : bool, default False 
        Re-create FASTQ even if it already written 

    Returns 
    list of Path 
        Paths to all FASTQ files written. 

    """
    input_dir = Path(input_dir) 
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True) 

    out_files: list[Path] = []
    key_map = build_ab1_output_key_map(input_dir)

    for ab1 in sorted(input_dir.rglob("*.ab1")):
        out_fq = output_dir / f"{key_map[ab1]}.fastq"
        if out_fq.exists() and not overwrite:
            out_files.append(out_fq)
            continue 

        # SeqIO.convert handles Phred extract from ABIF tags. 
        SeqIO.convert(ab1, "abi", out_fq, "fastq")
        out_files.append(out_fq) 

    return out_files 


