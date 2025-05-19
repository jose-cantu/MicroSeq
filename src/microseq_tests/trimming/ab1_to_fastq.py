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
from Bio import SeqIO 

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

    for ab1 in sorted(input_dir.rglob("*.ab1")):
        out_fq = output_dir / f"{ab1.stem}.fastq"
        if out_fq.exists() and not overwrite:
            out_files.append(out_fq)
            continue 

        # SeqIO.convert handles Phred extract from ABIF tags. 
        SeqIO.convert(ab1, "abi", out_fq, "fastq")
        out_files.append(out_fq) 

    return out_files 


