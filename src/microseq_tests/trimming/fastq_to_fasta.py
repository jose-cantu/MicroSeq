# src/microseq_tests/trimming/fastq_to_fasta.py 
from __future__ import annotations 
from pathlib import Path 
from typing import Iterable 
from Bio import SeqIO 

def fastq_folder_to_fasta(input_dir: str | Path, 
                          out_fa: str | Path) -> Path:
    """Merges ``*.fastq`` files in ``input_dir`` into one FASTA file ``out_fa``.

    ``input_dir`` is searched recursively so FASTQ files in nested folders are
    also included.
    """
    input_dir = Path(input_dir)
    out_fa = Path(out_fa)
    records: list = []
    for fq in sorted(input_dir.rglob("*.fastq")):
        records.extend(SeqIO.parse(fq, "fastq"))
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, out_fa, "fasta") 
    return out_fa 

