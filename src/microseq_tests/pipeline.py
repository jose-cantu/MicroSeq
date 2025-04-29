"""
microseq_tests.pipeline
Thin wrappers around the five CLI stages.
Return an int exit-code (0 = success) & raise on fatal errors.
"""

from __future__ import annotations
from pathlib import Path
from typing import Union, Sequence
import logging

# pull the existing implementation functions
from microseq_tests.trimming.quality_trim import quality_trim
from microseq_tests.trimming.ab1_to_fastq import ab1_folder_to_fastq as ab1_to_fastq 
from microseq_tests.trimming.biopy_trim import trim_folder as biopy_trim
from microseq_tests.trimming.fastq_to_fasta import fastq_folder_to_fasta as fastq_to_fasta
from microseq_tests.assembly.de_novo_assembly import de_novo_assembly
from microseq_tests.blast.run_blast import run_blast
from microseq_tests.utility.add_taxonomy import run_taxonomy_join
from microseq_tests.post_blast_analysis import run as postblast_run

__all__ = [
        "run_trim",
        "run_assembly",
        "run_blast_stage",
        "run_add_tax",
        "run_postblast",
        "run_ab1_to_fastq",
        "run_fastq_to_fasta",
        ] 

PathLike = Union[str, Path]
L = logging.getLogger(__name__)

# ───────────────────────────────────────────────────────── trimming
def run_trim(input_path: PathLike,
             workdir: PathLike,
             sanger: bool = False) -> int:
    """
    QC-trim FASTQ or convert+trim an AB1 folder.

    Returns 0 on success.
    """
    work = Path(workdir)
    work.mkdir(parents=True, exist_ok=True)

    if sanger:
        fastq_dir = work / "raw_fastq"
        ab1_folder_to_fastq(Path(input_path), fastq_dir)
        biopy_trim(fastq_dir, work / "qc")
    else:
        out_fq = work / "qc" / "trimmed.fastq"
        quality_trim(input_path, out_fq)

    fastq_folder_to_fasta(work / "qc", work / "qc" / "trimmed.fasta")
    L.info("Trim finished → %s", work / "qc" / "trimmed.fasta")
    return 0


# ───────────────────────────────────────────────────────── assembly
def run_assembly(fasta_in: PathLike,
                 out_dir: PathLike) -> int:
    de_novo_assembly(Path(fasta_in), Path(out_dir))
    return 0


# ───────────────────────────────────────────────────────── BLAST
def run_blast_stage(fasta_in: PathLike,
                    db_key: str,
                    out_tsv: PathLike,
                    identity: float = 97.0,
                    qcov: float = 80.0,
                    max_target_seqs: int = 5,
                    threads: int = 1) -> int:
    run_blast(fasta_in, db_key, out_tsv,
              pct_id=identity, qcov=qcov,
              max_target_seqs=max_target_seqs,
              threads=threads)
    return 0


# ───────────────────────────────────────────────────────── add-taxonomy
def run_add_tax(hits: PathLike,
                taxonomy_tsv: PathLike,
                out_tsv: PathLike) -> int:
    run_taxonomy_join(hits, taxonomy_tsv, out_tsv)
    return 0


# ───────────────────────────────────────────────────────── post-BLAST
def run_postblast(blast_hits: PathLike,
                  metadata: PathLike,
                  out_biom: PathLike,
                  sample_col: str | None = None) -> int:
    postblast_run(blast_hits, metadata, out_biom,
                  write_csv=True, sample_col=sample_col)
    return 0

# ───────────────────────────────────────────────────────── AB1 → FASTQ
def run_ab1_to_fastq(
    input_dir: PathLike,
    output_dir: PathLike,
    *,
    overwrite: bool = False,
) -> int:
    """
    Convert every *.ab1 in *input_dir* to FASTQ files in *output_dir*.
    Returns 0 on success, 1 on failure (and logs the error).
    """
    try:
        written: Sequence[Path] = ab1_to_fastq(input_dir, output_dir, overwrite=overwrite)
        L.info("AB1→FASTQ wrote %d files to %s", len(written), output_dir)
        return 0
    except Exception:
        L.exception("AB1→FASTQ failed")
        return 1


# ───────────────────────────────────────────────────────── FASTQ → FASTA
def run_fastq_to_fasta(
    input_dir: PathLike,
    output_fasta: PathLike,
) -> int:
    """
    Merge all *.fastq in *input_dir* into a single FASTA *output_fasta*.
    """
    try:
        out = fastq_to_fasta(input_dir, output_fasta)
        L.info("FASTQ→FASTA wrote %s", out)
        return 0
    except Exception:
        L.exception("FASTQ→FASTA failed")
        return 1

