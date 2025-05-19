"""
microseq_tests.pipeline
Thin wrappers around the five CLI stages.
Return an int exit-code (0 = success) & raise on fatal errors.
"""

from __future__ import annotations
from pathlib import Path
from typing import Union, Sequence
import logging
import shutil

# pull the existing implementation functions
from microseq_tests.trimming.quality_trim import quality_trim
from microseq_tests.trimming.ab1_to_fastq import ab1_folder_to_fastq as ab1_to_fastq
from microseq_tests.trimming.biopy_trim import trim_folder as biopy_trim
from microseq_tests.trimming.fastq_to_fasta import (
    fastq_folder_to_fasta as fastq_to_fasta,
)
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
    "run_full_pipeline",
]

PathLike = Union[str, Path]
L = logging.getLogger(__name__)


# ───────────────────────────────────────────────────────── trimming
def run_trim(input_path: PathLike, workdir: PathLike, sanger: bool = False) -> int:
    """
    QC-trim FASTQ or convert+trim an AB1 folder.

    Returns 0 on success.
    """
    work = Path(workdir)
    work.mkdir(parents=True, exist_ok=True)
    (work / "qc").mkdir(parents=True, exist_ok=True)

    if sanger:
        fastq_dir = work / "raw_fastq"
        ab1_source = Path(input_path)
        if ab1_source.is_file():
            raw_ab1 = work / "raw_ab1"
            raw_ab1.mkdir(parents=True, exist_ok=True)
            shutil.copy(ab1_source, raw_ab1 / ab1_source.name)
            ab1_source = raw_ab1
        ab1_to_fastq(ab1_source, fastq_dir)
        biopy_trim(fastq_dir, work / "qc")
        trim_dir = work / "passed_qc_fastq"
    else:
        out_fq = work / "qc" / "trimmed.fastq"
        quality_trim(input_path, out_fq)
        trim_dir = work / "qc"

    fastq_to_fasta(trim_dir, work / "qc" / "trimmed.fasta")
    L.info("Trim finished → %s", work / "qc" / "trimmed.fasta")
    return 0


# ───────────────────────────────────────────────────────── assembly
def run_assembly(fasta_in: PathLike, out_dir: PathLike) -> int:
    de_novo_assembly(Path(fasta_in), Path(out_dir))
    return 0


# ───────────────────────────────────────────────────────── BLAST
def run_blast_stage(
    fasta_in: PathLike,
    db_key: str,
    out_tsv: PathLike,
    identity: float = 97.0,
    qcov: float = 80.0,
    max_target_seqs: int = 5,
    threads: int = 1,
    on_progress=None,
) -> int:
    run_blast(
        fasta_in,
        db_key,
        out_tsv,
        search_id=identity,
        search_qcov=qcov,
        max_target_seqs=max_target_seqs,
        threads=threads,
        on_progress=on_progress,
    )
    return 0


# ───────────────────────────────────────────────────────── add-taxonomy
def run_add_tax(hits: PathLike, taxonomy_tsv: PathLike, out_tsv: PathLike) -> int:
    run_taxonomy_join(hits, taxonomy_tsv, out_tsv)
    return 0


# ───────────────────────────────────────────────────────── post-BLAST
def run_postblast(
    blast_hits: PathLike,
    metadata: PathLike,
    out_biom: PathLike,
    sample_col: str | None = None,
    identity_th: float = 97.0,
    *,  # forces everything after passed by name, prevent position mistakes
    id_normaliser: str = "none",
    taxonomy_col: str = "auto",
    taxonomy_format: str = "auto",
    duplicate_policy: str = "error",
    **kwargs,
) -> int:
    postblast_run(
        blast_hits,
        metadata,
        out_biom,
        write_csv=True,
        sample_col=sample_col,
        identity_th=identity_th,
        id_normaliser=id_normaliser,
        taxonomy_col=taxonomy_col,
        taxonomy_format=taxonomy_format,
        duplicate_policy=duplicate_policy,
        **kwargs,
    )
    return 0


# # ────────────────────────────────────────────────── full workflow here =)
from pathlib import Path
from microseq_tests.utility.utils import load_config, expand_db_path


def run_full_pipeline(
    infile: Path,
    db_key: str,
    out_dir: Path | None = None,
    *,
    postblast: bool = False,
    identity: int = 97,
    qcov: int = 80,
    max_target_seqs: int = 5,
    threads: int = 1,
    metadata: Path | None = None,
    on_stage=None,
    on_progress=None,
) -> dict[str, Path]:
    """Run trim → FASTA merge → BLAST → taxonomy (+ optional post‑BLAST).

    *infile* may be FASTA, FASTQ, a single ``.ab1`` trace, or a directory of
    ``.ab1`` files.  Sanger mode is triggered automatically when *infile* is a
    directory or ends with ``.ab1``.
    """

    on_stage = on_stage or (lambda *_: None)
    on_progress = on_progress or (lambda *_: None)

    if out_dir is None:
        stem = infile.with_suffix("").name
        out_dir = infile.parent / f"{stem}_microseq"
    out_dir.mkdir(parents=True, exist_ok=True)
    on_progress(0)

    is_fasta = infile.suffix.lower() in {".fasta", ".fa", ".fna", ".fas"}

    paths = {
        "trimmed_fastq": None if is_fasta else out_dir / "qc" / "trimmed.fastq",
        "trimmed_fasta": infile if is_fasta else out_dir / "qc" / "trimmed.fasta",
        "fasta": infile if is_fasta else out_dir / "reads.fasta",
        "hits": out_dir / "hits.tsv",
        "tax": out_dir / "hits_tax.tsv",
        "biom": out_dir / "table.biom",
    }

    n_stages = (2 if is_fasta else 4) + int(postblast)
    step = 100 // n_stages
    pct = 0

    def subprog(off):
        return lambda p: on_progress(off + p * step // 100)

    if not is_fasta:
        # 1 – Trim
        on_stage("Trim")
        sanger = infile.is_dir() or infile.suffix.lower() == ".ab1"
        run_trim(infile, out_dir, sanger=sanger)
        pct += step
        on_progress(pct)

        # 2 – Merge FASTQs to FASTA
        on_stage("Convert")
        fastq_dir = out_dir / "qc"
        if not any(fastq_dir.glob("*.fastq")):
            alt = out_dir / "passed_qc_fastq"
            if alt.exists() and any(alt.glob("*.fastq")):
                fastq_dir = alt
        run_fastq_to_fasta(fastq_dir, paths["fasta"])
        pct += step
        on_progress(pct)

    # 3 – BLAST
    on_stage("BLAST")
    run_blast_stage(
        paths["fasta"],
        db_key,
        paths["hits"],
        identity=identity,
        qcov=qcov,
        max_target_seqs=max_target_seqs,
        threads=threads,
        on_progress=subprog(pct),
    )
    pct += step
    on_progress(pct)

    # 4 – Add taxonomy
    cfg = load_config()
    tax_template = cfg["databases"][db_key]["taxonomy"]
    tax_fp = Path(expand_db_path(tax_template))
    on_stage("Taxonomy")
    run_add_tax(paths["hits"], tax_fp, paths["tax"])
    pct += step
    on_progress(pct)

    # 5 – Optional post-BLAST
    if postblast:
        if metadata is None:
            raise ValueError("metadata file required for postblast stage")
        on_stage("Post-BLAST")
        run_postblast(paths["tax"], metadata, paths["biom"])
        pct += step
        on_progress(pct)

    on_progress(100)
    return paths


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
        written: Sequence[Path] = ab1_to_fastq(
            input_dir, output_dir, overwrite=overwrite
        )
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
