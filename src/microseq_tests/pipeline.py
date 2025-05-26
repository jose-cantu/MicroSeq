"""
microseq_tests.pipeline
Thin wrappers around the five CLI stages.
Return an int exit-code (0 = success) & raise on fatal errors.
"""

from __future__ import annotations
from pathlib import Path
from typing import Union, Sequence
import pandas as pd 
import logging
try:
    from PySide6.QtCore import QThread
except Exception:  # allow use without PySide
    class QThread:
        @staticmethod
        def currentThread():
            return None
import shutil

# pull the existing implementation functions
from microseq_tests.trimming.quality_trim import quality_trim
from microseq_tests.trimming.ab1_to_fastq import ab1_folder_to_fastq as ab1_to_fastq
from microseq_tests.trimming.biopy_trim import trim_folder as biopy_trim
from microseq_tests.trimming.fastq_to_fasta import (
    fastq_folder_to_fasta as fastq_to_fasta,
)
from Bio import SeqIO
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
def run_trim(
    input_path: PathLike,
    workdir: PathLike,
    sanger: bool = False,
    *,
    summary_tsv: PathLike | None = None,
    link_raw: bool = False,

) -> int:
    """Trim reads and convert if needed.

    When ``sanger=True`` the input_path must point to an ``.ab1`` file or
    directory.  Those traces are converted to FASTQ before trimming.
    With ``sanger=False`` the function expects standard FASTQ input.
    When *summary_tsv* is given, write one-line stats per file to that path.
    Set link_raw to ``True`` to symlink AB1 traces into ``raw_ab1`` instead of
    copying them.


    Returns 0 on success.
    """
    work = Path(workdir)
    work.mkdir(parents=True, exist_ok=True)
    (work / "qc").mkdir(parents=True, exist_ok=True)

    if sanger:
        fastq_dir = work / "raw_fastq"
        ab1_source = Path(input_path)
        dst = work / "raw_ab1"
        if dst.exists():
            if dst.is_symlink() or dst.is_file():
                dst.unlink()
            else:
                shutil.rmtree(dst)
        if link_raw:
            dst.symlink_to(ab1_source.resolve(), target_is_directory=ab1_source.is_dir())
        else:
            if ab1_source.is_dir():
                shutil.copytree(ab1_source, dst, dirs_exist_ok=True)
            else:
                dst.mkdir(parents=True, exist_ok=True)
                shutil.copy2(ab1_source, dst / ab1_source.name)
        ab1_to_fastq(dst, fastq_dir)

        biopy_trim(fastq_dir, work / "qc", combined_tsv=summary_tsv)
        trim_dir = work / "passed_qc_fastq"
    else:
        out_fq = work / "qc" / "trimmed.fastq"
        quality_trim(input_path, out_fq)
        trim_dir = work / "qc"
        if summary_tsv:
            reads = bases = qsum = 0
            for rec in SeqIO.parse(out_fq, "fastq"):
                ph = rec.letter_annotations["phred_quality"]
                reads += 1
                bases += len(rec)
                qsum += sum(ph)
            avg_q = qsum / bases if bases else 0
            avg_len = bases / reads if reads else 0
            summary_fp = Path(summary_tsv)
            summary_fp.parent.mkdir(parents=True, exist_ok=True)
            write_header = not summary_fp.exists()
            with open(summary_fp, "a") as comb:
                if write_header:
                    comb.write("file\treads\tavg_len\tavg_q\n")
                comb.write(
                    f"{Path(input_path).name}\t{reads}\t{avg_len:.1f}\t{avg_q:.2f}\n"
                )

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
    blast_task: str = "megablast", 
) -> int:
    from microseq_tests.blast.run_blast import BlastOptions # local import keeps AB1 safe 
     

    thr = QThread.currentThread()
    if thr and thr.isInterruptionRequested():
        raise RuntimeError("Cancelled")


    # ------ first pass ---------------------------
    opts = BlastOptions(task=blast_task) 

    run_blast(
        fasta_in,
        db_key,
        out_tsv,
        options=opts, 
        search_id=identity,
        search_qcov=qcov,
        max_target_seqs=max_target_seqs,
        threads=threads,
        on_progress=on_progress,
    )
    # ------- Sensitive-mode fallback ----------------------
    # if the user chose the fast algorithm and the best hit is below 90 ID / 90 qcov   # then redo the search but with the more sensitive comprehesive blastn 
    if blast_task == "megablast":
        try:
            hits = pd.read_csv(out_tsv, sep="\t", comment="#", usecols=["pident", "qcovhsp"], nrows=1, dtype=float) 
            # empty -> no rows or best identity / qcov below 90 
            if hits.empty:                       # no hits at all
                best_id = best_qcov = 0.0

            else:                                # grab the first-row values
                first = hits.iloc[0]
                best_id   = float(first.pident)
                best_qcov = float(first.qcovhsp)

            needs_retry = (best_id < 90.0) or (best_qcov < 90.0)



        except Exception as e: # malformed table -> playing it safe 
            L.warning("Fallback check failed (%s) - rerunning with blastn", e)
            needs_retry = True 

        if needs_retry:
            L.info("Fast run found no close hit - so rerunning in sensitive mode using blastn" "(best %.1f %% ID / %.1f %% qcov)",
                   best_id, best_qcov)

            opts.task="blastn" # flip the dataclass 
            run_blast(
                fasta_in,
                db_key,
                out_tsv,
                options=opts,   # ← switch task in question 
                search_id=identity,
                search_qcov=qcov,
                max_target_seqs=max_target_seqs,
                threads=threads,
                on_progress=on_progress,
                ) 


    if thr and thr.isInterruptionRequested():
        raise RuntimeError("Cancelled")
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
    threads: int = 4,
    blast_task: str = "megablast", 
    metadata: Path | None = None,
    summary_tsv: Path | None = None,
    on_stage=None,
    on_progress=None,
) -> dict[str, Path]:
    """Run trim → FASTA merge → BLAST → taxonomy (+ optional post‑BLAST).

    infile may be FASTA, FASTQ, a single ``.ab1`` trace, or a directory of
    ``.ab1`` files.  Sanger mode is triggered automatically when *infile* is a
    directory or ends with ``.ab1``.  When *summary_tsv* is ``None`` the summary
    is written to ``qc/trim_summary.tsv`` inside *out_dir*.
    """

    on_stage = on_stage or (lambda *_: None)
    on_progress = on_progress or (lambda *_: None)
    thr = QThread.currentThread()
    if thr and thr.isInterruptionRequested():
        raise RuntimeError("Cancelled")

    if out_dir is None:
        stem = infile.with_suffix("").name
        out_dir = infile.parent / f"{stem}_microseq"
    out_dir.mkdir(parents=True, exist_ok=True)
    on_progress(0)

    is_fasta = infile.suffix.lower() in {".fasta", ".fa", ".fna", ".fas"}

    if summary_tsv is None and not is_fasta:
        summary_tsv = out_dir / "qc" / "trim_summary.tsv"

    paths = {
        "trimmed_fastq": None if is_fasta else out_dir / "qc" / "trimmed.fastq",
        "trimmed_fasta": infile if is_fasta else out_dir / "qc" / "trimmed.fasta",
        "fasta": infile if is_fasta else out_dir / "reads.fasta",
        "hits": out_dir / "hits.tsv",
        "tax": out_dir / "hits_tax.tsv",
        "biom": out_dir / "table.biom",
        "trim_summary": summary_tsv if not is_fasta else None,
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
        run_trim(infile, out_dir, sanger=sanger, summary_tsv=summary_tsv)

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
        if thr and thr.isInterruptionRequested():
            raise RuntimeError("Cancelled")
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
        blast_task=blast_task, 
    )
    if thr and thr.isInterruptionRequested():
        raise RuntimeError("Cancelled")
    pct += step
    on_progress(pct)

    # 4 – Add taxonomy
    cfg = load_config()
    tax_template = cfg["databases"][db_key]["taxonomy"]
    tax_fp = Path(expand_db_path(tax_template))
    on_stage("Taxonomy")
    run_add_tax(paths["hits"], tax_fp, paths["tax"])
    if thr and thr.isInterruptionRequested():
        raise RuntimeError("Cancelled")
    pct += step
    on_progress(pct)

    # 5 – Optional post-BLAST
    if postblast:
        if metadata is None:
            raise ValueError("metadata file required for postblast stage")
        on_stage("Post-BLAST")
        run_postblast(paths["tax"], metadata, paths["biom"])
        if thr and thr.isInterruptionRequested():
            raise RuntimeError("Cancelled")
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
    Convert every *.ab1 in input_dir to FASTQ files in output_dir.
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
    Merge all *.fastq in input_dir into a single FASTA output_fasta.
    """
    try:
        out = fastq_to_fasta(input_dir, output_fasta)
        L.info("FASTQ→FASTA wrote %s", out)
        return 0
    except Exception:
        L.exception("FASTQ→FASTA failed")
        return 1
