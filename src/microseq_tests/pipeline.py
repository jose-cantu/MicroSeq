"""
microseq_tests.pipeline
Thin wrappers around the five CLI stages.
Return an int exit-code (0 = success) & raise on fatal errors.
"""

from __future__ import annotations
from pathlib import Path
from collections import Counter, defaultdict
import re 
from typing import Union, Sequence
from Bio import SeqIO 
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
from microseq_tests.trimming.quality_trim import  trim_fastq_inputs
from microseq_tests.trimming.ab1_to_fastq import ab1_folder_to_fastq as ab1_to_fastq
from microseq_tests.trimming.biopy_trim import trim_folder as biopy_trim
from microseq_tests.trimming.fastq_to_fasta import (
    fastq_folder_to_fasta as fastq_to_fasta,
)

from microseq_tests.assembly.de_novo_assembly import de_novo_assembly
from microseq_tests.assembly.paired_assembly import assemble_pairs
from microseq_tests.assembly.pairing import  (
        DETECTORS, DupPolicy, _extract_well, extract_sid_orientation, make_pattern_detector) 
from microseq_tests.blast.run_blast import run_blast
from microseq_tests.utility.add_taxonomy import run_taxonomy_join
from microseq_tests.utility.utils import load_config, expand_db_path 
from microseq_tests.post_blast_analysis import run as postblast_run

__all__ = [
    "run_trim",
    "run_assembly",
    "run_paired_assembly",
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
    mee_max: float | None = None,
    mee_min_len: int | None = None, 

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
        input_path = Path(input_path) 
        trim_dir = work / "qc"
        trim_fastq_inputs(input_path, trim_dir, summary_tsv=summary_tsv) 

    fastq_to_fasta(trim_dir, work / "qc" / "trimmed.fasta")
    L.info("Trim finished -> %s", work / "qc" / "trimmed.fasta")
    return 0


# ───────────────────────────────────────────────────────── assembly
def run_assembly(fasta_in: PathLike, out_dir: PathLike, *, threads: int | None = None, **kwargs) -> int:
    options = dict(kwargs)
    if threads is not None:
        options["threads"] = threads 
    de_novo_assembly(Path(fasta_in), Path(out_dir), **options)
    return 0

def run_paired_assembly(input_dir: PathLike, output_dir: PathLike, *, dup_policy: DupPolicy = DupPolicy.ERROR, cap3_options=None, fwd_pattern: str | None = None, rev_pattern: str | None = None, pairing_report: PathLike | None = None, enforce_same_well: bool = False, well_pattern: str | re.Pattern[str] | None = None) -> list[Path]:
    return assemble_pairs(Path(input_dir), Path(output_dir), dup_policy=dup_policy, cap3_options=cap3_options, fwd_pattern=fwd_pattern, rev_pattern=rev_pattern, pairing_report=pairing_report, enforce_same_well=enforce_same_well, well_pattern=well_pattern )

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

def _summarize_paired_candidates(
        directory: Path, fwd_pattern: str | None, rev_pattern: str | None, *, enforce_same_well: bool = False, well_pattern: str | re.Pattern[str] | None = None
) -> str:
    """Return a human-readable summary of how many forward/reverse reads were detected.

    This is used to make the paired-assembly error actionable by reporting what the
    pairing logic was able to recognise inside the generated FASTA directory.
    """

    detectors = list(DETECTORS)
    if fwd_pattern and rev_pattern:
        detectors = [make_pattern_detector(fwd_pattern, rev_pattern), *detectors]

    counts: Counter[str] = Counter()
    sample_orients: dict[str, set[str]] = defaultdict(set)
    unknown: list[str] = [] 
    sid_wells: dict[str, set[str]] = defaultdict(set)
    missing_well: set[str] = set() 

    for suffix in ("*.fasta", "*.fastq", "*.ab1", "*.seq"):
        for fp in sorted(directory.rglob(suffix)):
            sid, orient = extract_sid_orientation(fp.name, detectors=detectors)
            well = _extract_well(fp.name, pattern=well_pattern) if enforce_same_well else None
            key = sid if not enforce_same_well else (well or sid)
            if orient in ("F", "R"):
                counts[orient] += 1
                sample_orients[key].add(orient)
                if enforce_same_well:
                    if well:
                        sid_wells[sid].add(orient)
                    else:
                        missing_well.add(fp.name)
                continue

            counts["unknown"] += 1
            if len(unknown) < 5:
                unknown.append(fp.name)

    missing = sorted(sid for sid, orients in sample_orients.items() if len(orients) == 1)

    parts = [
        f"detected {counts['F']} forward and {counts['R']} reverse reads",
        f"{counts['unknown']} files without F/R tokens",
    ]

    if missing:
        parts.append("samples missing a mate: " + ", ".join(missing[:5]))
        if len(missing) > 5:
            parts[-1] += ", ..."

    if unknown:
        parts.append("example unrecognised filenames: " + ", ".join(unknown))

    if fwd_pattern and rev_pattern:
        parts.append(
            f"custom patterns in use (forward: {fwd_pattern!r}, reverse: {rev_pattern!r})"
        )
    
    if enforce_same_well:
        if missing_well:
            preview = ", ".join(sorted(missing_well)[:3])
            suffix = "..." if len(missing_well) > 3 else "" 
            parts.append(f"{len(missing_well)} files missing plate well codes (e.g., {preview}{suffix}")
        multi_well = sorted(sid for sid, wells in sid_wells.items() if len(wells) > 1)
        if multi_well:
           preview = ", ".join(multi_well[:3])
           suffix = "..." if len(multi_well) > 3 else ""
           parts.append(f"sample IDs spanning multiple wells: {preview}{suffix}")

    return "; ".join(parts)

def _suggest_pairing_patterns(directory: Path) -> str: 
    """ Suggest custom --fwd-pattern/--rev-pattern examples based on filesnames."""
    token_rx = re.compile(r"([A-Za-z0-9]+[FR])", re.I)
    fwd_tokens: Counter[str] = Counter() 
    rev_tokens: Counter[str] = Counter() 
    names: list[str] = [] 

    for fp in sorted(directory.glob("*.fasta")):
        names.append(fp.name)
        for tok in token_rx.findall(fp.name):
            tok = tok.upper()
            if tok.endswith("F"):
                fwd_tokens[tok] += 1 
            elif tok.endswith("R"):
                rev_tokens[tok] += 1 

    def _common(counter: Counter[str]) -> str | None:
        return counter.most_common(1)[0][0] if counter else None 

    fwd = _common(fwd_tokens)
    rev = _common(rev_tokens) 

    suggestions: list[str] = []
    # Concrete regex examples based on the most common tokens MicroSeq detected.
    if fwd or rev:
        suggestions.append(
                "Example flag: --fwd-pattern "
            f"'{fwd or 'FORWARD_TOKEN'}' --rev-pattern '{rev or 'REVERSE_TOKEN'}'" 
        ) 
        suggestions.append(
            "Regex variant (case-insensitive): "
            f"--fwd-pattern '(?i){fwd or '27F'}' --rev-pattern '(?i){rev or '1492R'}'"
        )

    # Offer rename guidance using the detected tokens and the first few filenames as context.

    if names:
        rename_hint = (
            "Rename option: include primer tokens in filenames so MicroSeq can auto-detect "
            f"mates (e.g., sample_{fwd or '27F'}.fasta / sample_{rev or '1492R'}.fasta); "
            "examples seen: " + ", ".join(names[:3])
        )
        suggestions.append(rename_hint)

    if not suggestions:
        # Fallback message when no obvious tokens are present.
        suggestions.append(
            "No primer-like tokens detected; rename files to include forward/reverse labels "
            "(e.g., sample_27F / sample_1492R) or pass --fwd-pattern/--rev-pattern explicitly."
        )

    return " ".join(suggestions) 

def run_full_pipeline(
    infile: Path,
    db_key: str,
    out_dir: Path | None = None,
    *,
    mode: str = "single", 
    postblast: bool = False,
    identity: int = 97,
    qcov: int = 80,
    max_target_seqs: int = 5,
    threads: int = 4,
    blast_task: str = "megablast", 
    dup_policy: DupPolicy = DupPolicy.ERROR, 
    cap3_options=None, 
    fwd_pattern: str | None = None,
    rev_pattern: str | None = None,
    enforce_same_well: bool = False,
    well_pattern: str | re.Pattern[str] | None = None, 
    metadata: Path | None = None,
    summary_tsv: Path | None = None,
    mee_max: float | None = None,
    mee_min_len: int | None = None, 
    on_stage=None,
    on_progress=None,
) -> dict[str, Path]:
    """Run trim -> FASTA merge -> BLAST -> taxonomy (+ optional post‑BLAST).

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
    using_paired = mode == "paired" 

    if summary_tsv is None and not is_fasta:
        summary_tsv = out_dir / "qc" / "trim_summary.tsv"

    paths = {
        "trimmed_fastq": None if (is_fasta or using_paired) else out_dir / "qc" / "trimmed.fastq",
        "trimmed_fasta": infile if is_fasta else out_dir / "qc" / "trimmed.fasta",
        "fasta": infile if is_fasta else out_dir / "reads.fasta",
        "hits": out_dir / "hits.tsv",
        "tax": out_dir / "hits_tax.tsv",
        "biom": out_dir / "table.biom",
        "trim_summary": summary_tsv if (not is_fasta and not using_paired) else None,
    }

    if using_paired:
        # Assembly + BLAST + taxonomy + postblast (optional here if using) 
        n_stages = (3 if is_fasta else 5) + int(postblast)   
    else:
        n_stages = (2 if is_fasta else 4) + int(postblast)
    step = 100 // n_stages
    pct = 0

    def subprog(off):
        return lambda p: on_progress(off + p * step // 100)

    if using_paired:
        def _merge_cap3_contigs(files: Sequence[Path], destination: Path, map_tsv: Path | None = None) -> None:
            destination.parent.mkdir(parents=True, exist_ok=True)
            
            map_tsv = map_tsv or destination.with_suffix(".tsv")
            with destination.open("w", encoding="utf-8") as out, map_tsv.open(
                    "w", encoding="utf-8" 
            ) as manifest:
                manifest.write("qseqid\tsample\n")
                
                for fp in files:
                    sid = Path(fp).parent.name
                    for idx, rec in enumerate(SeqIO.parse(fp, "fasta"), 1):
                        rec.id = f"{sid}_c{idx}"
                        rec.name = rec.id
                        rec.description = ""
                        SeqIO.write(rec, out, "fasta")
                        manifest.write(f"{rec.id}\t{sid}\n")
                   
        def _fastq_to_paired_fastas(source_dir: Path, dest_dir: Path) -> list[Path]:
            dest_dir.mkdir(parents=True, exist_ok=True)
            written: list[Path] = []
            for fq in sorted(source_dir.rglob("*.fastq")):
                records = list(SeqIO.parse(fq, "fastq"))
                if not records:
                    continue
                out_fp = dest_dir / f"{fq.stem}.fasta"
                SeqIO.write(records, out_fp, "fasta")
                written.append(out_fp)
            return written

        assembly_input = infile
        if not is_fasta:
            on_stage("Trim")
            run_trim(infile, out_dir, sanger=True, summary_tsv=summary_tsv,)
            if thr and thr.isInterruptionRequested():
                raise RuntimeError("Cancelled")
            pct += step
            on_progress(pct)

            on_stage("Convert")
            fastq_dir = out_dir / "passed_qc_fastq"
            if not any(fastq_dir.glob("*.fastq")):
                fastq_dir = out_dir / "qc"
            run_fastq_to_fasta(fastq_dir, paths["trimmed_fasta"])
            if thr and thr.isInterruptionRequested():
                raise RuntimeError("Cancelled")
            pct += step
            on_progress(pct)

            fasta_dir = out_dir / "qc" / "paired_fasta"
            generated_fastas = _fastq_to_paired_fastas(fastq_dir, fasta_dir)
            if not generated_fastas:
                raise ValueError(f"No FASTA files generated from {fastq_dir}")
            assembly_input = fasta_dir
        

        pairing_report = out_dir / "qc" / "pairing_report.tsv"
        pairing_report.parent.mkdir(parents=True, exist_ok=True)

        on_stage("Paired assembly")
        contig_paths = assemble_pairs(
            assembly_input, 
            out_dir / "asm",
            dup_policy=dup_policy,
            cap3_options=cap3_options,
            fwd_pattern=fwd_pattern,
            rev_pattern=rev_pattern,
            pairing_report=pairing_report,
            enforce_same_well=enforce_same_well,
            well_pattern=well_pattern 
        )
        if not contig_paths:
           summary = _summarize_paired_candidates(assembly_input, fwd_pattern, rev_pattern, enforce_same_well=enforce_same_well, well_pattern=well_pattern)
           suggestions = _suggest_pairing_patterns(assembly_input)
           raise ValueError(
                f"No paired reads detected in {assembly_input}; {summary}. "
                "If your primer names differ, provide --fwd-pattern/--rev-pattern."
                f" {suggestions}"
            )  
        merged_contigs = out_dir / "asm" / "paired_contigs.fasta"
        contig_map = out_dir / "asm" / "contig_map.tsv"
        _merge_cap3_contigs(contig_paths, merged_contigs, contig_map)
        paths["fasta"] = merged_contigs
        paths["trimmed_fasta"] = merged_contigs

        pct += step
        on_progress(pct)

    elif not is_fasta:
        # 1 – Trim
        on_stage("Trim")

        sanger = infile.is_dir() or infile.suffix.lower() == ".ab1"
        run_trim(infile, out_dir, sanger=sanger, summary_tsv=summary_tsv,)

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

    if using_paired:
        contig_map = out_dir / "asm" / "contig_map.tsv"
        try:
            if contig_map.exists():
                tax_df = pd.read_csv(paths["tax"], sep="\t")
                map_df = pd.read_csv(contig_map, sep="\t")
                
                merge_key = None
                if "qseqid" in tax_df.columns:
                    merge_key = "qseqid"
                elif "sample_id" in tax_df.columns:
                    merge_key = "sample_id"
                    map_df = map_df.rename(columns={"qseqid": "sample_id"})

                if merge_key and "sample" not in tax_df.columns:
                    tax_df = tax_df.merge(map_df, on=merge_key, how="left") 

                    tax_df.to_csv(paths["tax"], sep="\t", index=False)
        except Exception as exc:
            L.warning("Failed to merge contig map into taxonomy table: %s", exc)

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


# ───────────────────────────────────────────────────────── AB1 -> FASTQ
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
        L.info("AB1->FASTQ wrote %d files to %s", len(written), output_dir)
        return 0
    except Exception:
        L.exception("AB1->FASTQ failed")
        return 1


# ───────────────────────────────────────────────────────── FASTQ -> FASTA
def run_fastq_to_fasta(
    input_dir: PathLike,
    output_fasta: PathLike,
) -> int:
    """
    Merge all *.fastq in input_dir into a single FASTA output_fasta.
    """
    try:
        out = fastq_to_fasta(input_dir, output_fasta)
        L.info("FASTQ->FASTA wrote %s", out)
        return 0
    except Exception:
        L.exception("FASTQ->FASTA failed")
        return 1
