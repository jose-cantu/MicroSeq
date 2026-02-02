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
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
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
from microseq_tests.assembly.paired_assembly import assemble_pairs, _build_keep_separate_pairs 
from microseq_tests.assembly.pairing import  (
        DETECTORS, DupPolicy, _extract_well, _strip_well_token, extract_sid_orientation, make_pattern_detector) 
from microseq_tests.assembly.cap3_report import parse_cap3_reports
from microseq_tests.assembly.cap3_profiles import resolve_cap3_profile
from microseq_tests.blast.run_blast import run_blast
from microseq_tests.utility.add_taxonomy import run_taxonomy_join
from microseq_tests.utility.utils import load_config, expand_db_path
from microseq_tests.vsearch_tools import ( 
    collapse_replicates_grouped, 
    chimera_check_ref, 
)
from microseq_tests.utility.id_normaliser import NORMALISERS
from microseq_tests.utility.io_utils import write_fasta_and_qual_from_fastq
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
    return assemble_pairs(Path(input_dir), Path(output_dir), dup_policy=dup_policy, cap3_options=cap3_options, fwd_pattern=fwd_pattern, rev_pattern=rev_pattern, pairing_report=pairing_report, enforce_same_well=enforce_same_well, well_pattern=well_pattern, use_qual=use_qual )

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

def _collect_pairing_catalog(
    directory: Path,
    *,
    fwd_pattern: str | None,
    rev_pattern: str | None,
    dup_policy: DupPolicy,
    enforce_same_well: bool,
    well_pattern: str | re.Pattern[str] | None,
) -> tuple[dict[str, dict[str, list[Path]]], set[str]]:
    """Return paired sample mappings and missing-mate sample IDs."""
    detectors = list(DETECTORS)
    if fwd_pattern and rev_pattern:
        detectors = [make_pattern_detector(fwd_pattern, rev_pattern), *detectors]

    pairs: dict[str, dict[str, list[Path]]] = defaultdict(lambda: {"F": [], "R": []})

    if directory.is_file():
        for record in SeqIO.parse(directory, "fasta"):
            sid, orient = extract_sid_orientation(record.id, detectors=detectors)
            if orient not in ("F", "R"):
                continue
            if not enforce_same_well:
                sid = _strip_well_token(sid)
            well = _extract_well(record.id, pattern=well_pattern) if enforce_same_well else None
            if enforce_same_well and not well:
                continue
            key = well if (enforce_same_well and well) else sid
            pairs[key][orient].append(directory)
    else:
        for suffix in ("*.fasta", "*.fa", "*.fna"):
            for fp in sorted(directory.glob(suffix)):
                sid, orient = extract_sid_orientation(fp.name, detectors=detectors)
                if orient not in ("F", "R"):
                    continue
                if not enforce_same_well:
                    sid = _strip_well_token(sid)
                well = _extract_well(fp.name, pattern=well_pattern) if enforce_same_well else None
                if enforce_same_well and not well:
                    continue
                key = well if (enforce_same_well and well) else sid
                pairs[key][orient].append(fp)

    paired_samples: dict[str, dict[str, list[Path]]] = {}
    missing_samples: set[str] = set()
    for sid, entries in pairs.items():
        f_sources = entries["F"]
        r_sources = entries["R"]
        if not f_sources or not r_sources:
            missing_samples.add(sid)
            continue
        if dup_policy == DupPolicy.KEEP_SEPARATE:
            primer_pairs = _build_keep_separate_pairs(f_sources, r_sources)
            for idx, (fwd, rev) in enumerate(primer_pairs, start=1):
                sample_key = sid if len(primer_pairs) == 1 else f"{sid}_{idx}"
                paired_samples[sample_key] = {"F": [fwd], "R": [rev]}
        else:
            paired_samples[sid] = {"F": f_sources, "R": r_sources}

    return paired_samples, missing_samples


def _read_first_record(path: Path, fmt: str) -> SeqRecord | None:
    for record in SeqIO.parse(path, fmt):
        return record
    return None


def _load_qual_record(qual_path: Path, record_id: str) -> list[int] | None:
    if not qual_path.exists():
        return None
    for qual_rec in SeqIO.parse(qual_path, "qual"):
        if qual_rec.id == record_id:
            return qual_rec.letter_annotations.get("phred_quality")
    return None


def _evaluate_overlap(
    fwd_seq: str,
    rev_seq: str,
    fwd_qual: list[int] | None,
    rev_qual: list[int] | None,
    *,
    min_overlap: int = 100,
    min_identity: float = 0.8,
    min_quality: float = 20.0,
) -> dict[str, float | str | None]:
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.5

    alignment = aligner.align(fwd_seq, rev_seq)[0]
    overlap_len = 0
    matches = 0
    quality_vals: list[int] = []

    for (start1, end1), (start2, end2) in zip(alignment.aligned[0], alignment.aligned[1]):
        segment_len = end1 - start1
        overlap_len += segment_len
        for idx in range(segment_len):
            base1 = fwd_seq[start1 + idx]
            base2 = rev_seq[start2 + idx]
            if base1 == base2:
                matches += 1
            if fwd_qual is not None and rev_qual is not None:
                quality_vals.append(min(fwd_qual[start1 + idx], rev_qual[start2 + idx]))

    identity = matches / overlap_len if overlap_len else 0.0
    overlap_quality = None
    if quality_vals:
        overlap_quality = sum(quality_vals) / len(quality_vals)

    if overlap_len < min_overlap:
        status = "overlap_too_short"
    elif identity < min_identity:
        status = "overlap_identity_low"
    elif overlap_quality is not None and overlap_quality < min_quality:
        status = "overlap_quality_low"
    else:
        status = "ok"

    return {
        "overlap_len": overlap_len,
        "identity": identity,
        "overlap_quality": overlap_quality,
        "status": status,
    }


def _write_overlap_audit(
    paired_samples: dict[str, dict[str, list[Path]]],
    output_tsv: Path,
) -> dict[str, str]:
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    status_map: dict[str, str] = {}
    with output_tsv.open("w", encoding="utf-8") as fh:
        fh.write("sample_id\toverlap_len\toverlap_identity\toverlap_quality\tstatus\n")
        for sample_id, entries in sorted(paired_samples.items()):
            if not entries["F"] or not entries["R"]:
                status = "overlap_missing_reads"
                status_map[sample_id] = status
                fh.write(f"{sample_id}\t0\t0\t\t{status}\n")
                continue
            fwd_path = entries["F"][0]
            rev_path = entries["R"][0]
            if not fwd_path.exists() or not rev_path.exists():
                status = "overlap_missing_reads"
                status_map[sample_id] = status
                fh.write(f"{sample_id}\t0\t0\t\t{status}\n")
                continue
            fwd_record = _read_first_record(fwd_path, "fasta")
            rev_record = _read_first_record(rev_path, "fasta")
            if fwd_record is None or rev_record is None:
                status = "overlap_missing_reads"
                status_map[sample_id] = status
                fh.write(f"{sample_id}\t0\t0\t\t{status}\n")
                continue

            fwd_qual = _load_qual_record(Path(f"{fwd_path}.qual"), fwd_record.id)
            rev_qual = _load_qual_record(Path(f"{rev_path}.qual"), rev_record.id)
            audit = _evaluate_overlap(
                str(fwd_record.seq),
                str(rev_record.seq),
                fwd_qual,
                rev_qual,
            )
            status = str(audit["status"])
            status_map[sample_id] = status
            overlap_quality = audit["overlap_quality"]
            fh.write(
                f"{sample_id}\t{audit['overlap_len']}\t{audit['identity']:.4f}\t"
                f"{'' if overlap_quality is None else f'{overlap_quality:.2f}'}\t{status}\n"
            )
    return status_map


def _build_blast_inputs(
    asm_dir: Path,
    paired_samples: dict[str, dict[str, list[Path]]],
    missing_samples: set[str],
    output_fasta: Path,
    output_tsv: Path,
) -> None:
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    records_to_write: list[SeqIO.SeqRecord] = []
    rows: list[dict[str, str]] = []

    for sample_id in sorted(set(paired_samples) | set(missing_samples)):
        if sample_id in missing_samples:
            rows.append(
                {
                    "sample_id": sample_id,
                    "blast_payload": "pair_missing",
                    "payload_ids": "",
                    "reason": "pair_missing",
                }
            )
            continue

        contigs_path = asm_dir / sample_id / f"{sample_id}_paired.fasta.cap.contigs"
        singlets_path = asm_dir / sample_id / f"{sample_id}_paired.fasta.cap.singlets"

        contigs = list(SeqIO.parse(contigs_path, "fasta")) if contigs_path.exists() else []
        singlets = list(SeqIO.parse(singlets_path, "fasta")) if singlets_path.exists() else []

        if contigs:
            payload = "contig"
            source_records = contigs
            reason = "contigs_present"
            label = "cap3_c"
        elif singlets:
            payload = "singlet"
            source_records = singlets
            reason = "singlets_only"
            label = "cap3_s"
        else:
            payload = "no_payload"
            source_records = []
            reason = "cap3_no_output"
            label = "cap3"

        payload_ids: list[str] = []
        for idx, record in enumerate(source_records, start=1):
            new_id = f"{sample_id}|{payload}|{label}{idx}"
            payload_ids.append(f"{new_id}={record.id}")
            record.id = new_id
            record.name = new_id
            record.description = ""
            records_to_write.append(record)

        rows.append(
            {
                "sample_id": sample_id,
                "blast_payload": payload,
                "payload_ids": ";".join(payload_ids),
                "reason": reason,
            }
        )

    with output_fasta.open("w", encoding="utf-8") as fh:
        SeqIO.write(records_to_write, fh, "fasta")

    with output_tsv.open("w", encoding="utf-8") as fh:
        fh.write("sample_id\tblast_payload\tpayload_ids\treason\n")
        for row in rows:
            fh.write("\t".join([row["sample_id"], row["blast_payload"], row["payload_ids"], row["reason"]]) + "\n")


def _resolve_cap3_options(
    cap3_profile: str,
    cap3_options: Sequence[str] | None,
    cap3_extra_args: Sequence[str] | None,
) -> list[str]:
    args = resolve_cap3_profile(cap3_profile)
    if cap3_options:
        args.extend(cap3_options)
    if cap3_extra_args:
        args.extend(cap3_extra_args)
    return args


def _merge_cap3_contigs(
    files: Sequence[Path],
    destination: Path,
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8") as out:
        for fp in files:
            sid = Path(fp).parent.name
            for idx, rec in enumerate(SeqIO.parse(fp, "fasta"), 1):
                rec.id = f"{sid}|contig|cap3_c{idx}"
                rec.name = rec.id
                rec.description = ""
                SeqIO.write(rec, out, "fasta")

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
    cap3_profile: str = "strict",
    cap3_extra_args: Sequence[str] | None = None,
    cap3_use_qual: bool = True,
    write_blast_inputs: bool = True,
    use_blast_inputs: bool = True,
    fwd_pattern: str | None = None,
    rev_pattern: str | None = None,
    enforce_same_well: bool = False,
    well_pattern: str | re.Pattern[str] | None = None, 
    metadata: Path | None = None,
    summary_tsv: Path | None = None,
    mee_max: float | None = None,
    mee_min_len: int | None = None,
    collapse_replicates: bool = False,
    replicate_id_th: float | None = None,
    min_replicate_size: int | None = None,
    chimera_mode: str = "off",
    chimera_db: PathLike | None = None,
    overlap_audit: bool = False, 
    on_stage=None,
    on_progress=None,
) -> dict[str, Path | None]:
    """
    Run trim -> FASTA merge -> BLAST -> taxonomy (+ optional post‑BLAST).

    infile may be FASTA, FASTQ, a single ``.ab1`` trace, or a directory of
    ``.ab1`` files.  Sanger mode is triggered automatically when *infile* is a
    directory or ends with ``.ab1``.  When *summary_tsv* is ``None`` the summary
    is written to ``qc/trim_summary.tsv`` inside *out_dir*. Set *overlap_audit* 
    to emit ``qc/overlap_audit.tsv`` and annotate paired assembly summaries with
    overlap status codes. CAP3 profiles and extra args can be supplied via
    *cap3_profile* and *cap3_extra_args* to tune paired assembly parameters.
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

    cfg = load_config()
    chimera_db_path: Path | None = None
    if chimera_mode != "off":
        if chimera_db:
            chimera_db_path = Path(chimera_db).expanduser().resolve()
        else:
            chimera_ref = cfg["databases"].get(db_key, {}).get("chimera_ref")
            if chimera_ref:
                chimera_db_path = Path(expand_db_path(chimera_ref))
        if chimera_db_path is None:
            raise ValueError(
                f"No chimera reference configured for '{db_key}'. "
                "Provide --chimera-db or add databases.<key>.chimera_ref."
            )

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
        "chimera_db": chimera_db_path,
        "collapsed_fasta": None,
        "nonchimera_fasta": None,
        "chimera_report": None,
        "replicate_map": None,
        "replicate_weights": None,
        "assembly_summary": out_dir / "asm" / "assembly_summary.tsv" if using_paired else None,
        "overlap_audit": out_dir / "qc" / "overlap_audit.tsv" if (using_paired and overlap_audit) else None,
    }

    extra_stages = int(collapse_replicates) + int(chimera_mode != "off")

    if using_paired:
        # Assembly + BLAST + taxonomy + postblast (optional here if using) 
        n_stages = (3 if is_fasta else 5) + extra_stages + int(postblast)   
    else:
        n_stages = (2 if is_fasta else 4) + extra_stages + int(postblast)
    step = 100 // n_stages
    pct = 0

    def subprog(off):
        return lambda p: on_progress(off + p * step // 100)

    if using_paired:
        def _fastq_to_paired_fastas(
            source_dir: Path,
            dest_dir: Path,
            *,
            use_qual: bool,
        ) -> list[Path]:
            dest_dir.mkdir(parents=True, exist_ok=True)
            written: list[Path] = []
            for fq in sorted(source_dir.rglob("*.fastq")): 
                out_fp = dest_dir / f"{fq.stem}.fasta"
                if use_qual:
                    result = write_fasta_and_qual_from_fastq(fq, out_fp)
                    if result is None:
                        continue
                    written.append(result[0])
                else:
                    records = list(SeqIO.parse(fq, "fastq"))
                    if not records:
                        continue
                    out_fp.parent.mkdir(parents=True, exist_ok=True)
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
            generated_fastas = _fastq_to_paired_fastas(fastq_dir, fasta_dir, use_qual=cap3_use_qual)
            if not generated_fastas:
                raise ValueError(f"No FASTA files generated from {fastq_dir}")
            assembly_input = fasta_dir
        

        pairing_report = out_dir / "qc" / "pairing_report.tsv"
        pairing_report.parent.mkdir(parents=True, exist_ok=True)

        on_stage("Paired assembly")
        cap3_args = _resolve_cap3_options(cap3_profile, cap3_options, cap3_extra_args)
        contig_paths = assemble_pairs(
            assembly_input, 
            out_dir / "asm",
            dup_policy=dup_policy,
            cap3_options=cap3_args,
            fwd_pattern=fwd_pattern,
            rev_pattern=rev_pattern,
            pairing_report=pairing_report,
            enforce_same_well=enforce_same_well,
            well_pattern=well_pattern,
            use_qual=cap3_use_qual
        )
        if not contig_paths:
           summary = _summarize_paired_candidates(assembly_input, fwd_pattern, rev_pattern, enforce_same_well=enforce_same_well, well_pattern=well_pattern)
           suggestions = _suggest_pairing_patterns(assembly_input)
           raise ValueError(
                f"No paired reads detected in {assembly_input}; {summary}. "
                "If your primer names differ, provide --fwd-pattern/--rev-pattern."
                f" {suggestions}"
            ) 
        paired_samples, missing_samples = _collect_pairing_catalog(
            assembly_input,
            fwd_pattern=fwd_pattern,
            rev_pattern=rev_pattern,
            dup_policy=dup_policy,
            enforce_same_well=enforce_same_well,
            well_pattern=well_pattern,
        )
        if pairing_report.exists():
            report_ids = set(
                line.split("\t", 1)[0]
                for line in pairing_report.read_text(encoding="utf-8").splitlines()[1:]
                if line.strip()
            )
            missing_samples.difference_update(report_ids)
            for sid in report_ids:
                paired_samples.setdefault(sid, {"F": [], "R": []})

        overlap_status = None
        if overlap_audit and paths["overlap_audit"]:
            if Path(assembly_input).is_dir():
                overlap_status = _write_overlap_audit(paired_samples, paths["overlap_audit"])
            else:
                L.warning("Overlap audit skipped; paired input is not a directory: %s", assembly_input)

        assembly_summary = paths["assembly_summary"]
        if assembly_summary:
            parse_cap3_reports(
                out_dir / "asm",
                sorted(set(paired_samples.keys()) | set(missing_samples)),
                output_tsv=assembly_summary,
                missing_samples=missing_samples,
                overlap_status=overlap_status,
            )
        merged_contigs = out_dir / "asm" / "paired_contigs.fasta"
        _merge_cap3_contigs(contig_paths, merged_contigs)

        if use_blast_inputs and not write_blast_inputs:
            L.warning("use_blast_inputs requested without write_blast_inputs; enabling blast input output.")
            write_blast_inputs = True

        blast_inputs_fasta = out_dir / "asm" / "blast_inputs.fasta"
        blast_inputs_tsv = out_dir / "asm" / "blast_inputs.tsv"
        if write_blast_inputs:
            _build_blast_inputs(
                out_dir / "asm",
                paired_samples,
                missing_samples,
                blast_inputs_fasta,
                blast_inputs_tsv,
            )

        if use_blast_inputs and write_blast_inputs:
            paths["fasta"] = blast_inputs_fasta
            paths["trimmed_fasta"] = blast_inputs_fasta
        else:
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

    # 3a - VSEARCH (opt in add on) Collapse and chimera stages 
    current_fasta = paths["fasta"]

    def _strip_sizes_and_filter_weights(
        source: Path,
        dest: Path,
        weights_in: Path | None,
        weights_out: Path,
    ) -> None:
        import re
        from Bio import SeqIO

        size_re = re.compile(r";size=(\d+);?")
        weight_map: dict[str, str] = {}
        if weights_in and weights_in.exists():
            for line in weights_in.read_text().splitlines()[1:]:
                if not line.strip():
                    continue
                qseqid, size = line.split("\t", 1)
                weight_map[qseqid] = size

        with dest.open("w", encoding="utf-8") as out_fh, weights_out.open(
            "w", encoding="utf-8"
        ) as weight_fh:
            weight_fh.write("qseqid\treplicate_size\n")
            for rec in SeqIO.parse(source, "fasta"):
                clean_id = size_re.sub("", rec.id)
                rec.id = clean_id
                rec.name = clean_id
                rec.description = ""
                SeqIO.write(rec, out_fh, "fasta")
                size = weight_map.get(clean_id)
                if size:
                    weight_fh.write(f"{clean_id}\t{size}\n")
                else:
                    weight_fh.write(f"{clean_id}\t1\n")

    if collapse_replicates:
        on_stage("Collapse replicates")
        collapsed = out_dir / "qc" / "replicates_collapsed.fasta"
        replicate_map = out_dir / "qc" / "replicate_map.tsv"
        replicate_weights = out_dir / "qc" / "replicate_weights.tsv"

        cfg_norm = cfg.get("metadata", {}).get("sample_id_normaliser", "strip_suffix")
        normaliser = NORMALISERS.get(cfg_norm, NORMALISERS["strip_suffix"])

        collapse_replicates_grouped(
            current_fasta,
            collapsed,
            group_fn=normaliser,
            min_size=min_replicate_size or 1,
            id_th=replicate_id_th,
            threads=threads,
            map_tsv=replicate_map,
            weights_tsv=replicate_weights,
        )
        paths["collapsed_fasta"] = collapsed
        paths["replicate_map"] = replicate_map
        paths["replicate_weights"] = replicate_weights
        current_fasta = collapsed
        pct += step
        on_progress(pct)

    if chimera_mode != "off":
        on_stage("Chimera check")
        nonchimera = out_dir / "qc" / "nonchimeras.fasta"
        report = out_dir / "qc" / "uchime_report.tsv"
        nonchimera, report = chimera_check_ref(
            current_fasta,
            nonchimera,
            reference=chimera_db_path,
            report_tsv=report,
            threads=threads,
            size_in=collapse_replicates,
        )
        paths["nonchimera_fasta"] = nonchimera
        paths["chimera_report"] = report
        current_fasta = nonchimera
        pct += step
        on_progress(pct)

    if collapse_replicates:
        cleaned = out_dir / "qc" / "replicates_clean.fasta"
        replicate_weights = out_dir / "qc" / "replicate_weights.tsv"
        weights_in = paths["replicate_weights"]
        _strip_sizes_and_filter_weights(
            current_fasta,
            cleaned,
            weights_in,
            replicate_weights,
        )
        paths["replicate_weights"] = replicate_weights
        paths["fasta"] = cleaned
    else:
        paths["fasta"] = current_fasta


    # 3b – BLAST
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
    tax_template = cfg["databases"][db_key]["taxonomy"]
    tax_fp = Path(expand_db_path(tax_template))
    on_stage("Taxonomy")
    run_add_tax(paths["hits"], tax_fp, paths["tax"])
    if thr and thr.isInterruptionRequested():
        raise RuntimeError("Cancelled")
    pct += step
    on_progress(pct)

    if using_paired:
        try:
            tax_df = pd.read_csv(paths["tax"], sep="\t")
            if "sample_id" not in tax_df.columns and "qseqid" in tax_df.columns:
                tax_df["sample_id"] = (
                    tax_df["qseqid"].astype(str).str.split("|").str[0]
                )
                tax_df.to_csv(paths["tax"], sep="\t", index=False) 
        except Exception as exc:
            L.warning("Failed to add paired sample_id column to taxonomy table: %s", exc)

    # 5 – Optional post-BLAST
    if postblast:
        if metadata is None:
            raise ValueError("metadata file required for postblast stage")
        on_stage("Post-BLAST")
        run_postblast(paths["tax"], metadata, paths["biom"], weights_tsv=paths["replicate_weights"],)
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
