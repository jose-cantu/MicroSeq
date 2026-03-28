"""
microseq_tests.pipeline
Thin wrappers around the five CLI stages.
Return an int exit-code (0 = success) & raise on fatal errors.
"""

from __future__ import annotations
from pathlib import Path
import os
from contextlib import nullcontext
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import re 
from typing import Callable, Union, Sequence
from dataclasses import dataclass
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
from numpy import identity
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
from microseq_tests.trimming.ab1_to_fastq import ab1_folder_to_fastq as ab1_to_fastq, build_ab1_output_key_map
from microseq_tests.trimming.biopy_trim import trim_folder as biopy_trim, build_trim_summary_row
from microseq_tests.trimming.ab1_qc import summarize_ab1_folder
from microseq_tests.trimming.fastq_to_fasta import (
    fastq_folder_to_fasta as fastq_to_fasta,
)
from microseq_tests.trimming.primer_trim import trim_primer_fastqs, update_trim_summary

from microseq_tests.assembly.de_novo_assembly import de_novo_assembly

from microseq_tests.assembly.paired_assembly import assemble_pairs, _build_keep_separate_pairs, _write_combined_fasta, write_process_logs, _write_pairing_report 

from microseq_tests.assembly.overlap_utils import (
    AlignedOverlapCandidate,
    iter_end_anchored_overlaps,
    select_best_overlap,
    pick_best_identity_candidate,
    rank_feasible_overlaps,
)
from microseq_tests.assembly.overlap_backends import (
    compute_overlap_candidates,
    resolve_overlap_engine,
    resolve_overlap_engine_order,
    resolve_overlap_engine_strategy,
)
from microseq_tests.assembly.pairing import  (
        DETECTORS, DupPolicy, _extract_well, _strip_well_token, extract_sid_orientation, group_pairs, iter_seq_files, make_pattern_detector,) 
from microseq_tests.assembly.cap3_report import parse_cap3_reports
from microseq_tests.assembly.cap3_profiles import resolve_cap3_profile
from microseq_tests.assembly.registry import list_assemblers, get_assembler_spec
from microseq_tests.assembly.two_read_merge import merge_two_reads
from microseq_tests.blast.run_blast import run_blast
from microseq_tests.utility.add_taxonomy import run_taxonomy_join
from microseq_tests.utility.utils import load_config, expand_db_path
from microseq_tests.primer_catalog import trim_presets
from microseq_tests.vsearch_tools import ( 
    collapse_replicates_grouped, 
    chimera_check_ref, 
    orient_reads as vsearch_orient_reads,
)
from microseq_tests.utility.id_normaliser import NORMALISERS, qseqid_to_sample_id 
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
    "stage_paired_fastas_from_fastq_dir",
    "run_full_pipeline",
    "run_compare_assemblers",
    "run_pairing_report",
    "run_assembly_summary",
]

PathLike = Union[str, Path]
L = logging.getLogger(__name__)


_INPUT_SEQ_SUFFIXES = {".ab1", ".fastq", ".fq", ".fasta", ".fa", ".fna", ".fas"}


def _default_run_output_dir(infile: Path) -> Path:
    """Return a sensible default ``*_microseq`` output path for the given input.

    For directory inputs, prefer anchoring output beside the deepest common
    directory containing discovered sequence files. This avoids creating the
    run folder one level above the actual data when users select a container
    folder that only nests a single sample directory.
    """

    if infile.is_file():
        return infile.parent / f"{infile.with_suffix('').name}_microseq"

    seq_files = [
        p for p in sorted(infile.rglob("*"))
        if p.is_file()
        and p.suffix.lower() in _INPUT_SEQ_SUFFIXES
        and not any(part.endswith("_microseq") for part in p.relative_to(infile).parts[:-1])
    ]
    if not seq_files:
        return infile.parent / f"{infile.name}_microseq"

    common = Path(os.path.commonpath([str(p.parent) for p in seq_files]))
    return common.parent / f"{common.name}_microseq"


def _check_cancel(stop_cb: Callable[[], bool] | None = None) -> None:
    if stop_cb and stop_cb():
        raise RuntimeError("Cancelled")
    thr = QThread.currentThread()
    if thr and thr.isInterruptionRequested():
        raise RuntimeError("Cancelled")


# Contract: qseqid is structural payload id (e.g., S1|contig|cap3_c1)
# sample_id is biological sample key (e.g., S1) 
BLAST_INPUTS_HEADERS: list[str] = [
    "sample_id", "blast_payload", "payload_ids", "reason",
    "selected_assembler_id", "selected_assembler_name", "selected_status",
    "payload_kind", "payload_n", "payload_entity_n", "payload_max_len",
    "ambiguity_flag", "safety_flag", "decision_source",
    "review_reason", "warning_flags",
    "structural_hypothesis_n", "hypotheses_with_hits_n", "missing_hits_n",
    "hypothesis_map", "source_id_map",
    "resolution_state", "resolved_hypothesis", "resolution_reason",
    "review_action", "advisory_reason",
    "trace_status", "trace_status_f", "trace_status_r", "trace_flags",
]

@dataclass
class ResolutionContractRow:
    sample_id: str
    review_action: str = "none"
    review_reason: str = ""
    advisory_reason: str = ""
    warning_flags: str = ""
    structural_hypothesis_n: int = 0
    hypotheses_with_hits_n: int = 0
    missing_hits_n: int = 0
    top_labels: str = ""
    resolution_state: str = "needs_review"
    resolved_hypothesis: str = ""
    resolution_reason: str = "taxonomy_missing"
    trace_status: str = "NA"
    trace_flags: str = ""

    def normalize(self) -> "ResolutionContractRow":
        allowed_states = {"unambiguous", "resolved_by_evidence", "needs_review"}
        if self.resolution_state not in allowed_states:
            self.resolution_state = "needs_review"
        allowed_actions = {"none", "queue"}
        if self.review_action not in allowed_actions:
            self.review_action = "queue" if self.resolution_state == "needs_review" else "none"

        flags = sorted({f.strip() for f in str(self.warning_flags or "").split(";") if f.strip()})
        self.warning_flags = ";".join(flags)
        if self.review_action != "queue":
            self.review_reason = ""
        if self.missing_hits_n < 0:
            self.missing_hits_n = 0
        return self

    @staticmethod
    def header() -> list[str]:
        return [
            "sample_id",
            "review_action",
            "review_reason",
            "advisory_reason",
            "warning_flags",
            "structural_hypothesis_n",
            "hypotheses_with_hits_n",
            "missing_hits_n",
            "top_labels",
            "resolution_state",
            "resolved_hypothesis",
            "resolution_reason",
            "trace_status",
            "trace_flags",
        ]

    def to_dict(self) -> dict[str, str]:
        self.normalize()
        return {
            "sample_id": self.sample_id,
            "review_action": self.review_action,
            "review_reason": self.review_reason,
            "advisory_reason": self.advisory_reason,
            "warning_flags": self.warning_flags,
            "structural_hypothesis_n": str(self.structural_hypothesis_n),
            "hypotheses_with_hits_n": str(self.hypotheses_with_hits_n),
            "missing_hits_n": str(self.missing_hits_n),
            "top_labels": self.top_labels,
            "resolution_state": self.resolution_state,
            "resolved_hypothesis": self.resolved_hypothesis,
            "resolution_reason": self.resolution_reason,
            "trace_status": self.trace_status,
            "trace_flags": self.trace_flags,
        }


_PRIMER_PRESETS: dict[str, dict[str, object]] = trim_presets()


def _resolve_trim_policy(config: dict | None = None) -> dict[str, object]:
    cfg = config if config is not None else load_config()
    trim_cfg = cfg.get("trim", {})
    sanger_cfg = trim_cfg.get("sanger", {})
    trimm_cfg = trim_cfg.get("trimmomatic", {})
    return {
        "min_len": int(trim_cfg.get("min_len", 200)),
        "phred": int(trim_cfg.get("phred", 33)),
        "sanger_method": str(sanger_cfg.get("method", "mott")).strip().lower(),
        "sanger_cutoff_q": int(sanger_cfg.get("cutoff_q", 20)),
        "sanger_window_size": int(sanger_cfg.get("window_size", 5)),
        "sanger_per_base_q": int(sanger_cfg.get("per_base_q", 20)),
        "sanger_file_q_threshold": float(sanger_cfg.get("file_q_threshold", 20.0)),
        "trimm_sliding_window_size": int(trimm_cfg.get("sliding_window_size", 5)),
        "trimm_sliding_window_q": int(trimm_cfg.get("sliding_window_q", 20)),
    }


def _normalize_primer_trim_cfg(cfg: dict) -> dict[str, object]:
    primer_cfg = dict(cfg.get("primer_trim", {}))
    mode = str(primer_cfg.get("mode", "")).strip().lower()
    if not mode:
        if "enabled" in primer_cfg:
            mode = "clip" if bool(primer_cfg.get("enabled")) else "off"
        else:
            mode = "off"
    if mode not in {"off", "detect", "clip"}:
        raise ValueError(f"Invalid primer_trim.mode '{mode}'. Expected one of: off, detect, clip.")

    stage = str(primer_cfg.get("stage", "post_quality")).strip().lower()
    if stage not in {"pre_quality", "post_quality"}:
        raise ValueError(
            f"Invalid primer_trim.stage '{stage}'. Expected one of: pre_quality, post_quality."
        )

    preset_configured = str(primer_cfg.get("preset", "")).strip()

    custom_fwd = list(primer_cfg.get("forward_primers", []) or [])
    custom_rev = list(primer_cfg.get("reverse_primers", []) or [])

    if preset_configured == "16S_27F_1492R":
        has_any_custom = bool(custom_fwd or custom_rev)
        if not has_any_custom:
            mode = "off"
        L.warning(
            "primer_trim.preset='16S_27F_1492R' was removed. "
            "Migrating to preset='' and %s. Configure custom synthetic flank sequences for Detect/Clip.",
            f"mode='{mode}'",
        )
        preset_configured = ""

    preset_cfg = _PRIMER_PRESETS.get(preset_configured) if preset_configured else None
    if preset_configured and preset_cfg is None:
        raise ValueError(
            f"Unknown primer_trim.preset '{preset_configured}'. Known presets: {', '.join(sorted(_PRIMER_PRESETS))}."
        )
    has_custom_fwd = bool(custom_fwd)
    has_custom_rev = bool(custom_rev)

    fwd = list(custom_fwd)
    rev = list(custom_rev)
    if preset_cfg and not fwd:
        fwd = list(preset_cfg.get("forward_primers", []) or [])
    if preset_cfg and not rev:
        rev = list(preset_cfg.get("reverse_primers", []) or [])

    if mode == "off":
        primer_source = "off"
        preset_effective = ""
    elif has_custom_fwd and has_custom_rev:
        primer_source = "custom"
        preset_effective = ""
    elif has_custom_fwd or has_custom_rev:
        primer_source = "mixed"
        preset_effective = ""
    elif preset_configured:
        primer_source = "preset"
        preset_effective = preset_configured
    else:
        primer_source = "custom"
        preset_effective = ""

    if mode != "off" and not (preset_effective or fwd or rev):
        raise ValueError(
            "Primer trimming mode is enabled but no primer sequences were provided. "
            "Set primer_trim.preset or forward_primers/reverse_primers."
        )

    return {
        "mode": mode,
        "stage": stage,
        "preset": preset_effective,
        "preset_configured": preset_configured,
        "primer_source": primer_source,
        "forward_primers": fwd,
        "reverse_primers": rev,
        "max_mismatch": int(primer_cfg.get("max_mismatch", preset_cfg.get("max_mismatch", 2) if preset_cfg else 2)),
        "max_search": int(primer_cfg.get("max_search", preset_cfg.get("max_search", 60) if preset_cfg else 60)),
        "max_primer_offset": int(primer_cfg.get("max_primer_offset", preset_cfg.get("max_primer_offset", 10) if preset_cfg else 10)),
        "iupac_mode": bool(primer_cfg.get("iupac_mode", preset_cfg.get("iupac_mode", True) if preset_cfg else True)),
        "post_quality_trim": primer_cfg.get("post_quality_trim", {}),
    }


def _reporting_flags(cfg: dict) -> dict[str, bool]:
    rep = cfg.get("reporting", {})
    return {
        "emit_engine_audit": bool(rep.get("emit_engine_audit", True)),
        "emit_per_sample_merge_report": bool(rep.get("emit_per_sample_merge_report", True)),
    }


def _load_tsv_rows_by_sample(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    rows: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if len(lines) < 2:
        return rows
    headers = lines[0].split("\t")
    try:
        sid_idx = headers.index("sample_id")
    except ValueError:
        return rows
    for line in lines[1:]:
        vals = line.split("\t")
        if len(vals) < len(headers):
            vals.extend([""] * (len(headers) - len(vals)))
        row = dict(zip(headers, vals, strict=False))
        sid = vals[sid_idx]
        if sid:
            rows[sid] = row
    return rows


def _write_tsv_rows_by_sample(
    path: Path,
    rows: dict[str, dict[str, str]],
    *,
    headers: list[str] | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = headers or BLAST_INPUTS_HEADERS
    with path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for sample_id in sorted(rows):
            row = rows[sample_id]
            fh.write("\t".join(str(row.get(c, "")) for c in cols) + "\n")


def _count_fasta_records(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(1 for _ in SeqIO.parse(path, "fasta"))


def _infer_compare_row_flags(status: str, warning: str) -> tuple[str, str, str]:
    status_norm = (status or "").strip().lower()
    warning_norm = (warning or "").strip().lower()
    ambiguity_flag = "1" if status_norm in {"ambiguous_topk", "ambiguous_overlap", "ambiguous_overlap_singlets", "merged_best_guess"} else "0"
    safety_flag = "none"
    if "high_conflict" in warning_norm or "high_conflict" in status_norm:
        safety_flag = "high_conflict"
    review_reason = ""
    if ambiguity_flag == "1":
        review_reason = "ambiguous_payload"
    return ambiguity_flag, safety_flag, review_reason


# ───────────────────────────────────────────────────────── trimming
def run_trim(
    input_path: PathLike,
    workdir: PathLike,
    sanger: bool = False,
    *,
    summary_tsv: PathLike | None = None,
    link_raw: bool = False,
    primer_cfg_override: dict | None = None,
    trace_qc_flags: str = "auto",

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
    cfg = load_config()
    if primer_cfg_override:
        merged = dict(cfg.get("primer_trim", {}))
        merged.update(primer_cfg_override)
        cfg["primer_trim"] = merged
    trim_policy = _resolve_trim_policy(cfg)
    policy_desc = (
        f"sanger={trim_policy['sanger_method']}(cutoff={trim_policy['sanger_cutoff_q']},window={trim_policy['sanger_window_size']},q={trim_policy['sanger_per_base_q']}); "
        f"fastq=trimmomatic(window={trim_policy['trimm_sliding_window_size']}:{trim_policy['trimm_sliding_window_q']},minlen={trim_policy['min_len']})"
    )
    L.info("Trim policy: %s", policy_desc)
    work = Path(workdir)
    work.mkdir(parents=True, exist_ok=True)
    (work / "qc").mkdir(parents=True, exist_ok=True)

    primer_cfg = _normalize_primer_trim_cfg(cfg)

    pre_quality_primer_results: dict[str, object] = {}

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

        quality_input_dir = fastq_dir
        if primer_cfg["mode"] != "off" and primer_cfg["stage"] == "pre_quality":
            pre_quality_primer_dir = work / "raw_fastq_primer_trim"
            pre_quality_primer_results = trim_primer_fastqs(
                fastq_dir,
                pre_quality_primer_dir,
                forward_primers=primer_cfg["forward_primers"],
                reverse_primers=primer_cfg["reverse_primers"],
                max_mismatch=int(primer_cfg["max_mismatch"]),
                max_search=int(primer_cfg["max_search"]),
                max_primer_offset=int(primer_cfg["max_primer_offset"]),
                iupac_mode=bool(primer_cfg["iupac_mode"]),
                report_path=(work / "qc" / "primer_trim_report.tsv") if primer_cfg["mode"] == "clip" else None,
                detect_report_path=work / "qc" / "primer_detect_report.tsv",
                mode=str(primer_cfg["mode"]),
            )
            if primer_cfg["mode"] == "clip":
                quality_input_dir = pre_quality_primer_dir

        apply_thresholds: bool | None
        if trace_qc_flags == "on":
            apply_thresholds = True
        elif trace_qc_flags == "off":
            apply_thresholds = False
        else:
            apply_thresholds = None

        trace_qc = summarize_ab1_folder(dst, apply_thresholds=apply_thresholds)
        trace_ab1_by_key = {key: ab1 for ab1, key in build_ab1_output_key_map(dst).items()}

        biopy_trim(
            quality_input_dir,
            work / "qc",
            combined_tsv=summary_tsv,
            trace_qc=trace_qc,
            trace_ab1_by_key=trace_ab1_by_key,
            trace_qc_apply_thresholds=apply_thresholds,
            method=str(trim_policy["sanger_method"]),
            cutoff_q=int(trim_policy["sanger_cutoff_q"]),
            window_size=int(trim_policy["sanger_window_size"]),
            per_base_q=int(trim_policy["sanger_per_base_q"]),
            file_q_threshold=float(trim_policy["sanger_file_q_threshold"]),
            min_len=int(trim_policy["min_len"]),
        )
        trim_dir = work / "passed_qc_fastq"
        if summary_tsv and primer_cfg["mode"] == "clip" and primer_cfg["stage"] == "pre_quality":
            update_trim_summary(Path(summary_tsv), pre_quality_primer_results)
    else:
        input_path = Path(input_path) 
        trim_dir = work / "qc"
        trim_fastq_inputs(input_path, trim_dir, summary_tsv=summary_tsv) 

    if primer_cfg["mode"] != "off" and primer_cfg["stage"] == "post_quality":
        primer_trim_dir = work / "passed_qc_fastq_primer_trim"
        primer_detect_report = work / "qc" / "primer_detect_report.tsv"
        primer_report = work / "qc" / "primer_trim_report.tsv"

        def _orientation_resolver(filename: str) -> str | None:
            _sid, orient = extract_sid_orientation(filename)
            return orient

        post_cfg = primer_cfg.get("post_quality_trim", {}) if isinstance(primer_cfg.get("post_quality_trim", {}), dict) else {}
        primer_results = trim_primer_fastqs(
            trim_dir,
            primer_trim_dir,
            forward_primers=primer_cfg["forward_primers"],
            reverse_primers=primer_cfg["reverse_primers"],
            max_mismatch=int(primer_cfg["max_mismatch"]),
            max_search=int(primer_cfg["max_search"]),
            max_primer_offset=int(primer_cfg["max_primer_offset"]),
            iupac_mode=bool(primer_cfg["iupac_mode"]),
            report_path=primer_report if primer_cfg["mode"] == "clip" else None,
            detect_report_path=primer_detect_report,
            orientation_resolver=_orientation_resolver,
            mode=str(primer_cfg["mode"]),
            post_quality_trim_enabled=bool(post_cfg.get("enabled", False)) and primer_cfg["mode"] == "clip",
            post_quality_method=str(post_cfg.get("method", trim_policy["sanger_method"])),
            post_quality_cutoff_q=int(post_cfg.get("cutoff_q", trim_policy["sanger_cutoff_q"])),
            post_quality_window_size=int(post_cfg.get("window_size", trim_policy["sanger_window_size"])),
            post_quality_per_base_q=int(post_cfg.get("per_base_q", trim_policy["sanger_per_base_q"])),
            post_quality_rescue_5prime_bases=int(post_cfg.get("rescue_5prime_bases", 0)),
            post_quality_min_len=int(cfg.get("trim", {}).get("min_len", 200)),
        )
        if summary_tsv and primer_cfg["mode"] == "clip":
            update_trim_summary(Path(summary_tsv), primer_results)
        if primer_cfg["mode"] == "clip":
            trim_dir = primer_trim_dir
        else:
            L.info("Primer detect mode active; reads were not clipped.")


    if summary_tsv:
        summary_path = Path(summary_tsv)
        if summary_path.exists():
            with summary_path.open("a", encoding="utf-8") as sfh:
                sfh.write("\t".join(build_trim_summary_row(
                    file_name="_policy",
                    reads=0,
                    avg_len=0.0,
                    avg_q=0.0,
                    qc_status="policy",
                    trim_policy=(
                        f"{policy_desc}; primer_mode={primer_cfg['mode']}; "
                        f"primer_stage={primer_cfg['stage']}; primer_preset={primer_cfg['preset'] or 'custom'}; primer_source={primer_cfg.get('primer_source','custom')}"
                    ),
                )) + "\n")

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

def run_paired_assembly(input_dir: PathLike, output_dir: PathLike, *, dup_policy: DupPolicy = DupPolicy.ERROR, cap3_options=None, fwd_pattern: str | None = None, rev_pattern: str | None = None, pairing_report: PathLike | None = None, enforce_same_well: bool = False, well_pattern: str | re.Pattern[str] | None = None, use_qual: bool = True, on_stage=None, on_progress=None) -> list[Path]:
    return assemble_pairs(Path(input_dir), Path(output_dir), dup_policy=dup_policy, cap3_options=cap3_options, fwd_pattern=fwd_pattern, rev_pattern=rev_pattern, pairing_report=pairing_report, enforce_same_well=enforce_same_well, well_pattern=well_pattern, use_qual=use_qual, on_stage=on_stage, on_progress=on_progress )

def stage_paired_fastas_from_fastq_dir(input_fastq_dir: PathLike, output_fasta_dir: PathLike, *, use_qual: bool = True,) -> list[Path]:
    """
    Converting each FASTQ in input_dir (recursively) into one FASTA in output_dir. 
    Returns list of written FASTA paths. Writes .qual alongside FASTA when 
    use_qual=True. 
    """
    input_fastq_dir = Path(input_fastq_dir)
    output_fasta_dir = Path(output_fasta_dir)
    output_fasta_dir.mkdir(parents=True, exist_ok=True)

    written_fasta_paths: list[Path] = [] 
    for input_fastq_path in sorted(input_fastq_dir.rglob("*.fastq")):
        output_fasta_path = output_fasta_dir / f"{input_fastq_path.stem}.fasta" 

        if use_qual:
            fasta_and_qual_paths = write_fasta_and_qual_from_fastq(
                input_fastq_path,
                output_fasta_path,
            )
            if fasta_and_qual_paths is None:
                continue 
            written_fasta_paths.append(fasta_and_qual_paths[0])
        else:
            records = list(SeqIO.parse(input_fastq_path, "fastq"))
            if not records:
                continue
            SeqIO.write(records, output_fasta_path, "fasta")
            written_fasta_paths.append(output_fasta_path)

    return written_fasta_paths  

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
    stop_cb: Callable[[], bool] | None = None,
) -> int:
    from microseq_tests.blast.run_blast import BlastOptions # local import keeps AB1 safe 
     

    _check_cancel(stop_cb)


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


    _check_cancel(stop_cb)
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
    stop_cb: Callable[[], bool] | None = None,
    **kwargs,
) -> int:
    _check_cancel(stop_cb)
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
    _check_cancel(stop_cb)
    return 0


# # ────────────────────────────────────────────────── full workflow here =)

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

    for fp in iter_seq_files(directory):
        sid, orient = extract_sid_orientation(fp.name, detectors=detectors)
        well = _extract_well(fp.name, pattern=well_pattern) if enforce_same_well else None 
        if orient in ("F", "R"):
            if not enforce_same_well:
                sid = _strip_well_token(sid)
            key = sid if not enforce_same_well else (well or sid)
            counts[orient] += 1 
            sample_orients[key].add(orient)
            if enforce_same_well:
                if well:
                    sid_wells[sid].add(orient)
                else:
                    missing_well.add(fp.name)
            continue 
        counts["unknown"]  += 1 
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

def _suggest_pairing_patterns_for_suffixes(
    directory: Path,
    *,
    suffixes: Sequence[str],
    context_label: str,
    fallback_message: str,
) -> str:
    """Suggest --fwd-pattern/--rev-pattern examples by scanning filename tokens."""

    token_rx = re.compile(r"([A-Za-z0-9]+[FR])", re.I)
    fwd_tokens: Counter[str] = Counter()
    rev_tokens: Counter[str] = Counter()
    names: list[str] = []

    allowed_suffixes = {s.lower() for s in suffixes}
    files: list[Path] = []
    if directory.is_file():
        if directory.suffix.lower() in allowed_suffixes:
            files = [directory]
    else:
        for suffix in allowed_suffixes:
            files.extend(sorted(directory.rglob(f"*{suffix}")))

    for fp in sorted(set(files)):
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
    if fwd or rev:
        suggestions.append(
            f"{context_label} Example flag: --fwd-pattern "
            f"'{fwd or 'FORWARD_TOKEN'}' --rev-pattern '{rev or 'REVERSE_TOKEN'}'"
        )
        suggestions.append(
            "Regex variant (case-insensitive): "
            f"--fwd-pattern '(?i){fwd or '27F'}' --rev-pattern '(?i){rev or '1492R'}'"
        )

    if names:
        rename_hint = (
            f"{context_label} Rename option: include forward/reverse primer labels in filenames "
            f"(e.g., sample_{fwd or '27F'}.fastq / sample_{rev or '1492R'}.ab1); "
            "examples seen: " + ", ".join(names[:3])
        )
        suggestions.append(rename_hint)

    if not suggestions:
        suggestions.append(fallback_message)

    return " ".join(suggestions)


def _suggest_pairing_patterns(directory: Path) -> str:
    """Suggest pairing labels for raw GUI inputs (AB1/FASTQ)."""

    return _suggest_pairing_patterns_for_suffixes(
        directory,
        suffixes=(".fastq", ".ab1"),
        context_label="Raw input naming suggestions:",
        fallback_message=(
            "Raw input naming suggestions: No forward/reverse primer labels detected in AB1/FASTQ filenames; "
            "rename files to include labels (e.g., sample_27F / sample_1492R) or pass --fwd-pattern/--rev-pattern explicitly."
        ),
    )


def _suggest_pairing_patterns_staged(directory: Path) -> str:
    """Provide staged/internal pairing guidance without raw-input fallback wording."""

    return _suggest_pairing_patterns_for_suffixes(
        directory,
        suffixes=(".fasta", ".fa", ".fna", ".fas"),
        context_label="Staged/internal pairing status:",
        fallback_message=(
            "Staged/internal pairing status: Pairing suggestions are based on raw input filenames; "
            "staged FASTA appears downstream of conversion. Use --fwd-pattern/--rev-pattern if pairing labels are non-standard."
        ),
    )

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


def _resolve_overlap_thresholds(cfg: dict | None = None) -> tuple[int, float, float]:
    config = cfg if cfg is not None else load_config()
    overlap_cfg = config.get("overlap_eval", {})
    merge_cfg = config.get("merge_two_reads", {})
    configured_overlap_engine = str(merge_cfg.get("overlap_engine", "ungapped")).strip().lower()
    overlap_engine = resolve_overlap_engine(configured_overlap_engine)
    anchor_tolerance_bases = int(merge_cfg.get("anchor_tolerance_bases", 30))
    min_overlap = int(overlap_cfg.get("min_overlap", merge_cfg.get("min_overlap", 100)))
    min_identity = float(overlap_cfg.get("min_identity", merge_cfg.get("min_identity", 0.8)))
    min_quality = float(overlap_cfg.get("min_quality", 20.0))
    return min_overlap, min_identity, min_quality




def _resolve_ambiguity_thresholds(cfg: dict | None = None) -> tuple[float, float]:
    config = cfg if cfg is not None else load_config()
    overlap_cfg = config.get("overlap_eval", {})
    merge_cfg = config.get("merge_two_reads", {})
    identity_delta = float(overlap_cfg.get("ambiguity_identity_delta", merge_cfg.get("ambiguity_identity_delta", 0.0025)))
    quality_epsilon = float(overlap_cfg.get("ambiguity_quality_epsilon", merge_cfg.get("ambiguity_quality_epsilon", 0.1)))
    return max(0.0, identity_delta), max(0.0, quality_epsilon)
def _classify_overlap_status(
    overlap_len: int,
    identity: float,
    overlap_quality: float | None,
    *,
    min_overlap: int,
    min_identity: float,
    min_quality: float,
    quality_mode: str = "warning",
) -> str:
    if overlap_len < min_overlap:
        return "overlap_too_short"
    if identity < min_identity:
        return "overlap_identity_low"
    if quality_mode == "blocking" and (overlap_quality is None or overlap_quality < min_quality):
        return "overlap_quality_low"
    return "ok"

def _collect_primer_trim_bases_by_sample(report_path: Path | None) -> dict[str, dict[str, int]]:
    """Aggregate primer-trimmed bases from primer_trim_report.tsv by sample/orientation."""
    if report_path is None or not report_path.exists():
        return {}
    rows = pd.read_csv(report_path, sep="\t")
    out: dict[str, dict[str, int]] = defaultdict(lambda: {"F": 0, "R": 0})
    for _, row in rows.iterrows():
        raw_name = row.get("file", "")
        if pd.isna(raw_name):
            continue
        name = str(raw_name).strip()
        if not name:
            continue
        parsed = extract_sid_orientation(Path(name).name)
        if not parsed:
            continue
        sample_id, orient = parsed
        if orient not in {"F", "R"}:
            continue
        try:
            bases_trimmed = int(float(row.get("bases_trimmed", 0)))
        except (TypeError, ValueError):
            bases_trimmed = 0
        out[sample_id][orient] += max(0, bases_trimmed)
    return dict(out)


def _write_overlap_audit(
    paired_samples: dict[str, dict[str, list[Path]]],
    output_tsv: Path,
    *,
    min_overlap: int | None = None,
    min_identity: float | None = None,
    min_quality: float | None = None,
    quality_mode: str = "warning",
    emit_engine_audit: bool = True,
    pretrim_paired_samples: dict[str, dict[str, list[Path]]] | None = None,
    primer_trim_bases_by_sample: dict[str, dict[str, int]] | None = None,
) -> dict[str, str]:
    config = load_config()
    merge_cfg = config.get("merge_two_reads", {})
    configured_overlap_engine = str(merge_cfg.get("overlap_engine", "ungapped")).strip().lower()
    overlap_engine = resolve_overlap_engine(configured_overlap_engine)
    overlap_strategy = resolve_overlap_engine_strategy(str(merge_cfg.get("overlap_engine_strategy", "single")))
    raw_order = merge_cfg.get("overlap_engine_order", ["ungapped", "biopython", "edlib"])
    if isinstance(raw_order, str):
        raw_order = [tok.strip() for tok in raw_order.split(",") if tok.strip()]
    if not isinstance(raw_order, (list, tuple)):
        raw_order = ["ungapped", "biopython", "edlib"]
    engine_order = resolve_overlap_engine_order(raw_order)

    anchor_tolerance_bases = int(merge_cfg.get("anchor_tolerance_bases", 30))

    if min_overlap is None or min_identity is None or min_quality is None:
        cfg_overlap, cfg_identity, cfg_quality = _resolve_overlap_thresholds(config)
        min_overlap = cfg_overlap if min_overlap is None else min_overlap
        min_identity = cfg_identity if min_identity is None else min_identity
        min_quality = cfg_quality if min_quality is None else min_quality
    ambiguity_identity_delta, ambiguity_quality_epsilon = _resolve_ambiguity_thresholds(config)

    def _blank_diag() -> dict[str, str]:
        return {
            "fwd_best_identity": "",
            "revcomp_best_identity": "",
            "fwd_best_overlap_len": "",
            "revcomp_best_overlap_len": "",
            "fwd_anchor_feasible": "no",
            "revcomp_anchor_feasible": "no",
            "identity_delta_revcomp_minus_fwd": "",
            "selected_vs_best_identity_delta": "",
            "top_candidate_count": "",
            "top2_identity_delta": "",
            "top2_overlap_len_delta": "",
            "top2_quality_delta": "",
            "tie_reason_code": "",
            "pretrim_best_identity": "",
            "posttrim_best_identity": "",
            "fwd_best_identity_any": "",
            "revcomp_best_identity_any": "",
            "fwd_best_overlap_len_any": "",
            "revcomp_best_overlap_len_any": "",
            "pretrim_best_overlap_len": "",
            "posttrim_selected_overlap_len": "",
            "pretrim_status": "",
            "posttrim_status": "",
            "ambiguity_identity_delta_used": "",
            "ambiguity_quality_epsilon_used": "",
            "primer_trim_bases_fwd": "0",
            "primer_trim_bases_rev": "0",
        }

    def _format_missing(sample_id: str, status: str) -> str:
        d = _blank_diag()
        return (
            f"{sample_id}\t0\t0\t\t\t{status}\t0\t\tno\tno"
            f"\t{d['fwd_best_identity']}\t{d['revcomp_best_identity']}\t{d['fwd_best_overlap_len']}\t{d['revcomp_best_overlap_len']}"
            f"\t{d['fwd_anchor_feasible']}\t{d['revcomp_anchor_feasible']}\t{d['identity_delta_revcomp_minus_fwd']}"
            f"\t{d['selected_vs_best_identity_delta']}\t{d['top_candidate_count']}\t{d['top2_identity_delta']}"
            f"\t{d['top2_overlap_len_delta']}\t{d['top2_quality_delta']}\t{d['tie_reason_code']}"
            f"\t{d['fwd_best_identity_any']}\t{d['revcomp_best_identity_any']}\t{d['fwd_best_overlap_len_any']}\t{d['revcomp_best_overlap_len_any']}"
            f"\t{d['pretrim_best_identity']}\t{d['posttrim_best_identity']}\t{d['pretrim_best_overlap_len']}\t{d['posttrim_selected_overlap_len']}"
            f"\t{d['pretrim_status']}\t{d['posttrim_status']}\t{d['ambiguity_identity_delta_used']}\t{d['ambiguity_quality_epsilon_used']}"
            f"\t{d['primer_trim_bases_fwd']}\t{d['primer_trim_bases_rev']}"
            f"\t{overlap_engine}\tno\t{configured_overlap_engine}\n"
        )

    def _best_by_orientation(cands: list[AlignedOverlapCandidate], orientation: str) -> AlignedOverlapCandidate | None:
        oriented = [c for c in cands if c.orientation == orientation and c.overlap_len >= min_overlap]
        if not oriented:
            return None
        return max(
            oriented,
            key=lambda c: (c.identity, c.overlap_len, -c.mismatches, c.overlap_quality if c.overlap_quality is not None else float("-inf")),
        )

    def _best_any_by_orientation(cands: list[AlignedOverlapCandidate], orientation: str) -> AlignedOverlapCandidate | None:
        oriented = [c for c in cands if c.orientation == orientation]
        if not oriented:
            return None
        return max(
            oriented,
            key=lambda c: (c.identity, c.overlap_len, -c.mismatches, c.overlap_quality if c.overlap_quality is not None else float("-inf")),
        )

    def _orientation_anchor_feasible(cands: list[AlignedOverlapCandidate], orientation: str) -> bool:
        return any(
            c.orientation == orientation
            and getattr(c, "end_anchored", True)
            and c.overlap_len >= min_overlap
            and c.identity >= min_identity
            and (
                quality_mode != "blocking"
                or (c.overlap_quality is not None and c.overlap_quality >= min_quality)
            )
            for c in cands
        )

    def _compute_diag_fields(
        overlap_candidates: list[AlignedOverlapCandidate],
        selected_overlap,
        best_identity_val: float,
        selected_status: str,
    ) -> dict[str, str]:
        d = _blank_diag()
        if not overlap_candidates:
            return d
        fwd_best = _best_by_orientation(overlap_candidates, "forward")
        rev_best = _best_by_orientation(overlap_candidates, "revcomp")
        fwd_best_any = _best_any_by_orientation(overlap_candidates, "forward")
        rev_best_any = _best_any_by_orientation(overlap_candidates, "revcomp")
        fwd_best_identity = fwd_best.identity if fwd_best is not None else None
        rev_best_identity = rev_best.identity if rev_best is not None else None
        d["fwd_best_identity"] = "" if fwd_best_identity is None else f"{fwd_best_identity:.4f}"
        d["revcomp_best_identity"] = "" if rev_best_identity is None else f"{rev_best_identity:.4f}"
        d["fwd_best_overlap_len"] = "" if fwd_best is None else str(fwd_best.overlap_len)
        d["revcomp_best_overlap_len"] = "" if rev_best is None else str(rev_best.overlap_len)
        d["fwd_anchor_feasible"] = "yes" if _orientation_anchor_feasible(overlap_candidates, "forward") else "no"
        d["revcomp_anchor_feasible"] = "yes" if _orientation_anchor_feasible(overlap_candidates, "revcomp") else "no"
        d["fwd_best_identity_any"] = "" if fwd_best_any is None else f"{fwd_best_any.identity:.4f}"
        d["revcomp_best_identity_any"] = "" if rev_best_any is None else f"{rev_best_any.identity:.4f}"
        d["fwd_best_overlap_len_any"] = "" if fwd_best_any is None else str(fwd_best_any.overlap_len)
        d["revcomp_best_overlap_len_any"] = "" if rev_best_any is None else str(rev_best_any.overlap_len)
        if fwd_best_identity is not None and rev_best_identity is not None:
            d["identity_delta_revcomp_minus_fwd"] = f"{(rev_best_identity - fwd_best_identity):.4f}"
        d["selected_vs_best_identity_delta"] = f"{(best_identity_val - selected_overlap.identity):.4f}"

        ranked = rank_feasible_overlaps(
            overlap_candidates,
            min_overlap=min_overlap,
            min_identity=min_identity,
            min_quality=min_quality,
            quality_mode=quality_mode,
        )
        if len(ranked) > 1:
            top1, top2 = ranked[0], ranked[1]
            d["top2_identity_delta"] = f"{abs(top1.identity - top2.identity):.4f}"
            d["top2_overlap_len_delta"] = str(abs(top1.overlap_len - top2.overlap_len))
            if top1.overlap_quality is not None and top2.overlap_quality is not None:
                d["top2_quality_delta"] = f"{abs(top1.overlap_quality - top2.overlap_quality):.4f}"
            if selected_status == "ambiguous_overlap":
                reason_parts: list[str] = []
                if top1.overlap_len == top2.overlap_len:
                    reason_parts.append("len_equal")
                if top1.mismatches == top2.mismatches:
                    reason_parts.append("mismatch_equal")
                if abs(top1.identity - top2.identity) <= ambiguity_identity_delta:
                    reason_parts.append("identity_eps")
                if (
                    top1.overlap_quality is not None
                    and top2.overlap_quality is not None
                    and abs(top1.overlap_quality - top2.overlap_quality) <= ambiguity_quality_epsilon
                ):
                    reason_parts.append("quality_eps")
                d["tie_reason_code"] = ";".join(reason_parts)
                d["top_candidate_count"] = str(
                    sum(
                        1
                        for cand in ranked
                        if cand.overlap_len == top1.overlap_len
                        and cand.mismatches == top1.mismatches
                        and (
                            abs(top1.identity - cand.identity) <= ambiguity_identity_delta
                            or (
                                top1.overlap_quality is not None
                                and cand.overlap_quality is not None
                                and abs(top1.overlap_quality - cand.overlap_quality) <= ambiguity_quality_epsilon
                            )
                        )
                    )
                )
        d["posttrim_best_identity"] = f"{best_identity_val:.4f}"
        d["posttrim_selected_overlap_len"] = str(selected_overlap.overlap_len)
        d["posttrim_status"] = selected_status
        d["ambiguity_identity_delta_used"] = f"{ambiguity_identity_delta:.4f}"
        d["ambiguity_quality_epsilon_used"] = f"{ambiguity_quality_epsilon:.4f}"
        return d

    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    engines_tsv = output_tsv.with_name("overlap_audit_engines.tsv")
    status_map: dict[str, str] = {}
    with output_tsv.open("w", encoding="utf-8") as fh, (
        engines_tsv.open("w", encoding="utf-8") if emit_engine_audit else nullcontext()
    ) as efh:
        fh.write(
            "sample_id	overlap_len	overlap_identity	overlap_quality	orientation	status"
            "	best_identity	best_identity_orientation	anchoring_feasible	end_anchored_possible"
            "	fwd_best_identity	revcomp_best_identity	fwd_best_overlap_len	revcomp_best_overlap_len"
            "	fwd_anchor_feasible	revcomp_anchor_feasible	identity_delta_revcomp_minus_fwd"
            "	selected_vs_best_identity_delta	top_candidate_count	top2_identity_delta"
            "	top2_overlap_len_delta	top2_quality_delta	tie_reason_code"
            "	fwd_best_identity_any	revcomp_best_identity_any	fwd_best_overlap_len_any	revcomp_best_overlap_len_any"
            "	pretrim_best_identity	posttrim_best_identity	pretrim_best_overlap_len	posttrim_selected_overlap_len"
            "	pretrim_status	posttrim_status	ambiguity_identity_delta_used	ambiguity_quality_epsilon_used"
            "	primer_trim_bases_fwd	primer_trim_bases_rev"
            "	selected_engine	fallback_used	configured_engine\n"
        )
        if efh:
            efh.write(
                "sample_id	engine	orientation	status	overlap_len	identity	overlap_quality"
                "	end_anchored	anchoring_feasible	selected_for_merge\n"
            )
        for sample_id, entries in sorted(paired_samples.items()):
            if not entries["F"] or not entries["R"]:
                status = "overlap_missing_reads"
                status_map[sample_id] = status
                fh.write(_format_missing(sample_id, status))
                continue
            fwd_path = entries["F"][0]
            rev_path = entries["R"][0]
            if not fwd_path.exists() or not rev_path.exists():
                status = "overlap_missing_reads"
                status_map[sample_id] = status
                fh.write(_format_missing(sample_id, status))
                continue
            fwd_record = _read_first_record(fwd_path, "fasta")
            rev_record = _read_first_record(rev_path, "fasta")
            if fwd_record is None or rev_record is None:
                status = "overlap_missing_reads"
                status_map[sample_id] = status
                fh.write(_format_missing(sample_id, status))
                continue

            fwd_qual = _load_qual_record(Path(f"{fwd_path}.qual"), fwd_record.id)
            rev_qual = _load_qual_record(Path(f"{rev_path}.qual"), rev_record.id)

            engines_to_run = [overlap_engine] if overlap_strategy == "single" else engine_order
            selected_engine = engines_to_run[0]
            selected_overlap = None
            best_identity_val = 0.0
            best_identity_orientation = ""
            per_engine_rows: list[tuple[str, object, str, bool, bool, list[AlignedOverlapCandidate]]] = []

            for idx, engine_name in enumerate(engines_to_run):
                if engine_name == "ungapped":
                    overlap_candidates = list(iter_end_anchored_overlaps(
                        str(fwd_record.seq),
                        str(rev_record.seq),
                        fwd_qual,
                        rev_qual,
                    ))
                else:
                    backend = compute_overlap_candidates(
                        str(fwd_record.seq),
                        str(rev_record.seq),
                        fwd_qual,
                        rev_qual,
                        engine=engine_name,
                        end_anchor_tolerance=anchor_tolerance_bases,
                    )
                    overlap_candidates = [
                        AlignedOverlapCandidate(
                            orientation=r.orientation,
                            overlap_len=r.overlap_len,
                            mismatches=r.mismatches,
                            identity=r.identity,
                            overlap_quality=r.overlap_quality,
                            aligned_fwd=r.aligned_fwd,
                            aligned_rev=r.aligned_rev,
                            indels=r.indels,
                            cigar=r.cigar,
                            overlap_span_cols=r.overlap_span_cols,
                            terminal_gap_cols=r.terminal_gap_cols,
                            internal_indels=r.internal_indels,
                            fwd_overlap_start=r.fwd_overlap_start,
                            fwd_overlap_end=r.fwd_overlap_end,
                            rev_overlap_start=r.rev_overlap_start,
                            rev_overlap_end=r.rev_overlap_end,
                            end_anchored=r.end_anchored,
                        )
                        for r in backend
                    ]
                overlap = select_best_overlap(
                    overlap_candidates,
                    min_overlap=min_overlap,
                    min_identity=min_identity,
                    min_quality=min_quality,
                    quality_mode=quality_mode,
                    ambiguity_identity_delta=ambiguity_identity_delta,
                    ambiguity_quality_epsilon=ambiguity_quality_epsilon,
                )
                best_identity = pick_best_identity_candidate(overlap_candidates, min_overlap=min_overlap)
                end_anchored_possible = any(getattr(c, "end_anchored", True) for c in overlap_candidates)
                anchoring_feasible = any(
                    getattr(c, "end_anchored", True)
                    and c.overlap_len >= min_overlap
                    and c.identity >= min_identity
                    and (
                        quality_mode != "blocking"
                        or (c.overlap_quality is not None and c.overlap_quality >= min_quality)
                    )
                    for c in overlap_candidates
                )
                if overlap.status in {"not_end_anchored", "ambiguous_overlap"}:
                    status = overlap.status
                else:
                    status = _classify_overlap_status(
                        overlap.overlap_len,
                        overlap.identity,
                        overlap.overlap_quality,
                        min_overlap=min_overlap,
                        min_identity=min_identity,
                        min_quality=min_quality,
                        quality_mode=quality_mode,
                    )
                per_engine_rows.append((engine_name, overlap, status, end_anchored_possible, anchoring_feasible, overlap_candidates))

                if selected_overlap is None:
                    selected_engine = engine_name
                    selected_overlap = overlap
                    best_identity_val = overlap.identity if best_identity is None else best_identity.identity
                    best_identity_orientation = overlap.orientation if best_identity is None else best_identity.orientation
                if overlap_strategy == "cascade" and status == "ok":
                    selected_engine = engine_name
                    selected_overlap = overlap
                    best_identity_val = overlap.identity if best_identity is None else best_identity.identity
                    best_identity_orientation = overlap.orientation if best_identity is None else best_identity.orientation
                    break
                if overlap_strategy == "all" and selected_overlap is not None:
                    prev_idx = engines_to_run.index(selected_engine)
                    prev_status = next(r[2] for r in per_engine_rows if r[0] == selected_engine)
                    prev_key = (1 if prev_status == "ok" else 0, selected_overlap.identity, selected_overlap.overlap_len, -prev_idx)
                    curr_key = (1 if status == "ok" else 0, overlap.identity, overlap.overlap_len, -idx)
                    if curr_key > prev_key:
                        selected_engine = engine_name
                        selected_overlap = overlap
                        best_identity_val = overlap.identity if best_identity is None else best_identity.identity
                        best_identity_orientation = overlap.orientation if best_identity is None else best_identity.orientation

            row = next((r for r in per_engine_rows if r[0] == selected_engine), None)
            selected_status = row[2] if row else "overlap_missing_reads"
            end_anchored_possible_selected = row[3] if row else False
            anchoring_feasible_selected = row[4] if row else False
            selected_candidates = row[5] if row else []
            status_map[sample_id] = selected_status

            if selected_overlap is None:
                selected_overlap = select_best_overlap([], min_overlap=min_overlap, min_identity=min_identity, min_quality=min_quality, quality_mode=quality_mode, ambiguity_identity_delta=ambiguity_identity_delta, ambiguity_quality_epsilon=ambiguity_quality_epsilon)

            fallback_used = "yes" if (overlap_strategy in {"cascade", "all"} and selected_engine != engines_to_run[0]) else "no"
            overlap_quality = selected_overlap.overlap_quality
            diag = _compute_diag_fields(selected_candidates, selected_overlap, best_identity_val, selected_status)

            pre_entries = pretrim_paired_samples.get(sample_id) if pretrim_paired_samples else None
            if pre_entries and pre_entries.get("F") and pre_entries.get("R"):
                pre_fwd = pre_entries["F"][0]
                pre_rev = pre_entries["R"][0]
                pre_f = _read_first_record(pre_fwd, "fasta") if pre_fwd.exists() else None
                pre_r = _read_first_record(pre_rev, "fasta") if pre_rev.exists() else None
                if pre_f is not None and pre_r is not None:
                    pre_candidates = list(iter_end_anchored_overlaps(str(pre_f.seq), str(pre_r.seq), None, None))
                    pre_best = pick_best_identity_candidate(pre_candidates, min_overlap=min_overlap)
                    pre_overlap = select_best_overlap(
                        pre_candidates,
                        min_overlap=min_overlap,
                        min_identity=min_identity,
                        min_quality=min_quality,
                        quality_mode=quality_mode,
                        ambiguity_identity_delta=ambiguity_identity_delta,
                        ambiguity_quality_epsilon=ambiguity_quality_epsilon,
                    )
                    pre_status = pre_overlap.status if pre_overlap.status in {"not_end_anchored", "ambiguous_overlap"} else _classify_overlap_status(
                        pre_overlap.overlap_len,
                        pre_overlap.identity,
                        pre_overlap.overlap_quality,
                        min_overlap=min_overlap,
                        min_identity=min_identity,
                        min_quality=min_quality,
                        quality_mode=quality_mode,
                    )
                    diag["pretrim_best_identity"] = f"{(pre_overlap.identity if pre_best is None else pre_best.identity):.4f}"
                    diag["pretrim_best_overlap_len"] = str(pre_overlap.overlap_len)
                    diag["pretrim_status"] = pre_status

            sample_primer = primer_trim_bases_by_sample.get(sample_id, {}) if primer_trim_bases_by_sample else {}
            diag["primer_trim_bases_fwd"] = str(sample_primer.get("F", 0))
            diag["primer_trim_bases_rev"] = str(sample_primer.get("R", 0))

            fh.write(
                f"{sample_id}\t{selected_overlap.overlap_len}\t{selected_overlap.identity:.4f}\t"
                f"{'' if overlap_quality is None else f'{overlap_quality:.2f}'}\t"
                f"{selected_overlap.orientation}\t{selected_status}\t{best_identity_val:.4f}\t{best_identity_orientation}"
                f"\t{'yes' if anchoring_feasible_selected else 'no'}\t{'yes' if end_anchored_possible_selected else 'no'}"
                f"\t{diag['fwd_best_identity']}\t{diag['revcomp_best_identity']}\t{diag['fwd_best_overlap_len']}\t{diag['revcomp_best_overlap_len']}"
                f"\t{diag['fwd_anchor_feasible']}\t{diag['revcomp_anchor_feasible']}\t{diag['identity_delta_revcomp_minus_fwd']}"
                f"\t{diag['selected_vs_best_identity_delta']}\t{diag['top_candidate_count']}\t{diag['top2_identity_delta']}"
                f"\t{diag['top2_overlap_len_delta']}\t{diag['top2_quality_delta']}\t{diag['tie_reason_code']}"
                f"\t{diag['fwd_best_identity_any']}\t{diag['revcomp_best_identity_any']}\t{diag['fwd_best_overlap_len_any']}\t{diag['revcomp_best_overlap_len_any']}"
                f"\t{diag['pretrim_best_identity']}\t{diag['posttrim_best_identity']}\t{diag['pretrim_best_overlap_len']}\t{diag['posttrim_selected_overlap_len']}"
                f"\t{diag['pretrim_status']}\t{diag['posttrim_status']}\t{diag['ambiguity_identity_delta_used']}\t{diag['ambiguity_quality_epsilon_used']}"
                f"\t{diag['primer_trim_bases_fwd']}\t{diag['primer_trim_bases_rev']}"
                f"\t{selected_engine}\t{fallback_used}\t{configured_overlap_engine}\n"
            )
            if efh:
                for engine_name, overlap, status, end_anchored_possible_e, anchoring_feasible_e, _cands in per_engine_rows:
                    efh.write(
                        f"{sample_id}\t{engine_name}\t{overlap.orientation}\t{status}\t{overlap.overlap_len}"
                        f"\t{overlap.identity:.4f}\t{'' if overlap.overlap_quality is None else f'{overlap.overlap_quality:.2f}'}"
                        f"\t{'yes' if end_anchored_possible_e else 'no'}\t{'yes' if anchoring_feasible_e else 'no'}"
                        f"\t{'yes' if engine_name == selected_engine else 'no'}\n"
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
                    "payload_kind": "none",
                    "payload_n": "0",
                    "payload_entity_n": "0",
                    "payload_max_len": "0",
                    "ambiguity_flag": "0",
                    "safety_flag": "none",
                    "decision_source": "auto",
                    "review_reason": "pair_missing",
                    "warning_flags": "",
                    "structural_hypothesis_n": "0",
                    "hypotheses_with_hits_n": "0",
                    "missing_hits_n": "0",
                    "hypothesis_map": "",
                    "source_id_map": "",
                }
            )
            continue

        contigs_path = asm_dir / sample_id / f"{sample_id}_paired.fasta.cap.contigs"
        singlets_path = asm_dir / sample_id / f"{sample_id}_paired.fasta.cap.singlets"

        contigs = list(SeqIO.parse(contigs_path, "fasta")) if contigs_path.exists() else []
        singlets = list(SeqIO.parse(singlets_path, "fasta")) if singlets_path.exists() else []

        # "contig or singlet" means CAP3 contigs are preferred and singlets are a fallback only when no contigs were produced for the sample. 
        source_group: tuple[str, str, list[SeqIO.SeqRecord]] | None = None 
        if contigs:
            source_group = ("contig", "cap3_c", contigs) 

        elif singlets:
            source_group = ("singlet", "cap3_s", singlets)

        payload_ids: list[str] = []
        hypothesis_map: list[str] = []
        source_id_map: list[str] = []

        payload_max_len_val = 0
        source_records_total = 0 

        if source_group:
            payload, label, source_records = source_group
            for idx, record in enumerate(source_records, start=1):
                new_id = f"{sample_id}|{payload}|{label}{idx}"
                source_id = str(record.id)
                structural_id = "struct1"
                payload_ids.append(f"{new_id}={source_id}")
                hypothesis_map.append(f"{new_id}={structural_id}")
                source_id_map.append(f"{new_id}={source_id}")
                record.id = new_id 
                record.name = new_id 
                record.description = ""
                records_to_write.append(record)
                source_records_total += 1
                payload_max_len_val = max(payload_max_len_val, len(record.seq))

        if source_group:
            payload = source_group[0]
            reason = "contig_present" if payload == "contig" else "singlets_only" 
        else:
            payload = "no_payload"
            reason = "cap3_no_output" 

        payload_kind = "none"
        payload_n = "0"
        payload_entity_n = "0"
        payload_max_len = "0"
        ambiguity_flag = "0"
        review_reason = ""
        warning_flags = ""
        structural_hypothesis_n = 0

        if payload in {"contig", "singlet"}:
            payload_kind = payload
            payload_n = str(source_records_total)
            payload_entity_n = payload_n
            payload_max_len = str(payload_max_len_val)  
            # CAP3 default emits payload entities; structural ambiguity is not implied by multi-contig/singlet.
            structural_hypothesis_n = 1
            if source_records_total > 1:
                warning_flags = "multi_payload"
        elif payload == "no_payload":
            review_reason = "no_payload"

        rows.append(
            {
                "sample_id": sample_id,
                "blast_payload": payload,
                "payload_ids": ";".join(payload_ids),
                "reason": reason,
                "payload_kind": payload_kind,
                "payload_n": payload_n,
                "payload_entity_n": payload_entity_n,
                "payload_max_len": payload_max_len,
                "ambiguity_flag": ambiguity_flag,
                "safety_flag": "none",
                "decision_source": "auto",
                "review_reason": review_reason,
                "warning_flags": warning_flags,
                "structural_hypothesis_n": str(structural_hypothesis_n),
                "hypotheses_with_hits_n": "0",
                "missing_hits_n": str(structural_hypothesis_n),
                "hypothesis_map": ";".join(hypothesis_map),
                "source_id_map": ";".join(source_id_map),
            }
        )

    with output_fasta.open("w", encoding="utf-8") as fh:
        SeqIO.write(records_to_write, fh, "fasta")

    headers = BLAST_INPUTS_HEADERS
    with output_tsv.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(headers) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(h, "") for h in headers) + "\n")


def _canonical_seq_key(name: str) -> str:
    key = str(name)
    for suffix in (".fastq.gz", ".fq.gz", ".fasta.gz", ".fa.gz", ".fna.gz", ".fas.gz", ".ab1.gz"):
        if key.lower().endswith(suffix):
            key = key[: -len(suffix)]
            break
    else:
        key = Path(key).stem
    if key.endswith("_trimmed"):
        key = key[: -len("_trimmed")]
    return key


def _extract_tax_label(value: object, rank: str = "species") -> str:
    text = str(value or "").strip()
    if not text:
        return ""

    rank_alias = {
        "kingdom": "k",
        "domain": "k",
        "phylum": "p",
        "class": "c",
        "order": "o",
        "family": "f",
        "genus": "g",
        "species": "s",
    }
    rank_key = rank_alias.get(str(rank).strip().lower(), "s")

    rank_map: dict[str, str] = {}
    for token in [tok.strip() for tok in text.split(";") if tok.strip()]:
        if "__" in token:
            prefix, val = token.split("__", 1)
            key = prefix.strip().lower()[:1]
            val_norm = val.strip()
            rank_map[key] = val_norm

    raw = rank_map.get(rank_key, "")
    if not raw:
        return ""

    norm = re.sub(r"\s+", " ", raw).strip()
    if norm.lower() in {"", "uncultured", "unclassified", "unknown", "na", "n/a"}:
        return ""
    return norm.lower()


def _trace_status_priority(status: str) -> int:
    order = {"na": 0, "pass": 1, "warn": 2, "fail": 3}
    return order.get(str(status or "").strip().lower(), 0)


def _worst_trace_status(statuses: Sequence[str]) -> str:
    cleaned = [str(s or "NA").strip().upper() for s in statuses if str(s or "").strip()]
    if not cleaned:
        return "NA"
    return max(cleaned, key=lambda s: _trace_status_priority(s))


def _load_trace_qc_by_sample(
    summary_tsv: Path | None,
    paired_samples: dict[str, dict[str, list[Path]]],
) -> dict[str, dict[str, str]]:
    if summary_tsv is None or not summary_tsv.exists():
        return {}

    try:
        df = pd.read_csv(summary_tsv, sep="\t", dtype=str, keep_default_na=False)
    except Exception:
        return {}

    required = {"file", "trace_qc_status", "trace_qc_flags"}
    if df.empty or not required.issubset(set(df.columns)):
        return {}

    trace_by_key: dict[str, dict[str, str]] = {}
    mixed_col = "trace_mixed_peak_frac" if "trace_mixed_peak_frac" in df.columns else None
    for _, row in df.iterrows():
        file_name = str(row.get("file", "") or "").strip()
        if not file_name or file_name.startswith("_"):
            continue
        key = _canonical_seq_key(file_name)
        status_raw = str(row.get("trace_qc_status", "") or "").strip().upper()
        flags_raw = str(row.get("trace_qc_flags", "") or "").strip()
        mixed_raw = str(row.get(mixed_col, "") if mixed_col else "").strip()
        trace_by_key[key] = {
            "trace_qc_status": status_raw if status_raw and status_raw != "NAN" else "NA",
            "trace_qc_flags": "" if flags_raw.upper() == "NAN" else flags_raw,
            "trace_mixed_peak_frac": "" if mixed_raw.upper() == "NAN" else mixed_raw,
        }

    out: dict[str, dict[str, str]] = {}
    for sample_id, orient_map in paired_samples.items():
        orient_status: dict[str, str] = {"F": "NA", "R": "NA"}
        flags: set[str] = set()
        mixed_values: list[float] = []
        for orient in ("F", "R"):
            files = orient_map.get(orient, [])
            file_statuses: list[str] = []
            for fp in files:
                key = _canonical_seq_key(fp.name)
                trace = trace_by_key.get(key)
                if not trace:
                    continue
                file_statuses.append(trace.get("trace_qc_status", "NA"))
                for flag in str(trace.get("trace_qc_flags", "")).split(";"):
                    flag = flag.strip()
                    if flag:
                        flags.add(flag)
                mixed_raw = str(trace.get("trace_mixed_peak_frac", "") or "").strip()
                if mixed_raw:
                    try:
                        mixed_values.append(float(mixed_raw))
                    except ValueError:
                        pass
            orient_status[orient] = _worst_trace_status(file_statuses)

        sample_status = _worst_trace_status([orient_status["F"], orient_status["R"]])
        out[sample_id] = {
            "trace_status": sample_status,
            "trace_status_f": orient_status["F"],
            "trace_status_r": orient_status["R"],
            "trace_flags": ";".join(sorted(flags)),
            "trace_mixed_peak_frac_max": f"{max(mixed_values):.4f}" if mixed_values else "",
        }
    return out


def _apply_trace_to_blast_payloads(
    blast_rows: dict[str, dict[str, str]],
    trace_rows: dict[str, dict[str, str]],
    *,
    trace_cfg: dict | None = None,
) -> dict[str, dict[str, str]]:
    if not blast_rows:
        return blast_rows

    cfg = trace_cfg if isinstance(trace_cfg, dict) else load_config().get("trace_qc", {})
    enable_mixture_inference = bool(cfg.get("enable_mixture_inference", False))
    try:
        mixture_thresh = float(cfg.get("mixture_suspect_threshold", cfg.get("mixed_frac_fail", 0.10)))
    except (TypeError, ValueError):
        mixture_thresh = 0.10

    out: dict[str, dict[str, str]] = {}
    for sample_id, row in blast_rows.items():
        merged = dict(row)
        trace = trace_rows.get(sample_id, {})
        trace_status = str(trace.get("trace_status", "NA") or "NA").upper()
        trace_flags = str(trace.get("trace_flags", "") or "")
        mixed_max_raw = str(trace.get("trace_mixed_peak_frac_max", "") or "")

        merged["trace_status"] = trace_status
        merged["trace_status_f"] = str(trace.get("trace_status_f", "NA") or "NA")
        merged["trace_status_r"] = str(trace.get("trace_status_r", "NA") or "NA")
        merged["trace_flags"] = trace_flags
        merged["trace_mixed_peak_frac_max"] = mixed_max_raw

        warning_set = {f for f in str(merged.get("warning_flags", "") or "").split(";") if f}
        safety_flag = str(merged.get("safety_flag", "none") or "none")
        review_reason = str(merged.get("review_reason", "") or "")

        if trace_status == "FAIL":
            merged["safety_flag"] = "trace_fail"
            if review_reason and review_reason != "trace_fail":
                warning_set.add(review_reason)
            merged["review_reason"] = "trace_fail"
        elif trace_status == "WARN":
            if safety_flag == "none":
                merged["safety_flag"] = "trace_warn"
            warning_set.add("trace_warn")

        if enable_mixture_inference and not merged.get("review_reason") and mixed_max_raw:
            try:
                if float(mixed_max_raw) >= mixture_thresh:
                    merged["review_reason"] = "mixture_suspected"
            except ValueError:
                pass

        merged["warning_flags"] = ";".join(sorted(warning_set))
        out[sample_id] = merged
    return out


def _apply_overlap_diagnostics_to_blast_payloads(
    blast_rows: dict[str, dict[str, str]],
    overlap_rows: dict[str, dict[str, str]] | None,
) -> dict[str, dict[str, str]]:
    if not blast_rows:
        return blast_rows
    overlap_rows = overlap_rows or {}
    out: dict[str, dict[str, str]] = {}

    for sample_id, row in blast_rows.items():
        merged = dict(row)
        warning_set = {f for f in str(merged.get("warning_flags", "") or "").split(";") if f}
        status = str((overlap_rows.get(sample_id) or {}).get("status", "") or "").strip().lower()

        if status in {"high_conflict"}:
            merged["safety_flag"] = "high_conflict"
            if not merged.get("review_reason"):
                merged["review_reason"] = "high_conflict"
        elif status in {"identity_low", "overlap_too_short", "quality_low", "ambiguous_overlap", "ambiguous_overlap_singlets", "merged_best_guess"}:
            warning_set.add(status)

        merged["warning_flags"] = ";".join(sorted(warning_set))
        out[sample_id] = merged

    return out


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



def run_compare_assemblers(
    input_dir: PathLike,
    out_dir: PathLike,
    *,
    assembler_ids: Sequence[str] | None = None, 
    fwd_pattern: str | None = None,
    rev_pattern: str | None = None,
    dup_policy: DupPolicy = DupPolicy.ERROR,
    enforce_same_well: bool = False,
    well_pattern: str | re.Pattern[str] | None = None,
    pairing_report: PathLike | None = None,
    use_qual: bool = True,
    threads: int = 1,
    on_stage=None,
    on_progress=None,
) -> Path:
    """Run all registered assemblers against staged paired FASTA inputs and emit a comparison table."""
    on_stage = on_stage or (lambda *_: None)
    on_progress = on_progress or (lambda *_: None)

    in_dir = Path(input_dir)
    if not in_dir.exists() or not in_dir.is_dir():
        raise ValueError(f"Assembler comparison requires a staged paired FASTA directory: {in_dir}")

    has_ab1_or_fastq = any(in_dir.rglob("*.ab1")) or any(in_dir.rglob("*.fastq")) or any(in_dir.rglob("*.fq"))
    has_fasta = any(in_dir.rglob("*.fasta")) or any(in_dir.rglob("*.fa")) or any(in_dir.rglob("*.fna")) or any(in_dir.rglob("*.fas"))
    if has_ab1_or_fastq and not has_fasta:
        raise ValueError(
            "Assembler comparison expects pre-staged paired FASTA inputs. "
            "Provide the paired FASTA staging directory (for example: qc/paired_fasta)."
        )

    out_root = Path(out_dir) / "asm" / "compare"
    out_root.mkdir(parents=True, exist_ok=True)

    paired_samples, missing_samples = _collect_pairing_catalog(
        in_dir,
        fwd_pattern=fwd_pattern,
        rev_pattern=rev_pattern,
        dup_policy=dup_policy,
        enforce_same_well=enforce_same_well,
        well_pattern=well_pattern,
    )
    if not paired_samples:
        raise ValueError("No paired samples available for assembler comparison")

    if pairing_report is not None:
        pairing_report_path = Path(pairing_report)
        pairing_report_path.parent.mkdir(parents=True, exist_ok=True)
        with pairing_report_path.open("w", encoding="utf-8") as fh:
            fh.write("sid\tF_path\tR_path\tdetector\n")
            for sid in sorted(set(paired_samples) | set(missing_samples)):
                entries = paired_samples.get(sid, {"F": [], "R": []})
                f_path = ";".join(str(p) for p in entries.get("F", []))
                r_path = ";".join(str(p) for p in entries.get("R", []))
                fh.write(f"{sid}\t{f_path}\t{r_path}\t\n")

    cfg = load_config()
    merge_cfg = cfg.get("merge_two_reads", {})
    cap3_exe = str(cfg.get("tools", {}).get("cap3", "cap3")).strip() or "cap3"
    overlap_min_overlap, overlap_min_identity, overlap_min_quality = _resolve_overlap_thresholds(cfg)
    quality_mode = str(cfg.get("overlap_eval", {}).get("quality_mode", "warning")).strip().lower()
    ambiguity_identity_delta = float(merge_cfg.get("ambiguity_identity_delta", 0.005))
    ambiguity_quality_epsilon = float(merge_cfg.get("ambiguity_quality_epsilon", 0.25))
    high_conflict_q_threshold = int(merge_cfg.get("high_conflict_q_threshold", 30))
    high_conflict_action = str(merge_cfg.get("high_conflict_action", "warn")).strip().lower()
    anchor_tolerance_bases = int(merge_cfg.get("anchor_tolerance_bases", 30))
    ambiguous_policy = str(merge_cfg.get("ambiguous_policy", "topk")).strip().lower()
    ambiguous_top_k = int(merge_cfg.get("ambiguous_top_k", 3))

    specs = list_assemblers()
    if assembler_ids:
        allowed = set(assembler_ids)
        specs = [spec for spec in specs if spec.id in allowed] 
        if not specs:
            known = ", ".join(sorted(spec.id for spec in list_assemblers()))
            asked = ", ".join(sorted(allowed)) 
            raise ValueError(f"No requested assemblers were found ({asked}). Known ids: {known}") 
    rows: list[dict[str, str]] = []
    sample_items = sorted(paired_samples.items())
    total = max(1, len(specs) * len(sample_items))
    done = 0
    heartbeat_interval = 5

    max_workers = max(1, int(threads or 1))

    spec_by_id = {spec.id: spec for spec in specs}

    def _run_compare_job(spec_id: str, sample_id: str, entries: dict[str, list[Path]]) -> tuple[str, str, dict[str, str]]:
        spec = spec_by_id[spec_id]
        spec_dir = out_root / spec.id.replace(":", "_")
        spec_dir.mkdir(parents=True, exist_ok=True)
        sample_dir = spec_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        status = "pair_missing"
        selected_engine = ""
        contig_len = ""
        warning = ""
        diag_code_for_machine = "pair_missing"
        diag_detail_for_human = "Missing foward/reverse pair for sample"
        payload_fasta = sample_dir / "payload.fasta"
        payload_written = False
        payload_kind = "none"
        payload_n = "0"
        payload_max_len = "0"
        decision_source = "auto"
        cap3_contigs_n = ""
        cap3_singlets_n = ""
        cap3_info_path = ""
        cap3_stdout_path = ""
        cap3_stderr_path = ""
        tool_name = ""
        tool_stdout_path = ""
        tool_stderr_path = ""

        if sample_id in missing_samples or not entries["F"] or not entries["R"]:
            status = "pair_missing"
        else:
            fwd_path = entries["F"][0]
            rev_path = entries["R"][0]
            try:
                if spec.kind == "merge_two_reads":
                    contig_path, report = merge_two_reads(
                        sample_id=sample_id,
                        fwd_path=fwd_path,
                        rev_path=rev_path,
                        output_dir=sample_dir,
                        min_overlap=overlap_min_overlap,
                        min_identity=overlap_min_identity,
                        min_quality=overlap_min_quality,
                        quality_mode=quality_mode,
                        ambiguity_identity_delta=ambiguity_identity_delta,
                        ambiguity_quality_epsilon=ambiguity_quality_epsilon,
                        high_conflict_q_threshold=high_conflict_q_threshold,
                        high_conflict_action=high_conflict_action,
                        overlap_engine=spec.overlap_engine or "biopython",
                        overlap_engine_strategy="single",
                        overlap_engine_order=None,
                        anchor_tolerance_bases=anchor_tolerance_bases,
                        ambiguous_policy=ambiguous_policy,
                        ambiguous_top_k=ambiguous_top_k,
                    )
                    status = str(getattr(report, "merge_status", "error") or "error")
                    report_engine = str(getattr(report, "overlap_engine", "") or "")
                    report_contig_len = getattr(report, "contig_len", 0)
                    report_warning = str(getattr(report, "merge_warning", "") or "")
                    report_overlap_len = getattr(report, "overlap_len", 0)
                    report_identity = float(getattr(report, "identity", 0.0) or 0.0)
                    report_orientation = str(getattr(report, "orientation", "") or "")

                    selected_engine = report_engine
                    contig_len = str(report_contig_len)
                    warning = report_warning
                    if status == "ambiguous_overlap" and not warning:
                        warning = "Overlap candidates were tied; no unique contig selected"
                    diag_code_for_machine = f"merge_{status}"
                    diag_detail_for_human = (
                        f"merge_status={status}; overlap_len={report_overlap_len};"
                        f"identity={report_identity:.4f}; engine={report_engine}"
                    )
                    tool_name = spec.display_name
                    merge_diag = [
                        f"sample_id={sample_id}",
                        f"assembler_id={spec.id}",
                        f"merge_status={status}",
                        f"overlap_engine={report_engine}",
                        f"orientation={report_orientation}",
                        f"overlap_len={report_overlap_len}",
                        f"identity={report_identity:.6f}",
                        f"warning={report_warning}",
                    ]
                    stdout_path, stderr_path = write_process_logs(
                        sample_dir,
                        sample_id,
                        spec.id,
                        ["merge_two_reads", f"--engine={spec.overlap_engine or 'biopython'}"],
                        "\n".join(merge_diag) + "\n",
                        "",
                    )
                    tool_stdout_path = str(stdout_path)
                    tool_stderr_path = str(stderr_path)
                    if contig_path and contig_path.exists():
                        payload_records = list(SeqIO.parse(contig_path, "fasta"))
                        if payload_records:
                            SeqIO.write(payload_records, payload_fasta, "fasta")
                            payload_written = True
                            payload_n = str(len(payload_records))
                            payload_kind = "contig_alt" if status == "ambiguous_topk" else "contig"
                            contig_len = str(max(len(rec.seq) for rec in payload_records))
                            payload_max_len = contig_len
                    else:
                        singlets_path = sample_dir / f"{sample_id}_paired.fasta.cap.singlets"
                        singlet_records = list(SeqIO.parse(singlets_path, "fasta")) if singlets_path.exists() else []
                        if singlet_records:
                            SeqIO.write(singlet_records, payload_fasta, "fasta")
                            payload_written = True
                            payload_kind = "singlet"
                            payload_n = str(len(singlet_records))
                            payload_max_len = str(max(len(rec.seq) for rec in singlet_records))
                else:
                    cap3_args = _resolve_cap3_options(spec.cap3_profile or "strict", None, None)
                    sample_fasta = sample_dir / f"{sample_id}_paired.fasta"
                    _write_combined_fasta([fwd_path, rev_path], sample_fasta, use_qual=use_qual)
                    import subprocess

                    cmd = [cap3_exe, sample_fasta.name, *cap3_args]
                    with tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as stdout_tmp, tempfile.TemporaryFile(mode="w+t", encoding="utf-8") as stderr_tmp:
                        cap3_run = subprocess.run(
                            cmd,
                            cwd=sample_dir,
                            stdout=stdout_tmp,
                            stderr=stderr_tmp,
                            text=True,
                            check=False,
                        )
                        stdout_tmp.seek(0)
                        stderr_tmp.seek(0)
                        stdout_text = stdout_tmp.read()
                        stderr_text = stderr_tmp.read()

                    stdout_path, stderr_path = write_process_logs(
                        sample_dir,
                        sample_id,
                        spec.id,
                        cmd,
                        stdout_text,
                        stderr_text,
                    )
                    cap3_stdout_path = str(stdout_path)
                    cap3_stderr_path = str(stderr_path)
                    tool_name = spec.display_name
                    tool_stdout_path = cap3_stdout_path
                    tool_stderr_path = cap3_stderr_path

                    cap_contigs = sample_dir / f"{sample_fasta.name}.cap.contigs"
                    cap_singlets = sample_dir / f"{sample_fasta.name}.cap.singlets"
                    cap_info = sample_dir / f"{sample_fasta.name}.cap.info"

                    payload_records = list(SeqIO.parse(cap_contigs, "fasta")) if cap_contigs.exists() else []
                    singlet_records = list(SeqIO.parse(cap_singlets, "fasta")) if cap_singlets.exists() else []
                    cap3_contigs_n = str(len(payload_records))
                    cap3_singlets_n = str(len(singlet_records))
                    cap3_info_path = str(cap_info) if cap_info.exists() else ""
                    selected_engine = "cap3"

                    if cap3_run.returncode == 0 and payload_records:
                        status = "assembled"
                        diag_code_for_machine = "cap3_contigs_present"
                    elif cap3_run.returncode == 0 and singlet_records:
                        status = "singlets_only"
                        diag_code_for_machine = "cap3_singlets_only"
                    elif cap3_run.returncode == 0:
                        status = "cap3_no_output"
                        diag_code_for_machine = "cap3_no_contigs_no_singlets"
                    elif payload_records:
                        status = "assembled"
                        diag_code_for_machine = "cap3_nonzero_exit_with_output"
                    elif singlet_records:
                        status = "singlets_only"
                        diag_code_for_machine = "cap3_nonzero_exit_with_output"
                    else:
                        status = "error"
                        diag_code_for_machine = "cap3_nonzero_exit_no_output"
                        warning = warning or (stderr_text or stdout_text or f"CAP3 exited {cap3_run.returncode}")

                    diag_detail_for_human = (
                        f"returncode={cap3_run.returncode}; profile={spec.cap3_profile or 'strict'}; "
                        f"input={sample_fasta.name}; contigs={len(payload_records)}; singlets={len(singlet_records)}"
                    )

                    payload_source = payload_records if payload_records else singlet_records
                    if payload_source:
                        SeqIO.write(payload_source, payload_fasta, "fasta")
                        payload_written = True
                        payload_n = str(len(payload_source))
                        payload_kind = "contig" if payload_records else "singlet"
                        payload_max_len = str(max(len(rec.seq) for rec in payload_source))
                        if payload_records:
                            contig_len = str(max(len(rec.seq) for rec in payload_records))
                            payload_max_len = contig_len
            except Exception as exc:
                status = "error"
                warning = str(exc)
                diag_code_for_machine = "exception"
                diag_detail_for_human = str(exc)

        ambiguity_flag, safety_flag, review_reason = _infer_compare_row_flags(status, warning)
        if not payload_written:
            payload_kind = "none"
            payload_n = "0"
            payload_max_len = "0"
        row = {
            "sample_id": sample_id,
            "assembler_id": spec.id,
            "assembler_name": spec.display_name,
            "dup_policy": str(dup_policy.value if isinstance(dup_policy, DupPolicy) else dup_policy),
            "status": status,
            "selected_engine": selected_engine,
            "contig_len": contig_len,
            "warnings": warning,
            "diag_code_for_machine": diag_code_for_machine,
            "diag_detail_for_human": diag_detail_for_human,
            "cap3_contigs_n": cap3_contigs_n,
            "cap3_singlets_n": cap3_singlets_n,
            "cap3_info_path": cap3_info_path,
            "cap3_stdout_path": cap3_stdout_path,
            "cap3_stderr_path": cap3_stderr_path,
            "tool_name": tool_name,
            "tool_stdout_path": tool_stdout_path,
            "tool_stderr_path": tool_stderr_path,
            "selection_trace_path": "",
            "winner_reason": "",
            "payload_fasta": str(payload_fasta) if payload_written else "",
            "payload_kind": payload_kind,
            "payload_n": payload_n,
            "payload_max_len": payload_max_len,
            "ambiguity_flag": ambiguity_flag,
            "safety_flag": safety_flag,
            "decision_source": decision_source,
            "review_reason": review_reason,
            "warning_flags": "",
            "structural_hypothesis_n": "0",
            "hypotheses_with_hits_n": "0",
            "missing_hits_n": "0",
            "hypothesis_map": "",
            "source_id_map": "",
        }

        return spec.id, sample_id, row

    job_items: list[tuple[str, str, dict[str, list[Path]]]] = [
        (spec.id, sample_id, entries)
        for spec in specs
        for sample_id, entries in sample_items
    ]

    completed_by_spec: dict[str, int] = defaultdict(int)
    on_stage("Compare assemblers")

    if max_workers == 1 or len(job_items) <= 1:
        for spec_id, sample_id, entries in job_items:
            done_spec_id, _sample_id, row = _run_compare_job(spec_id, sample_id, entries)
            rows.append(row)
            done += 1
            completed_by_spec[done_spec_id] += 1
            if completed_by_spec[done_spec_id] % heartbeat_interval == 0 or completed_by_spec[done_spec_id] == len(sample_items):
                spec_name = spec_by_id[done_spec_id].display_name
                heartbeat = f"Compare {spec_name}: {completed_by_spec[done_spec_id]}/{len(sample_items)} complete"
                L.info(heartbeat)
                on_stage(heartbeat)
            on_progress(int(done * 100 / total))
    else:
        future_map = {}
        with ThreadPoolExecutor(max_workers=min(max_workers, len(job_items))) as ex:
            for spec_id, sample_id, entries in job_items:
                future = ex.submit(_run_compare_job, spec_id, sample_id, entries)
                future_map[future] = (spec_id, sample_id)
            collected: dict[tuple[int, int], dict[str, str]] = {}
            spec_order = {spec.id: idx for idx, spec in enumerate(specs)}
            sample_order = {sample_id: idx for idx, (sample_id, _entries) in enumerate(sample_items)}
            for future in as_completed(future_map):
                spec_id, sample_id = future_map[future]
                done_spec_id, done_sample_id, row = future.result()
                collected[(spec_order[spec_id], sample_order[sample_id])] = row
                done += 1
                completed_by_spec[done_spec_id] += 1
                if completed_by_spec[done_spec_id] % heartbeat_interval == 0 or completed_by_spec[done_spec_id] == len(sample_items):
                    spec_name = spec_by_id[done_spec_id].display_name
                    heartbeat = f"Compare {spec_name}: {completed_by_spec[done_spec_id]}/{len(sample_items)} complete"
                    L.info(heartbeat)
                    on_stage(heartbeat)
                on_progress(int(done * 100 / total))
        for key in sorted(collected):
            rows.append(collected[key])

    out_tsv = Path(out_dir) / "asm" / "compare_assemblers.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    headers = [
        "sample_id", "assembler_id", "assembler_name", "dup_policy", "status", "selected_engine", "contig_len", "warnings",
        "diag_code_for_machine", "diag_detail_for_human", "cap3_contigs_n", "cap3_singlets_n", "cap3_info_path", "cap3_stdout_path", "cap3_stderr_path",
        "tool_name", "tool_stdout_path", "tool_stderr_path", "selection_trace_path", "winner_reason",
        "payload_fasta", "payload_kind", "payload_n", "payload_max_len", "ambiguity_flag", "safety_flag", "decision_source", "review_reason", "warning_flags",
        "structural_hypothesis_n", "hypotheses_with_hits_n", "missing_hits_n", "hypothesis_map", "source_id_map",
        "resolution_state", "resolved_hypothesis", "resolution_reason", "trace_status", "trace_status_f", "trace_status_r", "trace_flags",
    ]
    with out_tsv.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(headers) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(h, "") for h in headers) + "\n")
    return out_tsv

def _safe_trace_label(value: str) -> str:
    """Return a filesystem-safe identifier for selection-trace artifact names."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", (value or "").strip() or "sample")


def _write_selection_trace_tsv(
    trace_dir: Path,
    sample_id: str,
    mode: str,
    ranked_rows: list[dict[str, str]],
    winner_reason: str,
) -> Path:
    """Persist deterministic winner-selection diagnostics for one sample."""
    trace_dir.mkdir(parents=True, exist_ok=True)
    trace_path = trace_dir / f"{_safe_trace_label(sample_id)}.selection_trace.tsv"
    headers = [
        "sample_id",
        "mode",
        "rank_position",
        "is_winner",
        "winner_reason",
        "assembler_id",
        "assembler_name",
        "status",
        "success_rank",
        "payload_kind",
        "payload_rank",
        "ambiguity_flag",
        "amb_penalty",
        "normalized_len",
        "tie_break_key",
        "payload_fasta",
    ]
    with trace_path.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(headers) + "\n")
        for row in ranked_rows:
            fh.write("\t".join(str(row.get(h, "")) for h in headers) + "\n")
    return trace_path


def _select_best_compare_rows(
    compare_tsv: Path,
    *,
    assembler_mode: str,
    assembler_id: str | None = None,
) -> dict[str, dict[str, str]]:
    if not compare_tsv.exists():
        return {}
    df = pd.read_csv(compare_tsv, sep="\t", dtype=str).fillna("")
    if df.empty:
        return {}

    winners: dict[str, dict[str, str]] = {}
    mode = (assembler_mode or "").strip().lower()
    payload_rank = {"contig": 3, "contig_alt": 2, "singlet": 1, "none": 0}

    if "selection_trace_path" not in df.columns:
        df["selection_trace_path"] = ""
    if "winner_reason" not in df.columns:
        df["winner_reason"] = ""

    trace_root = compare_tsv.parent / "selection_trace"

    for sample_id, group in df.groupby("sample_id"):
        g = group.copy()
        g["_assembler_id"] = g.get("assembler_id", "").fillna("").astype(str)
        g["_assembler_name"] = g.get("assembler_name", "").fillna("").astype(str)
        g["_status"] = g.get("status", "").fillna("").astype(str)
        g["_payload_fasta"] = g.get("payload_fasta", "").fillna("").astype(str)
        g["_success_rank"] = g["_status"].map(lambda v: 2 if v in {"assembled", "merged"} else 1)

        if mode == "selected":
            if not assembler_id:
                raise ValueError("selected assembler mode requires assembler_id")
            g = g[g["_assembler_id"] == assembler_id].copy()
            if g.empty:
                continue
            g["_kind"] = g.get("payload_kind", "none").fillna("none").astype(str)
            g["_payload_rank"] = g["_kind"].map(lambda v: payload_rank.get(v, 0))
            g["_amb_flag"] = g.get("ambiguity_flag", "0").fillna("0").astype(str)
            g["_amb_penalty"] = g["_amb_flag"].map(lambda v: 1 if v == "1" else 0)
            g["_len"] = pd.to_numeric(g.get("contig_len", ""), errors="coerce").fillna(-1)
            g = g.sort_values(by=["_len", "_assembler_id"], ascending=[False, True], kind="mergesort")
            reason = (
                f"selected mode: restricted to assembler_id={assembler_id}; "
                "winner is longest contig_len, then assembler_id tie-break"
            )
        else:
            kind_raw = g["payload_kind"] if "payload_kind" in g.columns else "none"
            amb_raw = g["ambiguity_flag"] if "ambiguity_flag" in g.columns else "0"
            len_raw = g["payload_max_len"] if "payload_max_len" in g.columns else g.get("contig_len", "")

            g["_kind"] = pd.Series(kind_raw, index=g.index).fillna("none").astype(str)
            g["_amb_flag"] = pd.Series(amb_raw, index=g.index).fillna("0").astype(str)
            g["_len"] = pd.to_numeric(pd.Series(len_raw, index=g.index), errors="coerce").fillna(-1)
            g["_payload_rank"] = g["_kind"].map(lambda v: payload_rank.get(v, 0))
            g["_amb_penalty"] = g["_amb_flag"].map(lambda v: 1 if v == "1" else 0)
            g = g.sort_values(
                by=["_success_rank", "_payload_rank", "_amb_penalty", "_len", "_assembler_id"],
                ascending=[False, False, True, False, True],
                kind="mergesort",
            )
            reason = (
                "all mode: ranked by success tier (assembled and merged are equivalent), "
                "then payload_kind(contig>contig_alt>singlet>none), then ambiguity penalty (0 before 1), "
                "then normalized length desc, then assembler_id"
            )

        winner_idx = g.index[0]
        row = g.loc[winner_idx].to_dict()
        winner = {str(k): "" if pd.isna(v) else str(v) for k, v in row.items() if not str(k).startswith("_")}

        ranked_rows: list[dict[str, str]] = []
        for rank_position, (_, r) in enumerate(g.iterrows(), start=1):
            ranked_rows.append({
                "sample_id": str(sample_id),
                "mode": mode,
                "rank_position": str(rank_position),
                "is_winner": "1" if rank_position == 1 else "0",
                "winner_reason": reason if rank_position == 1 else "",
                "assembler_id": str(r.get("_assembler_id", "")),
                "assembler_name": str(r.get("_assembler_name", "")),
                "status": str(r.get("_status", "")),
                "success_rank": str(int(r.get("_success_rank", 1))),
                "payload_kind": str(r.get("_kind", "none")),
                "payload_rank": str(int(r.get("_payload_rank", 0))),
                "ambiguity_flag": str(r.get("_amb_flag", "0")),
                "amb_penalty": str(int(r.get("_amb_penalty", 0))),
                "normalized_len": str(int(r.get("_len", -1))),
                "tie_break_key": str(r.get("_assembler_id", "")),
                "payload_fasta": str(r.get("_payload_fasta", "")),
            })

        trace_path = _write_selection_trace_tsv(trace_root, str(sample_id), mode, ranked_rows, reason)
        df.loc[group.index, "selection_trace_path"] = str(trace_path)
        df.loc[group.index, "winner_reason"] = ""
        df.loc[winner_idx, "winner_reason"] = reason

        winner["selection_trace_path"] = str(trace_path)
        winner["winner_reason"] = reason
        winners[str(sample_id)] = winner

    base_cols = [
        "sample_id", "assembler_id", "assembler_name", "dup_policy", "status", "selected_engine", "contig_len", "warnings",
        "diag_code_for_machine", "diag_detail_for_human", "cap3_contigs_n", "cap3_singlets_n", "cap3_info_path", "cap3_stdout_path", "cap3_stderr_path",
        "tool_name", "tool_stdout_path", "tool_stderr_path", "selection_trace_path", "winner_reason",
        "payload_fasta", "payload_kind", "payload_n", "payload_max_len", "ambiguity_flag", "safety_flag", "decision_source", "review_reason", "warning_flags",
        "structural_hypothesis_n", "hypotheses_with_hits_n", "missing_hits_n", "hypothesis_map", "source_id_map",
        "resolution_state", "resolved_hypothesis", "resolution_reason", "trace_status", "trace_status_f", "trace_status_r", "trace_flags",
    ]
    ordered_cols = [c for c in base_cols if c in df.columns] + [c for c in df.columns if c not in base_cols]
    df = df[ordered_cols]
    df.to_csv(compare_tsv, sep="\t", index=False)
    return winners

def _build_selected_blast_inputs(
    selected_rows: dict[str, dict[str, str]],
    paired_samples: dict[str, dict[str, list[Path]]],
    missing_samples: set[str],
    output_fasta: Path,
    output_tsv: Path,
    *,
    no_payload_reason: str,
) -> dict[str, dict[str, str]]:
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    records: list[SeqRecord] = []
    rows: list[dict[str, str]] = []

    all_samples = sorted(set(paired_samples) | set(missing_samples))
    for sample_id in all_samples:
        chosen = selected_rows.get(sample_id, {})
        asm_id = chosen.get("assembler_id", "")
        payload_path = Path(chosen["payload_fasta"]) if chosen.get("payload_fasta") else None
        payload_records = list(SeqIO.parse(payload_path, "fasta")) if payload_path and payload_path.exists() else []

        payload = "pair_missing" if sample_id in missing_samples else "no_payload"
        reason = "pair_missing" if sample_id in missing_samples else no_payload_reason
        payload_kind = str(chosen.get("payload_kind", "none") or "none")
        payload_n = str(chosen.get("payload_n", "0") or "0")
        payload_max_len = str(chosen.get("payload_max_len", "0") or "0")
        ambiguity_flag = str(chosen.get("ambiguity_flag", "0") or "0")
        safety_flag = str(chosen.get("safety_flag", "none") or "none")
        decision_source = str(chosen.get("decision_source", "auto") or "auto")
        review_reason = str(chosen.get("review_reason", "") or "")
        warning_flags = str(chosen.get("warning_flags", "") or "")

        selected_status = str(chosen.get("status", "") or "").strip().lower()
        if payload_records and payload_kind == "none":
            if selected_status in {"ambiguous_topk"}:
                payload_kind = "contig_alt"
                ambiguity_flag = "1"
                if not review_reason:
                    review_reason = "ambiguous_payload"
            elif selected_status == "singlets_only":
                payload_kind = "singlet"
            else:
                payload_kind = "contig"

        claimed_payload = bool(chosen.get("payload_fasta")) or payload_kind != "none"

        payload_ids: list[str] = []
        hypothesis_map: list[str] = []
        source_id_map: list[str] = []
        if payload_records and payload_kind != "none":
            payload = "singlet" if payload_kind == "singlet" else "contig"
            reason = "selected_payload"
            label = asm_id.replace(":", "_") if asm_id else "selected"
            for idx, rec in enumerate(payload_records, start=1):
                source_id = str(rec.id)
                if payload_kind == "contig_alt":
                    suffix = f"hyp{idx}"
                    structural_id = suffix
                elif payload_kind == "singlet":
                    suffix = f"singlet{idx}"
                    structural_id = "struct1"
                else:
                    suffix = f"contig{idx}"
                    structural_id = "struct1"
                new_id = f"{sample_id}|{label}|{suffix}"
                payload_ids.append(f"{new_id}={source_id}")
                hypothesis_map.append(f"{new_id}={structural_id}")
                source_id_map.append(f"{new_id}={source_id}")
                rec.id = new_id
                rec.name = new_id
                rec.description = ""
                records.append(rec)
            payload_n = str(len(payload_records))
            payload_max_len = str(max(len(rec.seq) for rec in payload_records))
        elif sample_id in missing_samples:
            payload_kind = "none"
            payload_n = "0"
            payload_max_len = "0"
        elif claimed_payload:
            payload_kind = "none"
            payload_n = "0"
            payload_max_len = "0"
            ambiguity_flag = "0"
            review_reason = ""
            reason = "payload_missing_or_empty"
        else:
            payload_kind = "none"
            payload_n = "0"
            payload_max_len = "0"

        structural_hypothesis_n = 0
        if payload_records and payload_kind != "none":
            structural_hypothesis_n = len(payload_records) if payload_kind == "contig_alt" else 1
            if len(payload_records) > 1 and payload_kind != "contig_alt":
                warning_flags = ";".join(sorted(set([f for f in warning_flags.split(";") if f] + ["multi_payload"])))

        rows.append({
            "sample_id": sample_id,
            "blast_payload": payload,
            "payload_ids": ";".join(payload_ids),
            "reason": reason,
            "selected_assembler_id": asm_id,
            "selected_assembler_name": chosen.get("assembler_name", ""),
            "selected_status": chosen.get("status", ""),
            "payload_kind": payload_kind,
            "payload_n": payload_n,
            "payload_entity_n": payload_n,
            "payload_max_len": payload_max_len,
            "ambiguity_flag": ambiguity_flag,
            "safety_flag": safety_flag,
            "decision_source": decision_source,
            "review_reason": review_reason,
            "warning_flags": warning_flags,
            "structural_hypothesis_n": str(structural_hypothesis_n),
            "hypotheses_with_hits_n": "0",
            "missing_hits_n": str(structural_hypothesis_n),
            "hypothesis_map": ";".join(hypothesis_map),
            "source_id_map": ";".join(source_id_map),
        })

    if records:
        SeqIO.write(records, output_fasta, "fasta")
    else:
        output_fasta.write_text("", encoding="utf-8")

    headers = BLAST_INPUTS_HEADERS
    with output_tsv.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(headers) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(h, "") for h in headers) + "\n")

    return {row["sample_id"]: row for row in rows}

def _pick_advisory_reason(safety_flag: str, warning_flags: str) -> str:
    warnings = [f for f in str(warning_flags or "").split(";") if f]
    if safety_flag == "trace_warn":
        return "trace_warn"
    priority = ["trace_warn", "multi_payload", "identity_low", "overlap_too_short", "quality_low", "ambiguous_overlap", "ambiguous_overlap_singlets", "merged_best_guess"]
    present = set(warnings)
    for reason in priority:
        if reason in present:
            return reason
    return warnings[0] if warnings else ""


def _blocking_reason_priority(reason: str) -> int:
    order = {
        "trace_fail": 100,
        "high_conflict": 95,
        "hypothesis_map_invalid": 90,
        "mixture_suspected": 80,
        "pair_missing": 70,
        "no_payload": 65,
        "taxonomy_missing": 60,
        "no_hits": 55,
        "partial_hits": 50,
        "rank_missing": 45,
        "ambiguous_taxonomy": 40,
        "ambiguous_payload": 30,
    }
    return order.get(str(reason or "").strip(), 10)


def _resolve_sample_resolution(
    sample_id: str,
    per_hyp: pd.DataFrame,
    *,
    label_col: str,
    blast_meta: dict[str, str],
    tax_rank: str = "species",
    tax_table_ok: bool = True,
) -> dict[str, str]:
    per_hyp = per_hyp.copy() if per_hyp is not None else pd.DataFrame()

    qseq_to_structural: dict[str, str] = {}
    for pair in str(blast_meta.get("hypothesis_map", "") or "").split(";"):
        if not pair.strip():
            continue
        if "=" not in pair:
            raise ValueError(f"Invalid hypothesis_map entry (missing '='): {pair!r}")
        lhs, rhs = pair.split("=", 1)
        qseqid = lhs.strip()
        structural_id = rhs.strip()
        if not structural_id:
            raise ValueError(f"Invalid hypothesis_map entry (empty structural id): {pair!r}")
        if qseqid in qseq_to_structural and qseq_to_structural[qseqid] != structural_id:
            raise ValueError(f"Non-deterministic hypothesis_map for {qseqid}: {qseq_to_structural[qseqid]} vs {structural_id}")
        qseq_to_structural[qseqid] = structural_id

    hyp_rows: list[dict[str, str | float]] = []
    for _, hyp in per_hyp.iterrows():
        qseqid = str(hyp.get("qseqid", "") or "")
        structural = qseq_to_structural.get(qseqid, qseqid.split("|")[-1] if qseqid else "")
        label = _extract_tax_label(hyp.get(label_col, ""), rank=tax_rank)
        bitscore = float(hyp.get("bitscore")) if pd.notna(hyp.get("bitscore")) else float("-inf")
        pident = float(hyp.get("pident")) if pd.notna(hyp.get("pident")) else float("-inf")
        qcov = float(hyp.get("qcovhsp")) if pd.notna(hyp.get("qcovhsp")) else float("-inf")
        evalue = float(hyp.get("evalue")) if pd.notna(hyp.get("evalue")) else float("inf")
        hyp_rows.append(
            {
                "qseqid": qseqid,
                "structural": structural,
                "label": label,
                "bitscore": bitscore,
                "pident": pident,
                "qcovhsp": qcov,
                "evalue": evalue,
            }
        )

    hyp_rows = sorted(
        hyp_rows,
        key=lambda r: (-float(r["bitscore"]), -float(r["pident"]), -float(r["qcovhsp"]), float(r["evalue"]), str(r["qseqid"])),
    )

    labels = [str(h["label"]) for h in hyp_rows if str(h["label"])]
    unique_labels = sorted(set(labels))
    hypotheses_with_hits_n = len({str(h["structural"]) for h in hyp_rows if str(h["structural"])})

    try:
        structural_hypothesis_n = int(str(blast_meta.get("structural_hypothesis_n", "0") or "0"))
    except ValueError:
        structural_hypothesis_n = 0
    if structural_hypothesis_n <= 0:
        structural_hypothesis_n = hypotheses_with_hits_n

    missing_hits_n = max(0, structural_hypothesis_n - hypotheses_with_hits_n)

    default_reason = str(blast_meta.get("review_reason", "") or "")
    safety_flag = str(blast_meta.get("safety_flag", "none") or "none")
    trace_status = str(blast_meta.get("trace_status", "NA") or "NA")
    trace_flags = str(blast_meta.get("trace_flags", "") or "")
    warning_flags = str(blast_meta.get("warning_flags", "") or "")
    blast_payload = str(blast_meta.get("blast_payload", "") or "")

    resolved_hypothesis = str(hyp_rows[0]["structural"]) if hyp_rows else ""

    resolution_state = "needs_review"
    resolution_reason = "taxonomy_missing"
    if blast_payload == "pair_missing":
        resolution_reason = "pair_missing"
    elif blast_payload == "no_payload":
        resolution_reason = "no_payload"
    elif structural_hypothesis_n <= 0:
        resolution_reason = "no_structural_hypotheses"
    elif hypotheses_with_hits_n == 0:
        resolution_reason = "no_hits" if tax_table_ok else "taxonomy_missing"
    elif not labels and len(hyp_rows) > 0:
        resolution_reason = "rank_missing"
    elif missing_hits_n > 0:
        resolution_reason = "partial_hits"
    elif structural_hypothesis_n == 1:
        resolution_state = "unambiguous"
        resolution_reason = "single_hypothesis"
    elif len(unique_labels) == 1:
        resolution_state = "resolved_by_evidence"
        resolution_reason = f"hypotheses_agree_{tax_rank}"
    else:
        resolution_reason = "ambiguous_taxonomy"

    if len(hyp_rows) > 0 and not labels and resolution_reason not in {"pair_missing", "no_payload", "no_hits"}:
        resolution_state = "needs_review"
        resolution_reason = "rank_missing"

    if safety_flag in {"trace_fail", "high_conflict"}:
        resolution_state = "needs_review"
        resolution_reason = safety_flag

    if default_reason == "mixture_suspected":
        resolution_state = "needs_review"
        resolution_reason = "mixture_suspected"

    warning_set = {f for f in str(warning_flags or "").split(";") if f}

    review_action = "none"
    review_reason = ""
    if resolution_state == "needs_review":
        review_action = "queue"
        if default_reason:
            if _blocking_reason_priority(default_reason) >= _blocking_reason_priority(resolution_reason):
                review_reason = default_reason
                if default_reason != resolution_reason:
                    warning_set.add(resolution_reason)
            else:
                review_reason = resolution_reason
                warning_set.add(default_reason)
        else:
            review_reason = resolution_reason

    warning_flags = ";".join(sorted(warning_set))

    row = ResolutionContractRow(
        sample_id=str(sample_id),
        review_action=review_action,
        review_reason=review_reason,
        advisory_reason=_pick_advisory_reason(safety_flag, warning_flags),
        warning_flags=warning_flags,
        structural_hypothesis_n=structural_hypothesis_n,
        hypotheses_with_hits_n=hypotheses_with_hits_n,
        missing_hits_n=missing_hits_n,
        top_labels="|".join(unique_labels[:5]),
        resolution_state=resolution_state,
        resolved_hypothesis=resolved_hypothesis,
        resolution_reason=resolution_reason,
        trace_status=trace_status,
        trace_flags=trace_flags,
    ).to_dict()
    return row


def _annotate_resolution_from_tax(
    tax_tsv: Path,
    blast_rows: dict[str, dict[str, str]],
    *,
    tax_rank: str = "species",
) -> dict[str, dict[str, str]]:
    updated = {k: dict(v) for k, v in blast_rows.items()}

    tax_by_sample: dict[str, pd.DataFrame] = {}
    label_col: str | None = None
    tax_table_ok = False
    if tax_tsv.exists():
        try:
            df = pd.read_csv(tax_tsv, sep="\t")
        except Exception:
            df = pd.DataFrame()
        if "qseqid" in df.columns:
            label_col = "taxonomy" if "taxonomy" in df.columns else (
                "stitle" if "stitle" in df.columns else ("sseqid" if "sseqid" in df.columns else None)
            )
            tax_table_ok = label_col is not None
            if label_col is not None and not df.empty:
                temp = df.copy()
                temp["sample_id"] = temp["qseqid"].astype(str).str.split("|").str[0]
                for c in ("bitscore", "pident", "qcovhsp", "evalue"):
                    if c in temp.columns:
                        temp[c] = pd.to_numeric(temp[c], errors="coerce")
                sort_cols: list[str] = []
                ascending: list[bool] = []
                if "bitscore" in temp.columns:
                    sort_cols.append("bitscore")
                    ascending.append(False)
                if "pident" in temp.columns:
                    sort_cols.append("pident")
                    ascending.append(False)
                if "qcovhsp" in temp.columns:
                    sort_cols.append("qcovhsp")
                    ascending.append(False)
                if "evalue" in temp.columns:
                    sort_cols.append("evalue")
                    ascending.append(True)
                for sample_id, group in temp.groupby("sample_id"):
                    ranked = group.sort_values(by=sort_cols, ascending=ascending) if sort_cols else group
                    tax_by_sample[str(sample_id)] = ranked.groupby("qseqid", as_index=False).first()

    for sample_id in sorted(set(updated) | set(tax_by_sample)):
        per_hyp = tax_by_sample.get(sample_id, pd.DataFrame(columns=["qseqid", label_col or "taxonomy"]))
        try:
            resolved = _resolve_sample_resolution(
                sample_id,
                per_hyp,
                label_col=label_col or "taxonomy",
                blast_meta=updated.get(sample_id, {}),
                tax_rank=tax_rank,
                tax_table_ok=tax_table_ok,
            )
        except ValueError as exc:
            meta = dict(updated.get(sample_id, {}))
            warnings = {f for f in str(meta.get("warning_flags", "") or "").split(";") if f}
            warnings.add("hypothesis_map_invalid")
            fallback = ResolutionContractRow(
                sample_id=str(sample_id),
                review_action="queue",
                review_reason="hypothesis_map_invalid",
                advisory_reason="hypothesis_map_invalid",
                warning_flags=";".join(sorted(warnings)),
                structural_hypothesis_n=int(str(meta.get("structural_hypothesis_n", "0") or "0") or 0),
                hypotheses_with_hits_n=0,
                missing_hits_n=max(0, int(str(meta.get("structural_hypothesis_n", "0") or "0") or 0)),
                top_labels="",
                resolution_state="needs_review",
                resolved_hypothesis="",
                resolution_reason="hypothesis_map_invalid",
                trace_status=str(meta.get("trace_status", "NA") or "NA"),
                trace_flags=str(meta.get("trace_flags", "") or ""),
            ).to_dict()
            L.warning("Invalid hypothesis_map for sample %s: %s", sample_id, exc)
            resolved = fallback
        tgt = updated.setdefault(sample_id, {})
        tgt.update(resolved)
    return updated


def _write_review_queue_from_resolved(
    output_tsv: Path,
    resolved_rows: dict[str, dict[str, str]],
) -> None:
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    headers = ResolutionContractRow.header()
    with output_tsv.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(headers) + "\n")
        for sample_id in sorted(resolved_rows):
            row = ResolutionContractRow(sample_id=sample_id)
            data = dict(row.to_dict())
            data.update({k: str(v) for k, v in resolved_rows[sample_id].items()})
            # normalize again with provided values
            normalized = ResolutionContractRow(
                sample_id=data.get("sample_id", sample_id),
                review_action=data.get("review_action", "none"),
                review_reason=data.get("review_reason", ""),
                advisory_reason=data.get("advisory_reason", ""),
                warning_flags=data.get("warning_flags", ""),
                structural_hypothesis_n=int(data.get("structural_hypothesis_n", "0") or 0),
                hypotheses_with_hits_n=int(data.get("hypotheses_with_hits_n", "0") or 0),
                missing_hits_n=int(data.get("missing_hits_n", "0") or 0),
                top_labels=data.get("top_labels", ""),
                resolution_state=data.get("resolution_state", "needs_review"),
                resolved_hypothesis=data.get("resolved_hypothesis", ""),
                resolution_reason=data.get("resolution_reason", "taxonomy_missing"),
                trace_status=data.get("trace_status", "NA"),
                trace_flags=data.get("trace_flags", ""),
            ).to_dict()
            fh.write("\t".join(normalized.get(h, "") for h in headers) + "\n")


def _write_review_queue(
    tax_tsv: Path,
    output_tsv: Path,
    *,
    blast_rows: dict[str, dict[str, str]],
    tax_rank: str = "species",
) -> None:
    resolved_rows = _annotate_resolution_from_tax(tax_tsv, blast_rows, tax_rank=tax_rank)
    _write_review_queue_from_resolved(output_tsv, resolved_rows)

def _write_selected_assembly_summary(
    output_tsv: Path,
    paired_samples: dict[str, dict[str, list[Path]]],
    missing_samples: set[str],
    selected_rows: dict[str, dict[str, str]],
    blast_rows: dict[str, dict[str, str]],
    assembler_mode_label: str,
    *,
    primer_mode: str,
    primer_stage: str,
    primer_preset: str,
    primer_source: str,
    overlap_engine_strategy: str,
    overlap_engine_order: str,
    overlap_quality_mode: str,
    overlap_rows: dict[str, dict[str, str]] | None = None,
) -> None:
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    headers = [
        "sample_id", "status", "assembler", "contig_len", "blast_payload",
        "payload_kind", "payload_n", "payload_entity_n", "payload_max_len", "ambiguity_flag", "safety_flag", "decision_source", "review_reason", "warning_flags",
        "structural_hypothesis_n", "hypotheses_with_hits_n", "missing_hits_n", "hypothesis_map", "source_id_map",
        "resolution_state", "resolved_hypothesis", "resolution_reason", "trace_status", "trace_status_f", "trace_status_r", "trace_flags",
        "selected_engine", "configured_engine", "merge_status", "merge_overlap_len",
        "merge_identity", "overlaps_saved", "overlaps_removed", "primer_mode", "primer_stage",
        "primer_preset", "primer_source", "overlap_engine_strategy", "overlap_engine_order", "overlap_quality_mode",
        "selected_assembler_id", "selected_assembler_name", "assembler_mode", "selection_reason",
        "audit_status", "audit_overlap_identity", "audit_overlap_quality", "audit_overlap_orientation",
    ]
    with output_tsv.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(headers) + "\n")
        for sample_id in sorted(set(paired_samples) | set(missing_samples)):
            chosen = selected_rows.get(sample_id, {})
            blast = blast_rows.get(sample_id, {})
            overlap = (overlap_rows or {}).get(sample_id, {})
            row = {
                "sample_id": sample_id,
                "status": chosen.get("status", "pair_missing" if sample_id in missing_samples else "no_payload"),
                "assembler": chosen.get("assembler_id", "none"),
                "contig_len": chosen.get("contig_len", ""),
                "blast_payload": blast.get("blast_payload", ""),
                "payload_kind": blast.get("payload_kind", "none"),
                "payload_n": blast.get("payload_n", "0"),
                "payload_entity_n": blast.get("payload_entity_n", blast.get("payload_n", "0")),
                "payload_max_len": blast.get("payload_max_len", "0"),
                "ambiguity_flag": blast.get("ambiguity_flag", "0"),
                "safety_flag": blast.get("safety_flag", "none"),
                "decision_source": blast.get("decision_source", "auto"),
                "review_reason": blast.get("review_reason", ""),
                "warning_flags": blast.get("warning_flags", ""),
                "resolution_state": blast.get("resolution_state", ""),
                "review_action": blast.get("review_action", "none"),
                "advisory_reason": blast.get("advisory_reason", ""),
                "structural_hypothesis_n": blast.get("structural_hypothesis_n", "0"),
                "hypotheses_with_hits_n": blast.get("hypotheses_with_hits_n", "0"),
                "missing_hits_n": blast.get("missing_hits_n", "0"),
                "hypothesis_map": blast.get("hypothesis_map", ""),
                "source_id_map": blast.get("source_id_map", ""),
                "resolved_hypothesis": blast.get("resolved_hypothesis", ""),
                "resolution_reason": blast.get("resolution_reason", ""),
                "trace_status": blast.get("trace_status", ""),
                "trace_status_f": blast.get("trace_status_f", ""),
                "trace_status_r": blast.get("trace_status_r", ""),
                "trace_flags": blast.get("trace_flags", ""),
                "selected_engine": chosen.get("selected_engine", ""),
                "configured_engine": "",
                "merge_status": "",
                "merge_overlap_len": "",
                "merge_identity": "",
                "overlaps_saved": "",
                "overlaps_removed": "",
                "primer_mode": primer_mode,
                "primer_stage": primer_stage,
                "primer_preset": primer_preset,
                "primer_source": primer_source,
                "overlap_engine_strategy": overlap_engine_strategy,
                "overlap_engine_order": overlap_engine_order,
                "overlap_quality_mode": overlap_quality_mode,
                "selected_assembler_id": chosen.get("assembler_id", ""),
                "selected_assembler_name": chosen.get("assembler_name", ""),
                "assembler_mode": assembler_mode_label,
                "selection_reason": blast.get("reason", ""),
                "audit_status": overlap.get("status", ""),
                "audit_overlap_identity": overlap.get("overlap_identity", ""),
                "audit_overlap_quality": overlap.get("overlap_quality", ""),
                "audit_overlap_orientation": overlap.get("orientation", ""),
            }
            fh.write("\t".join(str(row.get(h, "")) for h in headers) + "\n")


def _apply_resolution_to_assembly_summary(
    assembly_summary: Path | None,
    blast_rows: dict[str, dict[str, str]],
) -> None:
    if assembly_summary is None or not assembly_summary.exists():
        return
    try:
        df = pd.read_csv(assembly_summary, sep="\t")
    except Exception:
        return
    if df.empty or "sample_id" not in df.columns:
        return

    for col in (
        "resolution_state",
        "resolved_hypothesis",
        "resolution_reason",
        "review_action",
        "advisory_reason",
        "warning_flags",
        "structural_hypothesis_n",
        "hypotheses_with_hits_n",
        "missing_hits_n",
        "hypothesis_map",
        "source_id_map",
        "trace_status",
        "trace_status_f",
        "trace_status_r",
        "trace_flags",
    ):
        if col not in df.columns:
            df[col] = ""

    for idx, row in df.iterrows():
        sid = str(row.get("sample_id", "") or "")
        payload = blast_rows.get(sid, {})
        for col in (
            "resolution_state",
            "resolved_hypothesis",
            "resolution_reason",
            "review_action",
            "advisory_reason",
            "warning_flags",
            "structural_hypothesis_n",
            "hypotheses_with_hits_n",
            "missing_hits_n",
            "hypothesis_map",
            "source_id_map",
            "trace_status",
            "trace_status_f",
            "trace_status_r",
            "trace_flags",
        ):
            if payload.get(col, ""):
                df.at[idx, col] = str(payload.get(col, ""))

    df.to_csv(assembly_summary, sep="\t", index=False)

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
    assembler_id: str | None = None,
    assembler_mode: str | None = None,
    write_blast_inputs: bool = True,
    use_blast_inputs: bool = True,
    fwd_pattern: str | None = None,
    rev_pattern: str | None = None,
    enforce_same_well: bool = False,
    well_pattern: str | re.Pattern[str] | None = None, 
    metadata: Path | None = None,
    summary_tsv: Path | None = None,
    collapse_replicates: bool = False,
    replicate_id_th: float | None = None,
    min_replicate_size: int | None = None,
    orient_reads: bool = False,
    orient_db: PathLike | None = None,
    chimera_mode: str = "off",
    chimera_db: PathLike | None = None,
    overlap_audit: bool = False, 
    primer_cfg_override: dict | None = None,
    on_stage=None,
    on_progress=None,
    stop_cb: Callable[[], bool] | None = None,
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
    _check_cancel(stop_cb)

    if out_dir is None:
        out_dir = _default_run_output_dir(infile)
    out_dir.mkdir(parents=True, exist_ok=True)
    on_progress(0)

    cfg = load_config()
    if primer_cfg_override:
        merged_primer = dict(cfg.get("primer_trim", {}))
        merged_primer.update(primer_cfg_override)
        cfg["primer_trim"] = merged_primer
    primer_cfg = _normalize_primer_trim_cfg(cfg)
    reporting_cfg = _reporting_flags(cfg)
    overlap_min_overlap, overlap_min_identity, overlap_min_quality = _resolve_overlap_thresholds(cfg)
    overlap_quality_mode = str(cfg.get("overlap_eval", {}).get("quality_mode", "warning")).strip().lower()
    if overlap_quality_mode not in {"blocking", "warning"}:
        overlap_quality_mode = "warning"
    merge_cfg = cfg.get("merge_two_reads", {})
    configured_overlap_engine = str(merge_cfg.get("overlap_engine", "auto")).strip().lower()
    overlap_strategy = str(merge_cfg.get("overlap_engine_strategy", "cascade")).strip().lower()
    raw_engine_order = merge_cfg.get("overlap_engine_order", ["biopython", "ungapped", "edlib"])
    if isinstance(raw_engine_order, str):
        raw_engine_order = [tok.strip() for tok in raw_engine_order.split(",") if tok.strip()]
    overlap_engine_order_str = ",".join(str(tok) for tok in raw_engine_order)
    orient_db_path: Path | None = None
    if orient_reads:
        if orient_db:
            orient_db_path = Path(orient_db).expanduser().resolve()
        else:
            orient_ref = cfg["databases"].get(db_key, {}).get("orient_ref")
            if not orient_ref:
                orient_ref = cfg["databases"].get(db_key, {}).get("chimera_ref")
            if orient_ref:
                orient_db_path = Path(expand_db_path(orient_ref))
        if orient_db_path is None:
            raise ValueError(
                f"No orient reference configured for '{db_key}'. "
                "Provide --orient-db or add databases.<key>.orient_ref."
            )
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
    resolved_assembler_mode = (assembler_mode or "").strip().lower() or None
    resolved_assembler_id = (assembler_id or "").strip() or None
    if resolved_assembler_mode not in {None, "cap3_default", "selected", "all"}:
        raise ValueError("assembler_mode must be one of: cap3_default, selected, all")
    if resolved_assembler_mode is None:
        if resolved_assembler_id in {None, "", "__cap3_default__"}:
            resolved_assembler_mode = "cap3_default"
        elif resolved_assembler_id == "__all__":
            resolved_assembler_mode = "all"
        else:
            resolved_assembler_mode = "selected"
    if resolved_assembler_mode == "all":
        resolved_assembler_id = None
    elif resolved_assembler_mode == "cap3_default":
        resolved_assembler_id = None

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
        "review_queue": out_dir / "qc" / "review_queue.tsv" if using_paired else None,
        "orient_fasta": None,
        "orient_notmatched": None,
        "orient_report": None,
    }

    extra_stages = (
        int(collapse_replicates) + int(chimera_mode != "off") + int(orient_reads)
    )

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


        assembly_input = infile
        pretrim_fasta_dir: Path | None = None
        if not is_fasta:
            on_stage("Trim")
            run_trim(infile, out_dir, sanger=True, summary_tsv=summary_tsv, primer_cfg_override=primer_cfg_override)
            _check_cancel(stop_cb)
            pct += step
            on_progress(pct)

            on_stage("Convert")
            fastq_dir = out_dir / "passed_qc_fastq_primer_trim"
            if not any(fastq_dir.glob("*.fastq")):
                fastq_dir = out_dir / "passed_qc_fastq"
            if not any(fastq_dir.glob("*.fastq")):
                fastq_dir = out_dir / "qc"
            run_fastq_to_fasta(fastq_dir, paths["trimmed_fasta"])
            _check_cancel(stop_cb)
            pct += step
            on_progress(pct)

            pretrim_fastq_dir = out_dir / "passed_qc_fastq"
            if any(pretrim_fastq_dir.glob("*.fastq")):
                pretrim_fasta_dir = out_dir / "qc" / "paired_fasta_pretrim"
                stage_paired_fastas_from_fastq_dir(pretrim_fastq_dir, pretrim_fasta_dir, use_qual=cap3_use_qual)

            fasta_dir = out_dir / "qc" / "paired_fasta"
            generated_fastas = stage_paired_fastas_from_fastq_dir(fastq_dir, fasta_dir, use_qual=cap3_use_qual)
            if not generated_fastas:
                raise ValueError(f"No FASTA files generated from {fastq_dir}")
            assembly_input = fasta_dir
        

        pairing_report = out_dir / "qc" / "pairing_report.tsv"
        pairing_report.parent.mkdir(parents=True, exist_ok=True)

        pretrim_paired_samples: dict[str, dict[str, list[Path]]] | None = None
        primer_trim_bases_by_sample = _collect_primer_trim_bases_by_sample(out_dir / "qc" / "primer_trim_report.tsv")

        on_stage("Paired assembly")
        paired_samples, missing_samples = _collect_pairing_catalog(
            assembly_input,
            fwd_pattern=fwd_pattern,
            rev_pattern=rev_pattern,
            dup_policy=dup_policy,
            enforce_same_well=enforce_same_well,
            well_pattern=well_pattern,
        )
        if not paired_samples:
           summary = _summarize_paired_candidates(assembly_input, fwd_pattern, rev_pattern, enforce_same_well=enforce_same_well, well_pattern=well_pattern)
           suggestions = _suggest_pairing_patterns_staged(assembly_input)
           raise ValueError(
                f"No paired reads detected in {assembly_input}; {summary}. "
                "If your primer names differ, provide --fwd-pattern/--rev-pattern."
                f" {suggestions}"
            )

        if pretrim_fasta_dir is not None and pretrim_fasta_dir.exists():
            pretrim_paired_samples, _pretrim_missing = _collect_pairing_catalog(
                pretrim_fasta_dir,
                fwd_pattern=fwd_pattern,
                rev_pattern=rev_pattern,
                dup_policy=dup_policy,
                enforce_same_well=enforce_same_well,
                well_pattern=well_pattern,
            )

        contig_paths: list[Path] = []
        if resolved_assembler_mode == "cap3_default":
            cap3_args = _resolve_cap3_options(cap3_profile, cap3_options, cap3_extra_args)
            assembly_progress = subprog(pct) 
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
                use_qual=cap3_use_qual,
                on_stage=on_stage,
                on_progress=assembly_progress,
            )
        else:
            if resolved_assembler_mode == "selected":
                get_assembler_spec(resolved_assembler_id or "")
            compare_progress = subprog(pct) 
            compare_assembler_ids = [resolved_assembler_id] if resolved_assembler_mode == "selected" and resolved_assembler_id else None 
            run_compare_assemblers(
                assembly_input,
                out_dir,
                assembler_ids=compare_assembler_ids, 
                fwd_pattern=fwd_pattern,
                rev_pattern=rev_pattern,
                dup_policy=dup_policy,
                enforce_same_well=enforce_same_well,
                well_pattern=well_pattern,
                pairing_report=pairing_report,
                use_qual=cap3_use_qual,
                threads=threads,
                on_stage=on_stage,
                on_progress=on_progress,
            )
        if resolved_assembler_mode == "cap3_default" and pairing_report.exists():
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
                overlap_status = _write_overlap_audit(
                    paired_samples,
                    paths["overlap_audit"],
                    min_overlap=overlap_min_overlap,
                    min_identity=overlap_min_identity,
                    min_quality=overlap_min_quality,
                    quality_mode=overlap_quality_mode,
                    emit_engine_audit=reporting_cfg["emit_engine_audit"],
                    pretrim_paired_samples=pretrim_paired_samples,
                    primer_trim_bases_by_sample=primer_trim_bases_by_sample,
                )
            else:
                L.warning("Overlap audit skipped; paired input is not a directory: %s", assembly_input)

        assembly_summary = paths["assembly_summary"]
        overlap_rows = _load_tsv_rows_by_sample(paths["overlap_audit"]) if paths.get("overlap_audit") else {}
        blast_payload_rows: dict[str, dict[str, str]] = {}
        trace_by_sample = _load_trace_qc_by_sample(summary_tsv, paired_samples)
        merged_contigs = out_dir / "asm" / "paired_contigs.fasta"

        if use_blast_inputs and not write_blast_inputs:
            L.warning("use_blast_inputs requested without write_blast_inputs; enabling blast input output.")
            write_blast_inputs = True

        blast_inputs_fasta = out_dir / "asm" / "blast_inputs.fasta"
        blast_inputs_tsv = out_dir / "asm" / "blast_inputs.tsv"

        if resolved_assembler_mode == "cap3_default":
            _merge_cap3_contigs(contig_paths, merged_contigs)
            if write_blast_inputs:
                _build_blast_inputs(
                    out_dir / "asm",
                    paired_samples,
                    missing_samples,
                    blast_inputs_fasta,
                    blast_inputs_tsv,
                )
                blast_payload_rows = _apply_trace_to_blast_payloads(
                    _apply_overlap_diagnostics_to_blast_payloads(
                        _load_tsv_rows_by_sample(blast_inputs_tsv),
                        overlap_rows,
                    ),
                    trace_by_sample,
                    trace_cfg=cfg.get("trace_qc", {}),
                )
                _write_tsv_rows_by_sample(blast_inputs_tsv, blast_payload_rows)
            if assembly_summary:
                parse_cap3_reports(
                    out_dir / "asm",
                    sorted(set(paired_samples.keys()) | set(missing_samples)),
                    output_tsv=assembly_summary,
                    missing_samples=missing_samples,
                    overlap_status=overlap_status,
                    overlap_rows=overlap_rows,
                    blast_payloads=blast_payload_rows,
                    primer_mode=str(primer_cfg["mode"]),
                    primer_stage=str(primer_cfg["stage"]),
                    primer_preset=str(primer_cfg.get("preset", "")),
                    primer_source=str(primer_cfg.get("primer_source", "custom")),
                    configured_engine=configured_overlap_engine,
                    overlap_engine_strategy=overlap_strategy,
                    overlap_engine_order=overlap_engine_order_str,
                    overlap_quality_mode=overlap_quality_mode,
                )

            if not reporting_cfg["emit_per_sample_merge_report"]:
                for sample_id in paired_samples.keys():
                    rp = (out_dir / "asm" / sample_id / f"{sample_id}_paired.merge_report.tsv")
                    if rp.exists():
                        rp.unlink()
        else:
            compare_tsv = out_dir / "asm" / "compare_assemblers.tsv"
            selected_rows = _select_best_compare_rows(
                compare_tsv,
                assembler_mode=resolved_assembler_mode,
                assembler_id=resolved_assembler_id,
            )
            no_payload_reason = "selected_backend_no_payload" if resolved_assembler_mode == "selected" else "winner_no_payload"
            blast_payload_rows = _apply_trace_to_blast_payloads(
                _apply_overlap_diagnostics_to_blast_payloads(
                    _build_selected_blast_inputs(
                selected_rows,
                paired_samples,
                missing_samples,
                blast_inputs_fasta,
                blast_inputs_tsv,
                no_payload_reason=no_payload_reason,
                    ),
                    overlap_rows,
                ),
                trace_by_sample,
                trace_cfg=cfg.get("trace_qc", {}),
            )
            _write_tsv_rows_by_sample(blast_inputs_tsv, blast_payload_rows)
            if assembly_summary:
                _write_selected_assembly_summary(
                    assembly_summary,
                    paired_samples,
                    missing_samples,
                    selected_rows,
                    blast_payload_rows,
                    resolved_assembler_mode,
                    primer_mode=str(primer_cfg["mode"]),
                    primer_stage=str(primer_cfg["stage"]),
                    primer_preset=str(primer_cfg.get("preset", "") or "custom"),
                    primer_source=str(primer_cfg.get("primer_source", "custom")),
                    overlap_engine_strategy=overlap_strategy,
                    overlap_engine_order=overlap_engine_order_str,
                    overlap_quality_mode=overlap_quality_mode,
                    overlap_rows=overlap_rows,
                )
            merged_contigs = blast_inputs_fasta

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
        run_trim(infile, out_dir, sanger=sanger, summary_tsv=summary_tsv, primer_cfg_override=primer_cfg_override)

        pct += step
        on_progress(pct)

        # 2 – Merge FASTQs to FASTA
        on_stage("Convert")
        fastq_dir = out_dir / "passed_qc_fastq_primer_trim"
        if not any(fastq_dir.glob("*.fastq")):
            alt = out_dir / "passed_qc_fastq"
            if alt.exists() and any(alt.glob("*.fastq")):
                fastq_dir = alt
            else:
                fastq_dir = out_dir / "qc"
        run_fastq_to_fasta(fastq_dir, paths["fasta"])
        _check_cancel(stop_cb)
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

    def _concat_fastas(source_a: Path, source_b: Path, dest: Path) -> None:
        from Bio import SeqIO

        dest.parent.mkdir(parents=True, exist_ok=True)
        with dest.open("w", encoding="utf-8") as out_fh:
            for fp in (source_a, source_b):
                if not fp.exists():
                    continue
                for rec in SeqIO.parse(fp, "fasta"):
                    SeqIO.write(rec, out_fh, "fasta")

    if orient_reads:
        on_stage("Orient reads")
        oriented = out_dir / "qc" / "oriented.fasta"
        notmatched = out_dir / "qc" / "orient_notmatched.fasta"
        orient_report = out_dir / "qc" / "orient_report.tsv"
        oriented, notmatched, orient_report = vsearch_orient_reads(
            current_fasta,
            oriented,
            reference=orient_db_path,
            notmatched_out=notmatched,
            tabbed_out=orient_report,
            threads=threads,
        )
        paths["orient_fasta"] = oriented
        paths["orient_notmatched"] = notmatched
        paths["orient_report"] = orient_report
        current_fasta = oriented
        pct += step
        on_progress(pct)

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

    if orient_reads and paths["orient_notmatched"]:
        merged = out_dir / "qc" / "oriented_combined.fasta"
        _concat_fastas(current_fasta, paths["orient_notmatched"], merged)
        current_fasta = merged

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
        stop_cb=stop_cb,
    )
    _check_cancel(stop_cb)
    pct += step
    on_progress(pct)

    # 4 – Add taxonomy
    tax_template = cfg["databases"][db_key]["taxonomy"]
    tax_fp = Path(expand_db_path(tax_template))
    on_stage("Taxonomy")
    run_add_tax(paths["hits"], tax_fp, paths["tax"])
    _check_cancel(stop_cb)
    pct += step
    on_progress(pct)

    if using_paired:
        try:
            tax_df = pd.read_csv(paths["tax"], sep="\t")
            if "qseqid" in tax_df.columns:
                normalized_sample_ids = tax_df["qseqid"].map(qseqid_to_sample_id)
                if "sample_id" not in tax_df.columns or not tax_df["sample_id"].astype(str).equals(normalized_sample_ids.astype(str)):
                    tax_df["sample_id"] = normalized_sample_ids
                    tax_df.to_csv(paths["tax"], sep="\t", index=False) 
        except Exception as exc:
            L.warning("Failed to add paired sample_id column to taxonomy table: %s", exc)

        try:
            resolved_rows = _annotate_resolution_from_tax(paths["tax"], blast_payload_rows)
            _apply_resolution_to_assembly_summary(paths.get("assembly_summary"), resolved_rows)
            _write_review_queue_from_resolved(paths["review_queue"], resolved_rows)
            _write_tsv_rows_by_sample(blast_inputs_tsv, resolved_rows)
            blast_payload_rows = resolved_rows
        except Exception as exc:
            L.warning("Failed to write review queue: %s", exc)

    # 5 – Optional post-BLAST
    if postblast:
        if metadata is None:
            raise ValueError("metadata file required for postblast stage")
        on_stage("Post-BLAST")
        run_postblast(paths["tax"], metadata, paths["biom"], weights_tsv=paths["replicate_weights"], stop_cb=stop_cb)
        _check_cancel(stop_cb)
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


def run_pairing_report(
    input_dir: PathLike,
    output_tsv: PathLike,
    *,
    fwd_pattern: str | None = None,
    rev_pattern: str | None = None,
    dup_policy: DupPolicy = DupPolicy.ERROR,
    enforce_same_well: bool = False,
    well_pattern: str | re.Pattern[str] | None = None,
) -> Path:
    """
    Write the canonical pairing report TSV for a staged paired-input directory.

    This reuses MicroSeq's pairing logic and returns the output TSV path.

    """
    # Normalize input directory as a path object 
    input_dir = Path(input_dir)

    # Normalize tsv output as a path object 
    output_tsv = Path(output_tsv) 

    # Ensure parent directory exists before writing 
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Reuse same pairing registry used by core codebase here 
    paired_samples, pairing_metadata = group_pairs(
        input_dir,
        dup_policy=dup_policy,
        fwd_pattern=fwd_pattern,
        rev_pattern=rev_pattern,
        enforce_same_well=enforce_same_well,
        well_pattern=well_pattern,
        return_metadata=True,
    )

    # Reuse the canonical TSV writer for pairing reports 
    _write_pairing_report(paired_samples, pairing_metadata, output_tsv)

    # Return path that was written 
    return output_tsv 


def run_assembly_summary(
    asm_dir: PathLike,
    pairing_input_dir: PathLike,
    output_tsv: PathLike,
    *,
    fwd_pattern: str | None = None,
    rev_pattern: str | None = None, 
    dup_policy: DupPolicy = DupPolicy.ERROR,
    enforce_same_well: bool = False,
    well_pattern: str | re.Pattern[str] | None = None,
) -> Path:
    """
    Write the canonical assembly summary TSV for a paired assembly run.

    This reuses MicroSeq's summary parser and returns the output TSV path.

    """
    # Normalize the assembly directory path.
    asm_dir = Path(asm_dir)

    # Normalize the staged pairing-input directory path.
    pairing_input_dir = Path(pairing_input_dir)

    # Normalize the output TSV path.
    output_tsv = Path(output_tsv)

    # Ensure the parent directory exists before writing.
    output_tsv.parent.mkdir(parents=True, exist_ok=True)

    # Reconstruct the same pairing catalog used to interpret per-sample outputs.
    paired_samples, missing_samples = _collect_pairing_catalog(
        pairing_input_dir,
        fwd_pattern=fwd_pattern,
        rev_pattern=rev_pattern,
        dup_policy=dup_policy,
        enforce_same_well=enforce_same_well,
        well_pattern=well_pattern,
    )

    # Keep both assembled and pair-missing sample IDs in the final summary.
    sample_ids = sorted(set(paired_samples) | set(missing_samples))

    # Reuse the canonical summary writer from cap3_report.py.
    parse_cap3_reports(
        asm_dir,
        sample_ids,
        output_tsv=output_tsv,
        missing_samples=missing_samples,
    )

    # Return the concrete path that was written.
    return output_tsv 
   
