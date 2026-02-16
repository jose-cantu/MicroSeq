# MicroSeq/src/microseq_tests/assembly/cap3_report.py 

from __future__ import annotations

""" 
Parse CAP3 output files into an assembly summary table. 

Rows are keyed by the CAP3 *sample_key* (the directory name that CAP3 writes
into). This means keep-separate mode keeps `sid_1`, `sid_2`, etc. as distinct
rows, which avoids ambiguity for contig/singlet counts and status values.
"""

from pathlib import Path 
import logging 
import re 
from typing import Iterable, Mapping 

from Bio import SeqIO 

L = logging.getLogger(__name__) 

_OVERLAP_SAVED_RX = re.compile(r"Number of overlaps saved:\s*(\d+)")
_OVERLAP_REMOVED_RX = re.compile(r"Number of overlaps removed:\s*(\d+)")
_CLIP_RX = re.compile(
    r"Clip\s+(?P<read>\S+)\s+left clip:\s*(?P<left>\d+),\s*right clip:\s*(?P<right>\d+)",
    re.I,
)


def _count_fasta_records(path: Path) -> int:
    if not path.exists():
        return 0
    return sum(1 for _ in SeqIO.parse(path, "fasta"))


def _contig_len_from_fasta(path: Path) -> int | None:
    if not path.exists():
        return None
    lengths = [len(rec.seq) for rec in SeqIO.parse(path, "fasta")]
    if not lengths:
        return None
    return max(lengths)


def parse_cap3_reports(
    asm_dir: Path,
    sample_keys: Iterable[str],
    *,
    output_tsv: Path | None = None,
    missing_samples: Iterable[str] | None = None,
    overlap_status: Mapping[str, str] | None = None,
    overlap_rows: Mapping[str, Mapping[str, str]] | None = None,
    blast_payloads: Mapping[str, Mapping[str, str]] | None = None,
    primer_mode: str = "off",
    primer_stage: str = "post_quality",
    primer_preset: str = "",
    primer_source: str = "custom",
    configured_engine: str = "auto",
    overlap_engine_strategy: str = "cascade",
    overlap_engine_order: str = "",
    overlap_quality_mode: str = "warning",
) -> list[dict[str, str | int | None]]:
    """Parse CAP3 info/contig/singlet outputs for the provided sample keys."""

    rows: list[dict[str, str | int | None]] = []
    missing_set = set(missing_samples or [])
    for sample_key in sample_keys:
        if sample_key in missing_set:
            rows.append(
                {
                    "sample_id": sample_key,
                    "overlaps_saved": 0,
                    "overlaps_removed": 0,
                    "contigs_count": 0,
                    "singlets_count": 0,
                    "clip_f_left": None,
                    "clip_f_right": None,
                    "clip_r_left": None,
                    "clip_r_right": None,
                    "status": "pair_missing",
                    "merge_status": None,
                    "merge_overlap_engine": None,
                    "merge_overlap_len": None,
                    "merge_identity": None,
                    "merge_qualities": None,
                    "merge_warning": None,
                    "merge_high_conflict_mismatches": None,
                    "cap3_validation": None,
                    "assembler": "none",
                    "selected_engine": None,
                    "configured_engine": configured_engine,
                    "fallback_used": "n/a",
                    "audit_status": None,
                    "audit_overlap_identity": None,
                    "audit_overlap_quality": None,
                    "audit_overlap_orientation": None,
                    "contig_len": None,
                    "blast_payload": (blast_payloads or {}).get(sample_key, {}).get("blast_payload"),
                    "primer_mode": primer_mode,
                    "primer_stage": primer_stage,
                    "primer_preset": primer_preset or "custom",
                    "primer_source": primer_source,
                    "overlap_engine_strategy": overlap_engine_strategy,
                    "overlap_engine_order": overlap_engine_order,
                    "overlap_quality_mode": overlap_quality_mode,
                }
            )
            continue
        sample_dir = asm_dir / sample_key
        info_path = sample_dir / f"{sample_key}_paired.fasta.cap.info"
        contigs_path = sample_dir / f"{sample_key}_paired.fasta.cap.contigs"
        singlets_path = sample_dir / f"{sample_key}_paired.fasta.cap.singlets"
        merge_report_path = sample_dir / f"{sample_key}_paired.merge_report.tsv"
        cap3_validation_path = sample_dir / f"{sample_key}_paired.cap3_validation.txt"

        overlaps_saved = 0
        overlaps_removed = 0
        clip_f_left = clip_f_right = None
        clip_r_left = clip_r_right = None

        if info_path.exists():
            for line in info_path.read_text(encoding="utf-8").splitlines():
                if match := _OVERLAP_SAVED_RX.search(line):
                    overlaps_saved = int(match.group(1))
                elif match := _OVERLAP_REMOVED_RX.search(line):
                    overlaps_removed = int(match.group(1))
                elif match := _CLIP_RX.search(line):
                    read_id = match.group("read")
                    orient = read_id[-1].upper() if read_id else ""
                    left = int(match.group("left"))
                    right = int(match.group("right"))
                    if orient == "F":
                        clip_f_left, clip_f_right = left, right
                    elif orient == "R":
                        clip_r_left, clip_r_right = left, right
        elif not merge_report_path.exists():
            L.warning("Missing CAP3 info file for %s", sample_key)

        contigs_count = _count_fasta_records(contigs_path)
        singlets_count = _count_fasta_records(singlets_path)

        if contigs_count > 0:
            status = "assembled"
        elif singlets_count > 0:
            status = "singlets_only"
        else:
            status = "cap3_no_output"

        merge_status = None
        merge_engine = None
        merge_overlap_len = None
        merge_identity = None
        merge_qualities = None
        merge_warning = None
        merge_high_conflict_mismatches = None
        cap3_validation = None
        contig_len = None
        if merge_report_path.exists():
            lines = merge_report_path.read_text(encoding="utf-8").splitlines()
            if len(lines) >= 2:
                headers = lines[0].split("\t")
                values = lines[1].split("\t")
                report = dict(zip(headers, values, strict=False))
                merge_status = report.get("merge_status")
                merge_engine = report.get("overlap_engine")
                merge_overlap_len = report.get("overlap_len")
                merge_identity = report.get("identity")
                merge_qualities = report.get("qualities")
                merge_warning = report.get("merge_warning")
                merge_high_conflict_mismatches = report.get("high_conflict_mismatches")
                contig_len = report.get("contig_len")

        overlap_row = (overlap_rows or {}).get(sample_key, {})
        blast_row = (blast_payloads or {}).get(sample_key, {})
        assembler = "merge_two_reads" if merge_status == "merged" else "cap3"
        selected_engine = merge_engine or overlap_row.get("selected_engine")
        fallback_used = "n/a"
        if selected_engine and overlap_engine_strategy in {"cascade", "all"}:
            first_engine = (overlap_engine_order.split(",", 1)[0].strip() if overlap_engine_order else "")
            if first_engine:
                fallback_used = "yes" if selected_engine != first_engine else "no"
        elif overlap_row.get("fallback_used"):
            raw_fb = str(overlap_row.get("fallback_used")).strip().lower()
            fallback_used = raw_fb if raw_fb in {"yes", "no", "n/a"} else "n/a"
        contig_len = _contig_len_from_fasta(contigs_path) or contig_len

        if cap3_validation_path.exists():
            cap3_validation = cap3_validation_path.read_text(encoding="utf-8").strip() or None
            if cap3_validation in {"failed", "rejected"}:
                status = "cap3_unverified"

        if overlap_status:
            audit_status = overlap_status.get(sample_key)
            if audit_status and audit_status != "ok":
                status = audit_status

        rows.append(
            {
                "sample_id": sample_key,
                "overlaps_saved": overlaps_saved,
                "overlaps_removed": overlaps_removed,
                "contigs_count": contigs_count,
                "singlets_count": singlets_count,
                "clip_f_left": clip_f_left,
                "clip_f_right": clip_f_right,
                "clip_r_left": clip_r_left,
                "clip_r_right": clip_r_right,
                "status": status,
                "merge_status": merge_status,
                "merge_overlap_engine": merge_engine,
                "merge_overlap_len": merge_overlap_len,
                "merge_identity": merge_identity,
                "merge_qualities": merge_qualities,
                "merge_warning": merge_warning,
                "merge_high_conflict_mismatches": merge_high_conflict_mismatches,
                "cap3_validation": cap3_validation,
                "assembler": assembler,
                "selected_engine": selected_engine,
                "configured_engine": configured_engine,
                "fallback_used": fallback_used,
                "audit_status": overlap_row.get("status"),
                "audit_overlap_identity": overlap_row.get("overlap_identity"),
                "audit_overlap_quality": overlap_row.get("overlap_quality"),
                "audit_overlap_orientation": overlap_row.get("orientation"),
                "contig_len": contig_len,
                "blast_payload": blast_row.get("blast_payload"),
                "primer_mode": primer_mode,
                "primer_stage": primer_stage,
                "primer_preset": primer_preset or "custom",
                "primer_source": primer_source,
                "overlap_engine_strategy": overlap_engine_strategy,
                "overlap_engine_order": overlap_engine_order,
                "overlap_quality_mode": overlap_quality_mode,
            }
        )

    if output_tsv:
        output_tsv.parent.mkdir(parents=True, exist_ok=True)
        headers = [
            "sample_id",
            "overlaps_saved",
            "overlaps_removed",
            "contigs_count",
            "singlets_count",
            "clip_f_left",
            "clip_f_right",
            "clip_r_left",
            "clip_r_right",
            "status",
            "merge_status",
            "merge_overlap_engine",
            "merge_overlap_len",
            "merge_identity",
            "merge_qualities",
            "merge_warning",
            "merge_high_conflict_mismatches",
            "cap3_validation",
            "assembler",
            "selected_engine",
            "configured_engine",
            "fallback_used",
            "audit_status",
            "audit_overlap_identity",
            "audit_overlap_quality",
            "audit_overlap_orientation",
            "contig_len",
            "blast_payload",
            "primer_mode",
            "primer_stage",
            "primer_preset",
            "primer_source",
            "overlap_engine_strategy",
            "overlap_engine_order",
            "overlap_quality_mode",
        ]
        with output_tsv.open("w", encoding="utf-8") as fh:
            fh.write("\t".join(headers) + "\n")
            for row in rows:
                fh.write("\t".join("" if row[h] is None else str(row[h]) for h in headers) + "\n")

    return rows
