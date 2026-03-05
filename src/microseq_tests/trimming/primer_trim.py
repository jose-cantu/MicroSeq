# microseq_tests/src/microseq_tests/trimming/primer_trim.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

from Bio import SeqIO

from microseq_tests.trimming.biopy_trim import trim_record_quality


IUPAC_BASES: dict[str, set[str]] = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "U": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}


@dataclass(frozen=True)
class PrimerTrimResult:
    file_name: str
    reads: int
    reads_trimmed: int
    bases_trimmed: int
    avg_trim_len: float
    avg_mismatches: float
    primer_hits: str
    avg_offset: float
    orientation_used: str
    orientation_source: str
    iupac_mode: str
    primer_bases_trimmed_total: int
    post_trim_bases_removed_total: int = 0
    post_trim_final_len_avg: float = 0.0


def _best_primer_match(
    seq: str,
    primers: Iterable[str],
    *,
    max_mismatch: int,
    max_search: int,
    max_primer_offset: int,
    iupac_mode: bool,
) -> tuple[int, int, int, str] | None:
    if not primers:
        return None
    seq = seq.upper()
    search_len = min(len(seq), max_search)
    best: tuple[int, int, int, str] | None = None
    for primer in primers:
        primer = primer.strip().upper()
        if not primer:
            continue
        p_len = len(primer)
        if p_len > search_len:
            continue
        max_offset = min(search_len - p_len, max_primer_offset)
        for offset in range(max_offset + 1):
            window = seq[offset : offset + p_len]
            if iupac_mode:
                mismatches = sum(
                    1
                    for base, p_base in zip(window, primer)
                    if base.upper() not in IUPAC_BASES.get(p_base.upper(), {p_base.upper()})
                )
            else:
                mismatches = sum(1 for a, b in zip(window, primer) if a != b)
            if mismatches > max_mismatch:
                continue
            candidate = (mismatches, offset, p_len, primer)
            if best is None or candidate < best:
                best = candidate
    return best


def trim_primer_fastqs(
    input_dir: Path,
    output_dir: Path,
    *,
    forward_primers: Iterable[str],
    reverse_primers: Iterable[str],
    max_mismatch: int = 2,
    max_search: int = 60,
    max_primer_offset: int = 10,
    iupac_mode: bool = True,
    report_path: Path | None = None,
    detect_report_path: Path | None = None,
    orientation_resolver: Callable[[str], str | None] | None = None,
    mode: str = "clip",
    post_quality_trim_enabled: bool = False,
    post_quality_method: str = "mott",
    post_quality_cutoff_q: int = 22,
    post_quality_window_size: int = 5,
    post_quality_per_base_q: int = 22,
    post_quality_rescue_5prime_bases: int = 0,
    post_quality_min_len: int = 50,
) -> dict[str, PrimerTrimResult]:
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, PrimerTrimResult] = {}
    fwd_primers = list(forward_primers)
    rev_primers = list(reverse_primers)

    mode_norm = mode.strip().lower()
    if mode_norm not in {"clip", "detect"}:
        raise ValueError(f"Invalid primer trim mode '{mode}'. Expected clip or detect.")

    report_lines = [
        "file\treads\treads_trimmed\tbases_trimmed\tavg_trim_len\tavg_mismatches\tprimer_hit\tavg_offset\torientation_used\torientation_source\tiupac_mode"
    ]
    detect_lines = [
        "file\treads_scanned\treads_with_primer_hit\thit_rate\tavg_mismatches\tavg_offset\tmatched_primers"
    ]
    for fastq_path in sorted(input_dir.glob("*.fastq")):
        reads = 0
        reads_trimmed = 0
        bases_trimmed = 0
        mismatches_total = 0
        offset_total = 0
        primer_hits: set[str] = set()
        trimmed_records = []
        post_trim_bases_removed = 0
        post_trim_final_len_sum = 0
        orientation_used = "combined"
        orientation_source = "combined"

        resolved_orient = orientation_resolver(fastq_path.name) if orientation_resolver else None
        if resolved_orient in {"F", "R"}:
            if resolved_orient == "F":
                candidate_primers = fwd_primers
                orientation_used = "forward"
            else:
                candidate_primers = rev_primers
                orientation_used = "reverse"
            orientation_source = "detector"
        else:
            upper_name = fastq_path.stem.upper()
            if upper_name.endswith("F") or "_F" in upper_name or "27F" in upper_name:
                candidate_primers = fwd_primers
                orientation_used = "forward"
                orientation_source = "filename"
            elif upper_name.endswith("R") or "_R" in upper_name or "1492R" in upper_name:
                candidate_primers = rev_primers
                orientation_used = "reverse"
                orientation_source = "filename"
            else:
                candidate_primers = fwd_primers + rev_primers

        for record in SeqIO.parse(fastq_path, "fastq"):
            reads += 1
            seq = str(record.seq)
            match = _best_primer_match(
                seq,
                candidate_primers,
                max_mismatch=max_mismatch,
                max_search=max_search,
                max_primer_offset=max_primer_offset,
                iupac_mode=iupac_mode,
            )
            if match is None:
                rec_for_post = record
            else:
                mismatches, offset, p_len, primer = match
                trim_len = offset + p_len
                if trim_len >= len(record) and mode_norm == "clip":
                    continue
                rec_for_post = record
                reads_trimmed += 1
                if mode_norm == "clip":
                    rec_for_post = record[trim_len:]
                    bases_trimmed += trim_len
                mismatches_total += mismatches
                offset_total += offset
                primer_hits.add(primer)

            if post_quality_trim_enabled:
                post_trimmed, left_trim_post, right_trim_post = trim_record_quality(
                    rec_for_post,
                    method=post_quality_method,
                    cutoff_q=post_quality_cutoff_q,
                    window_size=post_quality_window_size,
                    per_base_q=post_quality_per_base_q,
                    rescue_5prime_bases=post_quality_rescue_5prime_bases,
                    min_len=post_quality_min_len,
                )
                if post_trimmed is None:
                    continue
                post_trim_bases_removed += left_trim_post + right_trim_post
                rec_for_post = post_trimmed

            post_trim_final_len_sum += len(rec_for_post)
            trimmed_records.append(rec_for_post)
            continue

        if mode_norm == "clip":
            out_path = output_dir / fastq_path.name
            SeqIO.write(trimmed_records, out_path, "fastq")
        avg_trim_len = bases_trimmed / reads_trimmed if reads_trimmed else 0.0
        avg_mismatches = mismatches_total / reads_trimmed if reads_trimmed else 0.0
        avg_offset = offset_total / reads_trimmed if reads_trimmed else 0.0
        primer_hit_text = ";".join(sorted(primer_hits)) if primer_hits else ""
        post_trim_final_len_avg = (post_trim_final_len_sum / len(trimmed_records)) if trimmed_records else 0.0
        result = PrimerTrimResult(
            file_name=fastq_path.name,
            reads=reads,
            reads_trimmed=reads_trimmed,
            bases_trimmed=bases_trimmed,
            avg_trim_len=avg_trim_len,
            avg_mismatches=avg_mismatches,
            primer_hits=primer_hit_text,
            avg_offset=avg_offset,
            orientation_used=orientation_used,
            orientation_source=orientation_source,
            iupac_mode="on" if iupac_mode else "off",
            primer_bases_trimmed_total=bases_trimmed,
            post_trim_bases_removed_total=post_trim_bases_removed,
            post_trim_final_len_avg=post_trim_final_len_avg,
        )
        results[fastq_path.name] = result
        report_lines.append(
            "\t".join(
                [
                    fastq_path.name,
                    str(reads),
                    str(reads_trimmed),
                    str(bases_trimmed),
                    f"{avg_trim_len:.2f}",
                    f"{avg_mismatches:.2f}",
                    primer_hit_text or "-",
                    f"{avg_offset:.2f}",
                    orientation_used,
                    orientation_source,
                    "on" if iupac_mode else "off",
                ]
            )
        )
        hit_rate = (reads_trimmed / reads) if reads else 0.0
        detect_lines.append(
            "\t".join(
                [
                    fastq_path.name,
                    str(reads),
                    str(reads_trimmed),
                    f"{hit_rate:.4f}",
                    f"{avg_mismatches:.2f}",
                    f"{avg_offset:.2f}",
                    primer_hit_text or "-",
                ]
            )
        )

    if report_path:
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    if detect_report_path:
        detect_report_path.parent.mkdir(parents=True, exist_ok=True)
        detect_report_path.write_text("\n".join(detect_lines) + "\n", encoding="utf-8")

    return results


def update_trim_summary(
    summary_tsv: Path,
    primer_results: dict[str, PrimerTrimResult],
) -> None:
    summary_path = Path(summary_tsv)
    if not summary_path.exists():
        return
    lines = summary_path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return
    header = lines[0].split("\t")
    target_cols = [
        "post_trim_bases_removed",
        "post_trim_final_len_avg",
        "primer_trimmed",
        "primer_bases_trimmed",
        "primer_trim_len_avg",
        "primer_bases_trimmed_total",
        "primer_mismatch_avg",
        "primer_hit",
        "primer_offset_avg",
        "primer_orientation",
        "primer_orientation_source",
        "primer_iupac_mode",
    ]
    for col in target_cols:
        if col not in header:
            header.append(col)
    col_idx = {name: idx for idx, name in enumerate(header)}
    primer_lookup: dict[str, PrimerTrimResult] = {}
    for name, result in primer_results.items():
        primer_lookup[name] = result
        if name.endswith("_trimmed.fastq"):
            primer_lookup[name.replace("_trimmed.fastq", ".fastq")] = result
        if name.endswith("_trimmed.fq"):
            primer_lookup[name.replace("_trimmed.fq", ".fq")] = result

    updated_lines = ["\t".join(header)]

    for line in lines[1:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < len(header):
            parts.extend([""] * (len(header) - len(parts)))
        file_name = parts[0]
        result = primer_lookup.get(file_name)
        if result:
            parts[col_idx["post_trim_bases_removed"]] = str(result.post_trim_bases_removed_total)
            parts[col_idx["post_trim_final_len_avg"]] = f"{result.post_trim_final_len_avg:.2f}"
            parts[col_idx["primer_trimmed"]] = "yes" if result.reads_trimmed else "no"
            parts[col_idx["primer_bases_trimmed"]] = f"{result.avg_trim_len:.2f}"
            parts[col_idx["primer_trim_len_avg"]] = f"{result.avg_trim_len:.2f}"
            parts[col_idx["primer_bases_trimmed_total"]] = str(result.primer_bases_trimmed_total)
            parts[col_idx["primer_mismatch_avg"]] = f"{result.avg_mismatches:.2f}"
            parts[col_idx["primer_hit"]] = result.primer_hits or ""
            parts[col_idx["primer_offset_avg"]] = f"{result.avg_offset:.2f}"
            parts[col_idx["primer_orientation"]] = result.orientation_used
            parts[col_idx["primer_orientation_source"]] = result.orientation_source
            parts[col_idx["primer_iupac_mode"]] = result.iupac_mode
        updated_lines.append("\t".join(parts))

    summary_path.write_text("\n".join(updated_lines) + "\n", encoding="utf-8")
