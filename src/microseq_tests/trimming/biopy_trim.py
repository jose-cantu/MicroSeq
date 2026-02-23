from __future__ import annotations

from importlib import import_module
import logging
from pathlib import Path
from typing import Iterable, List, Optional

from Bio import SeqIO

from microseq_tests.trimming.ab1_qc import TRACE_QC_COLUMNS, TraceQcSummary

L = logging.getLogger(__name__)


def _trace_lookup_keys(file_stem: str) -> list[str]:
    """Return ordered candidate keys for trace-QC lookup."""
    stem = str(file_stem or "").strip()
    if not stem:
        return []
    out: list[str] = []

    def _add(value: str) -> None:
        if value and value not in out:
            out.append(value)

    _add(stem)
    if stem.endswith("_trimmed"):
        _add(stem[: -len("_trimmed")])

    # primer-trim suffix variants used by stage labels in some runs
    for suf in ("_primer", "_primer_trimmed", "_pre_quality", "_post_quality"):
        if stem.endswith(suf):
            _add(stem[: -len(suf)])

    return out

TRIM_SUMMARY_COLUMNS = [
    "file",
    "reads",
    "avg_len",
    "avg_q",
    "qc_status",
    "trim_policy",
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


def build_trim_summary_row(
    *,
    file_name: str,
    reads: int,
    avg_len: float,
    avg_q: float,
    qc_status: str,
    trim_policy: str = "",
) -> list[str]:
    row = ["" for _ in TRIM_SUMMARY_COLUMNS]
    row[0] = file_name
    row[1] = str(reads)
    row[2] = f"{avg_len:.1f}"
    row[3] = f"{avg_q:.2f}"
    row[4] = qc_status
    row[5] = trim_policy
    return row


def _tick_safe(bar, inc: int = 1) -> None:
    if hasattr(bar, "update"):
        bar.update(inc)


def _dynamic_trim_bounds(quals: list[int], win: int = 5, q: int = 20) -> tuple[int, int] | None:
    n = len(quals)
    if n < win:
        return None

    left = next((i for i in range(n - win + 1) if sum(quals[i : i + win]) / win >= q), None)
    if left is None:
        return None

    right = next((j for j in range(n, win - 1, -1) if sum(quals[j - win : j]) / win >= q), None)
    if right is None or right <= left:
        return None
    return left, right


def _mott_trim_bounds(quals: list[int], cutoff_q: int = 20) -> tuple[int, int] | None:
    if not quals:
        return None

    def _best_interval(scores: list[float]) -> tuple[int, int] | None:
        best_sum = float("-inf")
        best_start = best_end = -1
        cur_sum = 0.0
        cur_start = 0
        for i, val in enumerate(scores):
            if cur_sum <= 0:
                cur_start = i
                cur_sum = val
            else:
                cur_sum += val
            if cur_sum > best_sum:
                best_sum = cur_sum
                best_start = cur_start
                best_end = i + 1
        if best_sum <= 0 or best_start < 0:
            return None
        return best_start, best_end

    scores = [q - cutoff_q for q in quals]
    best = _best_interval(scores)
    if best is None:
        return None
    return best


def trim_record_quality(
    record,
    *,
    method: str = "mott",
    cutoff_q: int = 20,
    window_size: int = 5,
    per_base_q: int = 20,
    min_len: int = 1,
    rescue_5prime_bases: int = 0,
):
    quals: List[int] = record.letter_annotations["phred_quality"]
    mode = method.strip().lower()
    if mode == "legacy_window":
        bounds = _dynamic_trim_bounds(quals, window_size, per_base_q)
    else:
        bounds = _mott_trim_bounds(quals, cutoff_q)
    if bounds is None:
        return None, 0, 0
    left, right = bounds
    if rescue_5prime_bases > 0:
        left = max(0, left - rescue_5prime_bases)
    if right - left < min_len:
        return None, 0, 0
    trimmed = record[left:right]
    return trimmed, left, len(record) - right


def dynamic_trim(record, win: int = 5, q: int = 20):
    trimmed, _, _ = trim_record_quality(
        record,
        method="legacy_window",
        window_size=win,
        per_base_q=q,
        min_len=1,
    )
    return trimmed


def _trim_all(
    fastqs: Iterable[Path],
    *,
    method: str,
    cutoff_q: int,
    window_size: int,
    per_base_q: int,
    file_q_threshold: float,
    min_len: int,
    mee_max: float | None,
    mee_min_len: int | None,
    min_reads_kept: int | None,
    max_drop_fraction: float | None,
    passed_dir: Path,
    failed_dir: Path,
    stats_root: Path,
    comb: Optional[open],
    trace_qc: dict[str, TraceQcSummary] | None,
    bar,
) -> None:
    _ = (mee_max, mee_min_len, min_reads_kept, max_drop_fraction)
    policy_desc = (
        f"{method}:cutoff_q={cutoff_q}" if method != "legacy_window" else f"legacy_window:window={window_size},q={per_base_q}"
    )
    for fq in fastqs:
        base = fq.stem
        stats_path = stats_root / f"{base}_avg_qual.txt"
        trimmed_path = passed_dir / f"{base}_trimmed.fastq"

        reads = bases = qsum = 0
        trimmed_recs = []

        with stats_path.open("w", encoding="utf-8") as fh:
            fh.write("read_id\tleft_trim_bases\tright_trim_bases\tfinal_len\tmean_q\n")
            for rec in SeqIO.parse(fq, "fastq"):
                trimmed, left_trim, right_trim = trim_record_quality(
                    rec,
                    method=method,
                    cutoff_q=cutoff_q,
                    window_size=window_size,
                    per_base_q=per_base_q,
                    min_len=min_len,
                )
                if not trimmed:
                    continue
                ph = trimmed.letter_annotations["phred_quality"]
                mean_q = sum(ph) / len(ph) if ph else 0.0
                fh.write(f"{rec.id}\t{left_trim}\t{right_trim}\t{len(trimmed)}\t{mean_q:.2f}\n")
                reads += 1
                bases += len(trimmed)
                qsum += sum(ph)
                trimmed_recs.append(trimmed)

        avg_q = qsum / bases if bases else 0
        avg_len = bases / reads if reads else 0

        if avg_q < file_q_threshold:
            (failed_dir / fq.name).write_bytes(fq.read_bytes())
            (failed_dir / stats_path.name).write_bytes(stats_path.read_bytes())
            L.info("[FAIL] %s  (avgQ %.2f)", fq.name, avg_q)
            qc_status = "fail"
        else:
            SeqIO.write(trimmed_recs, trimmed_path, "fastq")
            L.info("[PASS] %s -> %s (avgQ %.2f)", fq.name, trimmed_path, avg_q)
            qc_status = "pass"

        if comb:
            row = build_trim_summary_row(
                file_name=fq.name,
                reads=reads,
                avg_len=avg_len,
                avg_q=avg_q,
                qc_status=qc_status,
                trim_policy=policy_desc,
            )
            trace_row = []
            if trace_qc:
                match = next((k for k in _trace_lookup_keys(fq.stem) if k in trace_qc), None)
                if match is not None:
                    trace_row = trace_qc[match].to_row()
            if not trace_row:
                trace_row = ["" for _ in TRACE_QC_COLUMNS]
            comb.write("\t".join(row + trace_row) + "\n")

        _tick_safe(bar)


def trim_folder(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    method: str = "mott",
    cutoff_q: int = 20,
    window_size: int = 5,
    per_base_q: int = 20,
    file_q_threshold: float = 20.0,
    min_len: int = 200,
    mee_max: float | None = None,
    mee_min_len: int | None = None,
    min_reads_kept: int | None = None,
    max_drop_fraction: float | None = None,
    combined_tsv: str | Path | None = None,
    trace_qc: dict[str, TraceQcSummary] | None = None,
    threads: int = 1,
    **kwargs,
) -> None:
    _ = (threads, kwargs)
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    passed_dir = output_dir.parent / "passed_qc_fastq"
    failed_dir = output_dir.parent / "failed_qc_fastq"
    for p in (output_dir, passed_dir, failed_dir):
        p.mkdir(parents=True, exist_ok=True)

    comb: Optional[open] = None
    if combined_tsv:
        comb = open(combined_tsv, "a", encoding="utf-8")
        if comb.tell() == 0:
            comb.write("\t".join(TRIM_SUMMARY_COLUMNS + TRACE_QC_COLUMNS) + "\n")

    fastqs = sorted(input_dir.glob("*.fastq"))

    prog = import_module("microseq_tests.utility.progress")
    cm_or_gen = prog.stage_bar(len(fastqs), desc="trim", unit="file")

    def _run(bar):
        _trim_all(
            fastqs,
            method=method,
            cutoff_q=cutoff_q,
            window_size=window_size,
            per_base_q=per_base_q,
            file_q_threshold=file_q_threshold,
            min_len=min_len,
            mee_max=mee_max,
            mee_min_len=mee_min_len,
            min_reads_kept=min_reads_kept,
            max_drop_fraction=max_drop_fraction,
            passed_dir=passed_dir,
            failed_dir=failed_dir,
            stats_root=output_dir,
            comb=comb,
            trace_qc=trace_qc,
            bar=bar,
        )

    if hasattr(cm_or_gen, "__enter__"):
        with cm_or_gen as bar:
            _run(bar)
    else:
        try:
            bar = next(cm_or_gen)
        except StopIteration:
            bar = cm_or_gen
        _run(bar)

        if hasattr(cm_or_gen, "close"):
            cm_or_gen.close()
        if hasattr(bar, "close") and bar is not cm_or_gen:
            bar.close()

    if comb:
        comb.close()
