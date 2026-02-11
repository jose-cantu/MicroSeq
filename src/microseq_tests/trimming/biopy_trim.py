# src/microseq_tests/trimming/biopy_trim.py
from __future__ import annotations

from importlib import import_module
import logging
from pathlib import Path
from typing import Iterable, List, Optional

from Bio import SeqIO

L = logging.getLogger(__name__)

TRIM_SUMMARY_COLUMNS = [
    "file",
    "reads",
    "avg_len",
    "avg_q",
    "qc_status",
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
) -> list[str]:
    row = ["" for _ in TRIM_SUMMARY_COLUMNS]
    row[0] = file_name
    row[1] = str(reads)
    row[2] = f"{avg_len:.1f}"
    row[3] = f"{avg_q:.2f}"
    row[4] = qc_status
    return row


def _tick_safe(bar, inc: int = 1) -> None:
    if hasattr(bar, "update"):
        bar.update(inc)


def dynamic_trim(record, win: int = 5, q: int = 20):
    quals: List[int] = record.letter_annotations["phred_quality"]
    n = len(quals)
    if n < win:
        return None

    left = next((i for i in range(n - win + 1) if sum(quals[i : i + win]) / win >= q), None)
    if left is None:
        return None

    right = next((j for j in range(n, win - 1, -1) if sum(quals[j - win : j]) / win >= q), None)
    if right is None or right <= left:
        return None

    return record[left:right]


def _trim_all(
    fastqs: Iterable[Path],
    *,
    window_size: int,
    per_base_q: int,
    file_q_threshold: float,
    mee_max: float | None,
    mee_min_len: int | None,
    min_reads_kept: int | None,
    max_drop_fraction: float | None,
    passed_dir: Path,
    failed_dir: Path,
    stats_root: Path,
    comb: Optional[open],
    bar,
) -> None:
    _ = (mee_max, mee_min_len, min_reads_kept, max_drop_fraction)
    for fq in fastqs:
        base = fq.stem
        stats_path = stats_root / f"{base}_avg_qual.txt"
        trimmed_path = passed_dir / f"{base}_trimmed.fastq"

        reads = bases = qsum = 0
        trimmed_recs = []

        for rec in SeqIO.parse(fq, "fastq"):
            trimmed = dynamic_trim(rec, window_size, per_base_q)
            if not trimmed:
                continue
            ph = trimmed.letter_annotations["phred_quality"]
            reads += 1
            bases += len(trimmed)
            qsum += sum(ph)
            trimmed_recs.append(trimmed)

        avg_q = qsum / bases if bases else 0
        avg_len = bases / reads if reads else 0

        with stats_path.open("w", encoding="utf-8") as fh:
            for r in trimmed_recs:
                ph = r.letter_annotations["phred_quality"]
                fh.write(f"{r.id}\t{len(r)}\t{sum(ph)/len(ph):.2f}\n")

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
            )
            comb.write("\t".join(row) + "\n")

        _tick_safe(bar)


def trim_folder(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    window_size: int = 5,
    per_base_q: int = 20,
    file_q_threshold: float = 20.0,
    mee_max: float | None = None,
    mee_min_len: int | None = None,
    min_reads_kept: int | None = None,
    max_drop_fraction: float | None = None,
    combined_tsv: str | Path | None = None,
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
            comb.write("\t".join(TRIM_SUMMARY_COLUMNS) + "\n")

    fastqs = sorted(input_dir.glob("*.fastq"))

    prog = import_module("microseq_tests.utility.progress")
    cm_or_gen = prog.stage_bar(len(fastqs), desc="trim", unit="file")

    def _run(bar):
        _trim_all(
            fastqs,
            window_size=window_size,
            per_base_q=per_base_q,
            file_q_threshold=file_q_threshold,
            mee_max=mee_max,
            mee_min_len=mee_min_len,
            min_reads_kept=min_reads_kept,
            max_drop_fraction=max_drop_fraction,
            passed_dir=passed_dir,
            failed_dir=failed_dir,
            stats_root=output_dir,
            comb=comb,
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
