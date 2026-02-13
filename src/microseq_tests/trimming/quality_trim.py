from __future__ import annotations

import argparse
import logging
import shutil
import subprocess
from pathlib import Path
from typing import Iterable

from Bio import SeqIO

from microseq_tests.trimming.biopy_trim import TRIM_SUMMARY_COLUMNS, build_trim_summary_row
from microseq_tests.utility.utils import load_config, setup_logging

L = logging.getLogger(__name__)

PathLike = str | Path


def quality_trim(input_file: PathLike, output_file: PathLike, *, threads: int = 1, **kwargs) -> Path:
    cfg = load_config()
    trimm = cfg["tools"]["trimmomatic"]
    trim_cfg = cfg.get("trim", {})
    trimm_cfg = trim_cfg.get("trimmomatic", {})
    window_size = int(kwargs.get("window_size", trimm_cfg.get("sliding_window_size", 5)))
    window_q = int(kwargs.get("window_q", trimm_cfg.get("sliding_window_q", 20)))
    min_len = int(kwargs.get("min_len", trim_cfg.get("min_len", 200)))
    phred = int(kwargs.get("phred", trim_cfg.get("phred", 33)))

    cmd = [
        trimm,
        "SE",
        "-threads",
        str(threads),
        f"-phred{phred}",
        str(input_file),
        str(output_file),
        f"SLIDINGWINDOW:{window_size}:{window_q}",
        f"MINLEN:{min_len}",
    ]

    L.info("RUN Trimmomatic: %s", " ".join(cmd))
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if result.stdout:
            L.info("Trimmomatic stdout:\n%s", result.stdout)
        if result.stderr:
            L.warning("Trimmomatic stderr:\n%s", result.stderr)
    except subprocess.CalledProcessError as exc:
        L.error("Trimmomatic failed (exit %s):\n%s", exc.returncode, exc.stderr)
        if exc.stdout:
            L.error("Trimmomatic stdout:\n%s", exc.stdout)
        raise

    out_path = Path(output_file).resolve()
    if not out_path.exists():
        raise FileNotFoundError(out_path)
    return out_path


def _looks_like_fastq(fp: Path) -> bool:
    if not fp.is_file():
        return False
    suffixes = {s.lower() for s in fp.suffixes}
    return ".fastq" in suffixes or ".fq" in suffixes


def _fastq_stats(path: Path) -> tuple[int, float, float]:
    reads = bases = qsum = 0
    for rec in SeqIO.parse(path, "fastq"):
        ph = rec.letter_annotations["phred_quality"]
        reads += 1
        bases += len(rec)
        qsum += sum(ph)
    avg_q = qsum / bases if bases else 0
    avg_len = bases / reads if reads else 0
    return reads, avg_len, avg_q


def _iter_fastq_sources(path: Path) -> Iterable[Path]:
    if path.is_dir():
        for fp in sorted(path.iterdir()):
            if _looks_like_fastq(fp):
                yield fp
    elif path.exists():
        yield path


def trim_fastq_inputs(input_path: PathLike, trim_dir: PathLike, *, summary_tsv: PathLike | None = None) -> Path:
    cfg = load_config()
    trim_cfg = cfg.get("trim", {})
    trimm_cfg = trim_cfg.get("trimmomatic", {})
    policy_desc = (
        f"trimmomatic:window={int(trimm_cfg.get('sliding_window_size', 5))},"
        f"q={int(trimm_cfg.get('sliding_window_q', 20))},minlen={int(trim_cfg.get('min_len', 200))}"
    )

    input_path = Path(input_path)
    trim_dir = Path(trim_dir)
    trim_dir.mkdir(parents=True, exist_ok=True)

    out_fq = trim_dir / "trimmed.fastq"
    sources = list(_iter_fastq_sources(input_path))
    if not sources:
        raise FileNotFoundError(f"No FASTQ files found in {input_path}")

    if out_fq.exists():
        out_fq.unlink()

    summary_rows: list[tuple[str, int, float, float, str]] = []

    for src in sources:
        tmp_out = trim_dir / f"{src.name}.trimmed.fastq"
        trimmed = quality_trim(src, tmp_out)

        reads, avg_len, avg_q = _fastq_stats(trimmed)
        if summary_tsv:
            summary_rows.append((src.name, reads, avg_len, avg_q, "pass"))

        with out_fq.open("ab") as combined, trimmed.open("rb") as fh:
            shutil.copyfileobj(fh, combined)

        trimmed.unlink(missing_ok=True)

    if summary_tsv and summary_rows:
        summary_fp = Path(summary_tsv)
        summary_fp.parent.mkdir(parents=True, exist_ok=True)

        write_header = not summary_fp.exists()

        total_reads = sum(reads for _, reads, _, _, _ in summary_rows)
        total_bases = sum(reads * avg_len for _, reads, avg_len, _, _ in summary_rows)
        total_qsum = sum(avg_q * reads * avg_len for _, reads, avg_len, avg_q, _ in summary_rows)

        combined_avg_len = total_bases / total_reads if total_reads else 0
        combined_avg_q = total_qsum / total_bases if total_bases else 0

        with summary_fp.open("a", encoding="utf-8") as comb:
            if write_header:
                comb.write("\t".join(TRIM_SUMMARY_COLUMNS) + "\n")
            for name, reads, avg_len, avg_q, qc_status in summary_rows:
                row_fields = build_trim_summary_row(
                    file_name=name,
                    reads=reads,
                    avg_len=avg_len,
                    avg_q=avg_q,
                    qc_status=qc_status,
                    trim_policy=policy_desc,
                )
                comb.write("\t".join(row_fields) + "\n")
            combined_row = build_trim_summary_row(
                file_name="_combined",
                reads=total_reads,
                avg_len=combined_avg_len,
                avg_q=combined_avg_q,
                qc_status="combined",
                trim_policy=policy_desc,
            )
            comb.write("\t".join(combined_row) + "\n")
    return out_fq


def main():
    setup_logging()
    parser = argparse.ArgumentParser(description="Quality trim using Trimmomatic")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ")
    parser.add_argument("-o", "--output", required=True, help="Output FASTQ")
    args = parser.parse_args()

    quality_trim(args.input, args.output)


if __name__ == "__main__":
    main()
