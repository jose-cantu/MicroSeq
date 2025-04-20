from __future__ import annotations
from pathlib import Path
from typing import List, Optional
from Bio import SeqIO                         
import logging, shutil

# ----------------------------------------------------------------------
# Single‑read sliding‑window trim
# ----------------------------------------------------------------------
def dynamic_trim(record, win: int = 5, q: int = 20):
    quals: List[int] = record.letter_annotations["phred_quality"]
    n = len(quals)
    if n < win:
        return None                            # too short

    # 5′ → 3′
    left = next(
        (i for i in range(n - win + 1)
         if sum(quals[i:i + win]) / win >= q),
        None,
    )
    if left is None:
        return None

    # 3′ → 5′
    right = next(
        (j for j in range(n, win - 1, -1)
         if sum(quals[j - win:j]) / win >= q),
        None,
    )
    if right is None or right <= left:
        return None

    return record[left:right]                 # trimmed SeqRecord


# ----------------------------------------------------------------------
# Folder‑level driver
# ----------------------------------------------------------------------
def trim_folder(
    input_dir: str | Path,
    output_dir: str | Path,
    *,
    window_size: int = 5,
    per_base_q: int = 20,
    file_q_threshold: float = 20.0,
    combined_tsv: str | Path | None = None,
) -> None:
    """
    For every *.fastq in input_dir ....

    • keep reads whose sliding‑window core avg Q ≥ per_base_q
    • write per‑read stats to *_avg_qual.txt
    • if file avg Q ≥ file_q_threshold  → trimmed FASTQ in passed_qc_fastq/
      else                                  move FASTQ to failed_qc_fastq/
    • append one summary row to *combined_tsv* (if provided)
    """
    input_dir  = Path(input_dir)
    output_dir = Path(output_dir)
    passed_dir = output_dir.parent / "passed_qc_fastq"
    failed_dir = output_dir.parent / "failed_qc_fastq"
    for p in (output_dir, passed_dir, failed_dir):
        p.mkdir(parents=True, exist_ok=True)

    comb: Optional[open] = None
    if combined_tsv:
        comb = open(combined_tsv, "a")
        if comb.tell() == 0:
            comb.write("file\treads\tavg_len\tavg_q\n")

    # ------------------------------------------------------------------
    for fq in input_dir.glob("*.fastq"):
        base = fq.stem
        stats_path   = output_dir / f"{base}_avg_qual.txt"
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
            qsum  += sum(ph)
            trimmed_recs.append(trimmed)

        # avg per file on kept bases
        avg_q   = qsum  / bases if bases else 0
        avg_len = bases / reads if reads else 0

        # -- save per‑read stats --------------------------------------
        with open(stats_path, "w") as fh:
            for r in trimmed_recs:
                ph = r.letter_annotations["phred_quality"]
                fh.write(f"{r.id}\t{len(r)}\t{sum(ph)/len(ph):.2f}\n")

        # -- pass / fail ----------------------------------------------
        if avg_q < file_q_threshold:
            (failed_dir / fq.name).write_bytes(fq.read_bytes())
            (failed_dir / stats_path.name).write_bytes(stats_path.read_bytes())
            logging.info("[FAIL] %s  (avgQ %.2f)", fq.name, avg_q)
        else:
            SeqIO.write(trimmed_recs, trimmed_path, "fastq")
            logging.info("[PASS] %s → %s (avgQ %.2f)", fq.name, trimmed_path, avg_q)
            if comb:
                comb.write(f"{fq.name}\t{reads}\t{avg_len:.1f}\t{avg_q:.2f}\n")

    if comb:
        comb.close()

