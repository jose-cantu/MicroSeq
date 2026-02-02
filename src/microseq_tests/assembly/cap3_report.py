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


def parse_cap3_reports(
    asm_dir: Path,
    sample_keys: Iterable[str],
    *,
    output_tsv: Path | None = None,
    missing_samples: Iterable[str] | None = None,
    overlap_status: Mapping[str, str] | None = None,
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
                }
            )
            continue
        sample_dir = asm_dir / sample_key
        info_path = sample_dir / f"{sample_key}_paired.fasta.cap.info"
        contigs_path = sample_dir / f"{sample_key}_paired.fasta.cap.contigs"
        singlets_path = sample_dir / f"{sample_key}_paired.fasta.cap.singlets"

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
        else:
            L.warning("Missing CAP3 info file for %s", sample_key)

        contigs_count = _count_fasta_records(contigs_path)
        singlets_count = _count_fasta_records(singlets_path)

        if contigs_count > 0:
            status = "assembled"
        elif singlets_count > 0:
            status = "singlets_only"
        else:
            status = "cap3_no_output"

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
        ]
        with output_tsv.open("w", encoding="utf-8") as fh:
            fh.write("\t".join(headers) + "\n")
            for row in rows:
                fh.write("\t".join("" if row[h] is None else str(row[h]) for h in headers) + "\n")

    return rows
