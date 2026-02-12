# src/microseq_tests/assembly/two_read_merge.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from microseq_tests.assembly.overlap_utils import (
    AlignedOverlapCandidate,
    OverlapResult,
    iter_end_anchored_overlaps,
    select_best_overlap,
)
from microseq_tests.assembly.overlap_backends import compute_overlap_candidates, resolve_overlap_engine

_IUPAC_PAIR_MAP = {
    frozenset({"A", "G"}): "R",
    frozenset({"C", "T"}): "Y",
    frozenset({"G", "C"}): "S",
    frozenset({"A", "T"}): "W",
    frozenset({"G", "T"}): "K",
    frozenset({"A", "C"}): "M",
}


def _iupac_for_pair(base_a: str, base_b: str) -> str:
    return _IUPAC_PAIR_MAP.get(frozenset({base_a.upper(), base_b.upper()}), "N")


class MergeInputError(ValueError):
    """Raised when merge input files do not contain exactly one sequence record."""

    def __init__(self, message: str, *, f_count: int | None = None, r_count: int | None = None):
        super().__init__(message)
        self.f_count = f_count
        self.r_count = r_count


@dataclass(frozen=True)
class MergeReport:
    sample_id: str
    overlap_engine: str
    orientation: str
    overlap_len: int
    identity: float
    mismatches: int
    contig_len: int
    merge_status: str
    qualities: str
    merge_warning: str
    high_conflict_mismatches: int


def _count_fasta_records(path: Path) -> int:
    return sum(1 for _ in SeqIO.parse(path, "fasta"))


def _read_singleton_fasta(path: Path) -> SeqRecord:
    iterator = SeqIO.parse(path, "fasta")
    first = next(iterator, None)
    if first is None:
        raise MergeInputError(f"No FASTA records found in {path}")
    second = next(iterator, None)
    if second is not None:
        raise MergeInputError(f"Expected singleton FASTA in {path} but found multiple records")
    return first


def _read_qual(path: Path, record_id: str) -> list[int] | None:
    if not path.exists():
        return None
    for rec in SeqIO.parse(path, "qual"):
        if rec.id == record_id:
            return rec.letter_annotations.get("phred_quality")
    return None


def _build_consensus(
    aln_a: str,
    aln_b: str,
    qual_a: Iterable[int] | None,
    qual_b: Iterable[int] | None,
) -> str:
    consensus: list[str] = []
    idx_a = 0
    idx_b = 0
    qual_a_list = list(qual_a) if qual_a is not None else []
    qual_b_list = list(qual_b) if qual_b is not None else []

    for base_a, base_b in zip(aln_a, aln_b):
        if base_a == "-":
            consensus.append(base_b)
            idx_b += 1
            continue
        if base_b == "-":
            consensus.append(base_a)
            idx_a += 1
            continue
        qa = qual_a_list[idx_a] if idx_a < len(qual_a_list) else None
        qb = qual_b_list[idx_b] if idx_b < len(qual_b_list) else None
        if base_a == base_b:
            consensus.append(base_a)
        else:
            if qa is None or qb is None:
                consensus.append("N")
            elif qa > qb:
                consensus.append(base_a)
            elif qb > qa:
                consensus.append(base_b)
            else:
                # unresolved low-confidence tie: preserve ambiguity with IUPAC code
                if qa is not None and qb is not None and qa <= 20 and qb <= 20:
                    consensus.append(_iupac_for_pair(base_a, base_b))
                else:
                    consensus.append("N")
        idx_a += 1
        idx_b += 1
    return "".join(consensus)




def _count_high_conflict_mismatches(
    aln_a: str,
    aln_b: str,
    qual_a: Iterable[int] | None,
    qual_b: Iterable[int] | None,
    *,
    q_threshold: int,
) -> int:
    """Count mismatch sites where both bases are high confidence."""
    if qual_a is None or qual_b is None:
        return 0

    count = 0
    idx_a = 0
    idx_b = 0
    qa_list = list(qual_a)
    qb_list = list(qual_b)
    for base_a, base_b in zip(aln_a, aln_b):
        if base_a == "-":
            idx_b += 1
            continue
        if base_b == "-":
            idx_a += 1
            continue

        qa = qa_list[idx_a] if idx_a < len(qa_list) else None
        qb = qb_list[idx_b] if idx_b < len(qb_list) else None
        if qa is not None and qb is not None and base_a != base_b and min(qa, qb) >= q_threshold:
            count += 1
        idx_a += 1
        idx_b += 1
    return count

def _write_stub_cap_info(path: Path) -> None:
    path.write_text(
        "Number of overlaps saved: 0\nNumber of overlaps removed: 0\n",
        encoding="utf-8",
    )


def merge_two_reads(
    *,
    sample_id: str,
    fwd_path: Path,
    rev_path: Path,
    output_dir: Path,
    min_overlap: int = 100,
    min_identity: float = 0.8,
    min_quality: float = 20.0,
    quality_mode: str = "warning",
    ambiguity_identity_delta: float = 0.0,
    ambiguity_quality_epsilon: float = 0.1,
    high_conflict_q_threshold: int = 30,
    high_conflict_action: str = "flag",
    overlap_engine: str = "ungapped",
    anchor_tolerance_bases: int = 30,
) -> tuple[Path | None, MergeReport]:
    f_count = _count_fasta_records(fwd_path)
    r_count = _count_fasta_records(rev_path)
    if f_count != 1 or r_count != 1:
        raise MergeInputError(
            f"Expected singleton FASTA inputs for {sample_id} (f_n={f_count}, r_n={r_count})",
            f_count=f_count,
            r_count=r_count,
        )

    fwd_record = _read_singleton_fasta(fwd_path)
    rev_record = _read_singleton_fasta(rev_path)

    fwd_qual = _read_qual(Path(f"{fwd_path}.qual"), fwd_record.id)
    rev_qual = _read_qual(Path(f"{rev_path}.qual"), rev_record.id)
    qualities = "present" if fwd_qual is not None and rev_qual is not None else "absent"

    resolved_overlap_engine = resolve_overlap_engine(overlap_engine)
    if resolved_overlap_engine == "ungapped":
        candidates = iter_end_anchored_overlaps(
            str(fwd_record.seq),
            str(rev_record.seq),
            fwd_qual,
            rev_qual,
        )
    else:
        backend_results = compute_overlap_candidates(
            str(fwd_record.seq),
            str(rev_record.seq),
            fwd_qual,
            rev_qual,
            engine=resolved_overlap_engine,
            end_anchor_tolerance=anchor_tolerance_bases,
        )
        candidates = [
            AlignedOverlapCandidate(
                orientation=result.orientation,
                overlap_len=result.overlap_len,
                mismatches=result.mismatches,
                identity=result.identity,
                overlap_quality=result.overlap_quality,
                aligned_fwd=result.aligned_fwd,
                aligned_rev=result.aligned_rev,
                indels=result.indels,
                cigar=result.cigar,
                overlap_span_cols=result.overlap_span_cols,
                terminal_gap_cols=result.terminal_gap_cols,
                internal_indels=result.internal_indels,
                fwd_overlap_start=result.fwd_overlap_start,
                fwd_overlap_end=result.fwd_overlap_end,
                rev_overlap_start=result.rev_overlap_start,
                rev_overlap_end=result.rev_overlap_end,
                end_anchored=result.end_anchored,
            )
            for result in backend_results
        ]
    overlap: OverlapResult = select_best_overlap(
        candidates,
        min_overlap=min_overlap,
        min_identity=min_identity,
        min_quality=min_quality,
        quality_mode=quality_mode,
        ambiguity_identity_delta=ambiguity_identity_delta,
        ambiguity_quality_epsilon=ambiguity_quality_epsilon,
    )
    rev_qual_for_consensus = rev_qual
    if overlap.orientation == "revcomp" and rev_qual is not None:
        rev_qual_for_consensus = list(reversed(rev_qual))

    contig_path = output_dir / f"{sample_id}_paired.fasta.cap.contigs"
    singlets_path = output_dir / f"{sample_id}_paired.fasta.cap.singlets"
    report_path = output_dir / f"{sample_id}_paired.merge_report.tsv"
    cap_info_path = output_dir / f"{sample_id}_paired.fasta.cap.info"

    if overlap.status == "ambiguous_overlap":
        SeqIO.write([fwd_record, rev_record], singlets_path, "fasta")
        report = MergeReport(
            sample_id,
            resolved_overlap_engine,
            overlap.orientation,
            overlap.overlap_len,
            overlap.identity,
            overlap.mismatches,
            0,
            "ambiguous_overlap",
            qualities,
            "",
            0,
        )
        _write_stub_cap_info(cap_info_path)
        _write_merge_report(report_path, report)
        return None, report

    if overlap.overlap_len < min_overlap:
        SeqIO.write([fwd_record, rev_record], singlets_path, "fasta")
        report = MergeReport(
            sample_id,
            resolved_overlap_engine,
            overlap.orientation,
            overlap.overlap_len,
            overlap.identity,
            overlap.mismatches,
            0,
            "overlap_too_short",
            qualities,
            "",
            0,
        )
        _write_stub_cap_info(cap_info_path)
        _write_merge_report(report_path, report)
        return None, report

    if overlap.status == "not_end_anchored":
        SeqIO.write([fwd_record, rev_record], singlets_path, "fasta")
        report = MergeReport(
            sample_id,
            resolved_overlap_engine,
            overlap.orientation,
            overlap.overlap_len,
            overlap.identity,
            overlap.mismatches,
            0,
            "not_end_anchored",
            qualities,
            "",
            0,
        )
        _write_stub_cap_info(cap_info_path)
        _write_merge_report(report_path, report)
        return None, report

    if overlap.identity < min_identity:
        SeqIO.write([fwd_record, rev_record], singlets_path, "fasta")
        report = MergeReport(
            sample_id,
            resolved_overlap_engine,
            overlap.orientation,
            overlap.overlap_len,
            overlap.identity,
            overlap.mismatches,
            0,
            "identity_low",
            qualities,
            "",
            0,
        )
        _write_stub_cap_info(cap_info_path)
        _write_merge_report(report_path, report)
        return None, report

    if (
        quality_mode == "blocking"
        and (overlap.overlap_quality is None or overlap.overlap_quality < min_quality)
    ):
        SeqIO.write([fwd_record, rev_record], singlets_path, "fasta")
        report = MergeReport(
            sample_id,
            resolved_overlap_engine,
            overlap.orientation,
            overlap.overlap_len,
            overlap.identity,
            overlap.mismatches,
            0,
            "quality_low",
            qualities,
            "",
            0,
        )
        _write_stub_cap_info(cap_info_path)
        _write_merge_report(report_path, report)
        return None, report

    high_conflict_mismatches = _count_high_conflict_mismatches(
        overlap.aligned_fwd,
        overlap.aligned_rev,
        fwd_qual,
        rev_qual_for_consensus,
        q_threshold=high_conflict_q_threshold,
    )

    merge_warning_tokens: list[str] = []
    if quality_mode == "warning" and qualities == "absent":
        merge_warning_tokens.append("qualities_absent")
    if (
        quality_mode == "warning"
        and overlap.overlap_quality is not None
        and overlap.overlap_quality < min_quality
    ):
        merge_warning_tokens.append("overlap_quality_low")
    action = high_conflict_action.strip().lower()
    if action not in {"flag", "route_cap3"}:
        action = "flag"

    if high_conflict_mismatches > 0 and action == "route_cap3":
        SeqIO.write([fwd_record, rev_record], singlets_path, "fasta")
        merge_warning_tokens.append("high_conflict")
        report = MergeReport(
            sample_id,
            resolved_overlap_engine,
            overlap.orientation,
            overlap.overlap_len,
            overlap.identity,
            overlap.mismatches,
            0,
            "high_conflict",
            qualities,
            ";".join(merge_warning_tokens),
            high_conflict_mismatches,
        )
        _write_stub_cap_info(cap_info_path)
        _write_merge_report(report_path, report)
        return None, report

    if high_conflict_mismatches > 0:
        merge_warning_tokens.append("high_conflict")

    merge_warning = ";".join(merge_warning_tokens)

    consensus = _build_consensus(
        overlap.aligned_fwd,
        overlap.aligned_rev,
        fwd_qual,
        rev_qual_for_consensus,
    )
    contig_record = SeqRecord(Seq(consensus), id=f"{sample_id}_merged", description="")
    SeqIO.write([contig_record], contig_path, "fasta")

    report = MergeReport(
        sample_id,
        resolved_overlap_engine,
        overlap.orientation,
        overlap.overlap_len,
        overlap.identity,
        overlap.mismatches,
        len(consensus),
        "merged",
        qualities,
        merge_warning,
        high_conflict_mismatches,
    )
    _write_stub_cap_info(cap_info_path)
    _write_merge_report(report_path, report)
    return contig_path, report


def _write_merge_report(path: Path, report: MergeReport) -> None:
    path.write_text(
        "sample_id\toverlap_engine\torientation\toverlap_len\tidentity\tmismatches\tcontig_len\tmerge_status\tqualities\tmerge_warning\thigh_conflict_mismatches\n"
        f"{report.sample_id}\t{report.overlap_engine}\t{report.orientation}\t{report.overlap_len}\t{report.identity:.4f}\t"
        f"{report.mismatches}\t{report.contig_len}\t{report.merge_status}\t{report.qualities}\t{report.merge_warning}\t{report.high_conflict_mismatches}\n",
        encoding="utf-8",
    )
