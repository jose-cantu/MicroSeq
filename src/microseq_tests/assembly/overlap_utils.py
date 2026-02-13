# src/microseq_tests/assembly/overlap_utils.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

from Bio.Seq import Seq


@dataclass(frozen=True)
class OverlapResult:
    orientation: str
    overlap_len: int
    identity: float
    mismatches: int
    overlap_quality: float | None
    aligned_fwd: str
    aligned_rev: str
    status: str = "ok"


@dataclass(frozen=True)
class OverlapCandidate:
    orientation: str
    seq_fwd: str
    seq_rev: str
    offset: int
    overlap_len: int
    mismatches: int
    identity: float
    overlap_quality: float | None

    def as_result(self) -> OverlapResult:
        aligned_fwd, aligned_rev = _build_aligned_strings(self.seq_fwd, self.seq_rev, self.offset)
        return OverlapResult(
            orientation=self.orientation,
            overlap_len=self.overlap_len,
            identity=self.identity,
            mismatches=self.mismatches,
            overlap_quality=self.overlap_quality,
            aligned_fwd=aligned_fwd,
            aligned_rev=aligned_rev,
        )


@dataclass(frozen=True)
class AlignedOverlapCandidate:
    orientation: str
    overlap_len: int
    mismatches: int
    identity: float
    overlap_quality: float | None
    aligned_fwd: str
    aligned_rev: str
    indels: int = 0
    cigar: str = ""
    overlap_span_cols: int = 0
    terminal_gap_cols: int = 0
    internal_indels: int = 0
    fwd_overlap_start: int = 0
    fwd_overlap_end: int = 0
    rev_overlap_start: int = 0
    rev_overlap_end: int = 0
    end_anchored: bool = True

    def as_result(self) -> OverlapResult:
        return OverlapResult(
            orientation=self.orientation,
            overlap_len=self.overlap_len,
            identity=self.identity,
            mismatches=self.mismatches,
            overlap_quality=self.overlap_quality,
            aligned_fwd=self.aligned_fwd,
            aligned_rev=self.aligned_rev,
        )


def _build_aligned_strings(seq_fwd: str, seq_rev: str, offset: int) -> tuple[str, str]:
    if offset >= 0:
        left_pad_f = 0
        left_pad_r = offset
    else:
        left_pad_f = -offset
        left_pad_r = 0

    len_total = max(left_pad_f + len(seq_fwd), left_pad_r + len(seq_rev))

    f_chars = ["-"] * len_total
    r_chars = ["-"] * len_total
    for idx, base in enumerate(seq_fwd):
        f_chars[left_pad_f + idx] = base
    for idx, base in enumerate(seq_rev):
        r_chars[left_pad_r + idx] = base
    return "".join(f_chars), "".join(r_chars)


def _compute_overlap_metrics_from_offset(
    seq_fwd: str,
    seq_rev: str,
    offset: int,
    fwd_qual: list[int] | None,
    rev_qual: list[int] | None,
) -> tuple[int, int, float, float | None]:
    f_start = max(0, offset)
    r_start = max(0, -offset)
    overlap_len = min(len(seq_fwd) - f_start, len(seq_rev) - r_start)
    if overlap_len <= 0:
        return 0, 0, 0.0, None

    f_slice = seq_fwd[f_start : f_start + overlap_len]
    r_slice = seq_rev[r_start : r_start + overlap_len]

    matches = 0
    mismatches = 0
    quality_vals: list[int] = []
    for idx, (base_f, base_r) in enumerate(zip(f_slice, r_slice)):
        if base_f == base_r:
            matches += 1
        else:
            mismatches += 1

        if fwd_qual is not None and rev_qual is not None:
            qf_idx = f_start + idx
            qr_idx = r_start + idx
            if qf_idx < len(fwd_qual) and qr_idx < len(rev_qual):
                quality_vals.append(min(fwd_qual[qf_idx], rev_qual[qr_idx]))

    identity = matches / overlap_len if overlap_len else 0.0
    overlap_quality = None
    if quality_vals:
        overlap_quality = sum(quality_vals) / len(quality_vals)
    return overlap_len, mismatches, identity, overlap_quality


def _end_anchored_candidates_for_orientation(
    seq_fwd: str,
    seq_rev: str,
    fwd_qual: list[int] | None,
    rev_qual: list[int] | None,
    *,
    orientation: str,
) -> list[OverlapCandidate]:
    candidates: list[OverlapCandidate] = []
    min_offset = -len(seq_rev) + 1
    max_offset = len(seq_fwd) - 1

    for offset in range(min_offset, max_offset + 1):
        f_start = max(0, offset)
        r_start = max(0, -offset)
        overlap_len = min(len(seq_fwd) - f_start, len(seq_rev) - r_start)
        if overlap_len <= 0:
            continue

        uses_f_suffix = (f_start + overlap_len) == len(seq_fwd)
        uses_r_suffix = (r_start + overlap_len) == len(seq_rev)
        uses_f_prefix = f_start == 0
        uses_r_prefix = r_start == 0
        if not ((uses_f_suffix and uses_r_prefix) or (uses_r_suffix and uses_f_prefix)):
            continue

        ov_len, mismatches, identity, ov_quality = _compute_overlap_metrics_from_offset(
            seq_fwd,
            seq_rev,
            offset,
            fwd_qual,
            rev_qual,
        )
        candidates.append(
            OverlapCandidate(
                orientation=orientation,
                seq_fwd=seq_fwd,
                seq_rev=seq_rev,
                offset=offset,
                overlap_len=ov_len,
                mismatches=mismatches,
                identity=identity,
                overlap_quality=ov_quality,
            )
        )
    return candidates


def iter_end_anchored_overlaps(
    seq_fwd: str,
    seq_rev: str,
    fwd_qual: list[int] | None = None,
    rev_qual: list[int] | None = None,
) -> list[OverlapCandidate]:
    forward = _end_anchored_candidates_for_orientation(
        seq_fwd,
        seq_rev,
        fwd_qual,
        rev_qual,
        orientation="forward",
    )

    revcomp_seq = str(Seq(seq_rev).reverse_complement())
    revcomp_qual = list(reversed(rev_qual)) if rev_qual is not None else None
    revcomp = _end_anchored_candidates_for_orientation(
        seq_fwd,
        revcomp_seq,
        fwd_qual,
        revcomp_qual,
        orientation="revcomp",
    )
    return [*forward, *revcomp]


def _quality_sort_value(overlap_quality: float | None) -> tuple[int, float]:
    if overlap_quality is None:
        return (0, float("-inf"))
    return (1, overlap_quality)


def _candidate_sort_key(candidate: OverlapCandidate | AlignedOverlapCandidate) -> tuple[int, int, tuple[int, float], float]:
    """Sort overlaps by longest overlap, then fewest mismatches, then quality, then identity."""
    return (
        candidate.overlap_len,
        -candidate.mismatches,
        _quality_sort_value(candidate.overlap_quality),
        candidate.identity,
    )


def _pick_best_candidate(candidates: Sequence[OverlapCandidate | AlignedOverlapCandidate]) -> OverlapCandidate | AlignedOverlapCandidate | None:
    if not candidates:
        return None
    return max(candidates, key=_candidate_sort_key)


def pick_best_identity_candidate(
    candidates: Sequence[OverlapCandidate | AlignedOverlapCandidate],
    *,
    min_overlap: int = 0,
) -> OverlapCandidate | AlignedOverlapCandidate | None:
    """Return the candidate with highest identity (ties -> longer overlap, fewer mismatches)."""
    filtered = [c for c in candidates if c.overlap_len >= min_overlap]
    if not filtered:
        return None
    return max(filtered, key=lambda c: (c.identity, c.overlap_len, -c.mismatches, _quality_sort_value(c.overlap_quality)))


def _is_ambiguous_top_pair(
    first: OverlapCandidate | AlignedOverlapCandidate,
    second: OverlapCandidate | AlignedOverlapCandidate,
    *,
    identity_delta: float,
    quality_epsilon: float,
) -> bool:
    """Return True when top candidates are not uniquely best by overlap metrics."""
    if first.overlap_len != second.overlap_len:
        return False
    if first.mismatches != second.mismatches:
        return False

    first_q = first.overlap_quality
    second_q = second.overlap_quality
    if first_q is not None or second_q is not None:
        if first_q is None or second_q is None:
            return False
        return abs(first_q - second_q) <= quality_epsilon

    return abs(first.identity - second.identity) <= identity_delta


def _is_feasible_candidate(
    candidate: OverlapCandidate | AlignedOverlapCandidate,
    *,
    min_overlap: int,
    min_identity: float,
    min_quality: float,
    quality_mode: str,
) -> bool:
    if hasattr(candidate, "end_anchored") and not getattr(candidate, "end_anchored"):
        return False
    if candidate.overlap_len < min_overlap:
        return False
    if candidate.identity < min_identity:
        return False
    if quality_mode == "blocking":
        if candidate.overlap_quality is None:
            return False
        if candidate.overlap_quality < min_quality:
            return False
    return True


def select_best_overlap(
    candidates: Sequence[OverlapCandidate | AlignedOverlapCandidate],
    *,
    min_overlap: int,
    min_identity: float,
    min_quality: float,
    quality_mode: str = "warning",
    ambiguity_identity_delta: float = 0.0,
    ambiguity_quality_epsilon: float = 0.1,
) -> OverlapResult:
    candidates = list(candidates)
    feasible = [
        c
        for c in candidates
        if _is_feasible_candidate(
            c,
            min_overlap=min_overlap,
            min_identity=min_identity,
            min_quality=min_quality,
            quality_mode=quality_mode,
        )
    ]

    if feasible:
        ranked = sorted(feasible, key=_candidate_sort_key, reverse=True)
        chosen = ranked[0]
        if len(ranked) > 1 and _is_ambiguous_top_pair(
            chosen,
            ranked[1],
            identity_delta=ambiguity_identity_delta,
            quality_epsilon=ambiguity_quality_epsilon,
        ):
            result = chosen.as_result()
            return OverlapResult(
                result.orientation,
                result.overlap_len,
                result.identity,
                result.mismatches,
                result.overlap_quality,
                result.aligned_fwd,
                result.aligned_rev,
                "ambiguous_overlap",
            )
        return chosen.as_result()

    best_overall = _pick_best_candidate(candidates)
    if best_overall is None:
        return OverlapResult("forward", 0, 0.0, 0, None, "", "")

    anchored_candidates = [c for c in candidates if not hasattr(c, "end_anchored") or getattr(c, "end_anchored")]
    if anchored_candidates:
        return _pick_best_candidate(anchored_candidates).as_result()

    if hasattr(best_overall, "end_anchored") and not getattr(best_overall, "end_anchored"):
        result = best_overall.as_result()
        return OverlapResult(
            result.orientation,
            result.overlap_len,
            result.identity,
            result.mismatches,
            result.overlap_quality,
            result.aligned_fwd,
            result.aligned_rev,
            "not_end_anchored",
        )
    return best_overall.as_result()


def best_pairwise_overlap(
    seq_fwd: str,
    seq_rev: str,
    fwd_qual: list[int] | None = None,
    rev_qual: list[int] | None = None,
) -> OverlapResult:
    candidates = iter_end_anchored_overlaps(seq_fwd, seq_rev, fwd_qual, rev_qual)
    chosen = _pick_best_candidate(candidates)
    if chosen is None:
        return OverlapResult("forward", 0, 0.0, 0, None, "", "")
    return chosen.as_result()
