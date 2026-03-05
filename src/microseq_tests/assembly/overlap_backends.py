from __future__ import annotations

from dataclasses import dataclass
import re

from Bio.Align import PairwiseAligner
from Bio.Seq import Seq


@dataclass(frozen=True)
class OverlapBackendResult:
    """Container for a single backend-produced overlap alignment and derived metrics."""
    orientation: str
    aligned_fwd: str
    aligned_rev: str
    overlap_len: int
    identity: float
    mismatches: int
    indels: int
    overlap_quality: float | None
    cigar: str
    overlap_span_cols: int
    terminal_gap_cols: int
    internal_indels: int
    fwd_overlap_start: int
    fwd_overlap_end: int
    rev_overlap_start: int
    rev_overlap_end: int
    end_anchored: bool


def _pairwise_aligner() -> PairwiseAligner:
    """Create and configure a semi-global PairwiseAligner for end-overlap detection."""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -2.0
    aligner.extend_gap_score = -0.5

    # End-free overlap alignment (semi-global behavior).
    compat_attrs = (
        "target_left_open_gap_score",
        "target_left_extend_gap_score",
        "target_right_open_gap_score",
        "target_right_extend_gap_score",
        "query_left_open_gap_score",
        "query_left_extend_gap_score",
        "query_right_open_gap_score",
        "query_right_extend_gap_score",
    )
    modern_attrs = (
        "open_left_insertion_score",
        "extend_left_insertion_score",
        "open_right_insertion_score",
        "extend_right_insertion_score",
        "open_left_deletion_score",
        "extend_left_deletion_score",
        "open_right_deletion_score",
        "extend_right_deletion_score",
    )
    for attr in (*compat_attrs, *modern_attrs):
        if hasattr(aligner, attr):
            setattr(aligner, attr, 0.0)
    return aligner


def _aligned_strings_from_biopython(alignment, seq_fwd: str, seq_rev: str) -> tuple[str, str]:
    """Reconstruct gapped forward/reverse strings from a Biopython alignment object."""
    t_blocks = alignment.aligned[0]
    q_blocks = alignment.aligned[1]
    f_parts: list[str] = []
    r_parts: list[str] = []
    t_pos = 0
    q_pos = 0
    for (t_start, t_end), (q_start, q_end) in zip(t_blocks, q_blocks):
        if t_start > t_pos:
            segment = seq_fwd[t_pos:t_start]
            f_parts.append(segment)
            r_parts.append("-" * len(segment))
        if q_start > q_pos:
            segment = seq_rev[q_pos:q_start]
            f_parts.append("-" * len(segment))
            r_parts.append(segment)

        f_parts.append(seq_fwd[t_start:t_end])
        r_parts.append(seq_rev[q_start:q_end])
        t_pos = t_end
        q_pos = q_end

    if t_pos < len(seq_fwd):
        segment = seq_fwd[t_pos:]
        f_parts.append(segment)
        r_parts.append("-" * len(segment))
    if q_pos < len(seq_rev):
        segment = seq_rev[q_pos:]
        f_parts.append("-" * len(segment))
        r_parts.append(segment)

    return "".join(f_parts), "".join(r_parts)


def _cigar_from_aligned(aligned_fwd: str, aligned_rev: str) -> str:
    """Build an extended CIGAR string from two already-gapped alignment strings."""
    ops: list[str] = []
    counts: list[int] = []

    for base_f, base_r in zip(aligned_fwd, aligned_rev):
        if base_f == "-":
            op = "I"
        elif base_r == "-":
            op = "D"
        elif base_f == base_r:
            op = "="
        else:
            op = "X"
        if ops and ops[-1] == op:
            counts[-1] += 1
        else:
            ops.append(op)
            counts.append(1)
    return "".join(f"{count}{op}" for count, op in zip(counts, ops))


def _metrics_from_alignment(
    aligned_fwd: str,
    aligned_rev: str,
    fwd_qual: list[int] | None,
    rev_qual: list[int] | None,
    *,
    end_anchor_tolerance: int,
) -> dict[str, int | float | bool | None]:
    """Compute overlap, identity, indel, geometry, and anchoring metrics from aligned strings."""
    both_cols = [idx for idx, (a, b) in enumerate(zip(aligned_fwd, aligned_rev)) if a != "-" and b != "-"]
    if not both_cols:
        return {
            "overlap_len": 0,
            "mismatches": 0,
            "indels": 0,
            "internal_indels": 0,
            "terminal_gap_cols": sum(1 for a, b in zip(aligned_fwd, aligned_rev) if (a == "-") ^ (b == "-")),
            "identity": 0.0,
            "overlap_quality": None,
            "overlap_span_cols": 0,
            "fwd_overlap_start": 0,
            "fwd_overlap_end": 0,
            "rev_overlap_start": 0,
            "rev_overlap_end": 0,
            "end_anchored": False,
        }

    start_col = both_cols[0]
    end_col = both_cols[-1]

    f_idx = 0
    r_idx = 0
    f_start = r_start = None
    f_end = r_end = None
    matches = 0
    mismatches = 0
    internal_indels = 0
    quality_vals: list[int] = []

    for col, (base_f, base_r) in enumerate(zip(aligned_fwd, aligned_rev)):
        if col < start_col or col > end_col:
            if base_f != "-":
                f_idx += 1
            if base_r != "-":
                r_idx += 1
            continue

        f_gap = base_f == "-"
        r_gap = base_r == "-"
        if not f_gap and not r_gap:
            if f_start is None:
                f_start = f_idx
                r_start = r_idx
            f_end = f_idx + 1
            r_end = r_idx + 1
            if base_f == base_r:
                matches += 1
            else:
                mismatches += 1
            if fwd_qual is not None and rev_qual is not None:
                if f_idx < len(fwd_qual) and r_idx < len(rev_qual):
                    quality_vals.append(min(fwd_qual[f_idx], rev_qual[r_idx]))
        else:
            internal_indels += 1

        if not f_gap:
            f_idx += 1
        if not r_gap:
            r_idx += 1

    f_start = 0 if f_start is None else f_start
    r_start = 0 if r_start is None else r_start
    f_end = f_start if f_end is None else f_end
    r_end = r_start if r_end is None else r_end

    terminal_gap_cols = sum(
        1
        for col, (a, b) in enumerate(zip(aligned_fwd, aligned_rev))
        if ((a == "-") ^ (b == "-")) and (col < start_col or col > end_col)
    )

    overlap_len = matches + mismatches
    denom = matches + mismatches + internal_indels
    identity = (matches / denom) if denom else 0.0
    overlap_quality = (sum(quality_vals) / len(quality_vals)) if quality_vals else None

    f_len = len(aligned_fwd.replace("-", ""))
    r_len = len(aligned_rev.replace("-", ""))
    anchored = (
        ((f_len - f_end) <= end_anchor_tolerance and r_start <= end_anchor_tolerance)
        or ((r_len - r_end) <= end_anchor_tolerance and f_start <= end_anchor_tolerance)
    )

    return {
        "overlap_len": overlap_len,
        "mismatches": mismatches,
        "indels": internal_indels + terminal_gap_cols,
        "internal_indels": internal_indels,
        "terminal_gap_cols": terminal_gap_cols,
        "identity": identity,
        "overlap_quality": overlap_quality,
        "overlap_span_cols": end_col - start_col + 1,
        "fwd_overlap_start": f_start,
        "fwd_overlap_end": f_end,
        "rev_overlap_start": r_start,
        "rev_overlap_end": r_end,
        "end_anchored": anchored,
    }


def _biopython_orientation(
    seq_fwd: str,
    seq_rev_oriented: str,
    fwd_qual: list[int] | None,
    rev_qual_oriented: list[int] | None,
    *,
    orientation: str,
    end_anchor_tolerance: int,
) -> OverlapBackendResult:
    """Compute one Biopython-backed overlap result for a single orientation."""
    aligner = _pairwise_aligner()
    alignment = aligner.align(seq_fwd, seq_rev_oriented)[0]
    aligned_fwd, aligned_rev = _aligned_strings_from_biopython(alignment, seq_fwd, seq_rev_oriented)
    m = _metrics_from_alignment(
        aligned_fwd,
        aligned_rev,
        fwd_qual,
        rev_qual_oriented,
        end_anchor_tolerance=end_anchor_tolerance,
    )
    return OverlapBackendResult(
        orientation=orientation,
        aligned_fwd=aligned_fwd,
        aligned_rev=aligned_rev,
        overlap_len=int(m["overlap_len"]),
        identity=float(m["identity"]),
        mismatches=int(m["mismatches"]),
        indels=int(m["indels"]),
        overlap_quality=m["overlap_quality"],
        cigar=_cigar_from_aligned(aligned_fwd, aligned_rev),
        overlap_span_cols=int(m["overlap_span_cols"]),
        terminal_gap_cols=int(m["terminal_gap_cols"]),
        internal_indels=int(m["internal_indels"]),
        fwd_overlap_start=int(m["fwd_overlap_start"]),
        fwd_overlap_end=int(m["fwd_overlap_end"]),
        rev_overlap_start=int(m["rev_overlap_start"]),
        rev_overlap_end=int(m["rev_overlap_end"]),
        end_anchored=bool(m["end_anchored"]),
    )


def _biopython_orientation_candidates(
    seq_fwd: str,
    seq_rev_oriented: str,
    fwd_qual: list[int] | None,
    rev_qual_oriented: list[int] | None,
    *,
    orientation: str,
    end_anchor_tolerance: int,
    max_alignments: int,
) -> list[OverlapBackendResult]:
    """Generate up to max_alignments Biopython candidates for one orientation."""
    aligner = _pairwise_aligner()
    alignments = aligner.align(seq_fwd, seq_rev_oriented)
    out: list[OverlapBackendResult] = []
    for idx, alignment in enumerate(alignments):
        if idx >= max_alignments:
            break
        aligned_fwd, aligned_rev = _aligned_strings_from_biopython(alignment, seq_fwd, seq_rev_oriented)
        m = _metrics_from_alignment(
            aligned_fwd,
            aligned_rev,
            fwd_qual,
            rev_qual_oriented,
            end_anchor_tolerance=end_anchor_tolerance,
        )
        out.append(
            OverlapBackendResult(
                orientation=orientation,
                aligned_fwd=aligned_fwd,
                aligned_rev=aligned_rev,
                overlap_len=int(m["overlap_len"]),
                identity=float(m["identity"]),
                mismatches=int(m["mismatches"]),
                indels=int(m["indels"]),
                overlap_quality=m["overlap_quality"],
                cigar=_cigar_from_aligned(aligned_fwd, aligned_rev),
                overlap_span_cols=int(m["overlap_span_cols"]),
                terminal_gap_cols=int(m["terminal_gap_cols"]),
                internal_indels=int(m["internal_indels"]),
                fwd_overlap_start=int(m["fwd_overlap_start"]),
                fwd_overlap_end=int(m["fwd_overlap_end"]),
                rev_overlap_start=int(m["rev_overlap_start"]),
                rev_overlap_end=int(m["rev_overlap_end"]),
                end_anchored=bool(m["end_anchored"]),
            )
        )
    return out


def compute_biopython_candidates(
    seq_fwd: str,
    seq_rev: str,
    fwd_qual: list[int] | None,
    rev_qual: list[int] | None,
    *,
    end_anchor_tolerance: int,
    max_alignments: int = 5,
) -> list[OverlapBackendResult]:
    """Compute overlap candidates for both forward and reverse-complement orientations with Biopython."""
    revcomp_seq = str(Seq(seq_rev).reverse_complement())
    revcomp_qual = list(reversed(rev_qual)) if rev_qual is not None else None
    return [
        *_biopython_orientation_candidates(
            seq_fwd,
            seq_rev,
            fwd_qual,
            rev_qual,
            orientation="forward",
            end_anchor_tolerance=end_anchor_tolerance,
            max_alignments=max_alignments,
        ),
        *_biopython_orientation_candidates(
            seq_fwd,
            revcomp_seq,
            fwd_qual,
            revcomp_qual,
            orientation="revcomp",
            end_anchor_tolerance=end_anchor_tolerance,
            max_alignments=max_alignments,
        ),
    ]


def _parse_cigar(cigar: str) -> list[tuple[int, str]]:
    """Parse an extended CIGAR string into (count, op) tuples."""
    tokens = re.findall(r"(\d+)([=XIDM])", cigar)
    return [(int(n), op) for n, op in tokens]


def _aligned_strings_from_edlib(
    seq_fwd: str,
    seq_rev_oriented: str,
    *,
    cigar: str,
    target_start: int,
) -> tuple[str, str]:
    """Expand an Edlib CIGAR and location into full-length aligned forward/reverse strings."""
    ops = _parse_cigar(cigar)
    q_idx = 0
    t_idx = target_start
    overlap_f: list[str] = []
    overlap_r: list[str] = []
    for count, op in ops:
        for _ in range(count):
            if op in {"=", "X", "M"}:
                overlap_f.append(seq_fwd[q_idx])
                overlap_r.append(seq_rev_oriented[t_idx])
                q_idx += 1
                t_idx += 1
            elif op == "I":
                overlap_f.append(seq_fwd[q_idx])
                overlap_r.append("-")
                q_idx += 1
            elif op == "D":
                overlap_f.append("-")
                overlap_r.append(seq_rev_oriented[t_idx])
                t_idx += 1

    q_end = q_idx
    t_end = t_idx

    parts_f: list[str] = []
    parts_r: list[str] = []

    if target_start > 0:
        parts_f.append("-" * target_start)
        parts_r.append(seq_rev_oriented[:target_start])

    parts_f.append("".join(overlap_f))
    parts_r.append("".join(overlap_r))

    if q_end < len(seq_fwd):
        parts_f.append(seq_fwd[q_end:])
        parts_r.append("-" * (len(seq_fwd) - q_end))
    if t_end < len(seq_rev_oriented):
        parts_f.append("-" * (len(seq_rev_oriented) - t_end))
        parts_r.append(seq_rev_oriented[t_end:])

    return "".join(parts_f), "".join(parts_r)


def _edlib_orientation(
    seq_fwd: str,
    seq_rev_oriented: str,
    fwd_qual: list[int] | None,
    rev_qual_oriented: list[int] | None,
    *,
    orientation: str,
    end_anchor_tolerance: int,
) -> OverlapBackendResult:
    """Compute one Edlib-backed overlap result for a single orientation."""
    import edlib  # type: ignore

    aln = edlib.align(seq_fwd, seq_rev_oriented, mode="HW", task="path")
    loc = aln.get("locations") or []
    target_start = int(loc[0][0]) if loc else 0
    aligned_fwd, aligned_rev = _aligned_strings_from_edlib(
        seq_fwd,
        seq_rev_oriented,
        cigar=aln.get("cigar", ""),
        target_start=target_start,
    )
    m = _metrics_from_alignment(
        aligned_fwd,
        aligned_rev,
        fwd_qual,
        rev_qual_oriented,
        end_anchor_tolerance=end_anchor_tolerance,
    )
    return OverlapBackendResult(
        orientation=orientation,
        aligned_fwd=aligned_fwd,
        aligned_rev=aligned_rev,
        overlap_len=int(m["overlap_len"]),
        identity=float(m["identity"]),
        mismatches=int(m["mismatches"]),
        indels=int(m["indels"]),
        overlap_quality=m["overlap_quality"],
        cigar=aln.get("cigar", ""),
        overlap_span_cols=int(m["overlap_span_cols"]),
        terminal_gap_cols=int(m["terminal_gap_cols"]),
        internal_indels=int(m["internal_indels"]),
        fwd_overlap_start=int(m["fwd_overlap_start"]),
        fwd_overlap_end=int(m["fwd_overlap_end"]),
        rev_overlap_start=int(m["rev_overlap_start"]),
        rev_overlap_end=int(m["rev_overlap_end"]),
        end_anchored=bool(m["end_anchored"]),
    )


def _edlib_orientation_candidates(
    seq_fwd: str,
    seq_rev_oriented: str,
    fwd_qual: list[int] | None,
    rev_qual_oriented: list[int] | None,
    *,
    orientation: str,
    end_anchor_tolerance: int,
) -> list[OverlapBackendResult]:
    """Generate path-consistent Edlib overlap candidates for one orientation."""
    import edlib  # type: ignore

    aln = edlib.align(seq_fwd, seq_rev_oriented, mode="HW", task="path")
    # Edlib returns a single path/cigar corresponding to the first location pair.
    # Keep candidate generation path-consistent by using only locations[0] here.
    locs = aln.get("locations") or [(0, max(0, len(seq_rev_oriented) - 1))]
    locs = [locs[0]]
    out: list[OverlapBackendResult] = []
    for loc in locs:
        target_start = int(loc[0])
        aligned_fwd, aligned_rev = _aligned_strings_from_edlib(
            seq_fwd,
            seq_rev_oriented,
            cigar=aln.get("cigar", ""),
            target_start=target_start,
        )
        m = _metrics_from_alignment(
            aligned_fwd,
            aligned_rev,
            fwd_qual,
            rev_qual_oriented,
            end_anchor_tolerance=end_anchor_tolerance,
        )
        out.append(
            OverlapBackendResult(
                orientation=orientation,
                aligned_fwd=aligned_fwd,
                aligned_rev=aligned_rev,
                overlap_len=int(m["overlap_len"]),
                identity=float(m["identity"]),
                mismatches=int(m["mismatches"]),
                indels=int(m["indels"]),
                overlap_quality=m["overlap_quality"],
                cigar=aln.get("cigar", ""),
                overlap_span_cols=int(m["overlap_span_cols"]),
                terminal_gap_cols=int(m["terminal_gap_cols"]),
                internal_indels=int(m["internal_indels"]),
                fwd_overlap_start=int(m["fwd_overlap_start"]),
                fwd_overlap_end=int(m["fwd_overlap_end"]),
                rev_overlap_start=int(m["rev_overlap_start"]),
                rev_overlap_end=int(m["rev_overlap_end"]),
                end_anchored=bool(m["end_anchored"]),
            )
        )
    return out


def compute_edlib_candidates(
    seq_fwd: str,
    seq_rev: str,
    fwd_qual: list[int] | None,
    rev_qual: list[int] | None,
    *,
    end_anchor_tolerance: int,
) -> list[OverlapBackendResult]:
    """Compute overlap candidates for both orientations using Edlib."""
    revcomp_seq = str(Seq(seq_rev).reverse_complement())
    revcomp_qual = list(reversed(rev_qual)) if rev_qual is not None else None
    return [
        *_edlib_orientation_candidates(
            seq_fwd,
            seq_rev,
            fwd_qual,
            rev_qual,
            orientation="forward",
            end_anchor_tolerance=end_anchor_tolerance,
        ),
        *_edlib_orientation_candidates(
            seq_fwd,
            revcomp_seq,
            fwd_qual,
            revcomp_qual,
            orientation="revcomp",
            end_anchor_tolerance=end_anchor_tolerance,
        ),
    ]


def compute_overlap_candidates(
    seq_fwd: str,
    seq_rev: str,
    fwd_qual: list[int] | None,
    rev_qual: list[int] | None,
    *,
    engine: str,
    end_anchor_tolerance: int = 30,
) -> list[OverlapBackendResult]:
    """Dispatch overlap candidate generation to the configured backend engine."""
    mode = resolve_overlap_engine(engine)
    if mode == "biopython":
        return compute_biopython_candidates(
            seq_fwd,
            seq_rev,
            fwd_qual,
            rev_qual,
            end_anchor_tolerance=end_anchor_tolerance,
        )
    if mode == "edlib":
        return compute_edlib_candidates(
            seq_fwd,
            seq_rev,
            fwd_qual,
            rev_qual,
            end_anchor_tolerance=end_anchor_tolerance,
        )
    raise ValueError(f"Unsupported overlap engine: {engine}")


def resolve_overlap_engine(engine: str) -> str:
    """Resolve the requested overlap backend mode, including auto-detection behavior."""
    mode = engine.strip().lower()
    if mode == "auto":
        try:
            import edlib as _edlib  # type: ignore # noqa: F401
        except Exception:
            mode = "biopython"
        else:
            mode = "edlib"
    return mode


def resolve_overlap_engine_strategy(strategy: str) -> str:
    mode = str(strategy).strip().lower()
    if mode not in {"single", "cascade", "all"}:
        return "single"
    return mode


def resolve_overlap_engine_order(order: list[str] | tuple[str, ...] | None) -> list[str]:
    if not order:
        return ["ungapped", "biopython", "edlib"]
    resolved: list[str] = []
    for raw in order:
        mode = resolve_overlap_engine(str(raw))
        if mode not in {"ungapped", "biopython", "edlib"}:
            raise ValueError(f"Unsupported overlap engine in order: {raw}")
        if mode not in resolved:
            resolved.append(mode)
    if not resolved:
        return ["ungapped", "biopython", "edlib"]
    return resolved
