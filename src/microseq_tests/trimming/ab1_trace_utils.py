from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True)
class Ab1TraceBundle:
    """Normalized AB1 trace primitives shared by GUI and QC."""

    traces_by_base: dict[str, list[int]]
    channels: list[list[int]]
    data_tags: tuple[str, str, str, str] | None
    channel_bases: tuple[str, str, str, str] | None
    channel_map: dict[str, str]
    basecalls: str
    peak_positions: list[int]  # viewer-safe (clamped)
    qc_peak_positions: list[int]  # QC-safe (filtered in-range)
    qc_pairs: list[tuple[int, int]]  # (call_index, trace_index), QC-safe
    ploc_total: int
    qualities: list[int]
    trace_len: int
    ploc_tag: str | None


def decode_basecalls(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, (bytes, bytearray)):
        return bytes(value).decode("ascii", errors="ignore").replace("\x00", "").strip()
    if isinstance(value, str):
        return value.replace("\x00", "").strip()
    if isinstance(value, Iterable):
        try:
            return "".join(chr(int(v)) for v in value if int(v) > 0)
        except (TypeError, ValueError):
            return ""
    return str(value).replace("\x00", "").strip()


def coerce_ints(value: object) -> list[int]:
    if value is None:
        return []
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, (bytes, bytearray)):
        return [int(b) for b in value]
    if isinstance(value, (list, tuple)):
        out: list[int] = []
        for v in value:
            try:
                out.append(int(v))
            except (TypeError, ValueError):
                continue
        return out
    try:
        return [int(value)]
    except (TypeError, ValueError):
        return []


def _normalize_peak_positions_base(peak_positions: list[int], trace_len: int) -> list[int]:
    if not peak_positions or trace_len <= 0:
        return []

    def _in_range(vals: list[int], *, one_based: bool) -> int:
        lo, hi = (1, trace_len) if one_based else (0, trace_len - 1)
        return sum(lo <= v <= hi for v in vals)

    zero_score = _in_range(peak_positions, one_based=False)
    one_score = _in_range(peak_positions, one_based=True)

    use_one_based = one_score > zero_score or (
        one_score == zero_score and min(peak_positions) >= 1 and max(peak_positions) >= trace_len
    )
    return [p - 1 for p in peak_positions] if use_one_based else list(peak_positions)


def normalize_peak_positions(peak_positions: list[int], trace_len: int) -> list[int]:
    """Normalize 0/1-based PLOC coordinates and clamp to valid trace bounds (viewer-safe)."""
    idx = _normalize_peak_positions_base(peak_positions, trace_len)
    return [max(0, min(trace_len - 1, i)) for i in idx] if idx and trace_len > 0 else []


def normalize_peak_positions_qc(peak_positions: list[int], trace_len: int) -> list[int]:
    """Normalize 0/1-based PLOC coordinates and drop out-of-range values (QC-safe)."""
    idx = _normalize_peak_positions_base(peak_positions, trace_len)
    return [i for i in idx if 0 <= i < trace_len] if idx and trace_len > 0 else []


def normalize_peak_pairs_qc(peak_positions: list[int], trace_len: int) -> list[tuple[int, int]]:
    """Return (call_index, trace_index) for in-range QC positions after normalization."""
    idx = _normalize_peak_positions_base(peak_positions, trace_len)
    if not idx or trace_len <= 0:
        return []
    return [(call_i, pos) for call_i, pos in enumerate(idx) if 0 <= pos < trace_len]


def trim_called_arrays(basecalls: str, peak_positions: list[int], qualities: list[int]) -> tuple[str, list[int], list[int]]:
    """Trim PBAS/PLOC together; qualities are optional and trimmed independently."""
    n_bp = min(len(basecalls), len(peak_positions))
    if n_bp <= 0:
        return "", [], []

    b = basecalls[:n_bp]
    p = peak_positions[:n_bp]
    q = qualities[: min(n_bp, len(qualities))] if qualities else []
    return b, p, q


def extract_ab1_trace_bundle(raw: dict) -> Ab1TraceBundle:
    channel_sets = [
        ("DATA9", "DATA10", "DATA11", "DATA12"),
        ("DATA1", "DATA2", "DATA3", "DATA4"),
    ]
    tags = next((ks for ks in channel_sets if all(k in raw for k in ks)), None)

    traces: dict[str, list[int]] = {"A": [], "C": [], "G": [], "T": []}
    channels: list[list[int]] = []
    channel_bases: tuple[str, str, str, str] | None = None
    channel_map: dict[str, str] = {}

    if tags:
        order_raw = raw.get("FWO_1", "")
        if isinstance(order_raw, (bytes, bytearray)):
            order = order_raw.decode("ascii", errors="ignore")
        else:
            order = str(order_raw or "")
        order = order.strip().upper()
        if len(order) == 4 and set(order).issubset({"A", "C", "G", "T"}):
            base_by_idx = list(order)
        else:
            base_by_idx = ["G", "A", "T", "C"]
        channel_bases = tuple(base_by_idx)  # type: ignore[assignment]

        for i, tag in enumerate(tags):
            vals = [int(v) for v in raw[tag]]
            channels.append(vals)
            base = base_by_idx[i]
            traces[base] = vals
            channel_map[base] = tag

    trace_len = max((len(v) for v in traces.values()), default=0)
    ploc_tag = "PLOC2" if raw.get("PLOC2") is not None else ("PLOC1" if raw.get("PLOC1") is not None else None)
    ploc_raw = coerce_ints(raw.get("PLOC2") if raw.get("PLOC2") is not None else raw.get("PLOC1"))
    peak_positions = normalize_peak_positions(ploc_raw, trace_len)
    qc_pairs = normalize_peak_pairs_qc(ploc_raw, trace_len)
    qc_peak_positions = [pos for _, pos in qc_pairs]
    qualities = coerce_ints(raw.get("PCON2") if raw.get("PCON2") is not None else raw.get("PCON1"))
    basecalls = decode_basecalls(raw.get("PBAS2") if raw.get("PBAS2") is not None else raw.get("PBAS1"))

    return Ab1TraceBundle(
        traces_by_base=traces,
        channels=channels,
        data_tags=tags,
        channel_bases=channel_bases,
        channel_map=channel_map,
        basecalls=basecalls,
        peak_positions=peak_positions,
        qc_peak_positions=qc_peak_positions,
        qc_pairs=qc_pairs,
        ploc_total=len(ploc_raw),
        qualities=qualities,
        trace_len=trace_len,
        ploc_tag=ploc_tag,
    )
