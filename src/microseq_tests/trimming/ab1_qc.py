from __future__ import annotations

from dataclasses import dataclass
import math
import statistics
from pathlib import Path
from typing import Iterable
import logging

from Bio import SeqIO

from microseq_tests.utility.utils import load_config
from microseq_tests.trimming.ab1_to_fastq import build_ab1_output_key_map

L = logging.getLogger(__name__)

TRACE_QC_COLUMNS = [
    "trace_has_data",
    "trace_n_bases",
    "trace_median_peak",
    "trace_noise_mad",
    "trace_snr",
    "trace_mixed_peak_frac",
    "trace_mixed_peak_p95",
    "trace_mixed_base_frac",
    "trace_n_base_frac",
    "trace_sn_tag_median",
    "trace_noise_tag_median",
    "trace_spacing_tag_median",
    "trace_qc_flags",
    "trace_qc_status",
]

IUPAC_MIXED = set("RYSWKMBDHV")


@dataclass(frozen=True)
class TraceQcParams:
    mixed_ratio_thresh: float = 0.25
    mixed_frac_fail: float = 0.10
    min_peak_height_for_mixed: float = 100.0
    snr_warn: float = 100.0
    snr_fail: float = 25.0
    peak_window: int = 8


@dataclass(frozen=True)
class TraceQcSummary:
    trace_has_data: bool
    trace_n_bases: int | None
    trace_median_peak: float | None
    trace_noise_mad: float | None
    trace_snr: float | None
    trace_mixed_peak_frac: float | None
    trace_mixed_peak_p95: float | None
    trace_mixed_base_frac: float | None
    trace_n_base_frac: float | None
    trace_sn_tag_median: float | None
    trace_noise_tag_median: float | None
    trace_spacing_tag_median: float | None
    trace_qc_flags: tuple[str, ...]
    trace_qc_status: str | None

    def to_row(self) -> list[str]:
        def _fmt(val: object) -> str:
            if val is None:
                return ""
            if isinstance(val, float) and math.isnan(val):
                return ""
            if isinstance(val, bool):
                return "1" if val else "0"
            if isinstance(val, (int, str)):
                return str(val)
            return f"{val:.3f}"

        values = [
            self.trace_has_data,
            self.trace_n_bases,
            self.trace_median_peak,
            self.trace_noise_mad,
            self.trace_snr,
            self.trace_mixed_peak_frac,
            self.trace_mixed_peak_p95,
            self.trace_mixed_base_frac,
            self.trace_n_base_frac,
            self.trace_sn_tag_median,
            self.trace_noise_tag_median,
            self.trace_spacing_tag_median,
            ";".join(self.trace_qc_flags),
            self.trace_qc_status,
        ]
        return [_fmt(val) for val in values]



def trace_qc_header() -> str:
    return "\t".join(TRACE_QC_COLUMNS)



def _trace_qc_defaults() -> TraceQcParams:
    cfg = load_config()
    trace_cfg = cfg.get("trace_qc", {})
    return TraceQcParams(
        mixed_ratio_thresh=trace_cfg.get("mixed_ratio_thresh", TraceQcParams.mixed_ratio_thresh),
        mixed_frac_fail=trace_cfg.get("mixed_frac_fail", TraceQcParams.mixed_frac_fail),
        min_peak_height_for_mixed=trace_cfg.get(
            "min_peak_height_for_mixed", TraceQcParams.min_peak_height_for_mixed
        ),
        snr_warn=trace_cfg.get("snr_warn", TraceQcParams.snr_warn),
        snr_fail=trace_cfg.get("snr_fail", TraceQcParams.snr_fail),
        peak_window=trace_cfg.get("peak_window", TraceQcParams.peak_window),
    )



def _trace_qc_apply_thresholds() -> bool:
    cfg = load_config()
    trace_cfg = cfg.get("trace_qc", {})
    return bool(trace_cfg.get("enable_flags", False))



def _coerce_basecalls(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("ascii", errors="ignore")
    if isinstance(value, str):
        return value
    if isinstance(value, Iterable):
        try:
            return "".join(chr(int(v)) for v in value)
        except (ValueError, TypeError):
            return ""
    return ""



def _coerce_float_values(value: object) -> list[float]:
    if value is None:
        return []
    if hasattr(value, "tolist"):
        value = value.tolist()
    if isinstance(value, (list, tuple)):
        out: list[float] = []
        for item in value:
            try:
                out.append(float(item))
            except (TypeError, ValueError):
                continue
        return out
    if isinstance(value, (bytes, bytearray)):
        try:
            return [float(value.decode("ascii", errors="ignore"))]
        except ValueError:
            return []
    try:
        return [float(value)]
    except (TypeError, ValueError):
        return []



def _median_or_none(values: list[float]) -> float | None:
    if not values:
        return None
    return float(statistics.median(values))



def _mad(values: Iterable[float]) -> float:
    vals = list(values)
    if not vals:
        return float("nan")
    med = statistics.median(vals)
    return statistics.median([abs(v - med) for v in vals])



def _percentile(values: Iterable[float], p: float) -> float:
    vals = sorted(values)
    if not vals:
        return float("nan")
    k = (len(vals) - 1) * p
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return float(vals[int(k)])
    return float(vals[f] + (vals[c] - vals[f]) * (k - f))



def _to_none_if_nan(value: float | None) -> float | None:
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    return value



def _select_trace_channels(raw: dict) -> list[list[int]]:
    if all(k in raw for k in ("DATA9", "DATA10", "DATA11", "DATA12")):
        keys = ("DATA9", "DATA10", "DATA11", "DATA12")
    elif all(k in raw for k in ("DATA1", "DATA2", "DATA3", "DATA4")):
        keys = ("DATA1", "DATA2", "DATA3", "DATA4")
    else:
        return []
    channels: list[list[int]] = []
    for key in keys:
        channels.append([int(v) for v in raw[key]])
    return channels



def _resolve_ploc(raw: dict) -> list[int]:
    ploc = raw.get("PLOC2") or raw.get("PLOC1")
    if not ploc:
        return []
    return [int(v) for v in ploc]



def _normalize_ploc(ploc: list[int], trace_len: int) -> list[int]:
    if not ploc:
        return []
    if max(ploc) < trace_len and min(ploc) >= 0:
        return ploc
    if max(ploc) - 1 < trace_len and min(ploc) >= 1:
        return [idx - 1 for idx in ploc]
    return []



def _build_intervals(idx: list[int], trace_len: int) -> list[tuple[int, int]]:
    if not idx:
        return []
    intervals: list[tuple[int, int]] = []
    for i, pos in enumerate(idx):
        prev_pos = idx[i - 1] if i > 0 else 0
        next_pos = idx[i + 1] if i + 1 < len(idx) else trace_len - 1
        lo = max(0, (prev_pos + pos) // 2)
        hi = min(trace_len - 1, (pos + next_pos) // 2)
        if hi < lo:
            continue
        intervals.append((lo, hi))
    return intervals



def _interval_channel_max(channels: list[list[int]], lo: int, hi: int) -> list[float]:
    if hi < lo:
        return []
    vals: list[float] = []
    for ch in channels:
        vals.append(float(max(ch[lo : hi + 1])))
    return vals



def read_ab1_vendor_metrics(record) -> dict[str, float | None]:
    raw = record.annotations.get("abif_raw", {})

    basecalls = _coerce_basecalls(raw.get("PBAS2") or raw.get("PBAS1") or str(record.seq))
    mixed = sum(base in IUPAC_MIXED for base in basecalls)
    n_count = basecalls.count("N")

    sn_values = _coerce_float_values(raw.get("S/N%1") or raw.get("S/N%"))
    noise_values = _coerce_float_values(raw.get("NOIS1") or raw.get("NOIS"))
    spacing_values = _coerce_float_values(raw.get("SPAC1") or raw.get("SPAC"))

    return {
        "mixed_base_frac": (mixed / len(basecalls)) if basecalls else None,
        "n_base_frac": (n_count / len(basecalls)) if basecalls else None,
        "sn_tag_median": _median_or_none(sn_values),
        "noise_tag_median": _median_or_none(noise_values),
        "spacing_tag_median": _median_or_none(spacing_values),
    }



def compute_trace_qc(
    channels: list[list[int]],
    ploc: list[int],
    params: TraceQcParams,
) -> tuple[float | None, float | None, float | None, float | None, float | None]:
    trace_len = len(channels[0]) if channels else 0
    idx = _normalize_ploc(ploc, trace_len)
    if not channels or not idx:
        return None, None, None, None, None

    intervals = _build_intervals(idx, trace_len)
    use_interval = bool(intervals)

    top_peaks: list[float] = []
    ratios: list[float] = []

    if use_interval:
        for lo, hi in intervals:
            vals = _interval_channel_max(channels, lo, hi)
            if len(vals) != 4:
                continue
            vals_sorted = sorted(vals)
            top = float(vals_sorted[-1])
            second = float(vals_sorted[-2])
            top_peaks.append(top)
            if top >= params.min_peak_height_for_mixed:
                ratios.append(second / top if top > 0 else 0.0)
    else:
        for pos in idx:
            vals = [ch[pos] for ch in channels]
            vals_sorted = sorted(vals)
            top = float(vals_sorted[-1])
            second = float(vals_sorted[-2])
            top_peaks.append(top)
            if top >= params.min_peak_height_for_mixed:
                ratios.append(second / top if top > 0 else 0.0)

    median_peak = float(statistics.median(top_peaks)) if top_peaks else float("nan")

    mask = [True] * trace_len
    for pos in idx:
        lo = max(0, pos - params.peak_window)
        hi = min(trace_len, pos + params.peak_window + 1)
        for i in range(lo, hi):
            mask[i] = False

    baseline_values: list[float] = []
    for ch in channels:
        baseline_values.extend(float(v) for i, v in enumerate(ch) if mask[i])

    noise_mad = float(_mad(baseline_values))
    snr = float(median_peak / (noise_mad + 1e-9))
    mixed_frac = float(sum(r > params.mixed_ratio_thresh for r in ratios) / max(1, len(ratios)))
    mixed_p95 = float(_percentile(ratios, 0.95)) if ratios else float("nan")

    return (
        _to_none_if_nan(median_peak),
        _to_none_if_nan(noise_mad),
        _to_none_if_nan(snr),
        _to_none_if_nan(mixed_frac),
        _to_none_if_nan(mixed_p95),
    )



def summarize_ab1_qc(
    path: Path,
    params: TraceQcParams | None = None,
    *,
    apply_thresholds: bool | None = None,
) -> TraceQcSummary:
    params = params or _trace_qc_defaults()
    if apply_thresholds is None:
        apply_thresholds = _trace_qc_apply_thresholds()

    record = SeqIO.read(path, "abi")
    raw = record.annotations.get("abif_raw", {})
    channels = _select_trace_channels(raw)
    ploc = _resolve_ploc(raw)

    vendor = read_ab1_vendor_metrics(record)
    median_peak, noise_mad, snr, mixed_frac, mixed_p95 = compute_trace_qc(channels, ploc, params)

    missing_primitives = not channels or not ploc

    flags: list[str] = []
    status: str | None = None
    if apply_thresholds:
        if missing_primitives:
            flags.append("TRACE_DATA_MISSING")
            status = "NA"
        else:
            status = "PASS"
            if snr is None:
                flags.append("SNR_UNDEFINED")
                status = "NA"
            elif snr < params.snr_fail:
                flags.append("LOW_SNR")
            elif snr < params.snr_warn:
                flags.append("LOW_SNR_WARN")

            if mixed_frac is None:
                flags.append("MIXED_UNDEFINED")
                status = "NA"
            elif mixed_frac > params.mixed_frac_fail:
                flags.append("MIXED_PEAKS")

            if status != "NA":
                if any(flag in {"LOW_SNR", "MIXED_PEAKS"} for flag in flags):
                    status = "FAIL"
                elif flags:
                    status = "WARN"

    n_bases = len(ploc) if ploc else None

    return TraceQcSummary(
        trace_has_data=bool(channels and ploc),
        trace_n_bases=n_bases,
        trace_median_peak=median_peak,
        trace_noise_mad=noise_mad,
        trace_snr=snr,
        trace_mixed_peak_frac=mixed_frac,
        trace_mixed_peak_p95=mixed_p95,
        trace_mixed_base_frac=vendor["mixed_base_frac"],
        trace_n_base_frac=vendor["n_base_frac"],
        trace_sn_tag_median=vendor["sn_tag_median"],
        trace_noise_tag_median=vendor["noise_tag_median"],
        trace_spacing_tag_median=vendor["spacing_tag_median"],
        trace_qc_flags=tuple(flags),
        trace_qc_status=status,
    )



def summarize_ab1_folder(
    path: Path,
    params: TraceQcParams | None = None,
    *,
    apply_thresholds: bool | None = None,
) -> dict[str, TraceQcSummary]:
    summaries: dict[str, TraceQcSummary] = {}
    key_map = build_ab1_output_key_map(path)
    for ab1 in sorted(path.rglob("*.ab1")):
        key = key_map.get(ab1, ab1.stem)
        if key in summaries:
            L.warning("Duplicate AB1 trace QC key %s for %s", key, ab1)
        try:
            summaries[key] = summarize_ab1_qc(ab1, params=params, apply_thresholds=apply_thresholds)
        except Exception as exc:
            L.warning("Trace QC parse failed for %s: %s", ab1, exc)
            summaries[key] = TraceQcSummary(
                trace_has_data=False,
                trace_n_bases=None,
                trace_median_peak=None,
                trace_noise_mad=None,
                trace_snr=None,
                trace_mixed_peak_frac=None,
                trace_mixed_peak_p95=None,
                trace_mixed_base_frac=None,
                trace_n_base_frac=None,
                trace_sn_tag_median=None,
                trace_noise_tag_median=None,
                trace_spacing_tag_median=None,
                trace_qc_flags=("TRACE_PARSE_ERROR",),
                trace_qc_status="NA",
            )
    return summaries
