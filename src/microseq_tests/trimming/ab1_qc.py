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
from microseq_tests.trimming.ab1_trace_utils import Ab1TraceBundle, extract_ab1_trace_bundle

L = logging.getLogger(__name__)

TRACE_QC_COLUMNS = [
    "trace_has_data",
    "trace_n_bases",
    "trace_median_peak",
    "trace_noise_mad",
    "trace_snr",
    "trace_snr_global_mad",
    "trace_snr_basecall_aware",
    "trace_snr_mode",
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
    mixed_min_run: int = 3
    core_min_quality: int = 20
    core_start_frac: float = 0.10
    core_end_frac: float = 0.85
    snr_mode: str = "global_mad"
    snr_warn_global_mad: float = 12.0
    snr_fail_global_mad: float = 7.0
    snr_warn_basecall_aware: float = 2.0
    snr_fail_basecall_aware: float = 1.3
    low_snr_fail_core_pcon_min: float = 20.0
    low_snr_fail_avg_q_min: float = 20.0
    low_snr_missing_quality_policy: str = "warn_if_quality_missing"
    snr_fail_hard: float = 0.8
    peak_window: int = 8
    ploc_local_half_window: int = 1


@dataclass(frozen=True)
class TraceQcSummary:
    trace_has_data: bool
    trace_n_bases: int | None
    trace_median_peak: float | None
    trace_noise_mad: float | None
    trace_snr: float | None
    trace_snr_global_mad: float | None
    trace_snr_basecall_aware: float | None
    trace_snr_mode: str | None
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
            self.trace_snr_global_mad,
            self.trace_snr_basecall_aware,
            self.trace_snr_mode,
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
        mixed_min_run=trace_cfg.get("mixed_min_run", TraceQcParams.mixed_min_run),
        core_min_quality=trace_cfg.get("core_min_quality", TraceQcParams.core_min_quality),
        core_start_frac=trace_cfg.get("core_start_frac", TraceQcParams.core_start_frac),
        core_end_frac=trace_cfg.get("core_end_frac", TraceQcParams.core_end_frac),
        snr_mode=str(trace_cfg.get("snr_mode", TraceQcParams.snr_mode) or TraceQcParams.snr_mode).strip().lower(),
        snr_warn_global_mad=trace_cfg.get("snr_warn_global_mad", TraceQcParams.snr_warn_global_mad),
        snr_fail_global_mad=trace_cfg.get("snr_fail_global_mad", TraceQcParams.snr_fail_global_mad),
        snr_warn_basecall_aware=trace_cfg.get("snr_warn_basecall_aware", TraceQcParams.snr_warn_basecall_aware),
        snr_fail_basecall_aware=trace_cfg.get("snr_fail_basecall_aware", TraceQcParams.snr_fail_basecall_aware),
        low_snr_fail_core_pcon_min=trace_cfg.get(
            "low_snr_fail_core_pcon_min", TraceQcParams.low_snr_fail_core_pcon_min
        ),
        low_snr_fail_avg_q_min=trace_cfg.get("low_snr_fail_avg_q_min", TraceQcParams.low_snr_fail_avg_q_min),
        low_snr_missing_quality_policy=str(
            trace_cfg.get("low_snr_missing_quality_policy", TraceQcParams.low_snr_missing_quality_policy)
            or TraceQcParams.low_snr_missing_quality_policy
        ).strip().lower(),
        snr_fail_hard=trace_cfg.get("snr_fail_hard", TraceQcParams.snr_fail_hard),
        peak_window=trace_cfg.get("peak_window", TraceQcParams.peak_window),
        ploc_local_half_window=trace_cfg.get("ploc_local_half_window", TraceQcParams.ploc_local_half_window),
    )


def _trace_qc_apply_thresholds() -> bool:
    cfg = load_config()
    trace_cfg = cfg.get("trace_qc", {})
    return bool(trace_cfg.get("enable_flags", False))


def _trace_qc_diagnostic_logging_enabled() -> bool:
    cfg = load_config()
    trace_cfg = cfg.get("trace_qc", {})
    return bool(trace_cfg.get("diagnostic_logging", True))


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


def _longest_true_run(values: list[bool]) -> int:
    best = 0
    cur = 0
    for v in values:
        if v:
            cur += 1
            best = max(best, cur)
        else:
            cur = 0
    return best


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




def _snr_severity(value: float | None, warn: float, fail: float) -> int:
    """0=PASS, 1=WARN, 2=FAIL, -1=undefined."""
    if value is None:
        return -1
    if value < fail:
        return 2
    if value < warn:
        return 1
    return 0

def _snr_thresholds(params: TraceQcParams) -> tuple[float, float]:
    if params.snr_mode == "basecall_aware":
        return float(params.snr_warn_basecall_aware), float(params.snr_fail_basecall_aware)
    return float(params.snr_warn_global_mad), float(params.snr_fail_global_mad)


def compute_trace_qc(
    bundle: Ab1TraceBundle,
    params: TraceQcParams,
) -> tuple[
    float | None,
    float | None,
    float | None,
    float | None,
    float | None,
    int | None,
    int | None,
    float | None,
    float | None,
    float | None,
]:
    channels = bundle.channels
    qpairs = list(bundle.qc_pairs)
    qualities = bundle.qualities
    trace_len = len(channels[0]) if channels else 0
    if not channels or not qpairs:
        return None, None, None, None, None, None, None, None, None, None

    n_idx = len(qpairs)
    core_start_i = max(0, min(n_idx - 1, int(n_idx * params.core_start_frac)))
    core_end_i = max(core_start_i + 1, min(n_idx, int(n_idx * params.core_end_frac)))
    core_pairs = qpairs[core_start_i:core_end_i]

    if qualities:
        q_len = len(qualities)
        core_pairs = [pair for pair in core_pairs if pair[0] < q_len and qualities[pair[0]] >= params.core_min_quality]

    if not core_pairs:
        core_pairs = qpairs[core_start_i:core_end_i]

    call_indices = [ci for ci, _ in core_pairs]
    core_positions = [pos for _, pos in core_pairs]

    top_peaks: list[float] = []
    ratios: list[float] = []
    mixed_events: list[bool] = []
    snr_per_base: list[float] = []

    for call_i, pos in core_pairs:
        local_lo = max(0, pos - params.ploc_local_half_window)
        local_hi = min(trace_len, pos + params.ploc_local_half_window + 1)
        vals = [float(max(ch[local_lo:local_hi])) for ch in channels]
        vals_sorted = sorted(vals)
        top = float(vals_sorted[-1])
        second = float(vals_sorted[-2])
        top_peaks.append(top)
        ratio = second / top if top > 0 else 0.0
        if top >= params.min_peak_height_for_mixed:
            ratios.append(ratio)
            mixed_events.append(ratio > params.mixed_ratio_thresh)
        else:
            mixed_events.append(False)

        if call_i < len(bundle.basecalls):
            base = bundle.basecalls[call_i].upper()
            called_trace = bundle.traces_by_base.get(base, [])
            if called_trace and base in {"A", "C", "G", "T"}:
                sig = float(max(called_trace[local_lo:local_hi]))
                other = [
                    float(max(bundle.traces_by_base[b][local_lo:local_hi]))
                    for b in ("A", "C", "G", "T")
                    if b != base and bundle.traces_by_base.get(b)
                ]
                if other:
                    noise = float(statistics.median(other))
                    snr_per_base.append(sig / (noise + 1e-9))

    median_peak = float(statistics.median(top_peaks)) if top_peaks else float("nan")
    core_lo = min(core_positions) if core_positions else None
    core_hi = max(core_positions) if core_positions else None

    mask = [True] * trace_len
    for pos in core_positions:
        lo = max(0, pos - params.peak_window)
        hi = min(trace_len, pos + params.peak_window + 1)
        for i in range(lo, hi):
            mask[i] = False

    baseline_values: list[float] = []
    if core_lo is not None and core_hi is not None:
        for ch in channels:
            baseline_values.extend(float(ch[i]) for i in range(core_lo, core_hi + 1) if mask[i])
        if not baseline_values:
            for ch in channels:
                baseline_values.extend(float(ch[i]) for i in range(core_lo, core_hi + 1))
    else:
        for ch in channels:
            baseline_values.extend(float(v) for v in ch)

    noise_mad = float(_mad(baseline_values))
    snr_global_mad = float(median_peak / (noise_mad + 1e-9))
    snr_basecall_aware = float(statistics.median(snr_per_base)) if snr_per_base else float("nan")
    snr_selected = snr_basecall_aware if params.snr_mode == "basecall_aware" else snr_global_mad

    mixed_frac = float(sum(mixed_events) / max(1, len(mixed_events)))
    mixed_p95 = float(_percentile(ratios, 0.95)) if ratios else float("nan")

    if _longest_true_run(mixed_events) < params.mixed_min_run:
        mixed_frac = 0.0

    core_qualities = [float(qualities[i]) for i in call_indices if i < len(qualities)] if qualities else []
    core_quality_median = float(statistics.median(core_qualities)) if core_qualities else None

    return (
        _to_none_if_nan(median_peak),
        _to_none_if_nan(noise_mad),
        _to_none_if_nan(snr_selected),
        _to_none_if_nan(mixed_frac),
        _to_none_if_nan(mixed_p95),
        core_lo,
        core_hi,
        _to_none_if_nan(core_quality_median),
        _to_none_if_nan(snr_global_mad),
        _to_none_if_nan(snr_basecall_aware),
    )


def summarize_ab1_qc(
    path: Path,
    params: TraceQcParams | None = None,
    *,
    apply_thresholds: bool | None = None,
    file_avg_q: float | None = None,
) -> TraceQcSummary:
    params = params or _trace_qc_defaults()
    if apply_thresholds is None:
        apply_thresholds = _trace_qc_apply_thresholds()

    record = SeqIO.read(path, "abi")
    raw = record.annotations.get("abif_raw", {})
    bundle = extract_ab1_trace_bundle(raw)

    vendor = read_ab1_vendor_metrics(record)
    (
        median_peak,
        noise_mad,
        snr,
        mixed_frac,
        mixed_p95,
        core_lo,
        core_hi,
        core_quality_median,
        snr_global_mad,
        snr_basecall_aware,
    ) = compute_trace_qc(bundle, params)

    missing_primitives = not bundle.channels or not bundle.qc_pairs
    snr_warn, snr_fail = _snr_thresholds(params)

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
            elif snr < snr_fail:
                quality_support_low = False
                if core_quality_median is not None and core_quality_median < params.low_snr_fail_core_pcon_min:
                    quality_support_low = True
                if file_avg_q is not None and file_avg_q < params.low_snr_fail_avg_q_min:
                    quality_support_low = True

                if quality_support_low:
                    flags.append("LOW_SNR")
                else:
                    if params.low_snr_missing_quality_policy == "fail_if_snr_below_hard_floor" and snr < params.snr_fail_hard:
                        flags.append("LOW_SNR_HARD_FAIL")
                    else:
                        flags.append("LOW_SNR_NO_QUALITY_SUPPORT")
            elif snr < snr_warn:
                flags.append("LOW_SNR_WARN")

            if mixed_frac is None:
                flags.append("MIXED_UNDEFINED")
                status = "NA"
            elif mixed_frac > params.mixed_frac_fail:
                flags.append("MIXED_PEAKS")

            if status != "NA":
                if any(flag in {"LOW_SNR", "LOW_SNR_HARD_FAIL", "MIXED_PEAKS"} for flag in flags):
                    status = "FAIL"
                elif flags:
                    status = "WARN"

            # Canonical mode decides FAIL; disagreement escalates to non-blocking WARN.
            global_warn, global_fail = float(params.snr_warn_global_mad), float(params.snr_fail_global_mad)
            aware_warn, aware_fail = float(params.snr_warn_basecall_aware), float(params.snr_fail_basecall_aware)
            sev_global = _snr_severity(snr_global_mad, global_warn, global_fail)
            sev_aware = _snr_severity(snr_basecall_aware, aware_warn, aware_fail)
            if status != "NA" and sev_global >= 0 and sev_aware >= 0 and sev_global != sev_aware:
                flags.append("SNR_MODE_DISCORDANT_WARN")
                if status == "PASS":
                    status = "WARN"

    if _trace_qc_diagnostic_logging_enabled():
        L.info(
            "AB1_TRACE_DIAG file=%s data_tags=%s channel_map=%s snr_mode=%s median_primary_peak=%s baseline_MAD=%s SNR_selected=%s SNR_global_mad=%s SNR_basecall_aware=%s mixed_peak_fraction=%s qc_window_start=%s qc_window_end=%s n_called_bases=%s ploc_total=%s ploc_in_range=%s core_pcon_median=%s snr_warn=%s snr_fail=%s",
            path.name,
            ",".join(bundle.data_tags or ()),
            bundle.channel_map,
            params.snr_mode,
            median_peak,
            noise_mad,
            snr,
            snr_global_mad,
            snr_basecall_aware,
            mixed_frac,
            core_lo,
            core_hi,
            len(bundle.peak_positions),
            bundle.ploc_total,
            len(bundle.qc_pairs),
            core_quality_median,
            snr_warn,
            snr_fail,
        )

    n_bases = len(bundle.qc_pairs) if bundle.qc_pairs else None

    return TraceQcSummary(
        trace_has_data=bool(bundle.channels and bundle.qc_pairs),
        trace_n_bases=n_bases,
        trace_median_peak=median_peak,
        trace_noise_mad=noise_mad,
        trace_snr=snr,
        trace_snr_global_mad=snr_global_mad,
        trace_snr_basecall_aware=snr_basecall_aware,
        trace_snr_mode=params.snr_mode,
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
    file_avg_q_by_key: dict[str, float] | None = None,
) -> dict[str, TraceQcSummary]:
    summaries: dict[str, TraceQcSummary] = {}
    key_map = build_ab1_output_key_map(path)
    for ab1 in sorted(path.rglob("*.ab1")):
        key = key_map.get(ab1, ab1.stem)
        if key in summaries:
            L.warning("Duplicate AB1 trace QC key %s for %s", key, ab1)
        try:
            avg_q = file_avg_q_by_key.get(key) if file_avg_q_by_key else None
            summaries[key] = summarize_ab1_qc(
                ab1,
                params=params,
                apply_thresholds=apply_thresholds,
                file_avg_q=avg_q,
            )
        except Exception as exc:
            L.warning("Trace QC parse failed for %s: %s", ab1, exc)
            summaries[key] = TraceQcSummary(
                trace_has_data=False,
                trace_n_bases=None,
                trace_median_peak=None,
                trace_noise_mad=None,
                trace_snr=None,
                trace_snr_global_mad=None,
                trace_snr_basecall_aware=None,
                trace_snr_mode=None,
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
