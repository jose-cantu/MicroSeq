from pathlib import Path

import pytest
pytest.importorskip("pandas")
from Bio import SeqIO

from microseq_tests.trimming.ab1_qc import TRACE_QC_COLUMNS, TraceQcParams, summarize_ab1_qc, compute_trace_qc
from microseq_tests.trimming.ab1_trace_utils import (
    Ab1TraceBundle,
    extract_ab1_trace_bundle,
    normalize_peak_positions,
    normalize_peak_positions_qc,
    trim_called_arrays,
)


def test_normalize_peak_positions_handles_one_based_and_clamps() -> None:
    assert normalize_peak_positions([1, 10, 12], trace_len=10) == [0, 9, 9]
    assert normalize_peak_positions([0, 3, 99], trace_len=10) == [0, 3, 9]


def test_normalize_peak_positions_qc_filters_out_of_range() -> None:
    assert normalize_peak_positions_qc([1, 10, 12], trace_len=10) == [0, 9]
    assert normalize_peak_positions_qc([0, 3, 99], trace_len=10) == [0, 3]


def test_trim_called_arrays_keeps_basecalls_without_qualities() -> None:
    basecalls, peaks, qualities = trim_called_arrays("ACGT", [1, 2, 3], [])
    assert basecalls == "ACG"
    assert peaks == [1, 2, 3]
    assert qualities == []


def test_extract_bundle_prefers_data9_to_data12() -> None:
    raw = {
        "DATA1": [1, 2, 3],
        "DATA2": [4, 5, 6],
        "DATA3": [7, 8, 9],
        "DATA4": [10, 11, 12],
        "DATA9": [13, 14, 15],
        "DATA10": [16, 17, 18],
        "DATA11": [19, 20, 21],
        "DATA12": [22, 23, 24],
        "FWO_1": b"ACGT",
        "PLOC2": [1, 2, 3],
        "PBAS2": b"ACG",
        "PCON2": [40, 40, 40],
    }
    bundle = extract_ab1_trace_bundle(raw)
    assert bundle.data_tags == ("DATA9", "DATA10", "DATA11", "DATA12")
    assert bundle.channel_map == {"A": "DATA9", "C": "DATA10", "G": "DATA11", "T": "DATA12"}
    assert bundle.peak_positions == [0, 1, 2]
    assert bundle.qc_peak_positions == [0, 1, 2]
    assert bundle.qc_pairs == [(0, 0), (1, 1), (2, 2)]


def test_qc_pairs_preserve_call_index_alignment_for_basecall_aware_snr() -> None:
    # call index 1 is intentionally out of range in qc_pairs, preserving index mapping is required.
    trace_a = [1, 1, 90, 1, 1]
    trace_c = [1, 1, 2, 1, 1]
    trace_g = [1, 1, 2, 1, 1]
    trace_t = [1, 1, 2, 1, 1]
    bundle = Ab1TraceBundle(
        traces_by_base={"A": trace_a, "C": trace_c, "G": trace_g, "T": trace_t},
        channels=[trace_a, trace_c, trace_g, trace_t],
        data_tags=("DATA9", "DATA10", "DATA11", "DATA12"),
        channel_bases=("A", "C", "G", "T"),
        channel_map={"A": "DATA9", "C": "DATA10", "G": "DATA11", "T": "DATA12"},
        basecalls="TAT",
        peak_positions=[2, 4, 2],
        qc_peak_positions=[2, 2],
        qc_pairs=[(0, 2), (2, 2)],  # call 1 dropped, call 2 retained
        ploc_total=3,
        qualities=[40, 0, 40],
        trace_len=5,
        ploc_tag="PLOC2",
    )
    params = TraceQcParams(
        snr_mode="basecall_aware",
        core_start_frac=0.0,
        core_end_frac=1.0,
        core_min_quality=0,
        peak_window=0,
        ploc_local_half_window=0,
    )
    _, _, snr_selected, _, _, _, _, core_pcon_median, _, snr_base = compute_trace_qc(bundle, params)
    assert snr_selected == snr_base
    assert snr_base is not None and 0.5 < snr_base < 2.0
    assert core_pcon_median == 40.0


def test_qc_and_viewer_use_identical_peak_index_arrays() -> None:
    pytest.importorskip("PySide6")
    from microseq_tests.gui.main_window import ChromatogramDialog

    ab1_path = Path("tests/paired_single_pair_ab1_demo_run/A+_27F_2026-01-29_C07.ab1")
    rec = SeqIO.read(ab1_path, "abi")
    raw = rec.annotations.get("abif_raw", {})
    bundle = extract_ab1_trace_bundle(raw)

    _traces, _basecalls, viewer_positions, _qualities = ChromatogramDialog._load_traces(ab1_path)

    assert viewer_positions == bundle.peak_positions[: len(viewer_positions)]


def test_trace_qc_emits_diagnostic_line(caplog: pytest.LogCaptureFixture) -> None:
    caplog.set_level("INFO", logger="microseq_tests.trimming.ab1_qc")
    ab1_path = Path("tests/paired_single_pair_ab1_demo_run/A+_1492R_2026-01-29_E07.ab1")

    summary = summarize_ab1_qc(ab1_path, apply_thresholds=False)

    assert summary.trace_has_data
    assert "AB1_TRACE_DIAG" in caplog.text
    assert "SNR_selected=" in caplog.text
    assert "SNR_global_mad=" in caplog.text
    assert "SNR_basecall_aware=" in caplog.text
    assert "snr_mode=" in caplog.text
    assert "ploc_total=" in caplog.text


def test_chromatogram_view_accepts_basecalls_without_qualities() -> None:
    pytest.importorskip("PySide6")
    from microseq_tests.gui.main_window import ChromatogramView

    view = ChromatogramView(
        {"A": [0, 5, 0], "C": [0, 0, 1], "G": [1, 0, 0], "T": [0, 0, 0]},
        basecalls="AC",
        peak_positions=[1, 2],
        qualities=[],
    )
    assert len(view._basecalls) == 2 and len(view._peak_positions) == 2


def test_trace_qc_header_has_dual_snr_columns() -> None:
    assert "trace_snr_global_mad" in TRACE_QC_COLUMNS
    assert "trace_snr_basecall_aware" in TRACE_QC_COLUMNS
    assert "trace_snr_mode" in TRACE_QC_COLUMNS
