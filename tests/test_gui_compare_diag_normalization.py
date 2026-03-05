from __future__ import annotations

import pytest

pytest.importorskip("pandas")
pytest.importorskip("Bio")
pytest.importorskip("PySide6")

from microseq_tests.gui.main_window import _normalize_legacy_compare_diag_detail


def test_normalize_legacy_compare_diag_detail_rewrites_numeric_engine_token() -> None:
    text = "merge_status=ambiguous_overlap; overlap_len=435;identity=0.9507; engine=435"
    out = _normalize_legacy_compare_diag_detail(text, "biopython")
    assert out.endswith("engine=biopython")


def test_normalize_legacy_compare_diag_detail_keeps_new_format() -> None:
    text = "merge_status=ambiguous_overlap; overlap_len=435;identity=0.9507; engine=biopython"
    out = _normalize_legacy_compare_diag_detail(text, "biopython")
    assert out == text
