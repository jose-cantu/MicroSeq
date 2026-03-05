"""Deprecated compatibility module.

Expected-error based QC modeling was removed for Sanger-focused workflows.
These shims are kept only to avoid import failures in external integrations.
"""

from __future__ import annotations


def expected_errors(phred):
    _ = phred
    raise RuntimeError("expected_errors() is deprecated and not supported in Sanger mode")


def qeff_from_mee_per_kb(mee_per_kb):
    _ = mee_per_kb
    raise RuntimeError("qeff_from_mee_per_kb() is deprecated and not supported in Sanger mode")
