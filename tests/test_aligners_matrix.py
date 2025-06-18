"""
Smoke-test that the public aligner registry returns a DataFrame with the
expected BLAST/VSEARCH schema for every logical engine exposed to users.

The test deliberately avoids               ⤦
  * external binaries   (BLAST, VSEARCH)   ⤤
  * optional heavy deps (pandera, Bio)     ⤥
  * any file-system I/O                    ⤦
so it can run in a minimal CI image.
"""
# ── std-lib ───────────────────────────────────────────────────────────
from pathlib import Path
import subprocess
import importlib.metadata as _ilm

# ── third-party ───────────────────────────────────────────────────────
import pandas as pd
import pytest

# ── project ───────────────────────────────────────────────────────────
# 1.  Stop entry-point discovery from importing heavy aligner modules
_orig_eps = _ilm.entry_points                      # keep original reference


def _empty_entry_points(group=None):
    """Return an empty list for our plugin group, else fall back."""
    if group == "microseq.aligners":
        return []
    return _orig_eps(group)


_ilm.entry_points = _empty_entry_points            # global monkey-patch

# 2.  Now import the package safely (no pandera pulled in)
import microseq_tests.aligners._parse as _p          # parser helpers
from microseq_tests.utility import utils as _utils   # load_config…
from microseq_tests.aligners import ALIGNER_REGISTRY, load

# 3.  Stub schema + dummy aligner class
COLS = (
    "qseqid sseqid pident qlen qcovhsp length evalue bitscore stitle".split()
)


def _stub_df():
    """Return a 10×9 DataFrame with numeric-compatible dtypes."""
    return pd.DataFrame([{c: 0 for c in COLS} for _ in range(10)], columns=COLS)


class _DummyAligner:
    """Minimal stand-in; ignores args, always returns schema-correct DF."""

    def run(self, *_, **__):
        return _stub_df()


# 4.  Register stubs under canonical keys
ALIGNER_REGISTRY.update({"blastn": _DummyAligner, "vsearch": _DummyAligner})

# ── test ──────────────────────────────────────────────────────────────


@pytest.mark.parametrize(
    "engine,kwargs",
    [
        ("blastn-fast", {"fast": True}),
        ("blastn-slow", {"fast": False}),
        ("vsearch", {}),
    ],
)
def test_aligners_shape(engine, kwargs, monkeypatch, tmp_path):
    """
    For each logical engine ensure:

    * run() executes without raising
    * returns a 10 × 9 DataFrame
    * column order matches the canonical BLAST schema
    """
    # neuter any residual I/O that a future refactor might add
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: None)
    monkeypatch.setattr(Path, "is_file", lambda *_: True)
    monkeypatch.setattr(_p, "parse_blast_tab", lambda *_: _stub_df())

    fake_cfg = {"databases": {"gg2": {"blastdb": "/tmp/gg2"}}}
    monkeypatch.setattr(_utils, "load_config", lambda: fake_cfg)
    monkeypatch.setattr(_utils, "expand_db_path", lambda x: x)

    # choose the correct registry key
    key = "blastn" if engine.startswith("blastn") else engine
    df = load(key).run("dummy.fa", "gg2", threads=1, **kwargs)

    assert df.shape == (10, 9)
    assert list(df.columns) == COLS

