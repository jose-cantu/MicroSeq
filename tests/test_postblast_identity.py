# tests/test_postblast_identity.py
from __future__ import annotations
import json
from pathlib import Path

import numpy as np          # for a dense BIOM matrix
import pandas as pd
import pytest

biom = pytest.importorskip("biom")

from microseq_tests.post_blast_analysis import run as postblast_run
from microseq_tests.utility.add_taxonomy import run_taxonomy_join
from microseq_tests import microseq as microseq_cli


# ---------------------------------------------------------------------------
# 1. CLI → pipeline passthrough
# ---------------------------------------------------------------------------

def test_cli_passes_identity(monkeypatch, tmp_path):
    """Ensure --post_blast_identity reaches postblast_run()."""

    # -- minimal synthetic inputs -------------------------------------------
    blast_tsv = tmp_path / "hits.tsv"
    blast_tsv.write_text(
        "qseqid\tsseqid\tpident\tevalue\tbitscore\ttaxonomy\n"
        "S1\tGG2_0001\t99.0\t1e-10\t200\tk__Bacteria; p__Proteobacteria\n"
    )
    meta_tsv = tmp_path / "meta.tsv"
    meta_tsv.write_text("sample_id\nS1\n")

    # -- monkey-capture postblast_run ---------------------------------------
    captured: dict[str, float] = {}

    def fake_run(blast_file, metadata_file, out_biom, **kw):
        captured.update(kw)                       # record kwargs
        # write a valid 1×1 BIOM so caller code doesn’t error
        table = biom.Table(np.array([[1]]), ["k__Bacteria"], ["S1"])
        with biom.util.biom_open(out_biom, "w") as fh:
            table.to_hdf5(fh, "test")

    monkeypatch.setattr("microseq_tests.microseq.postblast_run", fake_run)

    # -- invoke CLI (global flags precede sub-command) -----------------------
    argv = [
        "microseq",            # prog name
        "postblast",
        "-b", str(blast_tsv),
        "-m", str(meta_tsv),
        "-o", "dummy.biom",
        "--post_blast_identity", "95",
    ]
    monkeypatch.setattr("sys.argv", argv)
    microseq_cli.main()

    assert captured["identity_th"] == 95.0


# ---------------------------------------------------------------------------
# 2. Tail-pipeline sanity check
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("identity", [97.0, 95.0])
def test_end_to_end_tail(identity, tmp_path):
    """
    add_taxonomy → postblast
    then validate BIOM / CSV / JSON outputs.
    """

    # --- synthetic BLAST & taxonomy tables ---------------------------------
    hits = tmp_path / "hits.tsv"
    hits.write_text(
        "qseqid\tsseqid\tpident\tevalue\tbitscore\n"
        "S1\tGG2_0001\t99.0\t1e-10\t200\n"
        "S1\tGG2_0002\t94.0\t1e-9\t180\n"
    )

    tax = tmp_path / "taxonomy.tsv"
    tax.write_text(
        "sseqid\ttaxonomy\n"
        "GG2_0001\tk__Bacteria; p__Firmicutes\n"
        "GG2_0002\tk__Bacteria; p__Proteobacteria\n"
    )

    meta = tmp_path / "meta.tsv"
    meta.write_text("sample_id\tSource\nS1\tEye\n")

    # --- taxonomy join ------------------------------------------------------
    hits_tax = tmp_path / "hits_tax.tsv"
    run_taxonomy_join(hits, tax, hits_tax)

    # --- postblast ----------------------------------------------------------
    biom_out = tmp_path / "out.biom"
    postblast_run(
        hits_tax,
        meta,
        biom_out,
        sample_col="sample_id",
        identity_th=identity,
        taxonomy_col="taxonomy",
    )

    csv_out = biom_out.with_suffix(".csv")
    assert biom_out.exists() and csv_out.exists()

    # --- validate BIOM dimensions ------------------------------------------
    with biom.util.biom_open(biom_out) as fh:
        table = biom.Table.from_hdf5(fh)

    # _choose_best_hit() always keeps *one* row per sample
    assert table.shape == (1, 1)

    df = pd.read_csv(csv_out, index_col=0)
    assert len(df) == 1

    # --- JSON export (via BIOM API) ----------------------------------------
    json_out = biom_out.with_suffix(".json")
    json_out.write_text(table.to_json("MicroSeq"))
    # confirm it’s valid JSON
    json.loads(json_out.read_text())

