# tests/test_postblast_taxonomy_format.py
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from microseq_tests.post_blast_analysis import run as postblast_run


# ---------------------------------------------------------------------------
# Each taxonomy format should be parsed into the standard seven ranks.
# ---------------------------------------------------------------------------


def _make_hits(path: Path, lineage: str) -> None:
    """Write a minimal BLAST TSV with a taxonomy column."""
    path.write_text(
        "#sample_id\tsseqid\tpident\tevalue\tbitscore\ttaxonomy\n"
        f"S1\tseq1\t99.0\t1e-10\t200\t{lineage}\n"
    )


def _simple_meta(path: Path) -> None:
    """Write a one-column metadata TSV."""
    path.write_text("sample_id\nS1\n")


# ---------------------------------------------------------------------------
# Parametrised test over the three supported taxonomy formats.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "fmt,lineage",
    [
        (
            "gg2",
            "k__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; "
            "f__Lactobacillaceae; g__Lactobacillus; s__Lactobacillus acidophilus",
        ),
        (
            "silva",
            "d__Bacteria; p__Firmicutes; c__Bacilli; o__Lactobacillales; "
            "f__Lactobacillaceae; g__Lactobacillus; s__Lactobacillus acidophilus",
        ),
        (
            "ncbi",
            "Bacteria; Firmicutes; Bacilli; Lactobacillales; Lactobacillaceae; "
            "Lactobacillus; Lactobacillus acidophilus",
        ),
    ],
)
def test_taxonomy_formats(fmt: str, lineage: str, tmp_path: Path) -> None:
    hits = tmp_path / f"{fmt}.tsv"
    meta = tmp_path / "meta.tsv"
    out_biom = tmp_path / f"{fmt}.biom"

    _make_hits(hits, lineage)
    _simple_meta(meta)

    postblast_run(hits, meta, out_biom, taxonomy_format=fmt, taxonomy_col="taxonomy")

    tax_csv = out_biom.with_name(out_biom.stem + "_taxonomy_only.csv")
    assert out_biom.exists() and tax_csv.exists()

    df = pd.read_csv(tax_csv)
    # OTU_ID + seven taxonomy columns
    assert df.shape[1] == 8
