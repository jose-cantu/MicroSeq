# tests/test_postblast_replicate_sizes.py
from __future__ import annotations

from pathlib import Path

import pytest

"""
This test suite ensures that post-BLAST analysis correctly handles sequence 
abundance weights when generating biological observation matrices (BIOM).

1. test_postblast_weights_replicate_sizes:
   Verifies that a single sequence hit is correctly multiplied by its 
   'replicate_size' weight (e.g., a weight of 3 results in an abundance of 3 
   in the final table).

2. test_postblast_weights_sum_multiple_rows:
   Confirms that when multiple distinct sequences (replicates) map to the 
   same taxonomy for the same sample, their individual weights are 
   summed correctly (e.g., weights of 2 and 4 result in an abundance of 6).

The tests utilize the [BIOM-Format](https://biom-format.org) library to 
inspect the HDF5 output and [Pandas](https://pandas.pydata.org) for data 
validation.
"""

pytest.importorskip("pandas")
pytest.importorskip("numpy")
biom = pytest.importorskip("biom")

from microseq_tests.post_blast_analysis import run as postblast_run


def test_postblast_weights_replicate_sizes(tmp_path):
    hits = tmp_path / "hits_tax.tsv"
    hits.write_text(
        "qseqid\tsample_id\tsseqid\tpident\tevalue\tbitscore\ttaxonomy\n"
        "S1\tS1\tGG2_0001\t99.0\t1e-10\t200\tk__Bacteria; p__Firmicutes\n"
    )
    meta = tmp_path / "meta.tsv"
    meta.write_text("sample_id\nS1\n")
    weights = tmp_path / "replicate_weights.tsv"
    weights.write_text("qseqid\treplicate_size\nS1\t3\n")

    biom_out = tmp_path / "out.biom"
    postblast_run(hits, meta, biom_out, weights_tsv=weights)

    with biom.util.biom_open(str(biom_out)) as fh:
        table = biom.Table.from_hdf5(fh)
    df = table.to_dataframe(dense=True)
    assert df.loc["k__Bacteria; p__Firmicutes", "S1"] == 3


def test_postblast_weights_sum_multiple_rows(tmp_path):
    hits = tmp_path / "hits_tax.tsv"
    hits.write_text(
        "qseqid\tsample_id\tsseqid\tpident\tevalue\tbitscore\ttaxonomy\n"
        "S1a\tS1\tGG2_0001\t99.0\t1e-10\t200\tk__Bacteria; p__Firmicutes\n"
        "S1b\tS1\tGG2_0001\t99.0\t1e-10\t200\tk__Bacteria; p__Firmicutes\n"
    )
    meta = tmp_path / "meta.tsv"
    meta.write_text("sample_id\nS1\n")
    weights = tmp_path / "replicate_weights.tsv"
    weights.write_text("qseqid\treplicate_size\nS1a\t2\nS1b\t4\n")

    biom_out = tmp_path / "out.biom"
    postblast_run(hits, meta, biom_out, weights_tsv=weights)

    with biom.util.biom_open(str(biom_out)) as fh:
        table = biom.Table.from_hdf5(fh)
    df = table.to_dataframe(dense=True)
    assert df.loc["k__Bacteria; p__Firmicutes", "S1"] == 6
