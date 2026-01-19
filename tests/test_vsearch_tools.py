# tests/test_vsearch_tools.py

"""
This test suite verifies the functionality of vsearch-based sequence tools:

1. test_collapse_replicates_grouped: 
   Ensures sequences are grouped by a sample prefix and identical reads are 
   collapsed, updating the ';size=N' metadata correctly.

2. test_cluster_preserves_sizes: 
   Checks that near-identical sequences (based on an identity threshold) 
   are clustered together while accurately summing their total abundance.

3. test_chimera_check_ref_outputs_nonchimeras: 
   Validates that the chimera detection wrapper correctly identifies 
   non-chimeric sequences against a reference database and generates 
   the required report files.
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

import pytest

from Bio import SeqIO

from microseq_tests.vsearch_tools import (
    chimera_check_ref,
    collapse_replicates_grouped,
)


def _write_fasta(path: Path, records: dict[str, str]) -> None:
    path.write_text(
        "".join(f">{name}\n{seq}\n" for name, seq in records.items()),
        encoding="utf-8",
    )


@pytest.mark.skipif(shutil.which("vsearch") is None, reason="vsearch not installed")
def test_collapse_replicates_grouped(tmp_path):
    fasta_in = tmp_path / "reads.fasta"
    _write_fasta(
        fasta_in,
        {
            "SampleA_1": "ACGTACGT",
            "SampleA_2": "ACGTACGT",
            "SampleB_1": "ACGTACGT",
        },
    )
    fasta_out = tmp_path / "collapsed.fasta"
    collapse_replicates_grouped(
        fasta_in,
        fasta_out,
        group_fn=lambda value: value.split("_")[0],
        min_size=1,
    )

    records = list(SeqIO.parse(fasta_out, "fasta"))
    assert len(records) == 2
    sizes = []
    for rec in records:
        match = re.search(r";size=(\d+);?", rec.id)
        assert match
        sizes.append(int(match.group(1)))
    assert sorted(sizes) == [1, 2]


@pytest.mark.skipif(shutil.which("vsearch") is None, reason="vsearch not installed")
def test_cluster_preserves_sizes(tmp_path):
    fasta_in = tmp_path / "reads.fasta"
    _write_fasta(
        fasta_in,
        {
            "SampleA_1": "ACGTACGT",
            "SampleA_2": "ACGTACGT",
            "SampleA_3": "ACGTACGA",
        },
    )
    fasta_out = tmp_path / "collapsed.fasta"
    collapse_replicates_grouped(
        fasta_in,
        fasta_out,
        group_fn=lambda value: value.split("_")[0],
        min_size=1,
        id_th=0.99,
    )
    records = list(SeqIO.parse(fasta_out, "fasta"))
    assert len(records) == 1
    match = re.search(r";size=(\d+);?", records[0].id)
    assert match
    assert int(match.group(1)) == 3


@pytest.mark.skipif(shutil.which("vsearch") is None, reason="vsearch not installed")
def test_chimera_check_ref_outputs_nonchimeras(tmp_path):
    ref = tmp_path / "ref.fasta"
    query = tmp_path / "query.fasta"
    _write_fasta(ref, {"Ref1": "ACGTACGTACGT"})
    _write_fasta(query, {"Q1": "ACGTACGTACGT"})

    out = tmp_path / "nonchim.fasta"
    report = tmp_path / "uchime.tsv"
    nonchim, report = chimera_check_ref(query, out, reference=ref, report_tsv=report)

    assert nonchim.exists()
    assert report.exists()
    records = list(SeqIO.parse(nonchim, "fasta"))
    assert [rec.id for rec in records] == ["Q1"]
