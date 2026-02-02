from __future__ import annotations

from pathlib import Path
import logging
import sys
import types

import pytest

"""
This test suite validates the 'microseq_tests' pipeline for Sanger sequencing 
assembly and quality control. Key functionalities tested include:

1.  **File Conversions**: Verifies `write_fasta_and_qual_from_fastq` correctly 
    splits FASTQ files into FASTA and QUAL pairs while enforcing length parity.
2.  **Assembly Merging**: Ensures `_write_combined_fasta` handles missing QUAL 
    files gracefully and logs appropriate warnings.
3.  **BLAST Input Generation**: Confirms `_build_blast_inputs` correctly labels 
    contigs, singlets, and missing pairs for downstream alignment.
4.  **Pairing Logic**: Validates `_collect_pairing_catalog` accurately identifies 
    samples with matching forward/reverse reads vs. those with missing mates.
5.  **Report Parsing**: Tests `parse_cap3_reports` for its ability to scrape 
    assembly metrics (overlaps, clipping, counts) from CAP3 output files.
6.  **Overlap Auditing**: Checks `_evaluate_overlap` for classification logic 
    based on overlap length and sequence identity.
"""

if "pandas" not in sys.modules:
    sys.modules["pandas"] = types.ModuleType("pandas")
if "biom" not in sys.modules:
    biom_stub = types.ModuleType("biom")
    biom_stub.Table = type("Table", (), {})
    sys.modules["biom"] = biom_stub
    biom_util_stub = types.ModuleType("biom.util")
    biom_util_stub.biom_open = lambda *args, **kwargs: None
    sys.modules["biom.util"] = biom_util_stub
Bio = pytest.importorskip("Bio")
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from microseq_tests.assembly.cap3_report import parse_cap3_reports
from microseq_tests.assembly.paired_assembly import _write_combined_fasta
from microseq_tests.pipeline import (
    _build_blast_inputs,
    _collect_pairing_catalog,
    _evaluate_overlap,
)
from microseq_tests.assembly.pairing import DupPolicy
from microseq_tests.utility import io_utils


def _write_fasta(path: Path, records: list[SeqRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, path, "fasta")


def _write_qual(path: Path, records: list[SeqRecord]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, path, "qual")


def test_write_fasta_and_qual_from_fastq(tmp_path: Path) -> None:
    fastq_path = tmp_path / "reads.fastq"
    fasta_path = tmp_path / "reads.fasta"
    records = [
        SeqRecord(Seq("ACGT"), id="r1", description=""),
        SeqRecord(Seq("TTGC"), id="r2", description=""),
    ]
    records[0].letter_annotations["phred_quality"] = [40, 40, 39, 38]
    records[1].letter_annotations["phred_quality"] = [30, 30, 30, 30]
    SeqIO.write(records, fastq_path, "fastq")

    result = io_utils.write_fasta_and_qual_from_fastq(fastq_path, fasta_path)
    assert result is not None

    fasta_records = list(SeqIO.parse(fasta_path, "fasta"))
    qual_records = list(SeqIO.parse(fasta_path.with_suffix(".fasta.qual"), "qual"))
    assert [rec.id for rec in fasta_records] == [rec.id for rec in qual_records]
    for fasta_rec, qual_rec in zip(fasta_records, qual_records):
        quals = qual_rec.letter_annotations["phred_quality"]
        assert len(fasta_rec.seq) == len(quals)


def test_write_fasta_and_qual_from_fastq_mismatch(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fastq_path = tmp_path / "reads.fastq"
    fasta_path = tmp_path / "reads.fasta"
    bad_record = SeqRecord(Seq("ACGT"), id="bad", description="")
    bad_record._per_letter_annotations = {"phred_quality": [30, 30]}

    def fake_parse(*_args, **_kwargs):
        return [bad_record]

    monkeypatch.setattr(io_utils.SeqIO, "parse", fake_parse)

    with pytest.raises(ValueError, match="Quality length mismatch"):
        io_utils.write_fasta_and_qual_from_fastq(fastq_path, fasta_path)


def test_write_combined_fasta_with_missing_qual(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    fasta_a = tmp_path / "a.fasta"
    fasta_b = tmp_path / "b.fasta"
    dest = tmp_path / "combined.fasta"

    rec_a = SeqRecord(Seq("ACGT"), id="a1", description="")
    rec_a.letter_annotations["phred_quality"] = [30, 30, 30, 30]
    rec_b = SeqRecord(Seq("TTTT"), id="b1", description="")
    rec_b.letter_annotations["phred_quality"] = [20, 20, 20, 20]

    _write_fasta(fasta_a, [rec_a])
    _write_qual(Path(f"{fasta_a}.qual"), [rec_a])
    _write_fasta(fasta_b, [rec_b])

    with caplog.at_level(logging.WARNING):
        _write_combined_fasta([fasta_a, fasta_b], dest, use_qual=True)

    assert "Missing QUAL file" in caplog.text

    combined_fasta = list(SeqIO.parse(dest, "fasta"))
    combined_qual = list(SeqIO.parse(dest.with_name(f"{dest.name}.qual"), "qual"))
    assert [rec.id for rec in combined_fasta] == ["a1", "b1"]
    assert [rec.id for rec in combined_qual] == ["a1"]


def test_build_blast_inputs_manifest(tmp_path: Path) -> None:
    asm_dir = tmp_path / "asm"
    output_fasta = tmp_path / "blast_inputs.fasta"
    output_tsv = tmp_path / "blast_inputs.tsv"

    sample1_dir = asm_dir / "sample1"
    sample2_dir = asm_dir / "sample2"
    sample1_dir.mkdir(parents=True)
    sample2_dir.mkdir(parents=True)

    contig = SeqRecord(Seq("ACGT"), id="contigA", description="")
    singlet = SeqRecord(Seq("TTTT"), id="singletA", description="")
    _write_fasta(sample1_dir / "sample1_paired.fasta.cap.contigs", [contig])
    _write_fasta(sample2_dir / "sample2_paired.fasta.cap.singlets", [singlet])

    paired_samples = {
        "sample1": {"F": [Path("f1")], "R": [Path("r1")]},
        "sample2": {"F": [Path("f2")], "R": [Path("r2")]},
    }
    missing_samples = {"sample3"}

    _build_blast_inputs(asm_dir, paired_samples, missing_samples, output_fasta, output_tsv)

    rows = {
        line.split("\t", 1)[0]: line.split("\t")
        for line in output_tsv.read_text(encoding="utf-8").splitlines()[1:]
    }
    assert rows["sample1"][1] == "contig"
    assert rows["sample2"][1] == "singlet"
    assert rows["sample3"][1] == "pair_missing"

    payload_ids = rows["sample1"][2].split(";")[0]
    assert payload_ids.startswith("sample1|contig|cap3_c1=")

    fasta_ids = [rec.id for rec in SeqIO.parse(output_fasta, "fasta")]
    assert "sample1|contig|cap3_c1" in fasta_ids
    assert "sample2|singlet|cap3_s1" in fasta_ids


def test_collect_pairing_catalog_missing_samples(tmp_path: Path) -> None:
    rec = SeqRecord(Seq("ACGT"), id="seq1", description="")
    rec.letter_annotations["phred_quality"] = [30, 30, 30, 30]
    _write_fasta(tmp_path / "S1_27F.fasta", [rec])
    _write_fasta(tmp_path / "S1_1492R.fasta", [rec])
    _write_fasta(tmp_path / "S2_27F.fasta", [rec])

    paired_samples, missing_samples = _collect_pairing_catalog(
        tmp_path,
        fwd_pattern="27F",
        rev_pattern="1492R",
        dup_policy=DupPolicy.ERROR,
        enforce_same_well=False,
        well_pattern=None,
    )

    assert "S1" in paired_samples
    assert paired_samples["S1"]["F"]
    assert paired_samples["S1"]["R"]
    assert "S2" in missing_samples


def test_parse_cap3_reports_and_missing_samples(tmp_path: Path) -> None:
    asm_dir = tmp_path / "asm"
    sample_dir = asm_dir / "sample1"
    sample_dir.mkdir(parents=True)

    info_path = sample_dir / "sample1_paired.fasta.cap.info"
    info_path.write_text(
        "\n".join(
            [
                "Number of overlaps saved: 3",
                "Number of overlaps removed: 1",
                "Clip sampleF left clip: 10, right clip: 120",
                "Clip sampleR left clip: 5, right clip: 110",
            ]
        ),
        encoding="utf-8",
    )

    contigs = [
        SeqRecord(Seq("ACGT"), id="c1", description=""),
        SeqRecord(Seq("TGCA"), id="c2", description=""),
    ]
    singlets = [SeqRecord(Seq("AAAA"), id="s1", description="")]
    _write_fasta(sample_dir / "sample1_paired.fasta.cap.contigs", contigs)
    _write_fasta(sample_dir / "sample1_paired.fasta.cap.singlets", singlets)

    rows = parse_cap3_reports(
        asm_dir,
        ["sample1", "missing1"],
        missing_samples={"missing1"},
    )

    sample_row = next(row for row in rows if row["sample_id"] == "sample1")
    missing_row = next(row for row in rows if row["sample_id"] == "missing1")

    assert sample_row["overlaps_saved"] == 3
    assert sample_row["overlaps_removed"] == 1
    assert sample_row["contigs_count"] == 2
    assert sample_row["singlets_count"] == 1
    assert sample_row["clip_f_left"] == 10
    assert sample_row["clip_f_right"] == 120
    assert sample_row["clip_r_left"] == 5
    assert sample_row["clip_r_right"] == 110

    assert missing_row["status"] == "pair_missing"
    assert missing_row["contigs_count"] == 0
    assert missing_row["singlets_count"] == 0


def test_overlap_audit_classification() -> None:
    short = _evaluate_overlap("A" * 10, "A" * 10, None, None, min_overlap=20)
    assert short["status"] == "overlap_too_short"

    identity_low = _evaluate_overlap(
        "A" * 10,
        "A" * 8 + "T" * 2,
        None,
        None,
        min_overlap=5,
        min_identity=1.1,
    )
    assert identity_low["status"] == "overlap_identity_low"

