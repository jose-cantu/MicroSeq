from __future__ import annotations

from pathlib import Path
import logging
import sys
import types

import pytest
import importlib

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
from microseq_tests.assembly.overlap_utils import (
    AlignedOverlapCandidate,
    OverlapCandidate,
    OverlapResult,
    best_pairwise_overlap,
    iter_end_anchored_overlaps,
    select_best_overlap,
)
from microseq_tests.assembly.paired_assembly import _write_combined_fasta, assemble_pairs
from microseq_tests.assembly.two_read_merge import merge_two_reads
from microseq_tests.assembly.overlap_backends import compute_overlap_candidates, resolve_overlap_engine
from microseq_tests.pipeline import (
    _build_blast_inputs,
    _build_selected_blast_inputs,
    _classify_overlap_status,
    _collect_pairing_catalog,
    _evaluate_overlap,
    _write_overlap_audit,
    run_compare_assemblers,
)
from microseq_tests.assembly.pairing import DupPolicy
from microseq_tests.assembly.registry import AssemblerSpec
from microseq_tests.utility import io_utils


def _write_fasta(path: Path, records: list[SeqRecord]) -> None:
    """Write SeqRecord entries to a FASTA file, creating parent directories as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, path, "fasta")


def _write_qual(path: Path, records: list[SeqRecord]) -> None:
    """Write SeqRecord entries to a QUAL file, creating parent directories as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(records, path, "qual")


def test_write_fasta_and_qual_from_fastq(tmp_path: Path) -> None:
    """Verify FASTQ conversion writes paired FASTA and QUAL records with matching lengths."""
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
    """Ensure FASTQ conversion raises when sequence and quality lengths differ."""
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
    """Confirm combined FASTA writing logs warnings and omits missing QUAL inputs safely."""
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
    """Validate BLAST input FASTA/TSV manifest labeling for contig preferred fallback, singlets, and missing pairs."""
    asm_dir = tmp_path / "asm"
    output_fasta = tmp_path / "blast_inputs.fasta"
    output_tsv = tmp_path / "blast_inputs.tsv"

    sample1_dir = asm_dir / "sample1"
    sample2_dir = asm_dir / "sample2"
    sample4_dir = asm_dir / "sample4" 
    sample1_dir.mkdir(parents=True)
    sample2_dir.mkdir(parents=True)
    sample4_dir.mkdir(parents=True) 

    contig = SeqRecord(Seq("ACGT"), id="contigA", description="")
    singlet = SeqRecord(Seq("TTTT"), id="singletA", description="")
    _write_fasta(sample1_dir / "sample1_paired.fasta.cap.contigs", [contig])
    _write_fasta(sample1_dir / "sample1_paired.fasta.cap.singlets", [singlet])
    _write_fasta(sample2_dir / "sample2_paired.fasta.cap.singlets", [singlet])
    _write_fasta(sample4_dir / "sample4_paired.fasta.cap.contigs", [contig])

    paired_samples = {
        "sample1": {"F": [Path("f1")], "R": [Path("r1")]},
        "sample2": {"F": [Path("f2")], "R": [Path("r2")]},
        "sample4": {"F": [Path("f4")], "R": [Path("r4")]},
    }
    missing_samples = {"sample3"}

    _build_blast_inputs(asm_dir, paired_samples, missing_samples, output_fasta, output_tsv)

    rows = {
        line.split(" ", 1)[0]: line.split(" ")
        for line in output_tsv.read_text(encoding="utf-8").splitlines()[1:]
    }
    assert rows["sample1"][1] == "contig"
    assert rows["sample2"][1] == "singlet"
    assert rows["sample4"][1] == "contig"
    assert rows["sample3"][1] == "pair_missing"

    sample1_payload_ids = rows["sample1"][2]
    assert "sample1|contig|cap3_c1=" in sample1_payload_ids
    assert "sample1|singlet|cap3_s1=" not in sample1_payload_ids

    fasta_ids = [rec.id for rec in SeqIO.parse(output_fasta, "fasta")]
    assert "sample1|contig|cap3_c1" in fasta_ids
    assert "sample1|singlet|cap3_s1" not in fasta_ids
    assert "sample2|singlet|cap3_s1" in fasta_ids
    assert "sample4|contig|cap3_c1" in fasta_ids


def test_collect_pairing_catalog_missing_samples(tmp_path: Path) -> None:
    """Check pairing catalog groups complete pairs and marks incomplete samples as missing."""
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
    """Verify CAP3 report parsing captures metrics and emits missing-sample placeholder rows."""
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
    assert missing_row["configured_engine"] == "auto"
    assert missing_row["fallback_used"] == "n/a"




def test_parse_cap3_reports_rejected_validation_marks_unverified(tmp_path: Path) -> None:
    """Ensure rejected CAP3 validation marks assembled output as unverified."""
    asm_dir = tmp_path / "asm"
    sample_dir = asm_dir / "sample1"
    sample_dir.mkdir(parents=True)

    _write_fasta(sample_dir / "sample1_paired.fasta.cap.contigs", [SeqRecord(Seq("ACGT"), id="c1", description="")])
    (sample_dir / "sample1_paired.cap3_validation.txt").write_text("rejected\n", encoding="utf-8")

    rows = parse_cap3_reports(asm_dir, ["sample1"])
    row = rows[0]
    assert row["cap3_validation"] == "rejected"
    assert row["status"] == "cap3_unverified"


def test_parse_cap3_reports_unknown_validation_keeps_assembled_status(tmp_path: Path) -> None:
    """Ensure unknown CAP3 validation preserves assembled status."""
    asm_dir = tmp_path / "asm"
    sample_dir = asm_dir / "sample1"
    sample_dir.mkdir(parents=True)

    _write_fasta(sample_dir / "sample1_paired.fasta.cap.contigs", [SeqRecord(Seq("ACGT"), id="c1", description="")])
    (sample_dir / "sample1_paired.cap3_validation.txt").write_text("unknown\n", encoding="utf-8")

    rows = parse_cap3_reports(asm_dir, ["sample1"])
    row = rows[0]
    assert row["cap3_validation"] == "unknown"
    assert row["status"] == "assembled"

def test_overlap_audit_classification() -> None:
    """Confirm overlap audit classification handles short and low-identity cases."""
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




def test_overlap_audit_writes_best_identity_columns(tmp_path: Path) -> None:
    """Ensure overlap audit TSV includes and populates best-identity helper columns."""
    fwd = tmp_path / "S1_27F.fasta"
    rev = tmp_path / "S1_1492R.fasta"
    fwd.write_text(">f\n" + "A" * 120 + "\n", encoding="utf-8")
    rev.write_text(">r\n" + "A" * 119 + "T\n", encoding="utf-8")

    paired = {"S1": {"F": [fwd], "R": [rev]}}
    out = tmp_path / "overlap_audit.tsv"
    _write_overlap_audit(paired, out, min_overlap=50, min_identity=0.5, min_quality=20.0, quality_mode="warning")

    lines = out.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    values = lines[1].split("\t")
    assert "best_identity" in header
    assert "best_identity_orientation" in header
    assert "fwd_best_identity" in header
    assert "revcomp_best_identity" in header
    assert "top2_identity_delta" in header
    assert "pretrim_best_identity" in header
    assert "posttrim_selected_overlap_len" in header
    assert "fwd_best_identity_any" in header
    assert "ambiguity_identity_delta_used" in header
    row = dict(zip(header, values))
    assert float(row["best_identity"]) >= float(row["overlap_identity"])
    assert row["posttrim_status"] == row["status"]
    assert row["posttrim_selected_overlap_len"] == row["overlap_len"]


def test_overlap_audit_missing_reads_rows_have_expected_columns(tmp_path: Path) -> None:
    """Verify missing-read audit rows preserve the expected fixed column width."""
    out = tmp_path / "overlap_audit.tsv"
    paired = {"S1": {"F": [], "R": []}}
    _write_overlap_audit(paired, out)
    header = out.read_text(encoding="utf-8").splitlines()[0].split("\t")
    values = out.read_text(encoding="utf-8").splitlines()[1].split("\t")
    assert len(header) == 40
    assert len(values) == 40





def test_overlap_audit_exposes_any_orientation_metrics_for_short_high_identity_revcomp(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure *_any columns preserve short high-identity revcomp evidence below min_overlap."""
    from microseq_tests.assembly.overlap_backends import OverlapBackendResult

    fwd = tmp_path / "S2_27F.fasta"
    rev = tmp_path / "S2_1492R.fasta"
    fwd.write_text(">f\n" + "A" * 120 + "\n", encoding="utf-8")
    rev.write_text(">r\n" + "A" * 120 + "\n", encoding="utf-8")

    monkeypatch.setattr(
        "microseq_tests.pipeline.load_config",
        lambda: {
            "merge_two_reads": {
                "overlap_engine": "biopython",
                "overlap_engine_strategy": "single",
                "overlap_engine_order": ["biopython"],
                "anchor_tolerance_bases": 30,
                "min_overlap": 50,
                "min_identity": 0.5,
            },
            "overlap_eval": {
                "min_overlap": 50,
                "min_identity": 0.5,
                "min_quality": 20.0,
                "ambiguity_identity_delta": 0.003,
                "ambiguity_quality_epsilon": 0.2,
            },
        },
    )

    def _fake_backend(*_args, **_kwargs):
        return [
            OverlapBackendResult("forward", "A" * 80, "A" * 80, 80, 0.62, 30, 0, 30.0, "", 80, 0, 0, 0, 80, 0, 80, True),
            OverlapBackendResult("revcomp", "A" * 30, "A" * 30, 30, 0.98, 1, 0, 30.0, "", 30, 0, 0, 0, 30, 0, 30, True),
        ]

    monkeypatch.setattr("microseq_tests.pipeline.compute_overlap_candidates", _fake_backend)

    out = tmp_path / "overlap_audit.tsv"
    _write_overlap_audit({"S2": {"F": [fwd], "R": [rev]}}, out, min_overlap=50, min_identity=0.5, min_quality=20.0)

    header, values = [line.split("	") for line in out.read_text(encoding="utf-8").splitlines()]
    row = dict(zip(header, values))
    assert row["revcomp_best_identity"] == ""
    assert row["revcomp_best_identity_any"] == "0.9800"
    assert row["revcomp_best_overlap_len_any"] == "30"
    assert row["fwd_best_identity"] == "0.6200"
    assert row["ambiguity_identity_delta_used"] == "0.0030"
    assert row["ambiguity_quality_epsilon_used"] == "0.2000"

def test_overlap_audit_emit_engine_audit_false_has_no_engines_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fwd = tmp_path / "S1_27F.fasta"
    rev = tmp_path / "S1_1492R.fasta"
    fwd.write_text(">f\n" + "A" * 120 + "\n", encoding="utf-8")
    rev.write_text(">r\n" + "A" * 120 + "\n", encoding="utf-8")

    monkeypatch.setattr(
        "microseq_tests.pipeline.load_config",
        lambda: {
            "merge_two_reads": {
                "overlap_engine": "biopython",
                "overlap_engine_strategy": "single",
                "overlap_engine_order": ["biopython", "ungapped", "edlib"],
                "anchor_tolerance_bases": 30,
                "min_overlap": 50,
                "min_identity": 0.5,
            },
            "overlap_eval": {"min_overlap": 50, "min_identity": 0.5, "min_quality": 20.0},
        },
    )

    out = tmp_path / "overlap_audit.tsv"
    _write_overlap_audit(
        {"S1": {"F": [fwd], "R": [rev]}},
        out,
        min_overlap=50,
        min_identity=0.5,
        min_quality=20.0,
        emit_engine_audit=False,
    )
    assert out.exists()
    assert not out.with_name("overlap_audit_engines.tsv").exists()


def test_run_compare_assemblers_requires_staged_fasta(tmp_path: Path) -> None:
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "S1_27F.fastq").write_text("@r\nACGT\n+\n!!!!\n", encoding="utf-8")
    (raw / "S1_1492R.fastq").write_text("@r\nACGT\n+\n!!!!\n", encoding="utf-8")

    with pytest.raises(ValueError, match="pre-staged paired FASTA"):
        run_compare_assemblers(raw, tmp_path)


def test_run_compare_assemblers_merge_specs_use_single_engine(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    in_dir = tmp_path / "paired_fasta"
    in_dir.mkdir()
    (in_dir / "S1_27F.fasta").write_text(">f\n" + "A" * 120 + "\n", encoding="utf-8")
    (in_dir / "S1_1492R.fasta").write_text(">r\n" + "A" * 120 + "\n", encoding="utf-8")

    specs = (
        AssemblerSpec(id="merge_two_reads:biopython", display_name="b", kind="merge_two_reads", overlap_engine="biopython"),
        AssemblerSpec(id="merge_two_reads:ungapped", display_name="u", kind="merge_two_reads", overlap_engine="ungapped"),
        AssemblerSpec(id="merge_two_reads:edlib", display_name="e", kind="merge_two_reads", overlap_engine="edlib"),
    )
    monkeypatch.setattr("microseq_tests.pipeline.list_assemblers", lambda: specs)

    calls: list[tuple[str, str]] = []

    class FakeReport:
        def __init__(self, overlap_engine: str):
            self.merge_status = "merged"
            self.overlap_engine = overlap_engine
            self.merge_warning = ""
            self.overlap_len = 120
            self.identity = 1.0
            self.contig_len = 120
            self.mismatches = 0

    def fake_merge_two_reads(*, sample_id, fwd_path, rev_path, output_dir, overlap_engine, overlap_engine_strategy, **_kwargs):
        calls.append((overlap_engine, overlap_engine_strategy))
        contig = output_dir / f"{sample_id}.contig.fasta"
        contig.write_text(">c\n" + "A" * 120 + "\n", encoding="utf-8")
        return contig, FakeReport(overlap_engine)

    monkeypatch.setattr("microseq_tests.pipeline.merge_two_reads", fake_merge_two_reads)
    out_tsv = run_compare_assemblers(in_dir, tmp_path, fwd_pattern="27F", rev_pattern="1492R")

    assert out_tsv.exists()
    out_text = out_tsv.read_text(encoding="utf-8")
    assert "payload_kind" in out_text
    assert "payload_n" in out_text
    assert "decision_source" in out_text
    assert "merge_two_reads:biopython" in out_text
    assert "merge_two_reads:ungapped" in out_text
    assert "merge_two_reads:edlib" in out_text
    assert all(strategy == "single" for _engine, strategy in calls)
    assert {engine for engine, _strategy in calls} == {"biopython", "ungapped", "edlib"}

    compare_root = tmp_path / "asm" / "compare"
    assert (compare_root / "merge_two_reads_biopython" / "S1").exists()
    assert (compare_root / "merge_two_reads_ungapped" / "S1").exists()
    assert (compare_root / "merge_two_reads_edlib" / "S1").exists()

def test_run_compare_assemblers_selected_ids_only(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Selected assembler compare mode should emit rows only for requested assembler IDs."""
    in_dir = tmp_path / "paired_fasta"
    in_dir.mkdir(parents=True)

    (in_dir / "S1_27F.fasta").write_text(">S1_27F\nACGT\n", encoding="utf-8")
    (in_dir / "S1_1492R.fasta").write_text(">S1_1492R\nACGT\n", encoding="utf-8")

    class _DummyReport:
        status = "assembled"
        warning = ""
        selected_engine = "biopython"
        diag_detail_for_human = "ok"
        all_hypotheses = []

    def fake_merge_two_reads(*, sample_id, output_dir, overlap_engine, **_kwargs):
        contig = output_dir / f"{sample_id}.{overlap_engine}.fasta"
        contig.write_text(f">{sample_id}_{overlap_engine}\nACGT\n", encoding="utf-8")
        return contig, _DummyReport()

    monkeypatch.setattr("microseq_tests.pipeline.merge_two_reads", fake_merge_two_reads)

    out_tsv = run_compare_assemblers(
        in_dir,
        tmp_path,
        assembler_ids=["merge_two_reads:biopython"],
        fwd_pattern="27F",
        rev_pattern="1492R",
    )

    rows = out_tsv.read_text(encoding="utf-8").splitlines()
    assert len(rows) == 2
    assert "merge_two_reads:biopython" in rows[1]
    assert "merge_two_reads:edlib" not in "\n".join(rows)

def test_run_compare_assemblers_ambiguous_overlap_reports_zero_contig_len(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    in_dir = tmp_path / "paired_fasta"
    in_dir.mkdir()
    (in_dir / "S1_27F.fasta").write_text(">f\n" + "A" * 120 + "\n", encoding="utf-8")
    (in_dir / "S1_1492R.fasta").write_text(">r\n" + "A" * 120 + "\n", encoding="utf-8")

    specs = (
        AssemblerSpec(id="merge_two_reads:biopython", display_name="b", kind="merge_two_reads", overlap_engine="biopython"),
    )
    monkeypatch.setattr("microseq_tests.pipeline.list_assemblers", lambda: specs)

    class FakeReport:
        merge_status = "ambiguous_overlap"
        overlap_engine = "biopython"
        merge_warning = ""
        overlap_len = 80
        identity = 0.99
        contig_len = 0

    def fake_merge_two_reads(*, sample_id, output_dir, **_kwargs):
        return None, FakeReport()

    monkeypatch.setattr("microseq_tests.pipeline.merge_two_reads", fake_merge_two_reads)
    out_tsv = run_compare_assemblers(in_dir, tmp_path, fwd_pattern="27F", rev_pattern="1492R")

    text = out_tsv.read_text(encoding="utf-8")
    assert "ambiguous_overlap" in text
    assert "	biopython	0	" in text
    assert "Overlap candidates were tied; no unique contig selected" in text

def test_run_compare_assemblers_cap3_nonzero_exit_still_records_artifacts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    in_dir = tmp_path / "paired_fasta"
    in_dir.mkdir()
    (in_dir / "S1_27F.fasta").write_text(">f\n" + "A" * 120 + "\n", encoding="utf-8")
    (in_dir / "S1_1492R.fasta").write_text(">r\n" + "A" * 120 + "\n", encoding="utf-8")

    specs = (
        AssemblerSpec(id="cap3:strict", display_name="CAP3 strict", kind="cap3", cap3_profile="strict"),
    )
    monkeypatch.setattr("microseq_tests.pipeline.list_assemblers", lambda: specs)

    class FakeCompleted:
        def __init__(self) -> None:
            self.returncode = 1
            self.stdout = "CAP3 stdout"
            self.stderr = "CAP3 stderr"

    def fake_run(cmd, cwd=None, capture_output=None, text=None, check=None):
        assert check is False
        sample_fasta = Path(cwd) / cmd[1]
        contigs = sample_fasta.with_suffix(sample_fasta.suffix + ".cap.contigs")
        singlets = sample_fasta.with_suffix(sample_fasta.suffix + ".cap.singlets")
        info = sample_fasta.with_suffix(sample_fasta.suffix + ".cap.info")
        contigs.write_text(">c1\n" + "A" * 120 + "\n", encoding="utf-8")
        singlets.write_text(">s1\n" + "A" * 60 + "\n", encoding="utf-8")
        info.write_text("cap info\n", encoding="utf-8")
        return FakeCompleted()

    monkeypatch.setattr("subprocess.run", fake_run)

    out_tsv = run_compare_assemblers(in_dir, tmp_path, fwd_pattern="27F", rev_pattern="1492R")
    text = out_tsv.read_text(encoding="utf-8")

    assert "cap3_nonzero_exit_with_output" in text
    assert "returncode=1; profile=strict" in text
    assert "	assembled	cap3	120	" in text
    assert "logs/cap3/S1__cap3_strict.stdout.log" in text
    assert "logs/cap3/S1__cap3_strict.stderr.log" in text 

def test_assemble_pairs_writes_cap3_process_logs_on_success(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    in_dir = tmp_path / "paired_fasta"
    out_dir = tmp_path / "asm"
    in_dir.mkdir()
    (in_dir / "S1_27F.fasta").write_text(">f\n" + "A" * 120 + "\n", encoding="utf-8")
    (in_dir / "S1_1492R.fasta").write_text(">r\n" + "A" * 120 + "\n", encoding="utf-8")

    monkeypatch.setattr(
        "microseq_tests.assembly.paired_assembly.load_config",
        lambda: {
            "tools": {"cap3": "cap3"},
            "overlap_eval": {},
            "merge_two_reads": {"high_conflict_action": "route_cap3"},
        },
    )

    class FakeMergeReport:
        merge_status = "identity_low"
        overlap_engine = "biopython"
        orientation = "revcomp"
        overlap_len = 20
        identity = 0.5
        merge_warning = "identity_low"

    monkeypatch.setattr(
        "microseq_tests.assembly.paired_assembly.merge_two_reads",
        lambda **_kwargs: (None, FakeMergeReport()),
    )

    class FakeCompleted:
        def __init__(self, *, stdout: str = "", stderr: str = "", returncode: int = 0) -> None:
            self.stdout = stdout
            self.stderr = stderr
            self.returncode = returncode

    def fake_run(cmd, check=None, capture_output=None, text=None, cwd=None):
        if "--version" in cmd:
            return FakeCompleted(stdout="CAP3 1.0")
        sample_fasta = Path(cwd) / cmd[1]
        contigs = sample_fasta.with_suffix(sample_fasta.suffix + ".cap.contigs")
        contigs.write_text(">c1\n" + "A" * 120 + "\n", encoding="utf-8")
        return FakeCompleted(stdout="DETAILED DISPLAY OF CONTIGS", stderr="cap3 warning")

    monkeypatch.setattr("subprocess.run", fake_run)

    contigs = assemble_pairs(in_dir, out_dir, fwd_pattern="27F", rev_pattern="1492R", use_qual=False)
    assert contigs

    stdout_log = out_dir / "logs" / "cap3" / "S1__cap3_default.stdout.log"
    stderr_log = out_dir / "logs" / "cap3" / "S1__cap3_default.stderr.log"
    assert stdout_log.exists()
    assert stderr_log.exists()
    assert "# sample_id: S1" in stdout_log.read_text(encoding="utf-8")
    assert "DETAILED DISPLAY OF CONTIGS" in stdout_log.read_text(encoding="utf-8")


def test_assemble_pairs_writes_cap3_process_logs_on_failure(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    in_dir = tmp_path / "paired_fasta"
    out_dir = tmp_path / "asm"
    in_dir.mkdir()
    (in_dir / "S1_27F.fasta").write_text(">f\n" + "A" * 120 + "\n", encoding="utf-8")
    (in_dir / "S1_1492R.fasta").write_text(">r\n" + "A" * 120 + "\n", encoding="utf-8")

    monkeypatch.setattr(
        "microseq_tests.assembly.paired_assembly.load_config",
        lambda: {
            "tools": {"cap3": "cap3"},
            "overlap_eval": {},
            "merge_two_reads": {"high_conflict_action": "route_cap3"},
        },
    )

    class FakeMergeReport:
        merge_status = "identity_low"
        overlap_engine = "biopython"
        orientation = "revcomp"
        overlap_len = 20
        identity = 0.5
        merge_warning = "identity_low"

    monkeypatch.setattr(
        "microseq_tests.assembly.paired_assembly.merge_two_reads",
        lambda **_kwargs: (None, FakeMergeReport()),
    )

    class FakeCompleted:
        def __init__(self, *, stdout: str = "", stderr: str = "", returncode: int = 0) -> None:
            self.stdout = stdout
            self.stderr = stderr
            self.returncode = returncode

    def fake_run(cmd, check=None, capture_output=None, text=None, cwd=None):
        if "--version" in cmd:
            return FakeCompleted(stdout="CAP3 1.0")
        raise subprocess.CalledProcessError(returncode=2, cmd=cmd, output="CAP3 FAIL OUT", stderr="CAP3 FAIL ERR")

    import subprocess

    monkeypatch.setattr("subprocess.run", fake_run)

    with pytest.raises(subprocess.CalledProcessError):
        assemble_pairs(in_dir, out_dir, fwd_pattern="27F", rev_pattern="1492R", use_qual=False)

    stdout_log = out_dir / "logs" / "cap3" / "S1__cap3_default.stdout.log"
    stderr_log = out_dir / "logs" / "cap3" / "S1__cap3_default.stderr.log"
    assert stdout_log.exists()
    assert stderr_log.exists()
    assert "CAP3 FAIL OUT" in stdout_log.read_text(encoding="utf-8")
    assert "CAP3 FAIL ERR" in stderr_log.read_text(encoding="utf-8")





def test_overlap_helper_prefers_end_overlap_over_internal_match() -> None:
    """Ensure overlap helper does not accept internal motif matches as valid terminal overlaps."""
    fwd = "A" * 60 + "GATTACA" + "C" * 60
    rev = "N" * 40 + "GATTACA" + "N" * 40

    overlap = best_pairwise_overlap(fwd, rev)
    status = _classify_overlap_status(
        overlap.overlap_len,
        overlap.identity,
        overlap.overlap_quality,
        min_overlap=30,
        min_identity=0.9,
        min_quality=20.0,
        quality_mode="warning",
    )

    assert status in {"overlap_too_short", "overlap_identity_low"}


def test_gapped_backend_revcomp_indel_passes_where_ungapped_fails() -> None:
    # revcomp(seq_rev) has a 1bp internal deletion relative to fwd.
    """Check gapped backend can recover a revcomp overlap with an internal indel missed by ungapped mode."""
    fwd = "CAGATTTTCATATTATGCAGAAAATCTACT"
    rev = str(Seq("CAGATTTTCATATTATGCAGAAATCTACT").reverse_complement())

    ungapped = iter_end_anchored_overlaps(fwd, rev)
    chosen_ungapped = select_best_overlap(
        ungapped,
        min_overlap=15,
        min_identity=0.95,
        min_quality=20.0,
        quality_mode="warning",
    )
    assert chosen_ungapped.identity < 0.95

    gapped = compute_overlap_candidates(fwd, rev, None, None, engine="biopython")
    chosen_gapped = select_best_overlap(
        [
            AlignedOverlapCandidate(
                orientation=c.orientation,
                overlap_len=c.overlap_len,
                mismatches=c.mismatches,
                identity=c.identity,
                overlap_quality=c.overlap_quality,
                aligned_fwd=c.aligned_fwd,
                aligned_rev=c.aligned_rev,
                indels=c.indels,
                cigar=c.cigar,
            )
            for c in gapped
        ],
        min_overlap=15,
        min_identity=0.85,
        min_quality=20.0,
        quality_mode="warning",
    )
    assert chosen_gapped.orientation == "revcomp"
    assert chosen_gapped.identity >= 0.85


def test_backend_parity_biopython_vs_edlib_if_available() -> None:
    """Compare best-candidate orientation and approximate identity between Biopython and Edlib backends."""
    edlib = pytest.importorskip("edlib")
    assert edlib is not None
    fwd = "AAACCCGGGTTT"
    rev = str(Seq("AAACCCGGGTT").reverse_complement())
    b = compute_overlap_candidates(fwd, rev, None, None, engine="biopython")
    e = compute_overlap_candidates(fwd, rev, None, None, engine="edlib")
    b_best = max(b, key=lambda c: (c.identity, c.overlap_len))
    e_best = max(e, key=lambda c: (c.identity, c.overlap_len))
    assert b_best.orientation == e_best.orientation
    assert abs(b_best.identity - e_best.identity) < 0.15


def test_resolve_overlap_engine_auto_prefers_edlib_when_available() -> None:
    """Verify auto backend selection prefers Edlib when the dependency is importable."""
    edlib = pytest.importorskip("edlib")
    assert edlib is not None
    assert resolve_overlap_engine("auto") == "edlib"


def test_overlap_audit_uses_not_end_anchored_status(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure overlap audit marks non-terminal matches as not_end_anchored."""
    motif = "GATTACAGATTACA"
    fwd = tmp_path / "S4_27F.fasta"
    rev = tmp_path / "S4_1492R.fasta"
    fwd.write_text(f">f\nAAAA{motif}CCCC\n", encoding="utf-8")
    rev_oriented = f"TTTT{motif}GGGG"
    rev.write_text(f">r\n{str(Seq(rev_oriented).reverse_complement())}\n", encoding="utf-8")

    paired = {"S4": {"F": [fwd], "R": [rev]}}
    out = tmp_path / "overlap_audit.tsv"
    monkeypatch.setattr(
        "microseq_tests.pipeline.load_config",
        lambda: {
            "merge_two_reads": {"overlap_engine": "biopython", "anchor_tolerance_bases": 2, "min_overlap": 100, "min_identity": 0.8},
            "overlap_eval": {"min_overlap": 100, "min_identity": 0.8, "min_quality": 20.0},
        },
    )
    _write_overlap_audit(paired, out, min_overlap=10, min_identity=0.80, min_quality=20.0, quality_mode="warning")
    row = dict(zip(out.read_text(encoding="utf-8").splitlines()[0].split("\t"), out.read_text(encoding="utf-8").splitlines()[1].split("\t")))
    assert row["status"] in {"not_end_anchored", "overlap_too_short", "overlap_identity_low"}


def test_edlib_geometry_union_if_available() -> None:
    """Validate Edlib candidates retain full union geometry when flanking overhangs are present."""
    edlib = pytest.importorskip("edlib")
    assert edlib is not None
    overlap = "GGGTTTCCCAAAGGG"
    fwd = "AAAA" + overlap
    rev_oriented = overlap + "TTTT"
    rev = str(Seq(rev_oriented).reverse_complement())
    candidates = compute_overlap_candidates(fwd, rev, None, None, engine="edlib", end_anchor_tolerance=5)
    best = max(candidates, key=lambda c: (c.identity, c.overlap_len, c.end_anchored))
    assert best.orientation == "revcomp"
    assert best.overlap_len == len(overlap)
    assert best.identity >= 0.95
    assert "AAAA" in best.aligned_fwd.replace("-", "")
    assert "TTTT" in best.aligned_rev.replace("-", "")


def test_edlib_cigar_indel_semantics_via_stub(monkeypatch: pytest.MonkeyPatch) -> None:
    """Regression test for Edlib CIGAR I/D semantics using a stubbed Edlib module."""
    class FakeEdlib:
        @staticmethod
        def align(_query, _target, mode="HW", task="path"):
            assert mode == "HW"
            assert task == "path"
            # one match, insertion-to-target, one match, deletion-from-target, one match
            return {"cigar": "1=1I1=1D1=", "locations": [(0, 4), (1, 4)]}

    real_import = importlib.import_module

    def fake_import(name, package=None):
        if name == "edlib":
            return FakeEdlib
        return real_import(name, package)

    monkeypatch.setattr(importlib, "import_module", fake_import)
    monkeypatch.setattr("builtins.__import__", __import__)

    import sys
    sys.modules["edlib"] = FakeEdlib
    try:
        cands = compute_overlap_candidates("ABCDE", "ABCDE", None, None, engine="edlib", end_anchor_tolerance=2)
    finally:
        sys.modules.pop("edlib", None)

    # path-only safe mode: single location candidate per orientation
    assert len(cands) == 2
    fwd = next(c for c in cands if c.orientation == "forward")
    assert fwd.aligned_fwd.count("-") >= 1
    assert fwd.aligned_rev.count("-") >= 1


def test_merge_report_includes_overlap_engine_column(tmp_path: Path) -> None:
    """Ensure merge report output includes the overlap_engine metadata column."""
    fwd_path = tmp_path / "S5_27F.fasta"
    rev_path = tmp_path / "S5_1492R.fasta"
    fwd_path.write_text(">f\n" + "A" * 120 + "\n", encoding="utf-8")
    rev_path.write_text(">r\n" + "A" * 120 + "\n", encoding="utf-8")

    merge_two_reads(
        sample_id="S5",
        fwd_path=fwd_path,
        rev_path=rev_path,
        output_dir=tmp_path,
        min_overlap=100,
        min_identity=0.8,
        min_quality=20.0,
        quality_mode="warning",
        overlap_engine="auto",
    )

    report_path = tmp_path / "S5_paired.merge_report.tsv"
    header = report_path.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "overlap_engine" in header


def test_engine_matrix_invariants_for_simple_pair(tmp_path: Path) -> None:
    """Check merge invariants across overlap engines for a simple known-overlap read pair."""
    overlap = "GGGTTTCCCAAAGGG"
    fwd_seq = "AAAA" + overlap
    rev_oriented = overlap + "TTTT"
    rev_seq = str(Seq(rev_oriented).reverse_complement())

    expected = "AAAA" + overlap + "TTTT"
    engines = ["ungapped", "biopython", "auto"]
    if importlib.util.find_spec("edlib") is not None:
        engines.append("edlib")

    for eng in engines:
        d = tmp_path / eng
        d.mkdir()
        fwd_path = d / "S6_27F.fasta"
        rev_path = d / "S6_1492R.fasta"
        fwd_path.write_text(f">f\n{fwd_seq}\n", encoding="utf-8")
        rev_path.write_text(f">r\n{rev_seq}\n", encoding="utf-8")

        contig_path, report = merge_two_reads(
            sample_id="S6",
            fwd_path=fwd_path,
            rev_path=rev_path,
            output_dir=d,
            min_overlap=10,
            min_identity=0.8,
            min_quality=5.0,
            quality_mode="warning",
            overlap_engine=eng,
            anchor_tolerance_bases=5,
        )
        assert report.merge_status in {"merged", "ambiguous_overlap", "not_end_anchored", "identity_low", "overlap_too_short"}
        if report.merge_status == "merged":
            assert contig_path is not None
            merged = str(next(SeqIO.parse(contig_path, "fasta")).seq)
            assert merged.count("AAAA") == 1
            assert merged.count("TTTT") == 1
            assert len(merged) == len(fwd_seq) + len(rev_oriented) - report.overlap_len
            assert merged == expected


def test_gapped_overlap_identity_not_penalized_by_terminal_flanks() -> None:
    """Verify terminal flanks do not reduce identity for gapped overlap scoring."""
    overlap = "GGGTTTCCCAAAGGG"
    fwd = "AAAA" + overlap
    rev_oriented = overlap + "TTTT"
    rev = str(Seq(rev_oriented).reverse_complement())

    candidates = compute_overlap_candidates(fwd, rev, None, None, engine="biopython", end_anchor_tolerance=5)
    best = max(candidates, key=lambda c: (c.identity, c.overlap_len))

    assert best.orientation == "revcomp"
    assert best.overlap_len == len(overlap)
    assert best.identity == pytest.approx(1.0)
    assert best.terminal_gap_cols > 0
    assert best.end_anchored is True


def test_merge_two_reads_gapped_preserves_union_geometry(tmp_path: Path) -> None:
    """Ensure gapped merging preserves expected union sequence geometry."""
    overlap = "GGGTTTCCCAAAGGG"
    fwd_seq = "AAAA" + overlap
    rev_oriented = overlap + "TTTT"
    rev_seq = str(Seq(rev_oriented).reverse_complement())

    fwd_path = tmp_path / "S2_27F.fasta"
    rev_path = tmp_path / "S2_1492R.fasta"
    fwd_path.write_text(f">f\n{fwd_seq}\n", encoding="utf-8")
    rev_path.write_text(f">r\n{rev_seq}\n", encoding="utf-8")

    contig_path, report = merge_two_reads(
        sample_id="S2",
        fwd_path=fwd_path,
        rev_path=rev_path,
        output_dir=tmp_path,
        min_overlap=10,
        min_identity=0.95,
        min_quality=5.0,
        quality_mode="warning",
        overlap_engine="biopython",
        anchor_tolerance_bases=5,
    )

    assert contig_path is not None
    merged = next(SeqIO.parse(contig_path, "fasta"))
    assert str(merged.seq) == ("AAAA" + overlap + "TTTT")
    assert report.overlap_len == len(overlap)


def test_end_anchoring_guard_rejects_internal_repeat_match(tmp_path: Path) -> None:
    """Confirm end-anchoring guard rejects internal repeat-only alignments."""
    motif = "GATTACAGATTACA"
    fwd_seq = "AAAA" + motif + "CCCC"
    rev_oriented = "TTTT" + motif + "GGGG"
    rev_seq = str(Seq(rev_oriented).reverse_complement())

    fwd_path = tmp_path / "S3_27F.fasta"
    rev_path = tmp_path / "S3_1492R.fasta"
    fwd_path.write_text(f">f\n{fwd_seq}\n", encoding="utf-8")
    rev_path.write_text(f">r\n{rev_seq}\n", encoding="utf-8")

    contig_path, report = merge_two_reads(
        sample_id="S3",
        fwd_path=fwd_path,
        rev_path=rev_path,
        output_dir=tmp_path,
        min_overlap=10,
        min_identity=0.90,
        min_quality=5.0,
        quality_mode="warning",
        overlap_engine="biopython",
        anchor_tolerance_bases=2,
    )

    assert contig_path is None
    assert report.merge_status in {"not_end_anchored", "identity_low", "overlap_too_short"}


def test_merge_two_reads_ambiguous_topk_emits_alternative_payloads(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    fwd_path = tmp_path / "S4_27F.fasta"
    rev_path = tmp_path / "S4_1492R.fasta"
    fwd_path.write_text(">f\n" + "A" * 120 + "\n", encoding="utf-8")
    rev_path.write_text(">r\n" + "T" * 120 + "\n", encoding="utf-8")

    from microseq_tests.assembly.overlap_utils import AlignedOverlapCandidate, OverlapResult

    fake_candidates = [
        AlignedOverlapCandidate("revcomp", 60, 1, 0.98, 30.0, "A" * 120, "A" * 120),
        AlignedOverlapCandidate("revcomp", 60, 1, 0.98, 30.0, "A" * 110 + "C" * 10, "A" * 110 + "G" * 10),
        AlignedOverlapCandidate("revcomp", 59, 1, 0.97, 29.0, "A" * 100 + "C" * 20, "A" * 100 + "G" * 20),
    ]

    monkeypatch.setattr(
        "microseq_tests.assembly.two_read_merge.select_best_overlap",
        lambda *_a, **_k: OverlapResult("revcomp", 60, 0.98, 1, 30.0, "A" * 120, "A" * 120, "ambiguous_overlap"),
    )
    monkeypatch.setattr(
        "microseq_tests.assembly.two_read_merge.rank_feasible_overlaps",
        lambda *_a, **_k: fake_candidates,
    )

    contig_path, report = merge_two_reads(
        sample_id="S4",
        fwd_path=fwd_path,
        rev_path=rev_path,
        output_dir=tmp_path,
        min_overlap=40,
        min_identity=0.9,
        min_quality=20.0,
        quality_mode="warning",
        overlap_engine="biopython",
        ambiguous_policy="topk",
        ambiguous_top_k=3,
    )

    assert contig_path is not None
    assert report.merge_status == "ambiguous_topk"
    assert "policy=topk" in report.merge_warning
    alt_records = list(SeqIO.parse(contig_path, "fasta"))
    assert len(alt_records) == 3

def test_overlap_selection_prefers_long_feasible_candidate() -> None:
    """Ensure overlap selection prefers the longest feasible candidate over shorter alternatives."""
    fwd = "A" * 55 + "C" * 44 + "G"
    rev = "A" + "A" * 55 + "C" * 44

    candidates = iter_end_anchored_overlaps(fwd, rev)
    assert len(candidates) >= 2

    chosen = select_best_overlap(
        candidates,
        min_overlap=60,
        min_identity=0.9,
        min_quality=20.0,
        quality_mode="warning",
    )

    assert chosen.overlap_len == 100
    assert chosen.identity == pytest.approx(0.98)


def test_merge_two_reads_sets_qualities_absent_warning_and_cap_info(tmp_path: Path) -> None:
    """Verify missing QUAL inputs produce a warning state and CAP-style info output."""
    fwd_path = tmp_path / "S1_27F.fasta"
    rev_path = tmp_path / "S1_1492R.fasta"
    fwd_path.write_text(">f\n" + "A" * 120 + "\n", encoding="utf-8")
    rev_path.write_text(">r\n" + "A" * 120 + "\n", encoding="utf-8")

    contig_path, report = merge_two_reads(
        sample_id="S1",
        fwd_path=fwd_path,
        rev_path=rev_path,
        output_dir=tmp_path,
        min_overlap=100,
        min_identity=0.8,
        min_quality=20.0,
        quality_mode="warning",
    )

    assert contig_path is not None
    assert report.merge_warning == "qualities_absent"
    cap_info = tmp_path / "S1_paired.fasta.cap.info"
    assert cap_info.exists()


def test_build_selected_blast_inputs_legacy_ambiguous_topk_infers_contig_alt(tmp_path: Path) -> None:
    payload = tmp_path / "alt_payload.fasta"
    payload.write_text(">a1\nACGT\n>a2\nACGA\n", encoding="utf-8")

    rows = _build_selected_blast_inputs(
        selected_rows={
            "S1": {
                "assembler_id": "merge_two_reads:biopython",
                "assembler_name": "Merge two reads (Biopython overlap)",
                "status": "ambiguous_topk",
                "payload_fasta": str(payload),
            }
        },
        paired_samples={"S1": {"F": [tmp_path / "f.fasta"], "R": [tmp_path / "r.fasta"]}},
        missing_samples=set(),
        output_fasta=tmp_path / "blast_inputs.fasta",
        output_tsv=tmp_path / "blast_inputs.tsv",
        no_payload_reason="winner_no_payload",
    )

    row = rows["S1"]
    assert row["payload_kind"] == "contig_alt"
    assert row["ambiguity_flag"] == "1"
    assert "hyp1" in row["payload_ids"] and "hyp2" in row["payload_ids"]


def test_build_selected_blast_inputs_keeps_no_payload_reason_when_unclaimed(tmp_path: Path) -> None:
    rows = _build_selected_blast_inputs(
        selected_rows={"S1": {"assembler_id": "merge_two_reads:biopython", "status": "identity_low"}},
        paired_samples={"S1": {"F": [tmp_path / "f.fasta"], "R": [tmp_path / "r.fasta"]}},
        missing_samples=set(),
        output_fasta=tmp_path / "blast_inputs.fasta",
        output_tsv=tmp_path / "blast_inputs.tsv",
        no_payload_reason="winner_no_payload",
    )

    row = rows["S1"]
    assert row["payload_kind"] == "none"
    assert row["reason"] == "winner_no_payload"


def test_build_selected_blast_inputs_respects_singlet_payload_contract(tmp_path: Path) -> None:
    payload = tmp_path / "singlet_payload.fasta"
    payload.write_text(">S1_F\nAAAA\n>S1_R\nTTTT\n", encoding="utf-8")

    rows = _build_selected_blast_inputs(
        selected_rows={
            "S1": {
                "assembler_id": "cap3:relaxed",
                "assembler_name": "CAP3 (Relaxed profile)",
                "status": "singlets_only",
                "payload_fasta": str(payload),
                "payload_kind": "singlet",
                "payload_n": "2",
                "payload_max_len": "4",
                "ambiguity_flag": "0",
                "safety_flag": "none",
                "decision_source": "auto",
                "review_reason": "",
            }
        },
        paired_samples={"S1": {"F": [tmp_path / "f.fasta"], "R": [tmp_path / "r.fasta"]}},
        missing_samples=set(),
        output_fasta=tmp_path / "blast_inputs.fasta",
        output_tsv=tmp_path / "blast_inputs.tsv",
        no_payload_reason="winner_no_payload",
    )

    row = rows["S1"]
    assert row["blast_payload"] == "singlet"
    assert row["payload_kind"] == "singlet"
    assert row["payload_n"] == "2"
    assert "singlet1" in row["payload_ids"] and "singlet2" in row["payload_ids"]


def test_parse_cap3_reports_merge_only_no_missing_info_warning(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Ensure merge-only parsing path avoids CAP3-missing warnings while surfacing merge fields."""
    asm_dir = tmp_path / "asm"
    sample_dir = asm_dir / "sample1"
    sample_dir.mkdir(parents=True)

    (sample_dir / "sample1_paired.merge_report.tsv").write_text(
        "sample_id	orientation	overlap_len	identity	mismatches	contig_len	merge_status	qualities	merge_warning	high_conflict_mismatches\n"
        "sample1	forward	120	0.9900	1	125	merged	absent	qualities_absent	2\n",
        encoding="utf-8",
    )

    with caplog.at_level(logging.WARNING):
        rows = parse_cap3_reports(asm_dir, ["sample1"])

    assert "Missing CAP3 info file" not in caplog.text
    row = rows[0]
    assert row["merge_qualities"] == "absent"
    assert row["merge_warning"] == "qualities_absent"
    assert row["merge_high_conflict_mismatches"] == "2"


def test_merge_two_reads_blocking_rejects_missing_qualities(tmp_path: Path) -> None:
    """Ensure blocking quality mode rejects merges when QUAL data is absent."""
    fwd_path = tmp_path / "S1_27F.fasta"
    rev_path = tmp_path / "S1_1492R.fasta"
    fwd_path.write_text(">f\n" + "A" * 120 + "\n", encoding="utf-8")
    rev_path.write_text(">r\n" + "A" * 120 + "\n", encoding="utf-8")

    contig_path, report = merge_two_reads(
        sample_id="S1",
        fwd_path=fwd_path,
        rev_path=rev_path,
        output_dir=tmp_path,
        min_overlap=100,
        min_identity=0.8,
        min_quality=20.0,
        quality_mode="blocking",
    )

    assert contig_path is None
    assert report.merge_status == "quality_low"
    assert report.qualities == "absent"


def test_merge_two_reads_uses_iupac_for_low_quality_tie(tmp_path: Path) -> None:
    """Verify low-quality mismatch ties emit IUPAC ambiguity codes in merged consensus."""
    fwd_path = tmp_path / "S1_27F.fasta"
    rev_path = tmp_path / "S1_1492R.fasta"
    fwd_seq = "A" * 119 + "C"
    rev_seq = "A" * 119 + "T"
    fwd_path.write_text(f">f\n{fwd_seq}\n", encoding="utf-8")
    rev_path.write_text(f">r\n{rev_seq}\n", encoding="utf-8")

    f_qual = SeqRecord(Seq(fwd_seq), id="f", description="")
    r_qual = SeqRecord(Seq(rev_seq), id="r", description="")
    f_qual.letter_annotations["phred_quality"] = [10] * len(fwd_seq)
    r_qual.letter_annotations["phred_quality"] = [10] * len(rev_seq)
    SeqIO.write([f_qual], Path(f"{fwd_path}.qual"), "qual")
    SeqIO.write([r_qual], Path(f"{rev_path}.qual"), "qual")

    contig_path, report = merge_two_reads(
        sample_id="S1",
        fwd_path=fwd_path,
        rev_path=rev_path,
        output_dir=tmp_path,
        min_overlap=100,
        min_identity=0.8,
        min_quality=5.0,
        quality_mode="warning",
    )

    assert contig_path is not None
    merged = next(SeqIO.parse(contig_path, "fasta"))
    assert str(merged.seq).endswith("Y")
    assert report.merge_status == "merged"


def test_select_best_overlap_uses_quality_to_resolve_top_candidates() -> None:
    """Ensure candidate quality breaks identity ties during overlap selection."""
    candidates = [
        OverlapCandidate("forward", "A" * 120, "A" * 120, 0, 120, 1, 119 / 120, 30.0),
        OverlapCandidate("revcomp", "A" * 120, "A" * 120, 0, 120, 1, 119 / 120, 25.0),
    ]

    chosen = select_best_overlap(
        candidates,
        min_overlap=100,
        min_identity=0.8,
        min_quality=20.0,
        quality_mode="warning",
    )

    assert chosen.status == "ok"
    assert chosen.orientation == "forward"




def test_select_best_overlap_marks_ambiguous_when_quality_tied() -> None:
    """Ensure equal-quality top candidates are marked ambiguous."""
    candidates = [
        OverlapCandidate("forward", "A" * 120, "A" * 120, 0, 120, 1, 119 / 120, 30.0),
        OverlapCandidate("revcomp", "A" * 120, "A" * 120, 0, 120, 1, 119 / 120, 30.0),
    ]

    chosen = select_best_overlap(
        candidates,
        min_overlap=100,
        min_identity=0.8,
        min_quality=20.0,
        quality_mode="warning",
    )

    assert chosen.status == "ambiguous_overlap"


def test_select_best_overlap_marks_ambiguous_without_quality_when_identity_tied() -> None:
    """Ensure identity ties without quality data are marked ambiguous."""
    candidates = [
        OverlapCandidate("forward", "A" * 120, "A" * 120, 0, 120, 1, 119 / 120, None),
        OverlapCandidate("revcomp", "A" * 120, "A" * 120, 0, 120, 1, 119 / 120, None),
    ]

    chosen = select_best_overlap(
        candidates,
        min_overlap=100,
        min_identity=0.8,
        min_quality=20.0,
        quality_mode="warning",
    )

    assert chosen.status == "ambiguous_overlap"

def test_merge_two_reads_routes_high_conflict_to_cap3(tmp_path: Path) -> None:
    """Verify high-confidence conflicts can route merges to CAP3 instead of producing a contig."""
    fwd_path = tmp_path / "S1_27F.fasta"
    rev_path = tmp_path / "S1_1492R.fasta"
    fwd_seq = "A" * 119 + "C"
    rev_seq = "A" * 119 + "T"
    fwd_path.write_text(f">f\n{fwd_seq}\n", encoding="utf-8")
    rev_path.write_text(f">r\n{rev_seq}\n", encoding="utf-8")

    f_qual = SeqRecord(Seq(fwd_seq), id="f", description="")
    r_qual = SeqRecord(Seq(rev_seq), id="r", description="")
    f_qual.letter_annotations["phred_quality"] = [40] * len(fwd_seq)
    r_qual.letter_annotations["phred_quality"] = [40] * len(rev_seq)
    SeqIO.write([f_qual], Path(f"{fwd_path}.qual"), "qual")
    SeqIO.write([r_qual], Path(f"{rev_path}.qual"), "qual")

    contig_path, report = merge_two_reads(
        sample_id="S1",
        fwd_path=fwd_path,
        rev_path=rev_path,
        output_dir=tmp_path,
        min_overlap=100,
        min_identity=0.8,
        min_quality=5.0,
        quality_mode="warning",
        high_conflict_q_threshold=30,
        high_conflict_action="route_cap3",
    )

    assert contig_path is None
    assert report.merge_status == "high_conflict"
    assert report.high_conflict_mismatches == 1


def test_merge_two_reads_flags_high_conflict_on_merge(tmp_path: Path) -> None:
    """Verify high-confidence conflicts can be flagged while still allowing a merge."""
    fwd_path = tmp_path / "S1_27F.fasta"
    rev_path = tmp_path / "S1_1492R.fasta"
    fwd_seq = "A" * 119 + "C"
    rev_seq = "A" * 119 + "T"
    fwd_path.write_text(f">f\n{fwd_seq}\n", encoding="utf-8")
    rev_path.write_text(f">r\n{rev_seq}\n", encoding="utf-8")

    f_qual = SeqRecord(Seq(fwd_seq), id="f", description="")
    r_qual = SeqRecord(Seq(rev_seq), id="r", description="")
    f_qual.letter_annotations["phred_quality"] = [40] * len(fwd_seq)
    r_qual.letter_annotations["phred_quality"] = [40] * len(rev_seq)
    SeqIO.write([f_qual], Path(f"{fwd_path}.qual"), "qual")
    SeqIO.write([r_qual], Path(f"{rev_path}.qual"), "qual")

    contig_path, report = merge_two_reads(
        sample_id="S1",
        fwd_path=fwd_path,
        rev_path=rev_path,
        output_dir=tmp_path,
        min_overlap=100,
        min_identity=0.8,
        min_quality=5.0,
        quality_mode="warning",
        high_conflict_q_threshold=30,
        high_conflict_action="flag",
    )

    assert contig_path is not None
    assert report.merge_status == "merged"
    assert report.high_conflict_mismatches == 1
    assert "high_conflict" in report.merge_warning




def test_merge_two_reads_all_prefers_threshold_passing_engine(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Ensure strategy=all prefers an engine that truly passes thresholds over raw-status `ok` ties."""
    import microseq_tests.assembly.two_read_merge as m

    fwd_path = tmp_path / "S1_27F.fasta"
    rev_path = tmp_path / "S1_1492R.fasta"
    fwd_path.write_text(">f\n" + "A" * 120 + "\n", encoding="utf-8")
    rev_path.write_text(">r\n" + "A" * 120 + "\n", encoding="utf-8")

    monkeypatch.setattr(m, "iter_end_anchored_overlaps", lambda *args, **kwargs: [])
    monkeypatch.setattr(m, "compute_overlap_candidates", lambda *args, **kwargs: [])

    calls = {"n": 0}

    def fake_select(*_args, **_kwargs):
        calls["n"] += 1
        if calls["n"] == 1:
            return OverlapResult("forward", 110, 0.99, 1, None, "A" * 120, "A" * 120, "ok")
        return OverlapResult("forward", 105, 0.95, 2, 30.0, "A" * 120, "A" * 120, "ok")

    monkeypatch.setattr(m, "select_best_overlap", fake_select)

    contig_path, report = m.merge_two_reads(
        sample_id="S1",
        fwd_path=fwd_path,
        rev_path=rev_path,
        output_dir=tmp_path,
        min_overlap=100,
        min_identity=0.8,
        min_quality=20.0,
        quality_mode="blocking",
        overlap_engine_strategy="all",
        overlap_engine_order=["ungapped", "biopython"],
    )

    assert contig_path is not None
    assert report.merge_status == "merged"
    assert report.overlap_engine == "biopython"
def test_merge_two_reads_report_includes_high_conflict_column(tmp_path: Path) -> None:
    """Ensure merge report includes the high_conflict_mismatches output column."""
    fwd_path = tmp_path / "S1_27F.fasta"
    rev_path = tmp_path / "S1_1492R.fasta"
    fwd_path.write_text(">f\n" + "A" * 120 + "\n", encoding="utf-8")
    rev_path.write_text(">r\n" + "A" * 120 + "\n", encoding="utf-8")

    merge_two_reads(
        sample_id="S1",
        fwd_path=fwd_path,
        rev_path=rev_path,
        output_dir=tmp_path,
        min_overlap=100,
        min_identity=0.8,
        min_quality=20.0,
        quality_mode="warning",
    )

    report_path = tmp_path / "S1_paired.merge_report.tsv"
    lines = report_path.read_text(encoding="utf-8").splitlines()
    header_cols = lines[0].split("\t")
    value_cols = lines[1].split("\t")
    assert header_cols[-1] == "high_conflict_mismatches"
    assert len(header_cols) == len(value_cols)
