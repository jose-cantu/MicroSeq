from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("pandas")
pytest.importorskip("Bio")

from microseq_tests import pipeline


def test_paired_error_includes_suggestions(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Paired pipeline errors should surface the suggestion text."""

    infile = tmp_path / "input"
    infile.mkdir()

    fastq_dir = tmp_path / "out" / "passed_qc_fastq"
    fastq_dir.mkdir(parents=True)

    fastq_data = "@r1\nA\n+\n!\n"
    (fastq_dir / "S1_27F.fastq").write_text(fastq_data, encoding="utf-8")
    (fastq_dir / "S1_1492R.fastq").write_text(fastq_data, encoding="utf-8")

    def fake_run_trim(*_args, **_kwargs):
        return 0

    def fake_run_fastq_to_fasta(*_args, **_kwargs):
        return None

    def fake_fastq_to_paired_fastas(_fastq_dir: Path, fasta_dir: Path):
        fasta_dir.mkdir(parents=True, exist_ok=True)
        fwd = fasta_dir / "S1_27F.fasta"
        rev = fasta_dir / "S1_1492R.fasta"
        fwd.write_text(">fwd\nA\n", encoding="utf-8")
        rev.write_text(">rev\nA\n", encoding="utf-8")
        return [fwd, rev]

    monkeypatch.setattr(pipeline, "run_trim", fake_run_trim)
    monkeypatch.setattr(pipeline, "run_fastq_to_fasta", fake_run_fastq_to_fasta)
    monkeypatch.setattr(pipeline, "_fastq_to_paired_fastas", fake_fastq_to_paired_fastas)
    monkeypatch.setattr(pipeline, "assemble_pairs", lambda *_, **__: [])

    with pytest.raises(ValueError) as excinfo:
        pipeline.run_full_pipeline(
            infile,
            "nt",
            tmp_path / "out",
            mode="paired",
        )

    message = str(excinfo.value)
    assert "detected 1 forward and 1 reverse reads" in message
    assert "Example flag" in message


def test_selected_mode_preserves_pair_missing_when_compare_pairing_report_lists_all(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    """Compare-driven selected mode must keep pair_missing samples as pair_missing."""

    infile = tmp_path / "input"
    infile.mkdir()
    out_dir = tmp_path / "out"

    # stage minimal paired fasta inputs
    fastq_dir = out_dir / "passed_qc_fastq"
    fastq_dir.mkdir(parents=True)
    (fastq_dir / "S1_27F.fastq").write_text("@r1\nA\n+\n!\n", encoding="utf-8")
    (fastq_dir / "S1_1492R.fastq").write_text("@r2\nA\n+\n!\n", encoding="utf-8")
    (fastq_dir / "S2_27F.fastq").write_text("@r3\nA\n+\n!\n", encoding="utf-8")  # missing R

    def fake_run_trim(*_args, **_kwargs):
        return 0

    def fake_run_fastq_to_fasta(*_args, **_kwargs):
        return None

    def fake_fastq_to_paired_fastas(_fastq_dir: Path, fasta_dir: Path, *, use_qual: bool = True):
        fasta_dir.mkdir(parents=True, exist_ok=True)
        f1 = fasta_dir / "S1_27F.fasta"
        r1 = fasta_dir / "S1_1492R.fasta"
        f2 = fasta_dir / "S2_27F.fasta"
        f1.write_text(">S1_27F\nAC\n", encoding="utf-8")
        r1.write_text(">S1_1492R\nAC\n", encoding="utf-8")
        f2.write_text(">S2_27F\nAC\n", encoding="utf-8")
        return [f1, r1, f2]

    def fake_compare(
        _input_dir,
        output_dir,
        *,
        pairing_report=None,
        **_kwargs,
    ):
        asm = Path(output_dir) / "asm"
        asm.mkdir(parents=True, exist_ok=True)
        compare_tsv = asm / "compare_assemblers.tsv"
        payload = asm / "compare" / "cap3_strict" / "S1" / "payload.fasta"
        payload.parent.mkdir(parents=True, exist_ok=True)
        payload.write_text(">c1\nAC\n", encoding="utf-8")
        compare_tsv.write_text(
            "\t".join([
                "sample_id", "assembler_id", "assembler_name", "dup_policy", "status", "selected_engine", "contig_len", "warnings", "payload_fasta"
            ])
            + "\n"
            + "\t".join(["S1", "cap3:strict", "CAP3 (Strict profile)", "error", "assembled", "cap3", "2", "", str(payload)])
            + "\n",
            encoding="utf-8",
        )
        if pairing_report:
            Path(pairing_report).parent.mkdir(parents=True, exist_ok=True)
            # includes S2 as pair-missing line; selected-mode logic must not erase missing_samples from this
            Path(pairing_report).write_text(
                "sid\tF_path\tR_path\tdetector\n"
                "S1\t/path/f\t/path/r\t\n"
                "S2\t/path/f_only\t\t\n",
                encoding="utf-8",
            )
        return compare_tsv

    def fake_blast(*_args, **kwargs):
        Path(kwargs.get("out_tsv") if "out_tsv" in kwargs else _args[2]).write_text(
            "qseqid\tsseqid\tpident\tqcovhsp\n", encoding="utf-8"
        )
        return 0

    def fake_add_tax(_hits, _tax, out_tsv):
        Path(out_tsv).write_text("qseqid\ttaxonomy\n", encoding="utf-8")
        return 0

    monkeypatch.setattr(pipeline, "run_trim", fake_run_trim)
    monkeypatch.setattr(pipeline, "run_fastq_to_fasta", fake_run_fastq_to_fasta)
    monkeypatch.setattr(pipeline, "_fastq_to_paired_fastas", fake_fastq_to_paired_fastas)
    monkeypatch.setattr(pipeline, "run_compare_assemblers", fake_compare)
    monkeypatch.setattr(pipeline, "run_blast_stage", fake_blast)
    monkeypatch.setattr(pipeline, "run_add_tax", fake_add_tax)

    pipeline.run_full_pipeline(
        infile,
        "nt",
        out_dir,
        mode="paired",
        assembler_id="cap3:strict",
        write_blast_inputs=True,
        use_blast_inputs=True,
    )

    blast_inputs = (out_dir / "asm" / "blast_inputs.tsv").read_text(encoding="utf-8")
    assert "S2\tpair_missing" in blast_inputs


def test_selected_summary_includes_cap3_compatible_telemetry_columns(tmp_path: Path):
    summary = tmp_path / "assembly_summary.tsv"
    pipeline._write_selected_assembly_summary(
        summary,
        paired_samples={"S1": {"F": [tmp_path / "f.fasta"], "R": [tmp_path / "r.fasta"]}},
        missing_samples=set(),
        selected_rows={"S1": {"assembler_id": "cap3:strict", "assembler_name": "CAP3", "status": "assembled", "contig_len": "123", "selected_engine": "cap3"}},
        blast_rows={"S1": {"blast_payload": "contig", "reason": "selected_payload"}},
        assembler_mode_label="selected",
        primer_mode="clip",
        primer_stage="post_quality",
        primer_preset="16S",
        primer_source="preset",
        overlap_engine_strategy="cascade",
        overlap_engine_order="biopython,ungapped,edlib",
        overlap_quality_mode="warning",
        overlap_rows={"S1": {"status": "ok", "overlap_identity": "0.99", "overlap_quality": "30", "orientation": "forward"}},
    )

    header = summary.read_text(encoding="utf-8").splitlines()[0].split("\t")
    for col in [
        "primer_mode",
        "primer_stage",
        "primer_preset",
        "primer_source",
        "overlap_engine_strategy",
        "overlap_engine_order",
        "overlap_quality_mode",
    ]:
        assert col in header


def test_build_selected_blast_inputs_emits_payload_contract_and_hypothesis_ids(tmp_path: Path):
    payload = tmp_path / "payload.fasta"
    payload.write_text(">alt1\nACGT\n>alt2\nACGA\n", encoding="utf-8")

    rows = pipeline._build_selected_blast_inputs(
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
    assert row["payload_n"] == "2"
    assert row["ambiguity_flag"] == "1"
    assert "hyp1" in row["payload_ids"] and "hyp2" in row["payload_ids"]


def test_select_best_compare_rows_handles_missing_payload_contract_columns(tmp_path: Path):
    compare_tsv = tmp_path / "compare.tsv"
    compare_tsv.write_text(
        "sample_id\tassembler_id\tassembler_name\tstatus\tcontig_len\n"
        "S1\tcap3:relaxed\tCAP3\tassembled\t100\n",
        encoding="utf-8",
    )

    winners = pipeline._select_best_compare_rows(compare_tsv, assembler_mode="all")
    assert winners["S1"]["assembler_id"] == "cap3:relaxed"


def test_build_selected_blast_inputs_reconciles_missing_payload_file_to_none(tmp_path: Path):
    rows = pipeline._build_selected_blast_inputs(
        selected_rows={
            "S1": {
                "assembler_id": "merge_two_reads:biopython",
                "assembler_name": "Merge two reads (Biopython overlap)",
                "status": "merged",
                "payload_fasta": str(tmp_path / "missing_payload.fasta"),
                "payload_kind": "contig",
                "payload_n": "1",
                "payload_max_len": "1200",
            }
        },
        paired_samples={"S1": {"F": [tmp_path / "f.fasta"], "R": [tmp_path / "r.fasta"]}},
        missing_samples=set(),
        output_fasta=tmp_path / "blast_inputs.fasta",
        output_tsv=tmp_path / "blast_inputs.tsv",
        no_payload_reason="winner_no_payload",
    )

    row = rows["S1"]
    assert row["blast_payload"] == "no_payload"
    assert row["payload_kind"] == "none"
    assert row["reason"] == "payload_missing_or_empty"
