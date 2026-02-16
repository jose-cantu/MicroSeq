from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("pandas")
pytest.importorskip("Bio")

from microseq_tests.pipeline import _normalize_primer_trim_cfg, run_trim, run_full_pipeline, run_paired_assembly
from microseq_tests.primer_catalog import parse_primer_sequences, build_primer_cfg_override
from microseq_tests.assembly.cap3_report import parse_cap3_reports
from microseq_tests.trimming.biopy_trim import TRIM_SUMMARY_COLUMNS
import microseq_tests.pipeline as pipeline
import microseq_tests.primer_catalog as primer_catalog


def test_off_mode_ignores_invalid_override_sequences(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Off mode should not fail if override primer text is junk."""
    in_dir = tmp_path / "in"
    in_dir.mkdir()
    (in_dir / "a.fastq").write_text("@r\nACGT\n+\n!!!!\n", encoding="utf-8")

    def fake_trim_fastq_inputs(_input: Path, trim_dir: Path, summary_tsv=None):
        trim_dir.mkdir(parents=True, exist_ok=True)
        (trim_dir / "trimmed.fastq").write_text("@r\nACGT\n+\n!!!!\n", encoding="utf-8")
        return trim_dir / "trimmed.fastq"

    monkeypatch.setattr(pipeline, "trim_fastq_inputs", fake_trim_fastq_inputs)

    out_dir = tmp_path / "out"
    rc = run_trim(
        in_dir,
        out_dir,
        sanger=False,
        primer_cfg_override={
            "mode": "off",
            "forward_primers": ["NOT A DNA SEQ"],
            "reverse_primers": ["%%%%"],
        },
    )

    assert rc == 0
    assert (out_dir / "qc" / "trimmed.fasta").exists()


def test_blank_gui_sequences_keep_config_custom_primers() -> None:
    base_cfg = {
        "primer_trim": {
            "mode": "detect",
            "stage": "post_quality",
            "preset": "",
            "forward_primers": ["AGAGTTTGATCMTGGCTCAG"],
            "reverse_primers": ["TACGGYTACCTTGTTACGACTT"],
        }
    }
    # Mimic GUI override payload when custom text boxes are blank: no primer list keys.
    primer_override = {"mode": "detect", "stage": "post_quality"}
    merged = dict(base_cfg)
    merged_primer = dict(base_cfg["primer_trim"])
    merged_primer.update(primer_override)
    merged["primer_trim"] = merged_primer

    normalized = _normalize_primer_trim_cfg(merged)
    assert normalized["forward_primers"] == ["AGAGTTTGATCMTGGCTCAG"]
    assert normalized["reverse_primers"] == ["TACGGYTACCTTGTTACGACTT"]




def test_comment_header_only_custom_text_does_not_override_config_lists() -> None:
    override = build_primer_cfg_override(
        mode="detect",
        stage="post_quality",
        preset="",
        forward_raw=">27F\n#comment",
        reverse_raw="#just-note",
        for_preview=False,
    )
    assert "forward_primers" not in override
    assert "reverse_primers" not in override

    cfg = {
        "primer_trim": {
            "mode": "detect",
            "stage": "post_quality",
            "forward_primers": ["AGAGTTTGATCMTGGCTCAG"],
            "reverse_primers": ["TACGGYTACCTTGTTACGACTT"],
        }
    }
    merged = dict(cfg)
    primer = dict(cfg["primer_trim"])
    primer.update(override)
    merged["primer_trim"] = primer
    normalized = _normalize_primer_trim_cfg(merged)
    assert normalized["forward_primers"] == ["AGAGTTTGATCMTGGCTCAG"]
    assert normalized["reverse_primers"] == ["TACGGYTACCTTGTTACGACTT"]


def test_blank_gui_preset_does_not_override_config_preset() -> None:
    override = build_primer_cfg_override(
        mode="detect",
        stage="post_quality",
        preset="",
        forward_raw="",
        reverse_raw="",
        for_preview=False,
    )
    assert "preset" not in override

    cfg = {
        "primer_trim": {
            "mode": "detect",
            "stage": "post_quality",
            "preset": "16S_27F_1492R",
            "forward_primers": [],
            "reverse_primers": [],
        }
    }
    merged = dict(cfg)
    primer = dict(cfg["primer_trim"])
    primer.update(override)
    merged["primer_trim"] = primer
    normalized = _normalize_primer_trim_cfg(merged)
    assert normalized["preset"] == "16S_27F_1492R"
    assert normalized["forward_primers"] == ["AGAGTTTGATCMTGGCTCAG"]
    assert normalized["reverse_primers"] == ["TACGGYTACCTTGTTACGACTT"]


def test_normalize_primer_cfg_custom_both_sides_dominates_preset() -> None:
    cfg = {
        "primer_trim": {
            "mode": "detect",
            "stage": "post_quality",
            "preset": "16S_27F_1492R",
            "forward_primers": ["AACCGGTTAACC"],
            "reverse_primers": ["TTGGCCAATTGG"],
        }
    }
    normalized = _normalize_primer_trim_cfg(cfg)
    assert normalized["primer_source"] == "custom"
    assert normalized["preset"] == ""
    assert normalized["preset_configured"] == "16S_27F_1492R"
    assert normalized["forward_primers"] == ["AACCGGTTAACC"]
    assert normalized["reverse_primers"] == ["TTGGCCAATTGG"]


def test_normalize_primer_cfg_custom_one_side_is_mixed() -> None:
    cfg = {
        "primer_trim": {
            "mode": "detect",
            "stage": "post_quality",
            "preset": "16S_27F_1492R",
            "forward_primers": ["AACCGGTTAACC"],
            "reverse_primers": [],
        }
    }
    normalized = _normalize_primer_trim_cfg(cfg)
    assert normalized["primer_source"] == "mixed"
    assert normalized["preset"] == ""
    assert normalized["forward_primers"] == ["AACCGGTTAACC"]
    # reverse falls back to preset for mixed source
    assert normalized["reverse_primers"] == ["TACGGYTACCTTGTTACGACTT"]


def test_parse_cap3_reports_includes_primer_source_field() -> None:
    rows = parse_cap3_reports(
        Path("/tmp/nowhere"),
        ["missing1"],
        missing_samples={"missing1"},
        primer_mode="detect",
        primer_stage="post_quality",
        primer_preset="",
        primer_source="mixed",
    )
    assert rows[0]["primer_source"] == "mixed"


def test_load_primer_catalog_uses_packaged_resource_when_repo_config_missing(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    missing = tmp_path / "does-not-exist.yaml"
    packaged_root = tmp_path / "pkg"
    packaged_root.mkdir()
    packaged_file = packaged_root / "primer_catalog.yaml"
    packaged_file.write_text(
        """trim_presets:
  PACKAGED_ONLY:
    forward_primers: ['AAAAAAAACC']
    reverse_primers: ['TTTTTTTTGG']
""",
        encoding="utf-8",
    )

    monkeypatch.setattr(primer_catalog, "CATALOG_PATH", missing)
    monkeypatch.setattr(primer_catalog, "_DEFAULT_CATALOG", {"pairing_sets": {}, "trim_presets": {}})
    monkeypatch.setattr(primer_catalog.resources, "files", lambda _pkg: packaged_root)

    loaded = primer_catalog.load_primer_catalog(path=missing)
    assert "trim_presets" in loaded
    assert "PACKAGED_ONLY" in loaded["trim_presets"]
def test_primer_sequence_parser_accepts_label_and_comment_lines() -> None:
    text = """
    27F AGAGTTTGATCMTGGCTCAG # canonical forward
    reverse: TACGGYTACCTTGTTACGACTT
    """
    seqs = parse_primer_sequences(text)
    assert seqs == ["AGAGTTTGATCMTGGCTCAG", "TACGGYTACCTTGTTACGACTT"]




def test_primer_sequence_parser_accepts_fasta_headers_and_normalizes_u_to_t() -> None:
    text = """>27F
AGAGTTTGATCMTGGCTCAG
>1492R
UACGGYTACCTTGTTACGACTT
"""
    seqs = parse_primer_sequences(text)
    assert seqs == ["AGAGTTTGATCMTGGCTCAG", "TACGGYTACCTTGTTACGACTT"]

def test_pre_quality_clip_emits_report_and_summary_primer_columns(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    ab1_dir = tmp_path / "ab1"
    ab1_dir.mkdir()
    (ab1_dir / "s1.ab1").write_text("stub", encoding="utf-8")

    def fake_ab1_to_fastq(_src: Path, dst: Path) -> None:
        dst.mkdir(parents=True, exist_ok=True)
        (dst / "S1_27F.fastq").write_text(
            "@r1\nAGAGTTTGATCMTGGCTCAGAAAA\n+\n" + "I" * 24 + "\n",
            encoding="utf-8",
        )

    def fake_biopy_trim(_in: Path, out_dir: Path, combined_tsv: Path | None = None, **_kwargs) -> None:
        passed = out_dir.parent / "passed_qc_fastq"
        passed.mkdir(parents=True, exist_ok=True)
        (passed / "S1_27F.fastq").write_text("@r1\nAAAA\n+\nIIII\n", encoding="utf-8")
        if combined_tsv:
            combined_tsv.parent.mkdir(parents=True, exist_ok=True)
            header = "\t".join(TRIM_SUMMARY_COLUMNS)
            row = "S1_27F.fastq\t1\t4.0\t40.00\tpass"
            combined_tsv.write_text(f"{header}\n{row}\n", encoding="utf-8")

    monkeypatch.setattr(pipeline, "ab1_to_fastq", fake_ab1_to_fastq)
    monkeypatch.setattr(pipeline, "biopy_trim", fake_biopy_trim)

    summary = tmp_path / "out" / "qc" / "trim_summary.tsv"
    run_trim(
        ab1_dir,
        tmp_path / "out",
        sanger=True,
        summary_tsv=summary,
        primer_cfg_override={
            "mode": "clip",
            "stage": "pre_quality",
            "forward_primers": ["AGAGTTTGATCMTGGCTCAG"],
            "reverse_primers": [],
        },
    )

    assert (tmp_path / "out" / "qc" / "primer_trim_report.tsv").exists()
    text = summary.read_text(encoding="utf-8")
    assert "primer_trimmed" in text.splitlines()[0]
    assert "S1_27F.fastq" in text


def test_missing_pair_row_keeps_configured_invariants() -> None:
    rows = parse_cap3_reports(
        Path("/tmp/nowhere"),
        ["missing1"],
        missing_samples={"missing1"},
        primer_mode="detect",
        primer_stage="pre_quality",
        primer_preset="16S_27F_1492R",
        overlap_engine_strategy="cascade",
        overlap_engine_order="biopython,ungapped",
        overlap_quality_mode="warning",
    )
    row = rows[0]
    assert row["status"] == "pair_missing"
    assert row["primer_mode"] == "detect"
    assert row["primer_stage"] == "pre_quality"
    assert row["primer_preset"] == "16S_27F_1492R"
    assert row["overlap_engine_strategy"] == "cascade"
    assert row["configured_engine"] == "auto"
    assert row["fallback_used"] == "n/a"


def test_primer_override_kwargs_off_mode_preview_keeps_custom_sequences() -> None:
    normal = build_primer_cfg_override(
        mode="off",
        stage="post_quality",
        preset="",
        forward_raw="27F AGAGTTTGATCMTGGCTCAG",
        reverse_raw="1492R TACGGYTACCTTGTTACGACTT",
        for_preview=False,
    )
    preview = build_primer_cfg_override(
        mode="off",
        stage="post_quality",
        preset="",
        forward_raw="27F AGAGTTTGATCMTGGCTCAG",
        reverse_raw="1492R TACGGYTACCTTGTTACGACTT",
        for_preview=True,
    )

    assert normal["mode"] == "off"
    assert "forward_primers" not in normal
    assert "reverse_primers" not in normal
    assert preview["mode"] == "off"
    assert preview["forward_primers"] == ["AGAGTTTGATCMTGGCTCAG"]
    assert preview["reverse_primers"] == ["TACGGYTACCTTGTTACGACTT"]





def test_run_full_pipeline_passes_raw_primer_override_to_run_trim(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    infile = tmp_path / "reads.fastq"
    infile.write_text("@r\nACGT\n+\n!!!!\n", encoding="utf-8")

    seen: dict[str, object] = {}

    def fake_run_trim(_infile, _out_dir, *, sanger=False, summary_tsv=None, primer_cfg_override=None):
        seen["primer_cfg_override"] = primer_cfg_override
        raise RuntimeError("stop-after-trim")

    monkeypatch.setattr(pipeline, "run_trim", fake_run_trim)

    raw_override = {
        "mode": "detect",
        "stage": "post_quality",
        "preset": "16S_27F_1492R",
        "forward_primers": ["AACCGGTTAACC"],
    }

    with pytest.raises(RuntimeError, match="stop-after-trim"):
        run_full_pipeline(
            infile,
            "nt",
            out_dir=tmp_path / "out",
            mode="single",
            primer_cfg_override=raw_override,
        )

    assert seen["primer_cfg_override"] == raw_override


def test_run_paired_assembly_accepts_use_qual_and_forwards(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    captured: dict[str, object] = {}

    def fake_assemble_pairs(*args, **kwargs):
        captured.update(kwargs)
        return []

    monkeypatch.setattr(pipeline, "assemble_pairs", fake_assemble_pairs)

    run_paired_assembly(
        tmp_path,
        tmp_path / "asm",
        use_qual=False,
    )

    assert captured["use_qual"] is False
