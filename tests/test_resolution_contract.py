from pathlib import Path

import pytest

pd = pytest.importorskip("pandas")

import microseq_tests.pipeline as mp


def _write_tax(tmp_path: Path, rows: list[dict[str, object]]) -> Path:
    fp = tmp_path / "hits_tax.tsv"
    pd.DataFrame(rows).to_csv(fp, sep="\t", index=False)
    return fp


def _read_review(fp: Path) -> pd.DataFrame:
    return pd.read_csv(fp, sep="\t")


def test_review_queue_partial_hits_and_pair_missing(tmp_path: Path) -> None:
    tax = _write_tax(
        tmp_path,
        [
            {
                "qseqid": "S1|contig|cap3_c1",
                "taxonomy": "k__Bacteria;p__Actinobacteria;g__Cutibacterium;s__avidum",
                "bitscore": 100,
                "pident": 99.5,
                "qcovhsp": 100,
                "evalue": 0.0,
            }
        ],
    )
    blast_rows = {
        "S1": {
            "sample_id": "S1",
            "blast_payload": "contig",
            "structural_hypothesis_n": "3",
            "hypothesis_map": "S1|contig|cap3_c1=hyp1;S1|contig|cap3_c2=hyp2;S1|contig|cap3_c3=hyp3",
            "safety_flag": "none",
            "review_reason": "",
            "trace_status": "PASS",
            "trace_flags": "",
        },
        "S2": {
            "sample_id": "S2",
            "blast_payload": "pair_missing",
            "structural_hypothesis_n": "0",
            "hypothesis_map": "",
            "safety_flag": "none",
            "review_reason": "pair_missing",
            "trace_status": "NA",
            "trace_flags": "",
        },
    }
    out = tmp_path / "review_queue.tsv"
    mp._write_review_queue(tax, out, blast_rows=blast_rows)

    df = _read_review(out).set_index("sample_id")
    assert set(df.index) == {"S1", "S2"}

    s1 = df.loc["S1"]
    assert s1["resolution_state"] == "needs_review"
    assert s1["resolution_reason"] == "partial_hits"
    assert s1["review_action"] == "queue"
    assert s1["missing_hits_n"] == 2

    s2 = df.loc["S2"]
    assert s2["resolution_state"] == "needs_review"
    assert s2["resolution_reason"] == "pair_missing"
    assert s2["review_action"] == "queue"


def test_trace_fail_overrides_taxonomy_agreement(tmp_path: Path) -> None:
    tax = _write_tax(
        tmp_path,
        [
            {
                "qseqid": "S1|contig|cap3_c1",
                "taxonomy": "k__Bacteria;p__Firmicutes;g__Staphylococcus;s__epidermidis",
                "bitscore": 100,
                "pident": 99.9,
                "qcovhsp": 100,
                "evalue": 0.0,
            },
            {
                "qseqid": "S1|contig|cap3_c2",
                "taxonomy": "k__Bacteria;p__Firmicutes;g__Staphylococcus;s__epidermidis",
                "bitscore": 99,
                "pident": 99.8,
                "qcovhsp": 100,
                "evalue": 0.0,
            },
        ],
    )
    blast_rows = {
        "S1": {
            "sample_id": "S1",
            "blast_payload": "contig",
            "structural_hypothesis_n": "2",
            "hypothesis_map": "S1|contig|cap3_c1=hyp1;S1|contig|cap3_c2=hyp2",
            "safety_flag": "trace_fail",
            "review_reason": "trace_fail",
            "trace_status": "FAIL",
            "trace_flags": "LOW_SNR",
        }
    }
    out = tmp_path / "review_queue.tsv"
    mp._write_review_queue(tax, out, blast_rows=blast_rows)

    row = _read_review(out).iloc[0]
    assert row["resolution_state"] == "needs_review"
    assert row["resolution_reason"] == "trace_fail"
    assert row["review_action"] == "queue"


def test_tax_rank_parser_species_missing_genus_available() -> None:
    tax = "k__Bacteria;p__Actinobacteria;g__Cutibacterium"
    assert mp._extract_tax_label(tax, rank="species") == ""
    assert mp._extract_tax_label(tax, rank="genus") == "cutibacterium"


def test_resolved_hypothesis_deterministic_under_tie() -> None:
    per_hyp = pd.DataFrame(
        [
            {
                "qseqid": "S1|contig|cap3_c2",
                "taxonomy": "k__Bacteria;p__Firmicutes;g__Staphylococcus;s__epidermidis",
                "bitscore": 100,
                "pident": 99.9,
                "qcovhsp": 100,
                "evalue": 0.0,
            },
            {
                "qseqid": "S1|contig|cap3_c1",
                "taxonomy": "k__Bacteria;p__Firmicutes;g__Staphylococcus;s__epidermidis",
                "bitscore": 100,
                "pident": 99.9,
                "qcovhsp": 100,
                "evalue": 0.0,
            },
        ]
    )
    blast_meta = {
        "blast_payload": "contig",
        "structural_hypothesis_n": "2",
        "hypothesis_map": "S1|contig|cap3_c1=hyp1;S1|contig|cap3_c2=hyp2",
        "safety_flag": "none",
        "review_reason": "",
        "trace_status": "PASS",
        "trace_flags": "",
    }
    resolved = mp._resolve_sample_resolution(
        "S1", per_hyp, label_col="taxonomy", blast_meta=blast_meta, tax_rank="species"
    )
    assert resolved["resolved_hypothesis"] == "hyp1"


def test_taxonomy_missing_vs_no_hits_distinction(tmp_path: Path) -> None:
    blast_rows = {
        "S1": {
            "sample_id": "S1",
            "blast_payload": "contig",
            "structural_hypothesis_n": "1",
            "hypothesis_map": "S1|contig|cap3_c1=struct1",
            "safety_flag": "none",
            "review_reason": "",
            "trace_status": "PASS",
            "trace_flags": "",
        }
    }

    bad_tax = tmp_path / "bad_tax.tsv"
    bad_tax.write_text("not_qseqid\tfoo\n", encoding="utf-8")
    bad_resolved = mp._annotate_resolution_from_tax(bad_tax, blast_rows)
    assert bad_resolved["S1"]["resolution_reason"] == "taxonomy_missing"

    empty_tax = tmp_path / "empty_tax.tsv"
    pd.DataFrame(columns=["qseqid", "taxonomy", "bitscore", "pident", "qcovhsp", "evalue"]).to_csv(
        empty_tax, sep="\t", index=False
    )
    empty_resolved = mp._annotate_resolution_from_tax(empty_tax, blast_rows)
    assert empty_resolved["S1"]["resolution_reason"] == "no_hits"
