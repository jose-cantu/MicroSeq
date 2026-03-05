from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("pandas")
pytest.importorskip("Bio")

from microseq_tests.blast.run_blast import run_blast


def test_run_blast_empty_query_raises_with_selected_mode_hint(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    query = tmp_path / "blast_inputs.fasta"
    monkeypatch.setattr(
        "microseq_tests.blast.run_blast.load_config",
        lambda: {"databases": {"nt": {"blastdb": "dummy"}}},
    )
    query.write_text("", encoding="utf-8")
    (tmp_path / "blast_inputs.tsv").write_text(
        "sample_id\tblast_payload\tpayload_ids\treason\tselected_assembler_id\tselected_assembler_name\tselected_status\n"
        "S1\tno_payload\t\tselected_backend_no_payload\tmerge_two_reads:biopython\tMerge two reads (Biopython overlap)\tambiguous_overlap\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError) as excinfo:
        run_blast(query, "nt", tmp_path / "hits.tsv")

    message = str(excinfo.value)
    assert "nothing to BLAST" in message
    assert "selected_backend_no_payload=1" in message
    assert "ambiguous_overlap=1" in message

