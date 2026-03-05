from pathlib import Path

import pytest

pytest.importorskip("pandas")
pytest.importorskip("Bio")

import microseq_tests.pipeline as mp


def test_default_output_dir_for_nested_single_sample_dir(tmp_path: Path) -> None:
    selected = tmp_path / "container"
    deep = selected / "10292025_1080497"
    deep.mkdir(parents=True)
    (deep / "A01.ab1").write_bytes(b"x")

    out = mp._default_run_output_dir(selected)

    assert out == selected / "10292025_1080497_microseq"


def test_default_output_dir_for_input_root_with_files(tmp_path: Path) -> None:
    selected = tmp_path / "10292025_1080497"
    selected.mkdir(parents=True)
    (selected / "A01.ab1").write_bytes(b"x")
    (selected / "A02.fastq").write_text("@r\nACGT\n+\n!!!!\n", encoding="utf-8")

    out = mp._default_run_output_dir(selected)

    assert out == tmp_path / "10292025_1080497_microseq"
