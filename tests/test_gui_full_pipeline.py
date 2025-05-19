import pytest
from pathlib import Path

pytest.importorskip("PySide6")

from PySide6.QtCore import Qt
from microseq_tests.gui.main_window import MainWindow
import microseq_tests.pipeline as pl


def test_full_button_progress(qtbot, monkeypatch, tmp_path):
    # stub heavy pipeline stages to avoid running external tools
    for fn in ("run_trim", "run_fastq_to_fasta", "run_blast_stage", "run_add_tax", "run_postblast"):
        monkeypatch.setattr(pl, fn, lambda *a, **k: None)

    win = MainWindow()
    qtbot.addWidget(win)
    # set a dummy input file
    dummy = tmp_path / "dummy.fastq"
    dummy.write_text("@r1\nACGT\n+\n!!!!\n")
    win._infile = dummy

    qtbot.mouseClick(win.full_btn, Qt.LeftButton)
    qtbot.waitSignal(win._worker.finished, timeout=2000)

    assert win.progress.value() == 100
