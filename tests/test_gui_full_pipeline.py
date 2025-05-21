import os 
import pytest
from pathlib import Path

pytest.importorskip("PySide6")
pytest.importorskip("pytestqt")

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMessageBox 
from microseq_tests.gui.main_window import MainWindow
import microseq_tests.pipeline as pl


def test_full_button_progress(qtbot, monkeypatch, tmp_path):
    # stub heavy pipeline stages to avoid running external tools also modal dialogs to not block headless runs 
    monkeypatch.setattr(QMessageBox, "information", lambda *a, **k: QMessageBox.Ok)
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **k: QMessageBox.Ok)
    monkeypatch.setattr(QMessageBox, "critical", lambda *a, **k: QMessageBox.Ok) 

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

    # progress bar is reset to 0 in _done(); wait for that 
    qtbot.waitUntil(lambda: win.progress.value() == 0, timeout=1000)
    assert win.progress.value() == 0          # ← expect 0, not 100 
