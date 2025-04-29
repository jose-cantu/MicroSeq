# src/microseq_tests/gui

from __future__ import annotations 
import logging
L = logging.getLogger(__name__)

import sys, logging, traceback, subprocess, shlex  
from pathlib import Path 
from typing import Optional 

from PySide6.QtCore import Qt, QObject, QThread, Signal, Slot, QMetaObject, Q_ARG
from PySide6.QtWidgets import (
        QApplication, QMainWindow, QWidget, QFileDialog, QVBoxLayout,
        QHBoxLayout, QPushButton, QLabel, QTextEdit, QSpinBox, QMessageBox, QComboBox, QProgressBar   
        )

# ==== MicroSeq wrappers ---------- 
from microseq_tests.pipeline import (
        run_blast_stage, 
        run_trim,
        run_assembly,
        run_add_tax,
        run_postblast,
        )
from microseq_tests.utility.utils import setup_logging, load_config 

# Worker class ---------------- 
class Worker(QObject):
    """Background runner living in its own QThread.""" 
    finished = Signal(int) # exit-code 0 = success 
    log = Signal(str) # text lines for log 
    progress = Signal(int) 

    # Which stage here to run and its kwargs are injected at construction 
    def __init__(self, fn, *args, **kwargs):
        super().__init__()
        self._fn = fn 
        self._args = args 
        self._kwargs = kwargs 


    @Slot() # design to warn if any errors occur 
    def run(self):
        try:
            # --------- live banner so GUI shows immediate activity -- 
            L.info("▶ Starting BLAST…")
            self.progress.emit(0)
            rc = self._fn(*self._args, **self._kwargs)
            self.progress.emit(100) 
            L.info("✔ BLAST finished with rc=%s", rc)

        except Exception:
            self.log.emit(traceback.format_exc())
            rc = 1 
        self.finished.emit(rc) 


# ---- Main Window Constructor 
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("MicroSeq GUI α") # title for it 
        self.resize(800, 520) # size of app considering also log space here 

        # widgets --------------------------------------------------------
        self.fasta_lbl = QLabel("FASTA / AB1: —")
        browse_btn = QPushButton("Browse..")
        browse_btn.clicked.connect(self._choose_infile) 

        # label place holder and browse button wired to file picker 
        self.id_spin = QSpinBox()
        self.id_spin.setRange(50, 100) # for simplicity I set spinbox ID % 50-100 threshold 
        self.id_spin.setValue(97) # identity 97% default 
        self.id_spin.setSuffix(" % ID")

        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, 32)
        self.threads_spin.setValue(1)
        self.threads_spin.setSuffix(" CPU")

        # ----- alignment coverage (% of query aligned) 
        self.qcov_spin = QSpinBox()
        self.qcov_spin.setRange(10, 100)
        self.qcov_spin.setValue(80)
        self.qcov_spin.setSuffix(" % Q-cov")

        # -------- max target hits 
        self.hits_spin = QSpinBox()
        self.hits_spin.setRange(1, 500)
        self.hits_spin.setValue(5)
        self.hits_spin.setSuffix(" hits") 

        # ------- progress bar ----------
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)

        # here you just run blast ....
        self.run_btn = QPushButton("Run Blast")
        self.run_btn.clicked.connect(self._launch_blast)

        # lets you scroll log output and nowarp keeps long commadn lines intact 
        self.log_box = QTextEdit(readOnly=True, lineWrapMode=QTextEdit.NoWrap)

        # ---- Layout here will update UX -------------------------
        top = QHBoxLayout()
        top.addWidget(self.fasta_lbl)
        top.addWidget(browse_btn) 

        # horizontal row for file picker 
        mid = QHBoxLayout()
        
        # --- allows for database selection here in UX --------------------
        self.db_box = QComboBox()   
        self.db_box.addItems(load_config()["databases"].keys())


        mid.addWidget(QLabel("DB"))
        mid.addWidget(self.db_box)


        # --- concerning identity ---- just look at Qlabel it will tell you José 
        mid.addWidget(QLabel("Identity")) 
        mid.addWidget(self.id_spin)
        
        mid.addWidget(QLabel("Q-cov"))
        mid.addWidget(self.qcov_spin)

        mid.addWidget(QLabel("Max hits"))
        mid.addWidget(self.hits_spin) 

        mid.addWidget(QLabel("Threads"))
        mid.addWidget(self.threads_spin) 

        mid.addWidget(self.progress) 
        
        mid.addStretch()  # pushes Run button to the far right here  
        mid.addWidget(self.run_btn) 

        # vertical stack picker row, settings row, then logpane expands 
        outer = QVBoxLayout()
        outer.addLayout(top)
        outer.addLayout(mid)
        outer.addWidget(self.log_box) 

        # Embed composite layout as the window's central widget 
        root = QWidget(); root.setLayout(outer)
        self.setCentralWidget(root) 

        # state and logging connection 
        self._infile: Optional[Path] = None 
        setup_logging(level=logging.INFO) # file + console stderr 
        logging.getLogger().handlers.append(_QtHandler(self.log_box)) # appended to log_box 
    
    # ---- file picker --------------------------
    def _choose_infile(self):
# returns path string and intial dir = home filters removes other file types
        path, _ = QFileDialog.getOpenFileName(
                self, "Select FASTA/FASTQ/AB1", 
                str(Path.home()), "Seq files (*.fasta *.fa *.ab1)"
                )
        # stores path; updates label for user feedback 
        if path:
            self._infile = Path(path)
            self.fasta_lbl.setText(f"FASTA / AB1: {self._infile.name}") 
    
    # ---- Blast Stage Demo ------------------------------
    # guarding against accidental click 
    def _launch_blast(self):
        if not self._infile:
            QMessageBox.warning(self, "No input", "Choose a file first.")
            return

        self.progress.setValue(0)  # resets progress bar during each run 

        # derive output file beside input; disables button; logs starts 
        hits_path = self._infile.with_suffix(".hits.tsv")

        self.run_btn.setEnabled(False)
        self.log_box.append(f"\n▶ BLAST {self._infile.name} -> {hits_path.name}")

        # worker and thread wiring -------------------
        worker = Worker(
                run_blast_stage,
                self._infile,
                self.db_box.currentText(),  # default DB key selection 
                hits_path,
                identity=self.id_spin.value(),
                qcov=self.qcov_spin.value(),
                max_target_seqs=self.hits_spin.value(),
                threads=self.threads_spin.value(),
                )
        # injecting wrapper + args into Worker; moves into new thread 
        thread = QThread(self) # autodeleted with window 
        worker.moveToThread(thread) 

        # wire signal between log streaming, completion and thread shutdown

        worker.log.connect(self.log_box.append)
        worker.progress.connect(self.progress.setValue) 
        thread.started.connect(worker.run)
        worker.finished.connect(lambda rc: self._done(rc, hits_path))
        worker.finished.connect(thread.quit)

        # keep references so closeEvent can see them 
        self._worker = worker 
        self._thread = thread 

        print("DEBUG – about to start Worker with", self.db_box.currentText(),
        self.id_spin.value(), self.qcov_spin.value(),
        self.hits_spin.value(), self.threads_spin.value(), flush=True)
        thread.start() 

    # Callback -------------------------
    def _done(self, rc: int, out: Path):
        # appends outcome indicator to log 
        msg = "Success" if rc == 0 else f"Failed (exit {rc})"
        self.log_box.append(f"● {msg}\n") 
        # dialog only on success; always re-enable run button 
        if rc == 0:
            QMessageBox.information(self, "BLAST finished",
                                    f"Hits table written:\n{out}")
        self.run_btn.setEnabled(True)
        self._worker = None
        self._thread = None 
    # closeEvent ----------------------------
    def closeEvent(self, event):
        """
        Prevent the window from clsoing while BLAST thread is still running.
        """
        if getattr(self, "_thread", None) and self._thread.isRunning():
            # ask the thread to finish, then auto-close the window
            self._worker.finished.connect(lambda _: self.close())
            self._thread.requestInterruption()    # politely signal
            event.ignore()                        # keep window open
            self.log_box.append("⚠ Waiting for BLAST thread to finish…")
        else:
            event.accept()


# Custom QT Logging Handler --------------------------------
class _QtHandler(logging.Handler):
    """Relay root-logger records into the GUI log pane.""" 
    def __init__(self, box: QTextEdit):
        super().__init__(level=logging.INFO)
        self._box = box 
        self.setFormatter(logging.Formatter("%(asctime)s  %(levelname)s:  %(message)s"))

    def emit(self, record):
        # queued connection ensures thread-safety
        QMetaObject.invokeMethod(
            self._box, "append", Qt.QueuedConnection,
            Q_ARG(str, self.format(record))
        )

# Application entry point ------------------
def launch():
    app = QApplication(sys.argv) 
    win = MainWindow() 
    win.show() 
    sys.exit(app.exec())

# allow: python -m microseq_tests.gui 
if __name__ == "__main__":
    launch() 

        


