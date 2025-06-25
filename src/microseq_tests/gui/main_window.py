# src/microseq_tests/gui

from __future__ import annotations 
import logging
L = logging.getLogger(__name__)
import os
import inspect
import platform 
# Pick the first visible backend that has a server......... 
if "WSL_DISTRO_NAME" in os.environ: # running inside WSL 
    if os.environ.get("WAYLAND_DISPLAY"):  # WSLg comositor up for windows 11 
        # here let Qt auto-select 'wayland' which is best performance 
        pass 
    elif os.environ.get("DISPLAY"): # using external X-server/ Xwayland 
        os.environ.setdefault("QT_QPA_PLATFORM", "xcb") 
    else:                                            # headless CI, SSH 
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

elif platform.system() == "Linux":  # native Linux desktop 
    os.environ.setdefault("QT_QPA_PLATFORM", "xcb") # native X11 

    # macOS defaults to its own so no need for one here...... it would be 
    # 'Darwin' and defaults to 'cocoa' just incase for reference and i forget 

import sys, logging, traceback, subprocess, shlex, time
from pathlib import Path 
from typing import Optional


from PySide6.QtCore import (
        Qt, QObject, QThread, Signal, Slot, QMetaObject, Q_ARG, QSettings, QTimer)

from PySide6.QtWidgets import (
        QApplication, QMainWindow, QWidget, QFileDialog, QVBoxLayout,
        QHBoxLayout, QPushButton, QLabel, QTextEdit, QSpinBox, QMessageBox, QComboBox, QProgressBar, QCheckBox, QGroupBox, QRadioButton, 
        )
# log spam protection to avoid seg fault crash prevent worker thread emit log/progress faster than GUI thread can finish painting
# Qt-safe ring-buffer + logging 
from collections import deque 

# ==== MicroSeq wrappers ---------- 
from microseq_tests.pipeline import (
        run_blast_stage, 
        run_trim,
        run_assembly,
        run_add_tax,
        run_postblast,
        run_full_pipeline, 
        )
from microseq_tests.utility.utils import setup_logging, load_config 

# Logging class 
class LogBridge(QObject, logging.Handler):
    """Qt-aware logging.Handler: every record forwarded to a queued QT signal.
    The signal is auto-disconnected when the parent (MainWindow) vanishes."""
    sig = Signal(str)

    def __init__(self, parent):
        QObject.__init__(self, parent)
        logging.Handler.__init__(self, logging.INFO)
        self.setFormatter(logging.Formatter("%(asctime)s  %(levelname)s:  %(message)s")) 

    def emit(self, record):
        self.sig.emit(self.format(record)) 

# Worker class ---------------- 
class Worker(QObject):
    """Background runner living in its own QThread redirects every status message to the Python logging framework instead of a Qt signal this way a single log channel for whole application.""" 
    finished = Signal(object) # exit-code 0 = success
    log = Signal(str) 
    progress = Signal(int)
    status = Signal(str) # emits TRIM, BLASt etc... 

    # Which stage here to run and its kwargs are injected at construction 
    def __init__(self, fn, *args, **kwargs):
        super().__init__()
        self._fn = fn 
        self._args = args 
        self._kwargs = kwargs


    @Slot() # design to warn if any errors occur 
    def run(self):
        # Drop any duplicate that might have been supplied on mistake
        self._kwargs.pop("on_stage", None) 
        self._kwargs.pop("on_progress", None) # guard

        params = inspect.signature(self._fn).parameters 
        if "on_stage"    in params:
            self._kwargs["on_stage"]    = self.status.emit
        if "on_progress" in params:
            self._kwargs["on_progress"] = self.progress.emit


        try:
            self.log.emit(f"{self._fn.__name__} started")
            logging.info("%s started", self._fn.__name__)
            result = self._fn(
                *self._args,
                **self._kwargs,
            )
            self.log.emit(f"{self._fn.__name__} finished") 
            logging.info("%s finished", self._fn.__name__)
        except Exception as e:
            self.log.emit(traceback.format_exc()) 
            logging.exception("Worker crashed:")
            result = e
        self.finished.emit(result)   # dict | int | Exception
                 


# ---- Main Window Constructor 
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.settings = QSettings("MicroSeq", "MicroSeq") 
        self.setWindowTitle("MicroSeq GUI v1.0") # title for it 
        self.resize(800, 520) # size of app considering also log space here 

        # widgets --------------------------------------------------------
        self.fasta_lbl = QLabel("Input: —")
        browse_btn = QPushButton("Browse..")
        browse_btn.setToolTip("Select FASTA/FASTQ/AB1 file(s) or a folder")
        browse_btn.clicked.connect(self._choose_infile)

        # BLAST hits TSV picker
        self.hits_lbl = QLabel("Hits / post-BLAST TSV:")
        hits_btn = QPushButton("Select hits_tax.tsv for post-BLAST")
        hits_btn.setToolTip(
            "Pick an existing BLAST-result table hits.tsv or hits_tax.tsv to run Post-BLAST and/or make a BIOM table."
        )

        hits_btn.clicked.connect(self._choose_hits)

        # label place holder and browse button wired to file picker 
        self.id_spin = QSpinBox()
        self.id_spin.setRange(50, 100) # for simplicity I set spinbox ID % 50-100 threshold 
        self.id_spin.setValue(97) # identity 97% default 
        self.id_spin.setSuffix(" % ID")

        self.threads_spin = QSpinBox()
        self.threads_spin.setRange(1, 32)
        self.threads_spin.setValue(4)
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

        # ---- Alignmnet mode box ----------------
        mode_box = QGroupBox("Alignment mode")
        fast_rb = QRadioButton("Fast (megablast)")
        slow_rb = QRadioButton("Comprehensive - sensitivity (blastn)") 
        vbox = QVBoxLayout(mode_box)
        vbox.addWidget(fast_rb)
        vbox.addWidget(slow_rb)

        # restore previous choice (default = megablast)
        task = self.settings.value("blast_task", "megablast")
        (fast_rb if task == "megablast" else slow_rb).setChecked(True) 

        # save whever user toggles 
        fast_rb.toggled.connect(
            lambda on: on and self.settings.setValue("blast_task", "megablast"))
        slow_rb.toggled.connect(
            lambda on: on and self.settings.setValue("blast_task", "blastn")) 

        # ------- progress bar ----------
        self.progress = QProgressBar()
        self.progress.setRange(0, 100)
        self.progress.setValue(0)

        # here you just run blast ....
        self.run_btn = QPushButton("Run Blast")
        self.run_btn.clicked.connect(self._launch_blast)
        self.cancel_btn = QPushButton("Cancel")
        self.cancel_btn.setEnabled(False)
        self.cancel_btn.clicked.connect(self._cancel_run)

        # here you click to run QC or full pipeline 
        self.qc_btn = QPushButton("Run QC")
        self.full_btn = QPushButton("Full pipeline")
        self.postblast_btn = QPushButton("Post-BLAST")
        self.qc_btn.clicked.connect(self._launch_qc)
        self.full_btn.clicked.connect(self._launch_full)
        self.postblast_btn.clicked.connect(self._launch_postblast)
        

        # creating a checkbox here for post blast 
        self.biom_chk = QCheckBox("Make BIOM")
        self.biom_chk.setChecked(False)

        self.meta_btn = QPushButton("Browse metadata...")
        self.meta_btn.clicked.connect(self._choose_metadata) 

        # lets you scroll log output and nowarp keeps long commadn lines intact 
        self.log_box = QTextEdit(readOnly=True, lineWrapMode=QTextEdit.NoWrap)
        
        # Buffered Log output: avoid thousands of QTextEdit repaints 
        self._log_buf: deque[str] = deque(maxlen=200)
        self._flush_timer = QTimer(self, interval=80) # set at 80 fewer repaints 
        self._flush_timer.timeout.connect(self._flush_log)
        self._flush_timer.start() 

        # ---- Layout here will update UX -------------------------
        top = QHBoxLayout()
        top.addWidget(self.fasta_lbl)
        top.addWidget(browse_btn)
        top.addWidget(self.hits_lbl)
        top.addWidget(hits_btn)

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

        mid.addWidget(mode_box) # setting this here so its left aligned  

        mid.addWidget(self.progress)

        mid.addWidget(self.qc_btn)
        mid.addWidget(self.full_btn)
        mid.addWidget(self.postblast_btn)
        mid.addWidget(self.biom_chk)

        mid.addWidget(self.meta_btn) 
        
        mid.addStretch()  # pushes Run button to the far right here
        mid.addWidget(self.cancel_btn)
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
        self.hits_path: Optional[Path] = None
        self.meta_path: Optional[Path] = None
        self._current_stage: str = ""

        # placeholders to prevent AttributeErrors before first run
        self._thread: Optional[QThread] = None
        self._worker: Optional[QObject] = None 
        setup_logging(level=logging.INFO)  # file + console stderr + GUI output 
        root_logger = logging.getLogger() # grab singleton root logger
        self._log_handler = LogBridge(self) 
        self._log_handler.sig.connect(self._queue_log, Qt.QueuedConnection)  # GUI-safe handler 
        root_logger.addHandler(self._log_handler) # Plug into new root logger 
    
    # ---- file picker --------------------------
    def _choose_infile(self):

        """Select FASTA/FASTQ/AB1 file(s) or a folder of traces."""


        paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Select FASTA/FASTQ/AB1 file(s)",
            str(Path.home()),

            "Seq files (*.fasta *.fastq *.ab1);;All files (*)",
        )

        # If the user cancelled, offer a directory chooser instead

        if not paths:
            dir_path = QFileDialog.getExistingDirectory(
                self,
                "Select input folder",
                str(Path.home()),
            )
            if dir_path:
                self._infile = Path(dir_path)
        else:

            # multi-file selection is allowed but only the first file is used
            self._infile = Path(paths[0])


        if self._infile:
            label = (
                f"Input folder: {self._infile.name}"
                if self._infile.is_dir()
                else f"Input file: {self._infile.name}"
            )
            self.fasta_lbl.setText(label)
            # clean any previous metadata selection
            self.meta_path = None
            self.statusBar().clearMessage()

    # ------- choosing BLAST hits file --------------
    def _choose_hits(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select BLAST hits",
            str(Path.home()),
            "Hits TSV (*.tsv);;All files (*)",
        )
        if path:
            self.hits_path = Path(path)
            self.hits_lbl.setText(f"Hits: {self.hits_path.name}")
     
    # ------- chosing metadata file -------------- 
    def _choose_metadata(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select sample metadata",
            str(Path.home()),
            "Metadata files (*.csv *.tsv *.txt *.xlsx);;All files (*)"
        )
        if path:
            self.meta_path = Path(path)
            self.statusBar().showMessage(f"Metadata: {self.meta_path.name}")
    
    
    # ---- Blast Stage Demo ------------------------------
    # guarding against accidental click 
    def _launch_blast(self):
        if not self._infile:
            QMessageBox.warning(self, "No input", "Choose a file or folder first.")
            return

        self.progress.setValue(0)  # resets progress bar during each run

        task = self.settings.value("blast_task", "megablast") 

        # derive output file beside input; disables button; logs starts 
        hits_path = self._infile.with_suffix(".hits.tsv")

        self.run_btn.setEnabled(False)
        self.postblast_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
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
                blast_task=task,
                )
        # send Worker.log -> root logger -> _QtHandler -> GUI 
        worker.log.connect(lambda s: logging.info("%s", s), type=Qt.QueuedConnection) 
        t0 = time.time()
        _last_pct = -1 

        def _progress_with_eta(pct: int) -> None:
            nonlocal _last_pct 
            if pct == _last_pct:
                return 
            _last_pct = pct 

            self.progress.setValue(pct)
            if pct:
                eta = (time.time() - t0) * (100 - pct) / pct
                self.statusBar().showMessage(
                    f"BLAST {pct}% – ETA {int(eta//60)} m {int(eta%60)} s"
                )

        worker.progress.connect(_progress_with_eta, type=Qt.QueuedConnection) # UI repaint safe 
        # injecting wrapper + args into Worker; moves into new thread 
        thread = QThread()
        worker.moveToThread(thread)

        # wire signal between log streaming, completion and thread shutdown

        thread.started.connect(worker.run)
        # remember results so _finalise_ui can inspect them
        worker.finished.connect(self._remember_result, type=Qt.QueuedConnection)

        worker.finished.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater) # let Qt free QThread avoid manually doing it
        thread.finished.connect(self._on_job_done)

        # keep references so closeEvent can see them 
        self._worker = worker 
        self._thread = thread 

        print("DEBUG – about to start Worker with", self.db_box.currentText(),
        self.id_spin.value(), self.qcov_spin.value(),
        self.hits_spin.value(), self.threads_spin.value(), flush=True)
        thread.start() 
    
    # ---------- generic launcher for QC / full pipeline -------------------
    def _launch(self, fn, *args, **kw):
        """Start `fn` inside a Worker+QThread and wire its signals to the GUI."""
        self.progress.setValue(0)
        for b in (self.qc_btn, self.full_btn, self.run_btn, self.postblast_btn):
            b.setEnabled(False)
        self.cancel_btn.setEnabled(True)

        worker = Worker(fn, *args, **kw)
        worker.log.connect(lambda s: logging.info("%s", s), type=Qt.QueuedConnection)
        thread = QThread()
        worker.moveToThread(thread)

        t0 = time.time()
        _last_pct = -1 

        def _stage(msg: str) -> None:
            self._current_stage = msg
            self.statusBar().showMessage(msg)

        def _progress_with_eta(pct: int) -> None:
            nonlocal _last_pct 
            if pct == _last_pct:
                return 
            _last_pct = pct 

            self.progress.setValue(pct)
            if pct:
                eta = (time.time() - t0) * (100 - pct) / pct
                self.statusBar().showMessage(
                    f"{self._current_stage} {pct}% – ETA {int(eta//60)} m {int(eta%60)} s"
                )

        worker.progress.connect(_progress_with_eta, type=Qt.QueuedConnection)
        worker.status.connect(_stage, type=Qt.QueuedConnection)
        # remember the result object 
        worker.finished.connect(self._remember_result, type=Qt.QueuedConnection)
        
        worker.finished.connect(thread.quit) # ask thread to exit
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._on_job_done) # run _done when thread is gone

        self._worker = worker
        self._thread = thread
        thread.started.connect(worker.run)
        thread.start()
    

    # -------- Trim -> Convert -> BLAST -> Taxonomy ---------------
    def _launch_qc(self):
        if not self._infile:
            QMessageBox.warning(self, "No input", "Choose a file or folder first.")
            return

        task = self.settings.value("blast_task", "megablast") 

        self._launch(
            run_full_pipeline,
            self._infile,
            self.db_box.currentText(),
            threads=self.threads_spin.value(),
            postblast=self.biom_chk.isChecked(),
            metadata=None,        # Trim → Convert → BLAST → Tax
            blast_task=task, 
        )
    
    # ------ Run full pipeline with Post-Blast as well ------------- 
    def _launch_full(self):
        if not self._infile:
            QMessageBox.warning(self, "No input", "Choose a file or folder first.")
            return
        # if user asked for BIOM but hasn't chosen metadata, abort early 
        if self.biom_chk.isChecked() and not self.meta_path:
            QMessageBox.warning(
                self, "No metadata",
                "A metadata file is required to build the BIOM table.\n"
                "Click ‘Browse metadata…’ first or un-tick ‘Make BIOM’."
            )
            return
        # launch the pipeline here

        task = self.settings.value("blast_task", "megablast") 

        self._launch(
            run_full_pipeline,
            self._infile,
            self.db_box.currentText(),
            threads=self.threads_spin.value(),
            postblast=self.biom_chk.isChecked(), # source of truth decide via checkbox
            metadata=self.meta_path,         # None or Path run the Post-BLAST stage too
            blast_task=task, 
    )

    # ------ Run stand-alone Post-BLAST ---------------------------------
    def _launch_postblast(self):
        if not self.hits_path:
            QMessageBox.warning(self, "No hits", "Choose a BLAST hits file first.")
            return
        if not self.meta_path:
            QMessageBox.warning(self, "No metadata", "Choose a metadata file first.")
            return

        self.progress.setValue(0)
        out_biom = self.hits_path.with_suffix(".biom")

        self.postblast_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.log_box.append(
            f"\n▶ Post-BLAST {self.hits_path.name} -> {out_biom.name}"
        )

        worker = Worker(
            run_postblast,
            self.hits_path,
            self.meta_path,
            out_biom,
        )
        worker.log.connect(lambda s: logging.info("%s", s), type=Qt.QueuedConnection)
        thread = QThread()
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.finished.connect(self._remember_result, type=Qt.QueuedConnection)
        worker.finished.connect(thread.quit)
        worker.finished.connect(worker.deleteLater) # free QObject Worker 
        worker.finished.connect(thread.deleteLater) 
        thread.finished.connect(self._on_job_done)

        self._worker = worker
        self._thread = thread
        thread.start()

    # Helper method using here for batch log signals from worker threads
    @Slot(object)
    def _remember_result(self, result):
        """Runs in the GUI thread -> safe to touch self.* attributes."""
        self._last_result = result 

    @Slot(str)
    def _queue_log(self, line: str) -> None:
        """Enqueue a log line coming from a worker thread.""" 
        self._log_buf.append(line) 

    @Slot()
    def _flush_log(self) -> None:
        """Flush the ring-buffer into QTextEdit runs in GUI thread""" 
        if not self._log_buf:
            return 

        lines = "\n".join(self._log_buf)
        self._log_buf.clear()

        # NEVER touch a QWidget ddirectly schedule the append for the next
        # GUI-event turn so we cannot run inside its current paintEvent. 
        QMetaObject.invokeMethod(
            self.log_box,
            "append",
            Qt.QueuedConnection,
            Q_ARG(str, lines)

        ) 

    def _cancel_run(self):
        """Request interruption of the running worker thread."""
        if getattr(self, "_thread", None) and self._thread.isRunning():
            self._thread.requestInterruption()
            self.cancel_btn.setEnabled(False)
            self.log_box.append("Cancelling…") 


    @Slot()
    def _on_job_done(self):
        # let any pending paint events drains before closing widgets 
        QTimer.singleShot(0, self._safe_cleanup)

    def _safe_cleanup(self):
        if getattr(self, "_flush_timer", None) and self._flush_timer.isActive():
            self._flush_timer.stop()
        if getattr(self, "_thread", None):
            self._thread = None

        self._finalise_ui() # old _done body



    # Called when the background QThread has fully finished.
    # Read self._last_result (set in _launch / _launch_blast),
    # decide whether the run was a success, update the GUI, and
    # clean up Worker + QThread objects.
    def _finalise_ui(self):
        # interpret the saved result --------------------------------
        result = getattr(self, "_last_result", 1)     # default = error
        if isinstance(result, dict):                  # full‑pipeline success
            rc  = 0
            out = result.get("tax", Path())           # pick any key to show
        elif isinstance(result, RuntimeError) and str(result) == "Cancelled":
            rc  = None
            out = Path()
        elif isinstance(result, Exception):           # Worker caught an error
            rc  = 1
            out = Path()
            # friendlier message for the "no hits" situation 
            if "no blast hits" in str(result).lower():
                QMessageBox.information(
                    self, "Nothing to summarise",
                    ("BLAST finished, but no hits met the filters.\n"
                     "You can:\n"
                     " lower Identity / Q-cov, or\n"
                     "run again without the BIOM option.")
                )
            else:
                QMessageBox.warning(self, "Run failed", str(result)) 
        else:                                         # int from BLAST‑only path
            rc  = int(result)
            out = Path()

        # log + optional dialog -------------------------------------
        if rc is None:
            msg = "Cancelled"
        else:
            msg = "Success" if rc == 0 else f"Failed (exit {rc})"
        self.log_box.append(f"● {msg}\n")

        if rc == 0 and out:
            QMessageBox.information(
                self,
                "Pipeline finished",
                f"Last output file:\n{out}"
            )

        # re‑enable buttons -----------------------------------------
        for b in (self.qc_btn, self.full_btn, self.run_btn, self.postblast_btn):
            b.setEnabled(True)
        self.cancel_btn.setEnabled(False)

        # safe Qt clean‑up ------------------------------------------ 
        self._worker = None
        self._thread = None



    
    # closeEvent ----------------------------
    def closeEvent(self, event):
        """
        Prevent the window from clsoing while BLAST thread is still running.
        """
        if self._thread and self._thread.isRunning():
            if self._worker:

                # ask the thread to finish, then auto-close the window
                self._worker.finished.connect(lambda *_: self.close())
                self._thread.requestInterruption()    # politely signal
                event.ignore()                        # keep window open
                self.log_box.append("Waiting for BLAST thread to finish…")
                return 
            event.accept() 
            root = logging.getLogger() 
            if getattr(self, "_log_handler", None) in root.handlers:
                root.removeHandler(self._log_handler)
                self._log_handler = None 

# Application entry point ------------------
def launch():
    app = QApplication(sys.argv) 
    win = MainWindow() 
    win.show() 
    sys.exit(app.exec())

# allow: python -m microseq_tests.gui 
if __name__ == "__main__":
    launch() 

        


