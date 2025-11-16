# src/microseq_tests/gui

from __future__ import annotations
import logging
L = logging.getLogger(__name__)
import os
import inspect
import platform
import tracemalloc 
import atexit

# Pick the first visible backend that has a server.........
if "WSL_DISTRO_NAME" in os.environ:  # running inside WSL
    if os.environ.get("WAYLAND_DISPLAY"):  # WSLg comositor up for windows 11
        # here let Qt auto-select 'wayland' which is best performance
        pass
    elif os.environ.get("DISPLAY"):  # using external X-server/ Xwayland
        os.environ.setdefault("QT_QPA_PLATFORM", "xcb")
    else:  # headless CI, SSH
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

elif platform.system() == "Linux":  # native Linux desktop
    os.environ.setdefault("QT_QPA_PLATFORM", "xcb")  # native X11

    # macOS defaults to its own so no need for one here...... it would be
    # 'Darwin' and defaults to 'cocoa' just incase for reference and i forget

import sys, logging, traceback, subprocess, shlex, time
from pathlib import Path
from typing import Optional
import collections

from PySide6.QtCore import (
    QLine, Qt, QObject, QThread, Signal, Slot, QMetaObject, Q_ARG, QSettings, QTimer, QModelIndex
)
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QFileDialog, QVBoxLayout,
    QHBoxLayout, QPushButton, QLabel, QListView, QSpinBox, QMessageBox,
    QComboBox, QProgressBar, QCheckBox, QGroupBox, QRadioButton, QAbstractItemView, QLineEdit 
)
from PySide6 import QtCore

# ==== MicroSeq wrappers ----------
from microseq_tests.pipeline import (
    run_blast_stage,
    run_postblast,
    run_full_pipeline,
)
from microseq_tests.utility.utils import setup_logging, load_config


class LogBridge(QObject, logging.Handler):
    sig = Signal(str)

    def __init__(self, parent):
        super().__init__(parent)
        logging.Handler.__init__(self, logging.INFO)
        self.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s"))

    def emit(self, record):
        self.sig.emit(self.format(record))


# Logging class
class LogModel(QtCore.QAbstractListModel): # defining a new class
    MAX_ROWS = 200_000 # the tail cap

    def __init__(self, parent=None):
        super().__init__(parent) # calling constructor of parent class making sure Qt is setup correctly
        self._lines = collections.deque() # using deque here - 0(1) pops

    # mandetory overrides that QListView will ask of my model --- a flat list here vs tree model
    def rowCount(self, _parent=QtCore.QModelIndex()):
        return len(self._lines) # returning number from deque w/ _parent unused

    def data(self, index, role=QtCore.Qt.ItemDataRole.DisplayRole): # data to display specific item
        if role == QtCore.Qt.ItemDataRole.DisplayRole: # just display the text by first checking
            return self._lines[index.row()] # fetch line text
        return None

    # public API
    @QtCore.Slot(str)
    def append(self, line: str):
        if len(self._lines) >= self.MAX_ROWS: # overflow? drop oldest
            self.beginRemoveRows(QtCore.QModelIndex(), 0, 0)
            self._lines.popleft()
            self.endRemoveRows() # Finished removing can update again

        row = len(self._lines) # figure row end of current list place new row
        self.beginInsertRows(QtCore.QModelIndex(), row, row) # insert row end
        self._lines.append(line)
        self.endInsertRows() # Done inserting can update now


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
        if "on_stage" in params:
            self._kwargs["on_stage"] = self.status.emit
        if "on_progress" in params:
            self._kwargs["on_progress"] = self.progress.emit

        try:
            logging.info("%s started", self._fn.__name__)
            result = self._fn(*self._args, **self._kwargs)
            logging.info("%s finished", self._fn.__name__)
        except Exception as e:
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

        # -------- Assembly Mode ----------------------------------------
        # createa dropdown menu or combo box widget 
        self.mode_combo = QComboBox() 

        # add single option to the dropdown, storing single as internal data 
        self.mode_combo.addItem("Single", userData="single")

        # add the Paired option to dropdown, storing paired as internal data 
        self.mode_combo.addItem("Paired", userData="paired")

        # Create a single-line text input field for the forward read pattern 
        self.fwd_pattern_edit = QLineEdit() 

        # Add placeholder as a hint to uder that the field is optional its meant for regex that don't follow patterns I implemented 
        self.fwd_pattern_edit.setPlaceholderText("Forward regex (optional) use this for cases outside of automatic parsing of F/R")

        # Create a single-line text input field for the reverse read pattern 
        self.rev_pattern_edit = QLineEdit() 

        # Add placeholder text as a hint to user that the field is optional again for patterns not in the ones I set up via 'help' for common primers used 
        self.rev_pattern_edit.setPlaceholderText("Reverse regex (optional) use this for cases outside of automatic parsing of F/R")

        # ---- Restore previous choises from settings --------

        # Retrieve previous saved assembly mode from application settings, defaulting to "single" if not found 
        saved_mode = self.settings.value("assembly_mode", "single")

        # Find the infex of the saved mode in combo box data, ensuring valid index >= 0 is returned 
        idx = max(0, self.mode_combo.findData(saved_mode))

        # Set dropdown menu to display the saved/default choice 
        self.mode_combo.setCurrentIndex(idx) 

        # Immediately call function to update the UI based on intial mode selection (example: show/hide reverse pattern field) 
        self._on_mode_changed(idx) 

        # ------ Connect UI events to functions/settings storage -------- 

        # Connect signal emitted when dropdown index changes to the function handles the change 
        self.mode_combo.currentIndexChanged.connect(self._on_mode_changed)

        # Set the text of forward pattern field to the value saved in settings (defaulting to empty string) 
        self.fwd_pattern_edit.setText(self.settings.value("fwd_pattern", ""))

        # Set the text of reverse pattern field to the value saved in settings (defaulting to empty string) 
        self.rev_pattern_edit.setText(self.settings.value("rev_pattern", ""))

        # Connect the signal for text edits in the forward field to an anonymous function that saves the new text to settings. 
        self.fwd_pattern_edit.textEdited.connect(
            lambda txt: self.settings.setValue("fwd_pattern", txt)
        )

        # Connect the signal for text edits in the reverse field to an anonymous function that saves the new text to settings 
        self.rev_pattern_edit.textEdited.connect(
            lambda txt: self.settings.setValue("rev_pattern", txt)
        )

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

        # ----- Model/View Logging ------
        # The data model that holds the log lines
        self.log_model = LogModel(self)

        # The view widget that displays the data from the model
        self.log_view = QListView()
        self.log_view.setModel(self.log_model) # crucial link between the two

        # Performance tweak tells view all rows are the same height
        self.log_view.setUniformItemSizes(True)

        # connect to model's signal to auto-scroll to the bottom on new logs
        self.log_model.rowsInserted.connect(
            lambda _parent, _first, last: # _ to mark unused
                # Safely queue a call to the scrollTo method to run after the current            # event (the paint event) has finished. THis helps in avoiding 
                # and preventing re-entrancy crashes.
                QTimer.singleShot(
                    0,
                    lambda row=last: # capture the value 
                       self.log_view.scrollTo(
                           self.log_model.index(row, 0),
                           QAbstractItemView.ScrollHint.PositionAtBottom
                       ) 

                )
        )

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

        mid.addWidget(QLabel("Assembly"))
        mid.addWidget(self.mode_combo)
        mid.addWidget(self.fwd_pattern_edit)
        mid.addWidget(self.rev_pattern_edit) 

        mid.addWidget(QLabel("Threads"))
        mid.addWidget(self.threads_spin)

        mid.addWidget(mode_box) # setting this here so its left aligned

        mid.addWidget(self.progress)

        mid.addWidget(self.qc_btn)
        mid.addWidget(self.full_btn)
        mid.addWidget(self.postblast_btn)
        mid.addWidget(self.biom_chk)
        mid.addWidget(self.meta_btn)

        mid.addStretch() # pushes Run button to the far right here
        mid.addWidget(self.cancel_btn)
        mid.addWidget(self.run_btn)

        # vertical stack picker row, settings row, then logpane expands
        outer = QVBoxLayout()
        outer.addLayout(top)
        outer.addLayout(mid)
        outer.addWidget(self.log_view)

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
        self._worker: Optional[Worker] = None
        self._log_handler: Optional[LogBridge] = None
        setup_logging(level=logging.INFO) # file + console stderr + GUI output
        root_logger = logging.getLogger() # grab singleton root logger
        self._log_handler = LogBridge(self)
        self._log_handler.sig.connect(self._queue_log, QtCore.Qt.ConnectionType.QueuedConnection) # GUI-safe handler
        root_logger.addHandler(self._log_handler) # Plug into new root logger

    # ---- file picker --------------------------
    def _on_mode_changed(self, index: int):
        """
        This function is the event handler that runs every time the user changes the selection in the 'mode_combo' dropdown menu.
        """
        mode = self.mode_combo.itemData(index) 
        is_paired = mode == "paired" 
        for edit in (self.fwd_pattern_edit, self.rev_pattern_edit):
            edit.setVisible(is_paired)
            edit.setEnabled(is_paired)
        self.settings.setValue("assembly_mode", mode)

    def _assembly_kwargs(self) -> dict:
        """
        This function is meant to be used to gather all assembly configuration data into a dictionary.
        """ 
        mode = self.mode_combo.currentData()
        if mode == "paired":
            fwd = self.fwd_pattern_edit.text().strip() or None 
            rev = self.rev_pattern_edit.text().strip() or None 
        else:
            fwd = rev = None 
        return {"mode": mode, "fwd_pattern": fwd, "rev_pattern": rev} 


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
                self, "Select input folder", str(Path.home())
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
            self, "Select BLAST hits", str(Path.home()), "Hits TSV (*.tsv);;All files (*)"
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
            "Metadata files (*.csv *.tsv *.txt *.xlsx);;All files (*)",
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

        self.progress.setValue(0) # resets progress bar during each run
        task = self.settings.value("blast_task", "megablast")

        # derive output file beside input; disables button; logs starts
        hits_path = self._infile.with_suffix(".hits.tsv")

        self.run_btn.setEnabled(False)
        self.postblast_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.log_model.append(f"\n▶ BLAST {self._infile.name} -> {hits_path.name}")

        # worker and thread wiring -------------------
        worker = Worker(
            run_blast_stage,
            self._infile,
            self.db_box.currentText(), # default DB key selection
            hits_path,
            identity=self.id_spin.value(),
            qcov=self.qcov_spin.value(),
            max_target_seqs=self.hits_spin.value(),
            threads=self.threads_spin.value(),
            blast_task=task,
        )
        self._t0 = time.time()
        self._last_pct = -1   

        worker.progress.connect(self._progress_slot, type=QtCore.Qt.ConnectionType.QueuedConnection) # UI repaint safe

        worker.status.connect(self._stage_slot, type=QtCore.Qt.ConnectionType.QueuedConnection) 
        # injecting wrapper + args into Worker; moves into new thread
        thread = QThread()
        worker.moveToThread(thread)

        # wire signal between log streaming, completion and thread shutdown
        thread.started.connect(worker.run)
        # remember results so _finalise_ui can inspect them
        worker.finished.connect(self._remember_result, type=QtCore.Qt.ConnectionType.QueuedConnection)
        worker.finished.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater) # let Qt free QThread avoid manually doing it
        thread.finished.connect(self._on_job_done)
        # debugging thread cleanup 
        thread.finished.connect(lambda: logging.debug("QThread has finished")) 

        # keep references so closeEvent can see them
        self._worker = worker
        self._thread = thread
        thread.start()

    # ---------- generic launcher for QC / full pipeline -------------------
    def _launch(self, fn, *args, **kw):
        """Start `fn` inside a Worker+QThread and wire its signals to the GUI."""
        self.progress.setValue(0)
        for b in (self.qc_btn, self.full_btn, self.run_btn, self.postblast_btn):
            b.setEnabled(False)
        self.cancel_btn.setEnabled(True)

        worker = Worker(fn, *args, **kw)
        thread = QThread()
        # Worker here to the thread before connecting signals 
        worker.moveToThread(thread)


        # Local functions for the slots ---- 
        self._t0 = time.time() 
        self._last_pct = -1 

        # Connect signals to new, thread-safe GUI only @Slot methods 
        worker.progress.connect(self._progress_slot, type=QtCore.Qt.ConnectionType.QueuedConnection)
        worker.status.connect(self._stage_slot, type=QtCore.Qt.ConnectionType.QueuedConnection)
        # remember the result object
        worker.finished.connect(self._remember_result, type=QtCore.Qt.ConnectionType.QueuedConnection)
        worker.finished.connect(thread.quit) # ask thread to exit
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._on_job_done) # run _done when thread is gone 
        # for debugging thread cleanup
        thread.finished.connect(lambda: logging.debug("QThread has finished.")) 

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
            metadata=None,   # Trim → Convert → BLAST → Tax
            blast_task=task,
            **self._assembly_kwargs()
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
            metadata=self.meta_path,      # None or Path run the Post-BLAST stage too
            blast_task=task,
            **self._assembly_kwargs()
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
        self.log_model.append(
            f"\n▶ Post-BLAST {self.hits_path.name} -> {out_biom.name}"
        )

        worker = Worker(
            run_postblast, self.hits_path, self.meta_path, out_biom
        )
        thread = QThread()
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.finished.connect(self._remember_result, type=QtCore.Qt.ConnectionType.QueuedConnection)
        worker.finished.connect(thread.quit)
        thread.finished.connect(worker.deleteLater) # free QObject Worker
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._on_job_done)

        self._worker = worker
        self._thread = thread
        thread.start()

    # Helper method using here for batch log signals from worker threads
    @Slot(str) 
    def _stage_slot(self, msg: str) -> None: 
        """Update the status bar text - runs in the GUI thread."""
        self._current_stage = msg 
        self.statusBar().showMessage(msg) 

    @Slot(int)
    def _progress_slot(self, pct: int) -> None:
        """Update progress bar + ETA - runs in the GUI thread."""
        if not hasattr(self, "_t0"): # initialize timer on first progress signal
            self._t0 = time.time() 
        if pct == getattr(self, "_last_pct", -1):
            return 
        self._last_pct = pct 
        self.progress.setValue(pct) 
        if pct > 0:
            eta = (time.time() - self._t0) * (100 - pct) / pct 
            self.statusBar().showMessage(
                f"{self._current_stage} {pct}% - ETA "
                f"{int(eta//60)} m {int(eta%60)} s" 
            ) 


    @Slot(object)
    def _remember_result(self, result):
        """Runs in the GUI thread -> safe to touch self.* attributes."""
        self._last_result = result

    @Slot(str)
    def _queue_log(self, line: str):
        """Safely asks the GUI thread to call the Model's "append" method.
        This is the bridge from Python Logger to the GUI model."""
        QMetaObject.invokeMethod(
            self.log_model, "append",
            QtCore.Qt.ConnectionType.QueuedConnection,
            Q_ARG(str, line))

    def _cancel_run(self):
        """Request interruption of the running worker thread."""
        if self._thread and self._thread.isRunning():
            self._thread.requestInterruption()
            self.cancel_btn.setEnabled(False)
            self.log_model.append("Cancelling…")

    @Slot()
    def _on_job_done(self):
        # let any pending paint events drains before closing widgets
        QTimer.singleShot(0, self._safe_cleanup)

    def _safe_cleanup(self):
        self._thread = None
        self._finalise_ui() # old _done body

    # Called when the background QThread has fully finished.
    # Read self._last_result (set in _launch / _launch_blast),
    # decide whether the run was a success, update the GUI, and
    # clean up Worker + QThread objects.
    def _finalise_ui(self):
        """
        Runs after the worker thread has exited. 
        All UI mutations are queued with single Shot(0, ..) so they execute only
        after Qt has finished any pending paint events. 
        """
        # interpret the saved result --------------------------------
        result = getattr(self, "_last_result", 1)      # default = error

        # --- decide status/message -------------------------------- 
        if isinstance(result, dict):                   # full‑pipeline success
            rc, out = 0, result.get("tax", Path())        # pick any key to show
        elif isinstance(result, RuntimeError) and str(result) == "Cancelled":
            rc, out = None, Path()
        elif isinstance(result, Exception):            # Worker caught an error
            rc, out = 1, Path()
            err = str(result) 

            # friendlier message for the "no hits" situation
            if "no blast hits" in err.lower(): 
                dialog_fn = lambda: QMessageBox.information(
                    self, "Nothing to summarise",
                    ("BLAST finished, but no hits met the filters.\n"
                     "You can:\n"
                     " lower Identity / Q-cov, or\n"
                     "run again without the BIOM option.")
                )
            else:
                dialog_fn = lambda e=err: QMessageBox.warning(self, 
                    "Run Failed", e)
            QTimer.singleShot(0, dialog_fn)
        else:                                          # int from BLAST‑only path
            rc, out = int(result), Path()

        # log + optional dialog -------------------------------------
        msg = "Cancelled" if rc is None else ("Success" if rc == 0 else f"Failed (exit {rc})") 
            
        QTimer.singleShot(
            0, lambda m=msg: self.log_model.append(f"● {m}\n")
        )

        if rc == 0 and out:
            QTimer.singleShot(
                0,
                lambda p=out: QMessageBox.information(
                self,
                "Pipeline finished",
                f"Last output file:\n{p}") 
            )

        # re‑enable buttons -----------------------------------------
        def _reenable():

           for b in (self.qc_btn, self.full_btn, 
                     self.run_btn, self.postblast_btn):
               b.setEnabled(True)
           self.cancel_btn.setEnabled(False)

           # safe Qt clean‑up ------------------------------------------
           self._worker = None
           self._thread = None

        QTimer.singleShot(0, _reenable) 

    # closeEvent ----------------------------
    def closeEvent(self, event):
        """
        Prevent the window from clsoing while BLAST thread is still running.
        """
        if self._thread and self._thread.isRunning():
            if self._worker:
                # ask the thread to finish, then auto-close the window
                self._worker.finished.connect(lambda *_: self.close())
                self._thread.requestInterruption()     # politely signal
                event.ignore()                         # keep window open
                self.log_model.append("Waiting for BLAST thread to finish…")
                return
        event.accept()
        root = logging.getLogger()
        if self._log_handler and self._log_handler in root.handlers:
            root.removeHandler(self._log_handler)
            self._log_handler = None

# redirects all Qt internal C++ warnings here for debugging 
def qt_message_handler(mode, _context, message):
    """Redirect Qt messages to the Python logging module."""
    level = {
        QtCore.Qt.QtMsgType.QtDebugMsg:    logging.DEBUG,
        QtCore.Qt.QtMsgType.QtInfoMsg:     logging.INFO,
        QtCore.Qt.QtMsgType.QtWarningMsg:  logging.WARNING,
        QtCore.Qt.QtMsgType.QtCriticalMsg: logging.ERROR,
        QtCore.Qt.QtMsgType.QtFatalMsg:    logging.CRITICAL,
    }.get(mode, logging.CRITICAL)
    logging.log(level, "Qt: %s", message)


# Application entry point ------------------
def launch():
    tracemalloc.start()

    @atexit.register 
    def print_top10():
        snapshot = tracemalloc.take_snapshot()
        top = snapshot.statistics("lineno")[:10]
        logging.info("Top 10 memory allocations:\n%s", "\n".join(map(str, top)))

    app = QApplication(sys.argv)
    # So Qt Warnings show up such as QPainter inside the application log view 
    QtCore.qInstallMessageHandler(qt_message_handler) 
    win = MainWindow()
    win.show()
    sys.exit(app.exec())

# allow: python -m microseq_tests.gui
if __name__ == "__main__":
    launch()
