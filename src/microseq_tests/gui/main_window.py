# src/microseq_tests/gui

from __future__ import annotations
import logging

from numpy import partition
from pandas.core.generic import sample
L = logging.getLogger(__name__)
import os
import inspect
import platform
import tracemalloc 
import tempfile 
import shutil 
import atexit
import re
import pandas as pd 

import sys, logging, traceback, subprocess, shlex, time
from pathlib import Path
from typing import Optional
import collections

from microseq_tests.primer_catalog import pairing_label_sets, build_primer_cfg_override, trim_presets
from microseq_tests.assembly.registry import list_assemblers

PRIMER_SETS: dict[str, tuple[list[str], list[str]]] = pairing_label_sets()
TRIM_PRESETS: list[str] = sorted(trim_presets())

STATUS_COLOR_LEGEND = (
    "Status color legend:\n"
    "• assembled = green\n"
    "• singlets_only = yellow\n"
    "• cap3_no_output = orange\n"
    "• pair_missing = red\n"
    "• overlap_* = light blue"
)

TRACE_QC_COLOR_LEGEND = (
    "Trace QC color legend:\n"
    "• PASS = green\n"
    "• WARN = yellow\n"
    "• FAIL = red\n"
    "• NA = gray"
)


def _normalize_legacy_compare_diag_detail(detail: str, selected_engine: str) -> str:
    """Normalize legacy compare diagnostics where engine was printed as overlap_len."""
    text = (detail or "").strip()
    if not text:
        return text
    match = re.search(r"overlap_len=(\d+);identity=[^;]+;\s*engine=(\d+)$", text)
    if not match:
        return text
    overlap_len, engine_token = match.groups()
    if overlap_len != engine_token:
        return text
    engine_label = (selected_engine or "").strip()
    if not engine_label or engine_label.isdigit():
        return text
    return text[: match.start(2)] + engine_label

def _trace_flags_to_labels(flags: str, status: str, mixed_frac: str, review_reason: str = "") -> list[str]:
    labels: list[str] = []
    flag_set = {f.strip().upper() for f in str(flags or "").split(";") if f.strip()}
    status_norm = str(status or "NA").strip().upper()

    if "LOW_SNR" in flag_set:
        labels.append("Low signal / high background (fail)")
    elif "LOW_SNR_WARN" in flag_set:
        labels.append("Low signal warning")

    if "MIXED_PEAKS" in flag_set:
        labels.append("Mixed peaks detected")

    if flag_set.intersection({"TRACE_DATA_MISSING", "SNR_UNDEFINED", "MIXED_UNDEFINED", "TRACE_PARSE_ERROR"}):
        labels.append("Trace QC unavailable/NA")

    if str(review_reason or "").strip().lower() == "mixture_suspected":
        labels.append("Mixture suspected")

    if mixed_frac:
        labels.append(f"Max mixed-peak fraction: {mixed_frac}")

    if not labels:
        labels.append(f"Trace QC {status_norm}")
    return labels

from PySide6.QtCore import (
    QLine, Qt, QObject, QThread, Signal, Slot, QSettings, QTimer, QModelIndex, QProcess, QUrl
)
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QFileDialog, QVBoxLayout,
    QHBoxLayout, QPushButton, QLabel, QListView, QSpinBox, QMessageBox,
    QComboBox, QProgressBar, QCheckBox, QGroupBox, QRadioButton, QAbstractItemView, QLineEdit,
    QTabWidget, QTableWidget, QTableWidgetItem, QSplitter, QDialog, QPlainTextEdit, QHeaderView, QSizePolicy 
)
from PySide6 import QtCore
from PySide6.QtGui import QColor, QBrush, QDesktopServices 

# ==== MicroSeq wrappers ----------
from microseq_tests.pipeline import (
    run_blast_stage,
    run_postblast,
    run_full_pipeline,
    run_trim,
    run_compare_assemblers,
    _summarize_paired_candidates,
    _suggest_pairing_patterns
)
from microseq_tests.vsearch_tools import resolve_vsearch 
from microseq_tests.assembly.pairing import DupPolicy, PairingCandidate, analyze_pairing_candidates 
from microseq_tests.utility.utils import setup_logging, load_config


class LogBridge(QObject, logging.Handler):
    sig = Signal(str)

    def __init__(self, parent):
        super().__init__(parent)
        logging.Handler.__init__(self, logging.INFO)
        self.setFormatter(logging.Formatter("%(asctime)s %(levelname)s: %(message)s"))

    def emit(self, record):
        self.sig.emit(self.format(record))


def _set_header_tooltips(table: QTableWidget, tooltip_by_header: dict[str, str]) -> None:
    """Attach tooltips to table header labels by exact header text."""
    for col in range(table.columnCount()):
        item = table.horizontalHeaderItem(col)
        if not item:
            continue
        text = item.text().strip()
        tooltip = tooltip_by_header.get(text)
        if tooltip:
            item.setToolTip(tooltip)


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
        if role != QtCore.Qt.ItemDataRole.DisplayRole:
            return None
        if not index.isValid():
            return None
        row = index.row()
        if row < 0 or row >= len(self._lines):
            return None 
        return self._lines[row] # fetch line text 
    # public API
    @QtCore.Slot(str)
    def append(self, line: str):
        self.append_many([line]) 
   
    @QtCore.Slot(list)
    def append_many(self, lines: list[str]) -> None:
        if not lines:
            return 
        total = len(self._lines) + len(lines)
        overflow = max(0, total - self.MAX_ROWS)
        if overflow:
            self.beginRemoveRows(QtCore.QModelIndex(), 0, overflow -1)
            for _ in range(overflow):
                self._lines.popleft()
            self.endRemoveRows()

        row0 = len(self._lines)
        row1 = row0 + len(lines) - 1 
        self.beginInsertRows(QtCore.QModelIndex(), row0, row1)
        self._lines.extend(lines)
        self.endInsertRows() 

TEXT_FILE_EXTENSIONS = {
    ".txt",
    ".log",
    ".info",
    ".fasta",
    ".fa",
    ".fna",
    ".fastq",
    ".fq",
    ".tsv",
    ".csv",
}

def is_wsl() -> bool:
    """Checks if the code is executing inside Windows SubSystem for Linux."""
    return bool(os.environ.get("WSL_INTEROP") or os.environ.get("WSL_DISTRO_NAME"))

def has_gui_session() -> bool:
    """Checking if there an is an active windowing system (x11/Wayland/Win/macOS) to show a UI."""
    if platform.system() in ("Darwin", "Windows"):
        return True 
    return bool(os.environ.get("WAYLAND_DISPLAY") or os.environ.get("DISPLAY"))

def looks_texty(path: Path) -> bool:
    """ 
    Heuristic here I made to determine if a file is plain text. Due to the limit I'm using here of a max of 2MB inorder to deal with not overloading the machine and GUI by not opening huge byte files that are non text otherwise it risks being either not as responsive or worse freezing outright. Will test run this more robustly at a lter date for now the max is 2MB for opening files in the GUI which should be sufficient. 

    This function checks the extension first, then reads 4KB chunk of data to look for null bytes which is a binary marker if it does appear it will not load it in the GUI.    """ 
    if path.suffix.lower() in TEXT_FILE_EXTENSIONS:
        return True 
    try: 
       chunk = path.read_bytes()[:4096]
    except OSError:
        return False 
    return b"\x00" not in chunk

class TextViewerDialog(QDialog):
    def __init__(self, path: Path, parent=None, max_bytes: int = 2_000_000):
        super().__init__(parent)
        self.setWindowTitle(str(path))
        self.resize(900, 600)

        layout = QVBoxLayout(self)

        header = QLabel(f"{path} (showing up to {max_bytes:,} bytes)")
        header.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(header)

        editor = QPlainTextEdit(self)
        editor.setReadOnly(True)
        editor.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        layout.addWidget(editor)

        try:
            data = path.read_bytes()[:max_bytes]
            text = data.decode("utf-8", errors="replace")
        except OSError as exc:
            text = f"Failed to read file: {exc}"
        editor.setPlainText(text)

        btn_row = QHBoxLayout()
        copy_btn = QPushButton("Copy path", self)
        copy_btn.clicked.connect(lambda: QApplication.clipboard().setText(str(path)))
        btn_row.addWidget(copy_btn)
        btn_row.addStretch()

        close_btn = QPushButton("Close", self)
        close_btn.clicked.connect(self.close)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)

class TextBlobDialog(QDialog):
    def __init__(self, title: str, text: str, parent=None):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(800, 500)

        layout = QVBoxLayout(self)

        editor = QPlainTextEdit(self)
        editor.setReadOnly(True)
        editor.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        editor.setPlainText(text)
        layout.addWidget(editor)

        btn_row = QHBoxLayout()
        copy_btn = QPushButton("Copy", self)
        copy_btn.clicked.connect(lambda: QApplication.clipboard().setText(text))
        btn_row.addWidget(copy_btn)
        btn_row.addStretch()
        close_btn = QPushButton("Close", self)
        close_btn.clicked.connect(self.close)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)

class PairingPreviewDialog(QDialog):
    def __init__(
        self,
        *,
        summary: str,
        suggestions: str,
        candidates: list[PairingCandidate],
        enforce_same_well: bool,
        refresh_callback=None,
        open_callback=None,
        parent=None,
    ):
        super().__init__(parent)
        self.setWindowTitle("Pairing Preview")
        self.resize(1100, 650)
        self._candidates = candidates
        self._enforce_same_well = enforce_same_well
        self._refresh_callback = refresh_callback
        self._open_callback = open_callback
        self._renamed_paths: set[Path] = set()
        self._rename_map: list[tuple[Path, Path]] = []

        layout = QVBoxLayout(self)
        self._header = QLabel()
        self._header.setWordWrap(True)
        self._header.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(self._header)
        self._set_header(summary, suggestions, candidates)

        self.table = QTableWidget(0, 7, self)
        self.table.setHorizontalHeaderLabels(
            ["File", "Detected", "Sample ID", "Well", "Issues", "Suggested rename", "Status"]
        )
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table)

        self._populate_table()

        btn_row = QHBoxLayout()
        copy_btn = QPushButton("Copy rename map", self)
        copy_btn.clicked.connect(self._copy_rename_map)
        btn_row.addWidget(copy_btn)

        open_btn = QPushButton("Open containing folder", self)
        open_btn.clicked.connect(self._open_containing_folder)
        btn_row.addWidget(open_btn)

        copy_cmd_btn = QPushButton("Copy rename commands", self)
        copy_cmd_btn.clicked.connect(self._copy_rename_commands)
        btn_row.addWidget(copy_cmd_btn)

        export_btn = QPushButton("Export rename script…", self)
        export_btn.clicked.connect(self._export_rename_script)
        btn_row.addWidget(export_btn)

        apply_btn = QPushButton("Apply suggested renames…", self)
        apply_btn.clicked.connect(self._apply_renames)
        btn_row.addWidget(apply_btn)

        btn_row.addStretch()
        close_btn = QPushButton("Close", self)
        close_btn.clicked.connect(self.close)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)

    def _set_header(
        self,
        summary: str,
        suggestions: str,
        candidates: list[PairingCandidate],
    ) -> None:
        mismatch = any(
            "well mismatch" in issue
            for candidate in candidates
            for issue in candidate.issues
        )
        if self._enforce_same_well and mismatch:
            pairing_line = (
                "Pairs will not assemble until well codes match for each sample ID.\n"
            )
        else:
            pairing_line = (
                "Pairs will still assemble, but naming is non-canonical; suggested renames improve reproducibility.\n"
            )
        self._header.setText(
            f"{summary}\n"
            f"{pairing_line}"
            f"Suggestions: {suggestions}"
        )

    def _populate_table(self) -> None:
        self.table.setRowCount(0)
        self._rename_map = []
        for candidate in self._candidates:
            row = self.table.rowCount()
            self.table.insertRow(row)
            detected = candidate.orient or "-"
            issues = "; ".join(candidate.issues) if candidate.issues else "-"
            suggested = candidate.suggested_name or "-"
            status = "Issue corrected" if candidate.path in self._renamed_paths and not candidate.issues else "Issue" if candidate.issues else "OK"

            items = [
                candidate.path.name,
                detected,
                candidate.sid or "-",
                candidate.well or "-",
                issues,
                suggested,
                status,
            ]
            for col, value in enumerate(items):
                item = QTableWidgetItem(value)
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                self.table.setItem(row, col, item)

            if candidate.suggested_name and candidate.path.name != candidate.suggested_name:
                self._rename_map.append(
                    (candidate.path, candidate.path.with_name(candidate.suggested_name))
                )

    def _copy_rename_map(self) -> None:
        if not self._rename_map:
            QApplication.clipboard().setText("No rename suggestions available.")
            return
        lines = ["source\ttarget"]
        lines.extend(f"{src}\t{dest}" for src, dest in self._rename_map)
        QApplication.clipboard().setText("\n".join(lines))

    def _open_containing_folder(self) -> None:
        if not self._candidates:
            QMessageBox.information(self, "No files", "No files are available to open.")
            return
        folder = self._candidates[0].path.parent
        if self._open_callback is not None:
            ok = self._open_callback(folder)
            if ok:
                return
        QMessageBox.information(
            self,
            "Open manually",
            f"Unable to open folder via system handler:\n{folder}",
        )

    def _copy_rename_commands(self) -> None:
        if not self._rename_map:
            QApplication.clipboard().setText("No rename suggestions available.")
            return
        if os.name == "nt":
            lines = [
                f'Rename-Item -LiteralPath {shlex.quote(str(src))} -NewName {shlex.quote(dest.name)}'
                for src, dest in self._rename_map
            ]
        else:
            lines = [
                f"mv -n -- {shlex.quote(str(src))} {shlex.quote(str(dest))}"
                for src, dest in self._rename_map
            ]
        QApplication.clipboard().setText("\n".join(lines))

    def _rename_conflicts(self) -> tuple[set[Path], set[Path]]:
        targets = [dest for _, dest in self._rename_map]
        conflicts = {dest for dest in targets if dest.exists()}
        target_counts = collections.Counter(targets)
        duplicate_targets = {dest for dest, count in target_counts.items() if count > 1}
        return conflicts, duplicate_targets

    def _export_rename_script(self) -> None:
        if not self._rename_map:
            QMessageBox.information(self, "No renames", "No suggested renames to export.")
            return
        conflicts, duplicate_targets = self._rename_conflicts()
        if conflicts or duplicate_targets:
            details = []
            if conflicts:
                details.append(
                    "Targets already exist:\n" + "\n".join(str(p) for p in sorted(conflicts))
                )
            if duplicate_targets:
                details.append(
                    "Duplicate target names:\n"
                    + "\n".join(str(p) for p in sorted(duplicate_targets))
                )
            QMessageBox.warning(
                self,
                "Rename conflicts",
                "Unable to export renames due to conflicts.\n\n" + "\n\n".join(details),
            )
            return

        default_path = None
        if self._candidates:
            default_path = str(self._candidates[0].path.parent / "rename_map.tsv")

        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export rename script",
            default_path or "",
            "TSV (*.tsv);;Shell script (*.sh);;All files (*)",
        )
        if not path:
            return
        out_path = Path(path)
        if out_path.suffix.lower() == ".sh":
            if os.name == "nt":
                lines = [
                    f'Rename-Item -LiteralPath {shlex.quote(str(src))} -NewName {shlex.quote(dest.name)}'
                    for src, dest in self._rename_map
                ]
            else:
                lines = [
                    f"mv -n -- {shlex.quote(str(src))} {shlex.quote(str(dest))}"
                    for src, dest in self._rename_map
                ]
            out_path.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + "\n".join(lines) + "\n", encoding="utf-8")
        else:
            lines = ["source\ttarget"]
            lines.extend(f"{src}\t{dest}" for src, dest in self._rename_map)
            out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        QMessageBox.information(self, "Export complete", f"Rename script written to:\n{out_path}")

    def _apply_renames(self) -> None:
        if not self._rename_map:
            QMessageBox.information(self, "No renames", "No suggested renames to apply.")
            return
        confirm = QMessageBox.question(
            self,
            "Apply renames",
            "Apply the suggested renames to the files on disk? This cannot be undone.",
        )
        if confirm != QMessageBox.StandardButton.Yes:
            return

        conflicts, duplicate_targets = self._rename_conflicts()
        if conflicts or duplicate_targets:
            details = []
            if conflicts:
                details.append(
                    "Targets already exist:\n" + "\n".join(str(p) for p in sorted(conflicts))
                )
            if duplicate_targets:
                details.append(
                    "Duplicate target names:\n"
                    + "\n".join(str(p) for p in sorted(duplicate_targets))
                )
            QMessageBox.warning(
                self,
                "Rename conflicts",
                "Unable to apply renames due to conflicts.\n\n" + "\n\n".join(details),
            )
            return

        errors: list[str] = []
        for src, dest in self._rename_map:
            try:
                if src.exists():
                    src.rename(dest)
            except OSError as exc:
                errors.append(f"{src} → {dest}: {exc}")

        if errors:
            QMessageBox.warning(
                self,
                "Rename errors",
                "Some renames failed:\n" + "\n".join(errors),
            )
            return

        self._renamed_paths = {dest for _, dest in self._rename_map}
        if self._refresh_callback is not None:
            summary, suggestions, candidates = self._refresh_callback()
            self._candidates = candidates
            self._set_header(summary, suggestions, candidates)
            self._populate_table()
            QMessageBox.information(
                self,
                "Renames complete",
                "Suggested renames applied. Preview refreshed.",
            )
            return

        QMessageBox.information(
            self,
            "Renames complete",
            "Suggested renames applied. Re-run Preview Pairs to refresh results.",
        )
        self.close()


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
        if QThread.currentThread().isInterruptionRequested():
            self.finished.emit(RuntimeError("Cancelled"))
            return 
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
            if QThread.currentThread().isInterruptionRequested():
                raise RuntimeError("Cancelled")
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

        normal_geom = self.settings.value("window_normal_geometry", None)
        if isinstance(normal_geom, (list, tuple)) and len(normal_geom) == 4:
            try:
                x, y, w, h = [int(v) for v in normal_geom]
                self.setGeometry(x, y, w, h)
            except (TypeError, ValueError):
                pass
        else:
            # Legacy migration path from old saveGeometry blob key.
            geom = self.settings.value("window_geometry")
            if geom is not None:
                self.restoreGeometry(geom)

        self._start_maximized = self.settings.value("window_start_maximized", False, type=bool)

        # widgets --------------------------------------------------------
        self.fasta_lbl = QLabel("Input: —")
        browse_btn = QPushButton("Browse..")
        browse_btn.setToolTip("Select FASTA/FASTQ/AB1 file(s) or a folder {Hint: Press cancel to select a folder}")
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

        self.assembler_combo = QComboBox()
        self.assembler_combo.addItem("CAP3 default (legacy paired pipeline)", userData="__cap3_default__")
        self.assembler_combo.addItem("All assemblers (compare + pick best contig)", userData="__all__")
        for assembler in list_assemblers():
            self.assembler_combo.addItem(assembler.display_name, userData=assembler.id)
        self.assembler_combo.setToolTip(
                "Choose legacy CAP3 flow, a single assembler (faster: runs only that assembler), or all assemblers with auto-selection for BLAST payloads (in spirit of metawrap compare assebmlest and select best one)."
        )

        self.compare_assemblers_btn = QPushButton("Compare assemblers")
        self.compare_assemblers_btn.setToolTip("Compare all registered assemblers on paired FASTA inputs and write asm/compare_assemblers.tsv")
        self.assembler_combo.setVisible(False)
        self.assembler_combo.setEnabled(False)

        self.primer_set_combo = QComboBox()
        for label in PRIMER_SETS: 
            self.primer_set_combo.addItem(label)

        # Create a single-line text input field for the forward read tokens 
        self.fwd_pattern_edit = QLineEdit()
        self.fwd_pattern_edit.setPlaceholderText("Forward primer labels (comma-separated, e.g. 27F,8F,515F)")

        # Create a single-line text input field for the reverse read tokens 
        self.rev_pattern_edit = QLineEdit()
        self.rev_pattern_edit.setPlaceholderText("Reverse primer labels (comma-separated, e.g. 1492R, 806R)")

        self.detect_tokens_btn = QPushButton("Auto-detect Primers")

        self.detect_tokens_btn.setToolTip("Scan selected files/folders for forward/reverse primers")

        self.preview_pairs_btn = QPushButton("Preview Pairs")
        self.preview_pairs_btn.setToolTip("Summarize detected forward/reverse reads w/o running pipeline in a dry run")

        self.primer_preview_btn = QPushButton("Primer preview")
        self.primer_preview_btn.setToolTip("Run primer detect-only scan and show per-file hit metrics.")

        self.enforce_well_chk = QCheckBox("Enforce same plate well (A1-H12)")

        self.enforce_well_chk.setToolTip("Only pair forward/reverse reads when the same well plate code (e.g. B10)")

        self.cap3_profile_combo = QComboBox()
        self.cap3_profile_combo.addItem("Strict", userData="strict")
        self.cap3_profile_combo.addItem("Diagnostic", userData="diagnostic")
        self.cap3_profile_combo.addItem("Relaxed", userData="relaxed")
        self.cap3_profile_combo.setToolTip("Preset CAP3 overlap thresholds")

        self.cap3_extra_args_edit = QLineEdit()
        self.cap3_extra_args_edit.setPlaceholderText("Extra CAP3 args (e.g., -a 20 -b 15)")

        self.cap3_qual_chk = QCheckBox("Use per-base quality scores for assembly (its required for correct CAP3 scoring)")
        self.cap3_qual_chk.setToolTip("Use QUAL files to weight CAP3 assembly scoring it ensures its correct.")

        self.write_blast_inputs_chk = QCheckBox("Create BLAST input file (contigs or singlets)")
        self.write_blast_inputs_chk.setToolTip("Emit asm/blast_inputs.fasta + asm/blast_inputs.tsv.")

        self.use_blast_inputs_combo = QComboBox()
        self.use_blast_inputs_combo.addItem(
            "Paired contigs only (paired_contigs.fasta)", userData=False
        )
        self.use_blast_inputs_combo.addItem(
            "Contig or singlets (blast_inputs.fasta)", userData=True
        )

        self.overlap_audit_chk = QCheckBox("Overlap audit (diagnostic only)")
        self.overlap_audit_chk.setToolTip("Run overlap heuristic and write qc/overlap_audit.tsv.")

        self.primer_trim_mode_combo = QComboBox()
        self.primer_trim_mode_combo.addItems(["Off", "Detect", "Clip"])
        self.primer_trim_mode_combo.setCurrentText(self.settings.value("primer_trim_mode", "Off"))
        self.primer_trim_mode_combo.currentTextChanged.connect(
            lambda txt: self.settings.setValue("primer_trim_mode", txt)
        )

        self.primer_trim_stage_combo = QComboBox()
        self.primer_trim_stage_combo.addItem("Post-quality", userData="post_quality")
        self.primer_trim_stage_combo.addItem("Pre-quality", userData="pre_quality")
        stg = self.settings.value("primer_trim_stage", "post_quality")
        self.primer_trim_stage_combo.setCurrentIndex(max(0, self.primer_trim_stage_combo.findData(stg)))
        self.primer_trim_stage_combo.currentIndexChanged.connect(
            lambda _i: self.settings.setValue("primer_trim_stage", self.primer_trim_stage_combo.currentData())
        )

        self.primer_trim_preset_combo = QComboBox()
        self.primer_trim_preset_combo.addItem("Custom from config", userData="")
        for preset_name in TRIM_PRESETS:
            self.primer_trim_preset_combo.addItem(preset_name, userData=preset_name)
        pst = self.settings.value("primer_trim_preset", "")
        self.primer_trim_preset_combo.setCurrentIndex(max(0, self.primer_trim_preset_combo.findData(pst)))
        self.primer_trim_preset_combo.currentIndexChanged.connect(
            lambda _i: self.settings.setValue("primer_trim_preset", self.primer_trim_preset_combo.currentData())
        )

        self.primer_fwd_edit = QPlainTextEdit()
        self.primer_fwd_edit.setPlaceholderText("Forward synthetic flank sequences (one per line)")
        self.primer_fwd_edit.setFixedHeight(54)
        self.primer_fwd_edit.setPlainText(self.settings.value("primer_trim_fwd", ""))
        self.primer_fwd_edit.textChanged.connect(
            lambda: self.settings.setValue("primer_trim_fwd", self.primer_fwd_edit.toPlainText())
        )

        self.primer_rev_edit = QPlainTextEdit()
        self.primer_rev_edit.setPlaceholderText("Reverse synthetic flank sequences (one per line)")
        self.primer_rev_edit.setFixedHeight(54)
        self.primer_rev_edit.setPlainText(self.settings.value("primer_trim_rev", ""))
        self.primer_rev_edit.textChanged.connect(
            lambda: self.settings.setValue("primer_trim_rev", self.primer_rev_edit.toPlainText())
        )

        self.primer_save_btn = QPushButton("Use custom flanks for this run")
        self.primer_save_btn.setToolTip("Use current forward/reverse synthetic flank sequences for this session.")
        self.primer_save_btn.clicked.connect(lambda: self.primer_trim_preset_combo.setCurrentIndex(0))

        self.collapse_reps_chk = QCheckBox("Collapse replicates")
        self.collapse_reps_chk.setToolTip("Collapse technical replicates with vsearch (requires vsearch).")

        self.orient_reads_chk = QCheckBox("Orient reads")
        self.orient_reads_chk.setToolTip(
            "Orient reads against the selected DB reference before collapse/chimera. "
            "Recommended for mixed-orientation Sanger (forward+reverse primer runs)."
        )

        self.chimera_mode_combo = QComboBox()
        self.chimera_mode_combo.addItem("Off", userData="off")
        self.chimera_mode_combo.addItem("Reference (vsearch)", userData="reference")
        self.chimera_mode_combo.setToolTip(
            "Reference-based chimera filtering uses vsearch and the selected DB's chimera_ref by default."
        )
        
        self.dup_policy_lbl = QLabel("Duplicate Policy")
        self.dup_policy_combo = QComboBox() 
        for policy in DupPolicy:
            self.dup_policy_combo.addItem(policy.value, userData=policy) 

        self.advanced_regex_chk = QCheckBox("Show regex override")
        self.fwd_regex_edit = QLineEdit() 
        self.fwd_regex_edit.setPlaceholderText("Forward regex override (for advanced users who know regex)")
        self.rev_regex_edit = QLineEdit()
        self.rev_regex_edit.setPlaceholderText("Reverse regex override (for advanced users who know regex)")
        for edit in (self.fwd_regex_edit, self.rev_regex_edit):
            edit.setVisible(False)
            edit.setEnabled(False)
        self.advanced_regex_chk.setVisible(False) # Only shown in paired mode 

        # ---- Restore previous choises from settings --------

        # Retrieve previous saved assembly mode from application settings, defaulting to "single" if not found 
        saved_mode = self.settings.value("assembly_mode", "single")

        # Find the infex of the saved mode in combo box data, ensuring valid index >= 0 is returned 
        idx = max(0, self.mode_combo.findData(saved_mode))

        # Set dropdown menu to display the saved/default choice 
        self.mode_combo.setCurrentIndex(idx) 

        # Immediately call function to update the UI based on intial mode selection (example: show/hide reverse pattern field) 
        self._on_mode_changed(idx) 

        saved_assembler = self.settings.value("assembler_id", self.assembler_combo.itemData(0))
        asm_idx = self.assembler_combo.findData(saved_assembler)
        self.assembler_combo.setCurrentIndex(0 if asm_idx < 0 else asm_idx)
        self.assembler_combo.currentIndexChanged.connect(
            lambda _i: self.settings.setValue("assembler_id", self.assembler_combo.currentData())
        )

        # ------ Connect UI events to functions/settings storage -------- 

        # Connect signal emitted when dropdown index changes to the function handles the change 
        self.mode_combo.currentIndexChanged.connect(self._on_mode_changed)
        
        # Save the setting by default 27F and 1482R default to these if not set 
        saved_set = self.settings.value("primer_set", "16S (27F/1492R)")
        
        # Find position of saved primer set in dropdown menu `primer_set_combo` 
        idx = max(0, self.primer_set_combo.findText(saved_set))
        # updates visual display to show set that is saved from drop down menu 
        self.primer_set_combo.setCurrentIndex(idx)

        # Populate the token fields from settings or primer results from what saved 
        saved_fwd_tokens = self.settings.value("fwd_tokens", "")
        saved_rev_tokens = self.settings.value("rev_tokens", "") 
        default_fwd, default_rev = PRIMER_SETS.get(saved_set, PRIMER_SETS["16S (27F/1492R)"])

        # retrieve any custom regexes the user may have used from the previous run forward/reverse 
        self.fwd_pattern_edit.setText(saved_fwd_tokens or ", ".join(default_fwd))
        self.rev_pattern_edit.setText(saved_rev_tokens or ", ".join(default_rev))

        saved_fwd_regex = self.settings.value("fwd_regex", "")
        saved_rev_regex = self.settings.value("rev_regex", "") 
        self.fwd_regex_edit.setText(saved_fwd_regex)
        self.rev_regex_edit.setText(saved_rev_regex)
        
        # Set intial state of enfore_same_well checkbox default to False 
        self.enforce_well_chk.setChecked(self.settings.value("enforce_same_well", False, type=bool))

        self.collapse_reps_chk.setChecked(
            self.settings.value("collapse_replicates", False, type=bool)
        )
        self.orient_reads_chk.setChecked(
            self.settings.value("orient_reads", False, type=bool)
        )

        saved_chimera_mode = self.settings.value("chimera_mode", "off")
        chimera_idx = max(0, self.chimera_mode_combo.findData(saved_chimera_mode))
        self.chimera_mode_combo.setCurrentIndex(chimera_idx)
        
        saved_dup_policy = DupPolicy(self.settings.value("dup_policy", DupPolicy.ERROR.value))
        dup_idx = max(0, self.dup_policy_combo.findData(saved_dup_policy))
        self.dup_policy_combo.setCurrentIndex(dup_idx)

        saved_profile = self.settings.value("cap3_profile", "strict")
        profile_idx = max(0, self.cap3_profile_combo.findData(saved_profile))
        self.cap3_profile_combo.setCurrentIndex(profile_idx)

        saved_cap3_extra = self.settings.value("cap3_extra_args", "")
        self.cap3_extra_args_edit.setText(saved_cap3_extra)

        self.cap3_qual_chk.setChecked(
            self.settings.value("cap3_use_qual", True, type=bool)
        )
        self.write_blast_inputs_chk.setChecked(
            self.settings.value("write_blast_inputs", True, type=bool)
        )
        saved_use_blast_inputs = self.settings.value("use_blast_inputs", False, type=bool)
        use_blast_idx = max(0, self.use_blast_inputs_combo.findData(saved_use_blast_inputs))
        self.use_blast_inputs_combo.setCurrentIndex(use_blast_idx)

        self.overlap_audit_chk.setChecked(
            self.settings.value("overlap_audit", False, type=bool)
        )

        # ------ connecting user action event handlers -----------------
        # For when users change selection in drop down menu 
        self.primer_set_combo.currentIndexChanged.connect(self._on_primer_set_changed)
        self.detect_tokens_btn.clicked.connect(self._detect_tokens)
        self.preview_pairs_btn.clicked.connect(self._preview_pairs)
        self.primer_preview_btn.clicked.connect(self._preview_primers)
        self.compare_assemblers_btn.clicked.connect(self._compare_assemblers)
        self.advanced_regex_chk.toggled.connect(self._toggle_advanced_regex)
        self.enforce_well_chk.toggled.connect(
            lambda checked: self.settings.setValue("enforce_same_well", checked)
        )

        self.collapse_reps_chk.toggled.connect(
            lambda checked: self.settings.setValue("collapse_replicates", checked)
        )
        self.orient_reads_chk.toggled.connect(
            lambda checked: self.settings.setValue("orient_reads", checked)
        )

        self.chimera_mode_combo.currentIndexChanged.connect(
            lambda idx: self.settings.setValue(
                "chimera_mode", self.chimera_mode_combo.itemData(idx)
            )
        )
        
        self.dup_policy_combo.currentIndexChanged.connect(
            lambda idx: self.settings.setValue(
                "dup_policy", 
                getattr(data := self.dup_policy_combo.itemData(idx), "value", data),
            ) 
        )
        
        self.cap3_profile_combo.currentIndexChanged.connect(
            lambda idx: self.settings.setValue(
                "cap3_profile", self.cap3_profile_combo.itemData(idx)
            )
        )
        self.cap3_extra_args_edit.textEdited.connect(
            lambda txt: self.settings.setValue("cap3_extra_args", txt)
        )
        self.cap3_qual_chk.toggled.connect(
            lambda checked: self.settings.setValue("cap3_use_qual", checked)
        )
        self.write_blast_inputs_chk.toggled.connect(
            lambda checked: self.settings.setValue("write_blast_inputs", checked)
        )
        self.use_blast_inputs_combo.currentIndexChanged.connect(
            lambda idx: self.settings.setValue(
                "use_blast_inputs", self.use_blast_inputs_combo.itemData(idx)
            )
        )
        self.overlap_audit_chk.toggled.connect(
            lambda checked: self.settings.setValue("overlap_audit", checked)
        )

        # Connect the signal for text edits in the forward field to an anonymous function that saves the new text to settings. 
        self.fwd_pattern_edit.textEdited.connect(
            lambda txt: self._on_token_edited("fwd_tokens", txt)
        )

        # Connect the signal for text edits in the reverse field to an anonymous function that saves the new text to settings 
        self.rev_pattern_edit.textEdited.connect(
            lambda txt: self._on_token_edited("rev_tokens", txt)
        )

        self.fwd_regex_edit.textEdited.connect(
            lambda txt: self.settings.setValue("fwd_regex", txt)
        )

        self.rev_regex_edit.textEdited.connect(
            lambda txt: self.settings.setValue("rev_regex", txt)
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

        # ---- Output tables + tabs ----
        self.summary_filter = QComboBox()
        self.summary_filter.addItem("All statuses", userData=None)

        self.summary_table = QTableWidget(0, 15)
        self.summary_table.setHorizontalHeaderLabels(
            [
                "sample_id",
                "status",
                "assembler",
                "contig_len",
                "blast_payload",
                "selected_engine",
                "configured_engine",
                "merge_status",
                "merge_overlap_len",
                "merge_identity",
                "overlaps_saved",
                "overlaps_removed",
                "primer_mode",
                "primer_stage",
                "trace_qc",
            ]
        )
        _set_header_tooltips(self.summary_table, {"status": STATUS_COLOR_LEGEND})
        self.summary_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.summary_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.summary_table.itemSelectionChanged.connect(self._update_detail_panel)

        self.blast_table = QTableWidget(0, 5)
        self.blast_table.setHorizontalHeaderLabels(
            ["sample_id", "blast_payload", "reason", "payload_ids", "trace_qc"]
        )
        self.blast_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.blast_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.blast_table.itemSelectionChanged.connect(self._update_detail_panel)

        self.diagnostics_table = QTableWidget(0, 40)
        self.diagnostics_table.setHorizontalHeaderLabels(
            [
                "sample_id",
                "overlap_len",
                "overlap_identity",
                "overlap_quality",
                "orientation",
                "status",
                "best_identity",
                "best_identity_orientation",
                "anchoring_feasible",
                "end_anchored_possible",
                "fwd_best_identity",
                "revcomp_best_identity",
                "fwd_best_overlap_len",
                "revcomp_best_overlap_len",
                "fwd_anchor_feasible",
                "revcomp_anchor_feasible",
                "identity_delta_revcomp_minus_fwd",
                "selected_vs_best_identity_delta",
                "top_candidate_count",
                "top2_identity_delta",
                "top2_overlap_len_delta",
                "top2_quality_delta",
                "tie_reason_code",
                "fwd_best_identity_any",
                "revcomp_best_identity_any",
                "fwd_best_overlap_len_any",
                "revcomp_best_overlap_len_any",
                "pretrim_best_identity",
                "posttrim_best_identity",
                "pretrim_best_overlap_len",
                "posttrim_selected_overlap_len",
                "pretrim_status",
                "posttrim_status",
                "ambiguity_identity_delta_used",
                "ambiguity_quality_epsilon_used",
                "primer_trim_bases_fwd",
                "primer_trim_bases_rev",
                "selected_engine",
                "fallback_used",
                "overlap_engine",
            ]
        )
        _set_header_tooltips(
            self.diagnostics_table,
            {
                "status": STATUS_COLOR_LEGEND,
                "trace_qc": TRACE_QC_COLOR_LEGEND,
                "trace_qc_status": TRACE_QC_COLOR_LEGEND,
            },
        )
        self.diagnostics_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.diagnostics_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)

        self.compare_table = QTableWidget(0, 12)
        self.compare_table.setHorizontalHeaderLabels(
            [
                "sample_id",
                "assembler_id",
                "assembler_name",
                "dup_policy",
                "status",
                "selected_engine",
                "contig_len",
                "diag_code_for_machine",
                "diag_detail_for_human",
                "cap3_contigs_n",
                "cap3_singlets_n",
                "warnings",
            ]
        )
        self.compare_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.compare_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.compare_table.itemSelectionChanged.connect(self._update_detail_panel)

        self.open_blast_inputs_btn = QPushButton("Open blast_inputs.fasta")
        self.open_blast_inputs_btn.clicked.connect(lambda: self._open_path(self._blast_inputs_fasta, prefer_in_app=True))
        self.open_blast_inputs_btn.setEnabled(False)

        self.open_blast_inputs_system_btn = QPushButton("Open blast_inputs.fasta (system)")
        self.open_blast_inputs_system_btn.clicked.connect(
            lambda: self._open_path_system(self._blast_inputs_fasta)
        )
        self.open_blast_inputs_system_btn.setEnabled(False)


        self.open_asm_folder_btn = QPushButton("Open asm folder")
        self.open_asm_folder_btn.clicked.connect(lambda: self._open_path(self._asm_dir))
        self.open_asm_folder_btn.setEnabled(False)

        summary_tab = QWidget()
        summary_layout = QVBoxLayout(summary_tab)
        summary_filter_row = QHBoxLayout()
        summary_filter_row.addWidget(QLabel("Status filter"))
        summary_filter_row.addWidget(self.summary_filter)
        summary_filter_row.addStretch()
        summary_layout.addLayout(summary_filter_row)
        summary_layout.addWidget(self.summary_table)

        blast_tab = QWidget()
        blast_layout = QVBoxLayout(blast_tab)
        blast_buttons = QHBoxLayout()
        blast_buttons.addWidget(self.open_blast_inputs_btn)
        blast_buttons.addWidget(self.open_blast_inputs_system_btn) 
        blast_buttons.addWidget(self.open_asm_folder_btn)
        blast_buttons.addStretch()
        blast_layout.addLayout(blast_buttons)
        blast_layout.addWidget(self.blast_table)

        diagnostics_tab = QWidget()
        diagnostics_layout = QVBoxLayout(diagnostics_tab)
        diagnostics_layout.addWidget(self.diagnostics_table)

        compare_tab = QWidget()
        compare_layout = QVBoxLayout(compare_tab)
        compare_layout.addWidget(self.compare_table)

        logs_tab = QWidget()
        logs_layout = QVBoxLayout(logs_tab)
        logs_layout.addWidget(self.log_view)

        self.output_tabs = QTabWidget()
        self.output_tabs.addTab(logs_tab, "Logs")
        self.output_tabs.addTab(summary_tab, "Assembly Summary")
        self.output_tabs.addTab(blast_tab, "BLAST Inputs")
        self.output_tabs.addTab(diagnostics_tab, "Diagnostics")
        self.output_tabs.addTab(compare_tab, "Compare Assemblers")

        # ---- Detail panel ----
        self.detail_sample_lbl = QLabel("sample_id: ")
        self.detail_status_lbl = QLabel("status: ")
        self.detail_reason_lbl = QLabel("reason: ")
        self.detail_payload_lbl = QLabel("blast_payload: ")
        self.detail_payload_ids_lbl = QLabel("payload_ids: ")
        self.detail_overlap_lbl = QLabel("overlap: ")
        self.detail_assembly_engine_lbl = QLabel("assembly path: ")
        self.detail_trace_qc_lbl = QLabel("trace_qc: ")
        self.detail_trace_reason_lbl = QLabel("trace_qc_detail: ")
        self.detail_trace_qc_lbl.setToolTip(
            "Trace QC flags are chromatogram-signal indicators; contamination is not directly inferred from trace alone."
        )
        self.detail_trace_reason_lbl.setToolTip(
            "Trace QC flags are chromatogram-signal indicators; contamination is not directly inferred from trace alone."
        )
        for detail_lbl in (
            self.detail_reason_lbl,
            self.detail_payload_ids_lbl,
            self.detail_overlap_lbl,
            self.detail_trace_reason_lbl,
        ):
            detail_lbl.setWordWrap(True)
            detail_lbl.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
            detail_lbl.setMinimumWidth(0)

        self.detail_contigs_btn = QPushButton("Open contigs")
        self.detail_singlets_btn = QPushButton("Open singlets")
        self.detail_info_btn = QPushButton("Open cap.info")
        self.detail_contigs_system_btn = QPushButton("Open contigs (system)")
        self.detail_singlets_system_btn = QPushButton("Open singlets (system)")
        self.detail_info_system_btn = QPushButton("Open cap.info (system)")
        for btn in (
            self.detail_contigs_btn,
            self.detail_singlets_btn,
            self.detail_info_btn,
            self.detail_contigs_system_btn,
            self.detail_singlets_system_btn,
            self.detail_info_system_btn,
        ): 
            btn.setEnabled(False)
        self.detail_contigs_btn.clicked.connect(
            lambda: self._open_detail_paths("contigs", prefer_in_app=True)
        )
        self.detail_singlets_btn.clicked.connect(
            lambda: self._open_detail_paths("singlets", prefer_in_app=True)
        )
        self.detail_info_btn.clicked.connect(
            lambda: self._open_detail_paths("info", prefer_in_app=True)
        )
        self.detail_contigs_system_btn.clicked.connect(
            lambda: self._open_detail_paths("contigs", system=True)
        )
        self.detail_singlets_system_btn.clicked.connect(
            lambda: self._open_detail_paths("singlets", system=True)
        )
        self.detail_info_system_btn.clicked.connect(
            lambda: self._open_detail_paths("info", system=True) 
        )

        detail_box = QGroupBox("Details")
        detail_layout = QVBoxLayout(detail_box)
        detail_layout.addWidget(self.detail_sample_lbl)
        detail_layout.addWidget(self.detail_status_lbl)
        detail_layout.addWidget(self.detail_reason_lbl)
        detail_layout.addWidget(self.detail_payload_lbl)
        detail_layout.addWidget(self.detail_payload_ids_lbl)
        detail_layout.addWidget(self.detail_overlap_lbl)
        detail_layout.addWidget(self.detail_assembly_engine_lbl)
        detail_layout.addWidget(self.detail_trace_qc_lbl)
        detail_layout.addWidget(self.detail_trace_reason_lbl)
        detail_layout.addWidget(self.detail_contigs_btn)
        detail_layout.addWidget(self.detail_singlets_btn)
        detail_layout.addWidget(self.detail_info_btn)
        detail_layout.addWidget(self.detail_contigs_system_btn)
        detail_layout.addWidget(self.detail_singlets_system_btn)
        detail_layout.addWidget(self.detail_info_system_btn)
        detail_layout.addStretch()

        self.output_splitter = QSplitter()
        self.output_splitter.addWidget(self.output_tabs)
        self.output_splitter.addWidget(detail_box)
        self.output_splitter.setStretchFactor(0, 3)
        self.output_splitter.setStretchFactor(1, 1)

        self.summary_filter.currentIndexChanged.connect(self._apply_summary_filter)
        self.diagnostics_table.itemSelectionChanged.connect(self._update_detail_panel)
        self.summary_table.cellDoubleClicked.connect(
            lambda row, col: self._show_cell_text(self.summary_table, row, col)
        )
        self.blast_table.cellDoubleClicked.connect(
            lambda row, col: self._show_cell_text(self.blast_table, row, col)
        )
        self.diagnostics_table.cellDoubleClicked.connect(
            lambda row, col: self._show_cell_text(self.diagnostics_table, row, col)
        )
        self.compare_table.cellDoubleClicked.connect(
            lambda row, col: self._show_cell_text(self.compare_table, row, col)
        )


        self._blast_inputs_fasta: Optional[Path] = None
        self._asm_dir: Optional[Path] = None
        self._detail_paths: dict[str, list[Path]] = {"contigs": [], "singlets": [], "info": []} 
        self._summary_rows: dict[str, dict[str, str]] = {}
        self._blast_rows: dict[str, dict[str, str]] = {}
        self._audit_rows: dict[str, dict[str, str]] = {}
        self._compare_rows: dict[tuple[str, str], dict[str, str]] = {} 
        self._active_viewer: Optional[QDialog] = None
        self._active_pairing_preview: Optional[QDialog] = None

        # connect to model's signal to auto-scroll to the bottom on new logs if the user is already at the bottom or near it 
        sb = self.log_view.verticalScrollBar() 
        sb.valueChanged.connect(lambda v: setattr(self, "_autoscroll", (sb.maximum() - v) <= 2))

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

        pipeline_opts = QHBoxLayout()
        pipeline_opts.addWidget(QLabel("Pipeline options"))
        pipeline_opts.addWidget(self.collapse_reps_chk)
        pipeline_opts.addWidget(self.orient_reads_chk)
        pipeline_opts.addWidget(QLabel("Chimera"))
        pipeline_opts.addWidget(self.chimera_mode_combo)
        pipeline_opts.addStretch()

        cap3_opts = QHBoxLayout()
        cap3_opts.addWidget(QLabel("CAP3 profile"))
        cap3_opts.addWidget(self.cap3_profile_combo)
        cap3_opts.addWidget(self.cap3_extra_args_edit)
        cap3_opts.addWidget(self.cap3_qual_chk)
        cap3_opts.addWidget(self.write_blast_inputs_chk)
        cap3_opts.addWidget(self.use_blast_inputs_combo)
        cap3_opts.addWidget(self.overlap_audit_chk)
        cap3_opts.addWidget(QLabel("Primer trim"))
        cap3_opts.addWidget(self.primer_trim_mode_combo)
        cap3_opts.addWidget(QLabel("Primer stage"))
        cap3_opts.addWidget(self.primer_trim_stage_combo)
        cap3_opts.addWidget(QLabel("Known synthetic flank preset"))
        cap3_opts.addWidget(self.primer_trim_preset_combo)
        cap3_opts.addWidget(self.primer_save_btn)
        cap3_opts.addStretch()

        primer_seq_opts = QHBoxLayout()
        primer_seq_opts.addWidget(QLabel("Trim primers (F)"))
        primer_seq_opts.addWidget(self.primer_fwd_edit)
        primer_seq_opts.addWidget(QLabel("Trim primers (R)"))
        primer_seq_opts.addWidget(self.primer_rev_edit)

        pairing = QHBoxLayout()
        pairing.addWidget(QLabel("Primer set"))
        pairing.addWidget(self.primer_set_combo)
        pairing.addWidget(self.fwd_pattern_edit)
        pairing.addWidget(self.rev_pattern_edit)
        pairing.addWidget(self.detect_tokens_btn)
        pairing.addWidget(self.preview_pairs_btn)
        pairing.addWidget(self.primer_preview_btn)
        pairing.addWidget(self.compare_assemblers_btn)
        pairing.addWidget(QLabel("Assembler selection"))
        pairing.addWidget(self.assembler_combo)
        pairing.addWidget(self.dup_policy_lbl)
        pairing.addWidget(self.dup_policy_combo) 
        pairing.addWidget(self.enforce_well_chk)
        pairing.addWidget(self.advanced_regex_chk)
        pairing.addWidget(self.fwd_regex_edit)
        pairing.addWidget(self.rev_regex_edit) 

        # vertical stack picker row, settings row, then logpane expands
        outer = QVBoxLayout()
        outer.addLayout(top)
        outer.addLayout(mid)
        outer.addLayout(pipeline_opts)
        outer.addLayout(cap3_opts)
        outer.addLayout(primer_seq_opts)
        outer.addLayout(pairing) 
        outer.addWidget(self.output_splitter)

        # Embed composite layout as the window's central widget
        root = QWidget(); root.setLayout(outer)
        self.setCentralWidget(root)

        # state and logging connection
        self._infiles: list[Path] = []
        self._input_path: Optional[Path] = None 
        self._input_staging: Optional[tempfile.TemporaryDirectory[str]] = None 
        self.hits_path: Optional[Path] = None
        self.meta_path: Optional[Path] = None
        self._current_stage: str = ""

        # placeholders to prevent AttributeErrors before first run
        self._thread: Optional[QThread] = None
        self._worker: Optional[Worker] = None
        self._log_handler: Optional[LogBridge] = None
        self._pending_logs: list[str] = [] 
        self._log_flush_timer = QTimer(self)
        self._log_flush_timer.setInterval(33) # ~33Hz 
        self._log_flush_timer.setSingleShot(True) 
        self._log_flush_timer.timeout.connect(self._flush_logs)
        self._autoscroll = True
        self._active_box: Optional[QMessageBox] = None 
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
        for widget in (
            self.fwd_pattern_edit,
            self.rev_pattern_edit,
            self.primer_set_combo,
            self.detect_tokens_btn,
            self.preview_pairs_btn,
            self.primer_preview_btn,
            self.compare_assemblers_btn,
            self.assembler_combo,
            self.dup_policy_lbl,
            self.dup_policy_combo,
            self.advanced_regex_chk,
            self.fwd_regex_edit, 
            self.rev_regex_edit,
            self.enforce_well_chk,
            self.cap3_profile_combo,
            self.cap3_extra_args_edit,
            self.cap3_qual_chk,
            self.write_blast_inputs_chk,
            self.use_blast_inputs_combo,
            self.overlap_audit_chk,
            self.primer_trim_mode_combo,
            self.primer_trim_stage_combo,
            self.primer_trim_preset_combo,
            self.primer_fwd_edit,
            self.primer_rev_edit,
            self.primer_save_btn,
        ):
            widget.setVisible(is_paired)
            widget.setEnabled(is_paired)
        self._sync_primer_controls_visibility() 
        self.settings.setValue("assembly_mode", mode)

    def _assembly_kwargs(self) -> dict:
        """
        This function is meant to be used to gather all assembly configuration data into a dictionary.
        """ 
        mode = self.mode_combo.currentData()
        if mode == "paired":
            fwd, rev = self._current_patterns()
        else:
            fwd = rev = None
        raw_extra_args = self.cap3_extra_args_edit.text().strip()
        extra_args = raw_extra_args.split() if raw_extra_args else None
        write_blast_inputs = self.write_blast_inputs_chk.isChecked()
        use_blast_inputs = self.use_blast_inputs_combo.currentData()
        return {
            "mode": mode, 
            "assembler_id": self.assembler_combo.currentData(),
            "fwd_pattern": fwd,
            "rev_pattern": rev, 
            "enforce_same_well": self.enforce_well_chk.isChecked(),
            "dup_policy": self.dup_policy_combo.currentData(),
            "cap3_profile": self.cap3_profile_combo.currentData(),
            "cap3_extra_args": extra_args,
            "cap3_use_qual": self.cap3_qual_chk.isChecked(),
            "write_blast_inputs": write_blast_inputs,
            "use_blast_inputs": use_blast_inputs,
        }

    def _pipeline_kwargs(self) -> dict:
        return {
            "collapse_replicates": self.collapse_reps_chk.isChecked(),
            "orient_reads": self.orient_reads_chk.isChecked(),
            "chimera_mode": self.chimera_mode_combo.currentData(),
            "overlap_audit": self.overlap_audit_chk.isChecked(),
        }


    def _primer_override_kwargs(self, *, for_preview: bool = False) -> dict:
        mode_label = self.primer_trim_mode_combo.currentText().strip().lower()
        mode_map = {"off": "off", "detect": "detect", "clip": "clip"}
        mode = mode_map.get(mode_label, "off")
        stage = self.primer_trim_stage_combo.currentData()
        preset = self.primer_trim_preset_combo.currentData()
        forward_raw = self.primer_fwd_edit.toPlainText()
        reverse_raw = self.primer_rev_edit.toPlainText()
        override = build_primer_cfg_override(
            mode=mode,
            stage=stage,
            preset=preset,
            forward_raw=forward_raw,
            reverse_raw=reverse_raw,
            for_preview=for_preview,
        )
        return {"primer_cfg_override": override}
    
    def _tokens_to_regex(self, text: str) -> str | None:
        tokens = [t.strip() for t in text.split(",") if t.strip()]
        if not tokens:
            return None 
        escaped = [re.escape(tok) for tok in tokens] 
        return "|".join(escaped)

    def _current_patterns(self) -> tuple[str | None, str | None]:
        if self.advanced_regex_chk.isChecked():
            fwd = self.fwd_regex_edit.text().strip() or None 
            rev = self.rev_regex_edit.text().strip() or None 
        else:
            fwd = self._tokens_to_regex(self.fwd_pattern_edit.text())
            rev = self._tokens_to_regex(self.rev_pattern_edit.text())
        return fwd, rev

    def _current_token_list(self) -> list[str]:
        if self.advanced_regex_chk.isChecked():
            return []
        tokens = [
            token.strip()
            for raw in (self.fwd_pattern_edit.text(), self.rev_pattern_edit.text())
            for token in raw.split(",")
        ]
        cleaned = [tok for tok in tokens if tok]
        if cleaned:
            return cleaned
        label = self.primer_set_combo.currentText()
        fwd_tokens, rev_tokens = PRIMER_SETS.get(label, ([], []))
        return [*fwd_tokens, *rev_tokens]


    def _ensure_vsearch_available(self) -> bool:
        if not (
            self.collapse_reps_chk.isChecked()
            or self.orient_reads_chk.isChecked()
            or self.chimera_mode_combo.currentData() != "off"
        ):
            return True
        try:
            resolve_vsearch()
        except FileNotFoundError as exc:
            self._show_box(QMessageBox.Icon.Warning, "vsearch not found", str(exc))
            return False
        return True

    def _on_primer_set_changed(self, index: int):
        label = self.primer_set_combo.itemText(index)
        self.settings.setValue("primer_set", label)
        fwd_tokens, rev_tokens = PRIMER_SETS.get(label, ([], []))
        if fwd_tokens:
            self.fwd_pattern_edit.setText(", ".join(fwd_tokens))
        if rev_tokens:
            self.rev_pattern_edit.setText(", ".join(rev_tokens))

        self._sync_primer_controls_visibility() 

    def _sync_primer_controls_visibility(self):
        is_paired = self.mode_combo.currentData() == "paired"
        self.advanced_regex_chk.setVisible(is_paired)
        if not is_paired or not self.advanced_regex_chk.isChecked():
            for edit in (self.fwd_regex_edit, self.rev_regex_edit):
                edit.setVisible(False)
                edit.setEnabled(False) 

    def _toggle_advanced_regex(self, checked: bool):
        for edit in (self.fwd_regex_edit, self.rev_regex_edit):
            edit.setVisible(checked)
            edit.setEnabled(checked)

    def _on_token_edited(self, key: str, text: str):
        self.settings.setValue(key, text)
        if text:
            custom_idx = self.primer_set_combo.findText("Custom")
            if custom_idx >= 0:
                self.primer_set_combo.blockSignals(True)
                self.primer_set_combo.setCurrentIndex(custom_idx)
                self.settings.setValue("primer_set", "Custom")
                self.primer_set_combo.blockSignals(False) 

    def _detect_tokens(self):
        if not self._input_path:
            self._show_box(QMessageBox.Icon.Warning,"No input", "Choose a file or folder first please.")
            return 
        directory = self._input_path if self._input_path.is_dir() else self._input_path.parent

        fwd_tokens, rev_tokens = self._scan_tokens(directory)

        if not fwd_tokens and not rev_tokens:
            self._show_box(
                QMessageBox.Icon.Information,
                "No primer labels found",
                "No primer-like forward/reverse labels were detected in the filenames."
            )
            return 
        
        if fwd_tokens:
            joined = ", ".join(fwd_tokens)
            self.fwd_pattern_edit.setText(joined)
            self._on_token_edited("fwd_tokens", joined)
        if rev_tokens:
            joined = ", ".join(rev_tokens)
            self.rev_pattern_edit.setText(joined)
            self._on_token_edited("rev_tokens", joined)

        self._show_box(
            QMessageBox.Icon.Information,
            "Primers/primer labels detected",
            f"Forward primers/primer labels: {', '.join(fwd_tokens) or '-'}\n"
            f"Reverse primers/primer labels: {', '.join(rev_tokens) or '-'}\n"
        ) 

    def _scan_tokens(self, directory: Path) -> tuple[list[str], list[str]]:
        token_rx = re.compile(r"([A-Za-z0-9]+[FR])", re.I)
        fwd_counter: collections.Counter[str] = collections.Counter()
        rev_counter: collections.Counter[str] = collections.Counter() 

        for suffix in ("*.fasta", "*.fastq", "*.ab1", "*seq"):
            for fp in sorted(directory.rglob(suffix)):
                for tok in token_rx.findall(fp.name):
                    tok = tok.upper()
                    if tok.endswith("F"):
                        fwd_counter[tok] += 1 
                    elif tok.endswith("R"):
                        rev_counter[tok] += 1 

        def _top(counter: collections.Counter[str]) -> list[str]:
            return [tok for tok, _ in counter.most_common(5)] 

        return _top(fwd_counter), _top(rev_counter) 

    def _preview_pairs(self):
        if self.mode_combo.currentData() != "paired":
            self._show_box(QMessageBox.Icon.Information, "Not paired", "Preview is only available in paired mode.")
            return 

        if not self._input_path:
            self._show_box(QMessageBox.Icon.Warning, "No input", "Choose a file or folder first for set of inputs.")
            return 

        directory = self._input_path if self._input_path.is_dir() else self._input_path.parent 
        fwd, rev = self._current_patterns()
        def _refresh_preview():
            summary_text = _summarize_paired_candidates(
                directory,
                fwd,
                rev,
                enforce_same_well=self.enforce_well_chk.isChecked()
            )
            suggestions_text = _suggest_pairing_patterns(directory)
            candidates_list = analyze_pairing_candidates(
                directory,
                fwd,
                rev,
                known_tokens=self._current_token_list(),
                enforce_same_well=self.enforce_well_chk.isChecked(),
            )
            return summary_text, suggestions_text, candidates_list

        summary, suggestions, candidates = _refresh_preview() 
        dialog = PairingPreviewDialog(
            summary=summary,
            suggestions=suggestions,
            candidates=candidates,
            enforce_same_well=self.enforce_well_chk.isChecked(),
            refresh_callback=_refresh_preview,
            open_callback=self._open_external_detached,
            parent=self,
        )
        # Avoid nested event loops from exec();  keep this as show() to avoid the GUI being unresponsive from minimizing the window from preview_pairs 
        # aligns with the rest of the app's no-nested-event-loop stance that I will have moving foward given past issues.
        dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        dialog.setModal(False)
        dialog.setWindowModality(Qt.WindowModality.NonModal)
        self._active_pairing_preview = dialog
        dialog.destroyed.connect(lambda *_: setattr(self, "_active_pairing_preview", None))
        dialog.show()

    def _run_primer_preview_job(self, input_path: Path, primer_override: dict[str, object]) -> dict[str, object]:
        with tempfile.TemporaryDirectory(prefix="microseq_primer_preview_") as td:
            out_dir = Path(td)
            summary_tsv = out_dir / "qc" / "trim_summary.tsv"
            run_trim(
                input_path,
                out_dir,
                sanger=(input_path.suffix.lower() == ".ab1" or (input_path.is_dir() and any(input_path.glob("*.ab1")))),
                summary_tsv=summary_tsv,
                primer_cfg_override=primer_override,
            )
            report = out_dir / "qc" / "primer_detect_report.tsv"
            if not report.exists():
                return {
                    "primer_preview_notice": "No primer detect report was generated. Check primer_trim preset/sequences in config.",
                }
            return {
                "primer_preview_text": report.read_text(encoding="utf-8"),
            }

    def _preview_primers(self):
        if not self._input_path:
            self._show_box(QMessageBox.Icon.Warning, "No input", "Choose a file or folder first.")
            return

        try:
            primer_override = dict(self._primer_override_kwargs(for_preview=True).get("primer_cfg_override", {}))
        except ValueError as exc:
            self._show_box(QMessageBox.Icon.Warning, "Invalid primer sequence", str(exc))
            return
        primer_override["mode"] = "detect"
        self._launch(self._run_primer_preview_job, Path(self._input_path), primer_override)


    def _cleanup_input_staging(self) -> None:
        if self._input_staging:
            self._input_staging.cleanup()
            self._input_staging = None 
        self._input_path = None 

    def _stage_selected_files(self) -> None:
        if not self._infiles:
            self._cleanup_input_staging()
            return 

        if len(self._infiles) == 1:
            self._cleanup_input_staging()
            self._input_path = self._infiles[0]
            return 

        self._cleanup_input_staging()
        self._input_staging = tempfile.TemporaryDirectory(prefix="microseq_gui_inputs_")
        staging_dir = Path(self._input_staging.name)
        for idx, fp in enumerate(self._infiles, 1):
            dest = staging_dir / fp.name 
            if dest.exists():
                dest = staging_dir / f"{fp.stem}_{idx}{fp.suffix}"
            shutil.copy(fp, dest)
        self._input_path = staging_dir 

    def _choose_infile(self):
        """Select FASTA/FASTQ/AB1 file(s) or a folder of traces."""
        self._infiles = [] 
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
                self._infiles = [Path(dir_path)]
        else:
            # multi-file selection is allowed but only the first file is used
            self._infiles = [Path(p) for p in paths]
        
        self._stage_selected_files() 

        if self._input_path:
            label = (
                f"Input files: {len(self._infiles)} selected"
                if self._input_path.is_dir() and len(self._infiles) > 1 
                else f"Input folder: {self._input_path.name}"
                if self._input_path.is_dir()
                else f"Input file: {self._input_path.name}"
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
        if not self._input_path:
            self._show_box(QMessageBox.Icon.Warning, "No input", "Choose a file or folder first.")
            return

        self.progress.setValue(0) # resets progress bar during each run
        task = self.settings.value("blast_task", "megablast")

        # derive output file beside input; disables button; logs starts
        hits_path = self._input_path.with_suffix(".hits.tsv")

        self.run_btn.setEnabled(False)
        self.postblast_btn.setEnabled(False)
        self.cancel_btn.setEnabled(True)
        self.log_model.append(f"\n▶ BLAST {self._input_path.name} -> {hits_path.name}")

        # worker and thread wiring -------------------
        worker = Worker(
            run_blast_stage,
            self._input_path,
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

    def _compare_assemblers(self):
        if not self._input_path:
            self._show_box(QMessageBox.Icon.Warning, "No input", "Choose a file or folder first.")
            return
        if self.mode_combo.currentData() != "paired":
            self._show_box(QMessageBox.Icon.Information, "Paired mode required", "Assembler comparison requires paired mode.")
            return

        # Explicit contract: comparison runs on staged paired FASTA inputs.
        source_dir = self._input_path if self._input_path.is_dir() else self._input_path.parent
        has_fasta = any(source_dir.rglob("*.fasta")) or any(source_dir.rglob("*.fa")) or any(source_dir.rglob("*.fna")) or any(source_dir.rglob("*.fas"))
        has_raw = any(source_dir.rglob("*.ab1")) or any(source_dir.rglob("*.fastq")) or any(source_dir.rglob("*.fq"))
        if has_raw and not has_fasta:
            self._show_box(
                QMessageBox.Icon.Information,
                "Staged FASTA required",
                "Compare assemblers expects pre-staged paired FASTA inputs (e.g., qc/paired_fasta).",
            )
            return

        self.log_model.append("\n▶ Compare all assemblers")

        self._launch(
            run_compare_assemblers,
            source_dir,
            self.out_dir if self.out_dir else source_dir.parent,
            fwd_pattern=self._current_patterns()[0],
            rev_pattern=self._current_patterns()[1],
            enforce_same_well=self.enforce_well_chk.isChecked(),
        )


    # -------- Trim -> Convert -> BLAST -> Taxonomy ---------------
    def _launch_qc(self):
        if not self._input_path:
            self._show_box(QMessageBox.Icon.Warning, "No input", "Choose a file or folder first.")
            return
        if not self._ensure_vsearch_available():
            return
        primer_kw: dict = {}
        if self.mode_combo.currentData() == "paired":
            try:
                primer_kw = self._primer_override_kwargs()
            except ValueError as exc:
                self._show_box(QMessageBox.Icon.Warning, "Invalid primer sequence", str(exc))
                return
        task = self.settings.value("blast_task", "megablast")
        self._launch(
            run_full_pipeline,
            self._input_path,
            self.db_box.currentText(),
            threads=self.threads_spin.value(),
            postblast=self.biom_chk.isChecked(),
            metadata=None,   # Trim → Convert → BLAST → Tax
            blast_task=task,
            **self._assembly_kwargs(),
            **self._pipeline_kwargs(),
            **primer_kw,
        )

    # ------ Run full pipeline with Post-Blast as well -------------
    def _launch_full(self):
        if not self._input_path:
            self._show_box(QMessageBox.Icon.Warning, "No input", "Choose a file or folder first.")
            return
        if not self._ensure_vsearch_available():
            return
        primer_kw: dict = {}
        if self.mode_combo.currentData() == "paired":
            try:
                primer_kw = self._primer_override_kwargs()
            except ValueError as exc:
                self._show_box(QMessageBox.Icon.Warning, "Invalid primer sequence", str(exc))
                return
        # if user asked for BIOM but hasn't chosen metadata, abort early
        if self.biom_chk.isChecked() and not self.meta_path:
            self._show_box(QMessageBox.Icon.Warning,
                "No metadata",
                "A metadata file is required to build the BIOM table.\n"
                "Click ‘Browse metadata…’ first or un-tick ‘Make BIOM’."
            )
            return
        # launch the pipeline here
        task = self.settings.value("blast_task", "megablast")
        self._launch(
            run_full_pipeline,
            self._input_path,
            self.db_box.currentText(),
            threads=self.threads_spin.value(),
            postblast=self.biom_chk.isChecked(), # source of truth decide via checkbox
            metadata=self.meta_path,      # None or Path run the Post-BLAST stage too
            blast_task=task,
            **self._assembly_kwargs(),
            **self._pipeline_kwargs(),
            **primer_kw,
        )

    # ------ Run stand-alone Post-BLAST ---------------------------------
    def _launch_postblast(self):
        if not self.hits_path:
            self._show_box(QMessageBox.Icon.Warning, "No hits", "Choose a BLAST hits file first.")
            return
        if not self.meta_path:
            self._show_box(QMessageBox.Icon.Warning, "No metadata", "Choose a metadata file first.")
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
    def _queue_log(self, line: str):
        """Appends the log line directly to the model on the GUI thread.
        This is called from the Python Logger to update the GUI model."""
        self._pending_logs.append(line)
        if not self._log_flush_timer.isActive():
            self._log_flush_timer.start() 

    @Slot()
    def _flush_logs(self) -> None:
        if not self._pending_logs:
            return 
        lines = self._pending_logs
        self._pending_logs = [] 
        self.log_model.append_many(lines)

        if self._autoscroll:
            last = self.log_model.rowCount() - 1 
            if last >= 0:
                self.log_view.scrollTo(
                    self.log_model.index(last, 0),
                    QAbstractItemView.ScrollHint.PositionAtBottom
                ) 
    def _show_box(self, icon: QMessageBox.Icon, title: str, text: str) -> None:
        box = QMessageBox(self)
        box.setIcon(icon)
        box.setWindowTitle(title)
        box.setText(text)
        box.setStandardButtons(QMessageBox.Ok)
        box.setModal(True)
        box.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self._active_box = box
        box.finished.connect(lambda *_: setattr(self, "_active_box", None))
        box.open()

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
        self._flush_logs() 
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
        assert QThread.currentThread() == QApplication.instance().thread() 
        # interpret the saved result --------------------------------
        result = getattr(self, "_last_result", 1)      # default = error
        out: Optional[Path] = None 

        # --- decide status/message -------------------------------- 
        if isinstance(result, dict):                   # full‑pipeline success
            rc = 0 
            for key in ("tax", "hits_tax", "biom", "fasta", "assembly_summary"): 
                candidate = result.get(key)
                if candidate:
                    path = Path(candidate)
                    if path.exists():
                        out = path 
                        break 
        elif isinstance(result, RuntimeError) and str(result) == "Cancelled":
            rc = None 
        elif isinstance(result, Exception):            # Worker caught an error
            rc = 1 
            err = str(result) 

            # friendlier message for the "no hits" situation
            if "no blast hits" in err.lower(): 
                dialog_fn = lambda: self._show_box(
                    QMessageBox.Icon.Information, 
                    "Nothing to summarise",
                    ("BLAST finished, but no hits met the filters.\n"
                     "You can:\n"
                     " lower Identity / Q-cov, or\n"
                     "run again without the BIOM option.")
                )
            else:
                dialog_fn = lambda e=err: self._show_box(
                    QMessageBox.Icon.Warning,
                    "Run Failed", e)
            QTimer.singleShot(0, dialog_fn)
        elif isinstance(result, Path):
            rc, out = 0, result
        else:                                          # int from BLAST‑only path
            rc, out = int(result), Path()

        # log + optional dialog -------------------------------------
        msg = "Cancelled" if rc is None else ("Success" if rc == 0 else f"Failed (exit {rc})") 
            
        QTimer.singleShot(
            0, lambda m=msg: self.log_model.append(f"● {m}\n")
        )

        if rc == 0 and out is not None and not (isinstance(result, dict) and ("primer_preview_text" in result or "primer_preview_notice" in result)):
            QTimer.singleShot(
                0,
                lambda p=out: self._show_box(
                QMessageBox.Icon.Information,
                "Pipeline finished",
                f"Last output file:\n{p}") 
            )

        if rc == 0 and isinstance(result, dict) and "primer_preview_text" not in result and "primer_preview_notice" not in result:
            QTimer.singleShot(0, lambda r=result: self._load_output_tables(r))
        elif rc == 0 and isinstance(result, Path) and result.name == "compare_assemblers.tsv":
            QTimer.singleShot(0, lambda p=result: self._load_output_tables({"compare_assemblers": p}))

        if rc == 0 and isinstance(result, dict) and "primer_preview_text" in result:
            QTimer.singleShot(0, lambda t=result["primer_preview_text"]: TextBlobDialog("Primer preview (detect-only)", t, parent=self).exec())
        if rc == 0 and isinstance(result, dict) and "primer_preview_notice" in result:
            QTimer.singleShot(0, lambda m=result["primer_preview_notice"]: self._show_box(QMessageBox.Icon.Information, "Primer preview", m))

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

    def _open_path(self, path: Optional[Path], prefer_in_app: bool = False) -> None:
        if not path:
            self._show_box(QMessageBox.Icon.Information, "No file", "No file is available for this selection.")
            return
        if not path.exists():
            self._show_box(QMessageBox.Icon.Warning, "Missing file", f"File not found:\n{path}")
            return

        if prefer_in_app and path.is_file():
            self._open_text_viewer(path)
            return

        if self._open_external_detached(path):
            return

        self._open_path_fallback(path, reason="Unable to open via system handler.")

    def _open_path_system(self, path: Optional[Path]) -> None:
        if not path:
            self._show_box(QMessageBox.Icon.Information, "No file", "No file is available for this selection.")
            return
        if not path.exists():
            self._show_box(QMessageBox.Icon.Warning, "Missing file", f"File not found:\n{path}")
            return 
        if self._open_external_detached(path):
            return

        self._open_path_fallback(path, reason="Unable to open via system handler.")

    def _open_external_detached(self, path: Path) -> bool:
        if is_wsl():
            if shutil.which("wslview"):
                ok, _pid = QProcess.startDetached("wslview", [str(path)])
                return ok
            if shutil.which("explorer.exe") and shutil.which("wslpath"):
                try:
                    win_path = subprocess.check_output(
                        ["wslpath", "-w", str(path)],
                        text=True,
                    ).strip()
                except (OSError, subprocess.SubprocessError):
                    win_path = ""
                if win_path:
                    ok, _pid = QProcess.startDetached("explorer.exe", [win_path])
                    return ok

        if sys.platform == "darwin" and shutil.which("open"):
            ok, _pid = QProcess.startDetached("open", [str(path)])
            return ok

        if os.name == "nt" and shutil.which("explorer.exe"):
            ok, _pid = QProcess.startDetached("explorer.exe", [str(path)])
            return ok

        if shutil.which("xdg-open"):
            ok, _pid = QProcess.startDetached("xdg-open", [str(path)])
            return ok

        url = QUrl.fromLocalFile(str(path))
        return QDesktopServices.openUrl(url)

    def _open_text_viewer(self, path: Path) -> None:
        dialog = TextViewerDialog(path, self)
        dialog.setModal(True)
        dialog.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self._active_viewer = dialog
        dialog.finished.connect(lambda *_: setattr(self, "_active_viewer", None))
        dialog.open()

    def _open_path_fallback(self, path: Path, reason: str) -> None:
        if path.is_file() and looks_texty(path):
            self._open_text_viewer(path)
            return

        message = (
            f"{reason}\n\n"
            f"Path:\n{path}"
        )
        self._show_box(QMessageBox.Icon.Information, "Open manually", message)

    def _open_detail_paths(
        self,
        kind: str,
        prefer_in_app: bool = False,
        system: bool = False,
    ) -> None:
        paths = self._detail_paths.get(kind, [])
        if not paths:
            self._show_box(QMessageBox.Icon.Information, "No file", "No file is available for this selection.")
            return
        self._open_paths(paths, prefer_in_app=prefer_in_app, system=system)

    def _open_paths(
        self,
        paths: list[Path],
        prefer_in_app: bool = False,
        system: bool = False,
    ) -> None:
        if not paths:
            return
        if len(paths) > 5:
            reply = QMessageBox.question(
                self,
                "Open multiple files",
                f"Open {len(paths)} files?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if reply != QMessageBox.StandardButton.Yes:
                return
        for path in paths:
            if system:
                self._open_path_system(path)
            else:
                self._open_path(path, prefer_in_app=prefer_in_app)

    def _show_cell_text(self, table: QTableWidget, row: int, col: int) -> None:
        item = table.item(row, col)
        if not item:
            return
        text = item.text()
        if not text:
            return
        dialog = TextBlobDialog("Cell details", text, self)
        dialog.setModal(True)
        dialog.setAttribute(QtCore.Qt.WidgetAttribute.WA_DeleteOnClose, True)
        dialog.open()
 

    @Slot(object)
    def _close_after_worker_finished(self, _result=None) -> None: 
        self.close() 

    def _load_output_tables(self, paths: dict) -> None:
        assembly_summary = paths.get("assembly_summary")
        if assembly_summary:
            self._load_summary_table(Path(assembly_summary))

        asm_dir = None
        if assembly_summary:
            asm_dir = Path(assembly_summary).parent
        elif paths.get("fasta"):
            asm_dir = Path(paths.get("fasta")).parent

        blast_tsv = asm_dir / "blast_inputs.tsv" if asm_dir else None
        if blast_tsv and blast_tsv.exists():
            self._load_blast_inputs_table(blast_tsv)

        overlap_audit = paths.get("overlap_audit")
        if overlap_audit:
            self._load_overlap_audit_table(Path(overlap_audit))

        compare_tsv = paths.get("compare_assemblers")
        if compare_tsv:
            self._load_compare_assemblers_table(Path(compare_tsv))
        elif asm_dir:
            compare_candidate = asm_dir / "compare_assemblers.tsv"
            if compare_candidate.exists():
                self._load_compare_assemblers_table(compare_candidate)

        if asm_dir:
            self._asm_dir = asm_dir
            self._blast_inputs_fasta = self._asm_dir / "blast_inputs.fasta"
            self.open_blast_inputs_btn.setEnabled(self._blast_inputs_fasta.exists())
            self.open_blast_inputs_system_btn.setEnabled(self._blast_inputs_fasta.exists())
            self.open_asm_folder_btn.setEnabled(self._asm_dir.exists())

    @staticmethod
    def _fmt_table_value(value: object) -> str:
        if value is None:
            return ""
        text = str(value)
        if text.lower() == "nan":
            return ""
        return text

    def _load_summary_table(self, summary_path: Path) -> None:
        if not summary_path.exists():
            return
        df = pd.read_csv(summary_path, sep="\t")
        self.summary_table.setRowCount(0)
        self._summary_rows.clear()
        self.summary_filter.blockSignals(True)
        self.summary_filter.clear()
        self.summary_filter.addItem("All statuses", userData=None)

        status_colors = {
            "assembled": QColor("#C8E6C9"),
            "singlets_only": QColor("#FFF9C4"),
            "cap3_no_output": QColor("#FFE0B2"),
            "pair_missing": QColor("#FFCDD2"),
        }

        trace_colors = {
            "FAIL": QColor("#FFCDD2"),
            "WARN": QColor("#FFE0B2"),
            "PASS": QColor("#C8E6C9"),
            "NA": QColor("#E0E0E0"),
        }

        for row_idx, row in df.iterrows():
            self.summary_table.insertRow(row_idx)
            row_values = [
                self._fmt_table_value(row.get("sample_id", "")),
                self._fmt_table_value(row.get("status", "")),
                self._fmt_table_value(row.get("assembler", "")),
                self._fmt_table_value(row.get("contig_len", "")),
                self._fmt_table_value(row.get("blast_payload", "")),
                self._fmt_table_value(row.get("selected_engine", "")),
                self._fmt_table_value(row.get("configured_engine", "")),
                self._fmt_table_value(row.get("merge_status", "")),
                self._fmt_table_value(row.get("merge_overlap_len", "")),
                self._fmt_table_value(row.get("merge_identity", "")),
                self._fmt_table_value(row.get("overlaps_saved", "")),
                self._fmt_table_value(row.get("overlaps_removed", "")),
                self._fmt_table_value(row.get("primer_mode", "")),
                self._fmt_table_value(row.get("primer_stage", "")),
                self._fmt_table_value(row.get("trace_status", "NA")),
            ]
            for col, value in enumerate(row_values):
                item = QTableWidgetItem(value)
                item.setToolTip(value) 
                if col == 1:
                    status = value
                    color = status_colors.get(status, QColor("#E1F5FE") if status.startswith("overlap_") else None)
                    if color:
                        item.setBackground(QBrush(color))
                elif col == 14:
                    color = trace_colors.get(value.strip().upper(), QColor("#E0E0E0"))
                    item.setBackground(QBrush(color))
                self.summary_table.setItem(row_idx, col, item)
            self._summary_rows[row_values[0]] = {
                "status": row_values[1],
                "assembler": row_values[2],
                "contig_len": row_values[3],
                "blast_payload": row_values[4],
                "selected_engine": row_values[5],
                "configured_engine": row_values[6],
                "merge_status": row_values[7],
                "merge_overlap_len": row_values[8],
                "merge_identity": row_values[9],
                "overlaps_saved": row_values[10],
                "overlaps_removed": row_values[11],
                "primer_mode": row_values[12],
                "primer_stage": row_values[13],
                "trace_status": row_values[14],
                "trace_flags": self._fmt_table_value(row.get("trace_flags", "")),
                "trace_mixed_peak_frac_max": self._fmt_table_value(row.get("trace_mixed_peak_frac_max", "")),
                "review_reason": self._fmt_table_value(row.get("review_reason", "")),
            }
        for status in sorted(df["status"].dropna().unique()):
            self.summary_filter.addItem(status, userData=status)
        self.summary_filter.blockSignals(False)
        self._apply_summary_filter()
        self._configure_table_view(self.summary_table) 

    def _load_blast_inputs_table(self, blast_tsv: Path) -> None:
        df = pd.read_csv(blast_tsv, sep="\t")
        self.blast_table.setRowCount(0)
        self._blast_rows.clear()
        trace_colors = {
            "FAIL": QColor("#FFCDD2"),
            "WARN": QColor("#FFE0B2"),
            "PASS": QColor("#C8E6C9"),
            "NA": QColor("#E0E0E0"),
        }
        for row_idx, row in df.iterrows():
            self.blast_table.insertRow(row_idx)
            row_values = [
                self._fmt_table_value(row.get("sample_id", "")),
                self._fmt_table_value(row.get("blast_payload", "")),
                self._fmt_table_value(row.get("reason", "")),
                self._fmt_table_value(row.get("payload_ids", "")),
                self._fmt_table_value(row.get("trace_status", "NA")),
            ]
            for col, value in enumerate(row_values):
                item = QTableWidgetItem(value)
                item.setToolTip(value)
                if col == 4:
                    color = trace_colors.get(value.strip().upper(), QColor("#E0E0E0"))
                    item.setBackground(QBrush(color))
                self.blast_table.setItem(row_idx, col, item)
            self._blast_rows[row_values[0]] = {
                "blast_payload": row_values[1],
                "reason": row_values[2],
                "payload_ids": row_values[3],
                "trace_status": row_values[4],
                "trace_flags": self._fmt_table_value(row.get("trace_flags", "")),
                "trace_mixed_peak_frac_max": self._fmt_table_value(row.get("trace_mixed_peak_frac_max", "")),
                "review_reason": self._fmt_table_value(row.get("review_reason", "")),
            }
        self._configure_table_view(self.blast_table)

    def _load_overlap_audit_table(self, audit_path: Path) -> None:
        if not audit_path.exists():
            return
        df = pd.read_csv(audit_path, sep="\t")
        self.diagnostics_table.setRowCount(0)
        self._audit_rows.clear()
        for row_idx, row in df.iterrows():
            self.diagnostics_table.insertRow(row_idx)
            row_values = [
                self._fmt_table_value(row.get("sample_id", "")),
                self._fmt_table_value(row.get("overlap_len", "")),
                self._fmt_table_value(row.get("overlap_identity", "")),
                self._fmt_table_value(row.get("overlap_quality", "")),
                self._fmt_table_value(row.get("orientation", "")),
                self._fmt_table_value(row.get("status", "")),
                self._fmt_table_value(row.get("best_identity", "")),
                self._fmt_table_value(row.get("best_identity_orientation", "")),
                self._fmt_table_value(row.get("anchoring_feasible", "")),
                self._fmt_table_value(row.get("end_anchored_possible", "")),
                self._fmt_table_value(row.get("fwd_best_identity", "")),
                self._fmt_table_value(row.get("revcomp_best_identity", "")),
                self._fmt_table_value(row.get("fwd_best_overlap_len", "")),
                self._fmt_table_value(row.get("revcomp_best_overlap_len", "")),
                self._fmt_table_value(row.get("fwd_anchor_feasible", "")),
                self._fmt_table_value(row.get("revcomp_anchor_feasible", "")),
                self._fmt_table_value(row.get("identity_delta_revcomp_minus_fwd", "")),
                self._fmt_table_value(row.get("selected_vs_best_identity_delta", "")),
                self._fmt_table_value(row.get("top_candidate_count", "")),
                self._fmt_table_value(row.get("top2_identity_delta", "")),
                self._fmt_table_value(row.get("top2_overlap_len_delta", "")),
                self._fmt_table_value(row.get("top2_quality_delta", "")),
                self._fmt_table_value(row.get("tie_reason_code", "")),
                self._fmt_table_value(row.get("fwd_best_identity_any", "")),
                self._fmt_table_value(row.get("revcomp_best_identity_any", "")),
                self._fmt_table_value(row.get("fwd_best_overlap_len_any", "")),
                self._fmt_table_value(row.get("revcomp_best_overlap_len_any", "")),
                self._fmt_table_value(row.get("pretrim_best_identity", "")),
                self._fmt_table_value(row.get("posttrim_best_identity", "")),
                self._fmt_table_value(row.get("pretrim_best_overlap_len", "")),
                self._fmt_table_value(row.get("posttrim_selected_overlap_len", row.get("posttrim_best_overlap_len", ""))),
                self._fmt_table_value(row.get("pretrim_status", "")),
                self._fmt_table_value(row.get("posttrim_status", "")),
                self._fmt_table_value(row.get("ambiguity_identity_delta_used", "")),
                self._fmt_table_value(row.get("ambiguity_quality_epsilon_used", "")),
                self._fmt_table_value(row.get("primer_trim_bases_fwd", "")),
                self._fmt_table_value(row.get("primer_trim_bases_rev", "")),
                self._fmt_table_value(row.get("selected_engine", "")),
                self._fmt_table_value(row.get("fallback_used", "")),
                self._fmt_table_value(row.get("configured_engine", row.get("overlap_engine", ""))),
            ]
            for col, value in enumerate(row_values):
                item = QTableWidgetItem(value)
                item.setToolTip(value)
                self.diagnostics_table.setItem(row_idx, col, item)
            self._audit_rows[row_values[0]] = {
                "overlap_len": row_values[1],
                "overlap_identity": row_values[2],
                "overlap_quality": row_values[3],
                "orientation": row_values[4],
                "status": row_values[5],
                "best_identity": row_values[6],
                "best_identity_orientation": row_values[7],
                "anchoring_feasible": row_values[8],
                "end_anchored_possible": row_values[9],
                "fwd_best_identity": row_values[10],
                "revcomp_best_identity": row_values[11],
                "fwd_best_overlap_len": row_values[12],
                "revcomp_best_overlap_len": row_values[13],
                "fwd_anchor_feasible": row_values[14],
                "revcomp_anchor_feasible": row_values[15],
                "identity_delta_revcomp_minus_fwd": row_values[16],
                "selected_vs_best_identity_delta": row_values[17],
                "top_candidate_count": row_values[18],
                "top2_identity_delta": row_values[19],
                "top2_overlap_len_delta": row_values[20],
                "top2_quality_delta": row_values[21],
                "tie_reason_code": row_values[22],
                "fwd_best_identity_any": row_values[23],
                "revcomp_best_identity_any": row_values[24],
                "fwd_best_overlap_len_any": row_values[25],
                "revcomp_best_overlap_len_any": row_values[26],
                "pretrim_best_identity": row_values[27],
                "posttrim_best_identity": row_values[28],
                "pretrim_best_overlap_len": row_values[29],
                "posttrim_selected_overlap_len": row_values[30],
                "pretrim_status": row_values[31],
                "posttrim_status": row_values[32],
                "ambiguity_identity_delta_used": row_values[33],
                "ambiguity_quality_epsilon_used": row_values[34],
                "primer_trim_bases_fwd": row_values[35],
                "primer_trim_bases_rev": row_values[36],
                "selected_engine": row_values[37],
                "fallback_used": row_values[38],
                "overlap_engine": row_values[39],
                "configured_engine": row_values[39],
            }
        self._configure_table_view(self.diagnostics_table)

    def _load_compare_assemblers_table(self, compare_path: Path) -> None:
        if not compare_path.exists():
            return
        df = pd.read_csv(compare_path, sep="\t")
        self.compare_table.setRowCount(0)
        self._compare_rows.clear() 
        for row_idx, row in df.iterrows():
            self.compare_table.insertRow(row_idx)
            row_values = [
                self._fmt_table_value(row.get("sample_id", "")),
                self._fmt_table_value(row.get("assembler_id", "")),
                self._fmt_table_value(row.get("assembler_name", "")),
                self._fmt_table_value(row.get("dup_policy", "")),
                self._fmt_table_value(row.get("status", "")),
                self._fmt_table_value(row.get("selected_engine", "")),
                self._fmt_table_value(row.get("contig_len", "")),
                self._fmt_table_value(row.get("diag_code_for_machine", "")),
                self._fmt_table_value(row.get("diag_detail_for_human", "")),
                self._fmt_table_value(row.get("cap3_contigs_n", "")),
                self._fmt_table_value(row.get("cap3_singlets_n", "")),
                self._fmt_table_value(row.get("warnings", "")),
            ]
            for col, value in enumerate(row_values):
                item = QTableWidgetItem(value)
                item.setToolTip(value)
                self.compare_table.setItem(row_idx, col, item)
            sample_id = row_values[0]
            assembler_id = row_values[1]
            self._compare_rows[(sample_id, assembler_id)] = {
                "sample_id": sample_id,
                "assembler_id": assembler_id,
                "assembler_name": row_values[2],
                "dup_policy": row_values[3],
                "status": row_values[4],
                "selected_engine": row_values[5],
                "contig_len": row_values[6],
                "warnings": row_values[11],
                "diag_code_for_machine": self._fmt_table_value(row.get("diag_code_for_machine", "")),
                "diag_detail_for_human": self._normalize_compare_diag_detail(
                    self._fmt_table_value(row.get("diag_detail_for_human", "")),
                    row_values[5],
                ),
                "cap3_contigs_n": self._fmt_table_value(row.get("cap3_contigs_n", "")),
                "cap3_singlets_n": self._fmt_table_value(row.get("cap3_singlets_n", "")),
                "cap3_info_path": self._fmt_table_value(row.get("cap3_info_path", "")),
                "cap3_stdout_path": self._fmt_table_value(row.get("cap3_stdout_path", "")),
                "cap3_stderr_path": self._fmt_table_value(row.get("cap3_stderr_path", "")),
                "payload_fasta": self._fmt_table_value(row.get("payload_fasta", "")),
                "payload_kind": self._fmt_table_value(row.get("payload_kind", "none")),
                "payload_n": self._fmt_table_value(row.get("payload_n", "0")),
                "payload_max_len": self._fmt_table_value(row.get("payload_max_len", "0")),
                "ambiguity_flag": self._fmt_table_value(row.get("ambiguity_flag", "0")),
                "safety_flag": self._fmt_table_value(row.get("safety_flag", "none")),
                "review_action": self._fmt_table_value(row.get("review_action", row.get("status", "none"))),
                "review_reason": self._fmt_table_value(row.get("review_reason", "")),
                "advisory_reason": self._fmt_table_value(row.get("advisory_reason", "")),
                "warning_flags": self._fmt_table_value(row.get("warning_flags", "")),
            } 

        self._configure_table_view(self.compare_table)

    def _apply_summary_filter(self) -> None:
        selected = self.summary_filter.currentData()
        for row in range(self.summary_table.rowCount()):
            status_item = self.summary_table.item(row, 1)
            status = status_item.text() if status_item else ""
            self.summary_table.setRowHidden(row, selected is not None and status != selected)

    def _configure_table_view(self, table: QTableWidget) -> None:
        table.setWordWrap(False)
        header = table.horizontalHeader()
        if table.rowCount() <= 500:
            table.resizeColumnsToContents()
            header.setSectionResizeMode(QHeaderView.ResizeToContents)
        else:
            header.setSectionResizeMode(QHeaderView.ResizeToContents)
            header.setSectionResizeMode(0, QHeaderView.ResizeToContents)
            if table.columnCount() > 1:
                header.setSectionResizeMode(1, QHeaderView.ResizeToContents)
            if table.columnCount() > 2:
                header.setSectionResizeMode(2, QHeaderView.ResizeToContents)
            if table.columnCount() > 3:
                header.setSectionResizeMode(3, QHeaderView.Interactive)
                table.setColumnWidth(3, 260)
        _set_header_tooltips(
            table,
            {
                "status": STATUS_COLOR_LEGEND,
                "trace_qc": TRACE_QC_COLOR_LEGEND,
                "trace_qc_status": TRACE_QC_COLOR_LEGEND,
            },
        )
        header.setStretchLastSection(True)
        if table.rowCount() <= 200:
            table.resizeRowsToContents()

    def _selected_sample_ids(self, table: QTableWidget) -> list[str]:
        model = table.selectionModel()
        if not model:
            return []
        rows = model.selectedRows(0)
        if rows:
            return [table.item(idx.row(), 0).text() for idx in rows if table.item(idx.row(), 0)]
        current_row = table.currentRow()
        if current_row >= 0:
            item = table.item(current_row, 0)
            if item:
                return [item.text()]
        return []

    def _selected_compare_keys(self) -> list[tuple[str, str]]:
        """ This is to pair the sample id and asm id together to help with troubleshooting in the compare assembly tab.""" 
        model = self.compare_table.selectionModel() 
        if not model:
            return [] 
        rows = model.selectedRows(0)
        selected: list[tuple[str, str]] = [] 
        for idx in rows: 
            row = idx.row() 
            sid_item = self.compare_table.item(row, 0)
            asm_item = self.compare_table.item(row, 1)
            if sid_item and asm_item:
                selected.append((sid_item.text(), asm_item.text()))

            if selected:
                return selected 

            row = self.compare_table.currentRow()
            if row >= 0: 
                sid_item = self.compare_table.item(row, 0)
                asm_item = self.compare_table.item(row, 1)  
                if sid_item and asm_item:
                    return [(sid_item.text(), asm_item.text())]
            return [] 

    def _join_values(self, values: list[str]) -> str:
        cleaned = [v for v in values if v]
        if not cleaned:
            return "—"
        unique = list(dict.fromkeys(cleaned))
        if len(unique) <= 5:
            return ", ".join(unique)
        return f"{', '.join(unique[:5])}, … (+{len(unique) - 5} more)"

    def _normalize_compare_diag_detail(self, detail: str, selected_engine: str) -> str:
        return _normalize_legacy_compare_diag_detail(detail, selected_engine)

    def _set_detail_label(self, label: QLabel, prefix: str, summary: str, full: str | None = None) -> None:
        compact = (summary or "—").strip() or "—"
        label.setText(f"{prefix}: {compact}")
        tooltip = (full if full is not None else compact) or ""
        label.setToolTip(tooltip)

    def _update_detail_panel(self) -> None:
        table = None
        tab_index = self.output_tabs.currentIndex()
        if tab_index == 1:
            table = self.summary_table
        elif tab_index == 2:
            table = self.blast_table
        elif tab_index == 3:
            table = self.diagnostics_table
        elif tab_index == 4:
            table = self.compare_table
        else:
            table = self.summary_table if self.summary_table.selectedItems() else self.blast_table

        sample_ids = self._selected_sample_ids(table) if table else []
        if not sample_ids:
            return

        compare_keys = self._selected_compare_keys() if tab_index == 4 else [] 
        if compare_keys:
            rows = [self._compare_rows.get(k, {}) for k in compare_keys]
            self.detail_sample_lbl.setText(
                f"sample_id: {self._join_values([r.get('sample_id', '') for r in rows])}"
            )
            self.detail_status_lbl.setText(
                f"status: {self._join_values([r.get('status', '-') for r in rows])}"
            )
            reason_bits = [] 
            for r in rows:
                parts = [
                    r.get("diag_code_for_machine", ""),
                    r.get("diag_detail_for_human", ""),
                    r.get("warnings", ""),
                    r.get("safety_flag", ""),
                    r.get("review_action", ""),
                    r.get("review_reason", ""),
                    r.get("advisory_reason", ""),
                    r.get("warning_flags", ""),
                ]
                reason_bits.append(" | ".join(p for p in parts if p and p != "none") or "-")
            self.detail_reason_lbl.setText(f"reason: {self._join_values(reason_bits)}")
            payload_vals = [f"{r.get('payload_kind', 'none')} (n={r.get('payload_n', '0')})" for r in rows]
            self.detail_payload_lbl.setText(f"blast_payload: {self._join_values(payload_vals)}")
            cap3_meta = [] 
            for r in rows:
                if r.get("cap3_contigs_n") or r.get("cap3_singlets_n"):
                    cap3_meta.append(
                        f"{r.get('assembler_id', '')} contigs={r.get('cap3_contigs_n', '0')} singlets= {r.get('cap3_singlets_n', '0')}"
                    )
                else:
                    cap3_meta.append(r.get("assembler_id", ""))
            self.detail_payload_ids_lbl.setText(f"payload_ids: {self._join_values(cap3_meta)}")
            self.detail_assembly_engine_lbl.setText(
                f"assembly path: {self._join_values([r.get('assembler_name', '') for r in rows])}"
            )

            trace_status_vals = [str(r.get("trace_status", "NA") or "NA").upper() for r in rows]
            trace_labels = [
                "; ".join(
                    _trace_flags_to_labels(
                        r.get("trace_flags", ""),
                        r.get("trace_status", "NA"),
                        r.get("trace_mixed_peak_frac_max", ""),
                        r.get("review_reason", ""),
                    )
                )
                for r in rows
            ]
            self.detail_trace_qc_lbl.setText(f"trace_qc: {self._join_values(trace_status_vals)}")
            self.detail_trace_reason_lbl.setText(f"trace_qc_detail: {self._join_values(trace_labels)}")

            audit_values = [self._audit_rows.get(r.get("sample_id", ""), {}) for r in rows]
            statuses = [audit.get("status", "-") for audit in audit_values if audit]
            self.detail_overlap_lbl.setText(f"overlap: {self._join_values(statuses)}")

            for key in self._detail_paths:
                self._detail_paths[key] = []
            for r in rows:
                payload = r.get("payload_fasta", "")
                if payload:
                    p = Path(payload)
                    if p.exists():
                        self._detail_paths["contigs"].append(p)
                info_txt = r.get("cap3_info_path", "")
                if info_txt:
                    info_path = Path(info_txt)
                    if info_path.exists():
                        self._detail_paths["info"].append(info_path)
                        singlet_path = info_path.with_suffix(".singlets")
                        if singlet_path.exists():
                            self._detail_paths["singlets"].append(singlet_path)

            self.detail_contigs_btn.setEnabled(bool(self._detail_paths["contigs"]))
            self.detail_singlets_btn.setEnabled(bool(self._detail_paths["singlets"]))
            self.detail_info_btn.setEnabled(bool(self._detail_paths["info"]))
            self.detail_contigs_system_btn.setEnabled(bool(self._detail_paths["contigs"]))
            self.detail_singlets_system_btn.setEnabled(bool(self._detail_paths["singlets"]))
            self.detail_info_system_btn.setEnabled(bool(self._detail_paths["info"]))
            return

        summary_vals = [self._summary_rows.get(sid, {}).get("status", "—") for sid in sample_ids]
        reason_vals = [self._blast_rows.get(sid, {}).get("reason", "—") for sid in sample_ids]
        payload_vals = [self._blast_rows.get(sid, {}).get("blast_payload", "—") for sid in sample_ids]
        payload_ids_vals = [self._blast_rows.get(sid, {}).get("payload_ids", "—") for sid in sample_ids]

        self.detail_sample_lbl.setText(f"sample_id: {self._join_values(sample_ids)}")
        self.detail_status_lbl.setText(f"status: {self._join_values(summary_vals)}")
        self.detail_reason_lbl.setText(f"reason: {self._join_values(reason_vals)}")
        self.detail_payload_lbl.setText(f"blast_payload: {self._join_values(payload_vals)}")
        self.detail_payload_ids_lbl.setText(f"payload_ids: {self._join_values(payload_ids_vals)}")

        engine_vals: list[str] = []
        for sid in sample_ids:
            srow = self._summary_rows.get(sid, {})
            merge_status = srow.get("merge_status", "")
            status = srow.get("status", "")
            if merge_status == "merged":
                engine_vals.append("merged_two_reads")
            elif merge_status == "high_conflict":
                engine_vals.append("merged_two_reads -> CAP3 fallback (high_conflict)")
            elif merge_status in {"identity_low", "overlap_too_short", "quality_low", "ambiguous_overlap"}:
                engine_vals.append(f"merged_two_reads attempted ({merge_status}) -> CAP3 fallback")
            elif status == "assembled":
                engine_vals.append("CAP3 assembled")
            elif status == "singlets_only":
                engine_vals.append("CAP3 singlets_only")
            else:
                engine_vals.append("paired flow (see merge_status + status)")
        self.detail_assembly_engine_lbl.setText(f"assembly path: {self._join_values(engine_vals)}")
        
        trace_status_vals = [
            str(
                self._summary_rows.get(sid, {}).get("trace_status")
                or self._blast_rows.get(sid, {}).get("trace_status")
                or "NA"
            ).upper()
            for sid in sample_ids
        ]
        trace_labels = []
        for sid in sample_ids:
            src = self._summary_rows.get(sid, {})
            if not src:
                src = self._blast_rows.get(sid, {})
            trace_labels.append(
                "; ".join(
                    _trace_flags_to_labels(
                        src.get("trace_flags", ""),
                        src.get("trace_status", "NA"),
                        src.get("trace_mixed_peak_frac_max", ""),
                        src.get("review_reason", ""),
                    )
                )
            )

        self.detail_trace_qc_lbl.setText(f"trace_qc: {self._join_values(trace_status_vals)}")
        self.detail_trace_reason_lbl.setText(f"trace_qc_detail: {self._join_values(trace_labels)}")

        audit_values = [self._audit_rows.get(sid, {}) for sid in sample_ids]
        if len(sample_ids) == 1 and audit_values and audit_values[0]:
            audit = audit_values[0]
            overlap_summary = (
                f"{audit.get('status', '—')} (len {audit.get('overlap_len', '—')}, "
                f"id {audit.get('overlap_identity', '—')}, "
                f"q {audit.get('overlap_quality', '—')}, "
                f"orient {audit.get('orientation', '—')}, "
                f"engine {audit.get('selected_engine', '—')}, "
                f"fallback {audit.get('fallback_used', '—')})"
            )
            overlap_full = (
                f"{overlap_summary}; best-id {audit.get('best_identity', '—')} @ {audit.get('best_identity_orientation', '—')}; "
                f"fwd-best-id {audit.get('fwd_best_identity', '—')} (any {audit.get('fwd_best_identity_any', '—')}), "
                f"rev-best-id {audit.get('revcomp_best_identity', '—')} (any {audit.get('revcomp_best_identity_any', '—')}); "
                f"delta(rev-fwd) {audit.get('identity_delta_revcomp_minus_fwd', '—')}; "
                f"selected-vs-best-delta {audit.get('selected_vs_best_identity_delta', '—')}; "
                f"fwd-anchor-feasible {audit.get('fwd_anchor_feasible', '—')}, rev-anchor-feasible {audit.get('revcomp_anchor_feasible', '—')}; "
                f"top-cands {audit.get('top_candidate_count', '—')}, top2-id-delta {audit.get('top2_identity_delta', '—')}, "
                f"top2-len-delta {audit.get('top2_overlap_len_delta', '—')}, top2-q-delta {audit.get('top2_quality_delta', '—')}; "
                f"tie-reason {audit.get('tie_reason_code', '—')} (id_eps {audit.get('ambiguity_identity_delta_used', '—')}, "
                f"q_eps {audit.get('ambiguity_quality_epsilon_used', '—')}); "
                f"pretrim {audit.get('pretrim_status', '—')} id {audit.get('pretrim_best_identity', '—')} len {audit.get('pretrim_best_overlap_len', '—')}; "
                f"posttrim {audit.get('posttrim_status', '—')} id {audit.get('posttrim_best_identity', '—')} "
                f"len-selected {audit.get('posttrim_selected_overlap_len', '—')}; "
                f"primer-trim-bases F {audit.get('primer_trim_bases_fwd', '—')} R {audit.get('primer_trim_bases_rev', '—')}; "
                f"anchored-feasible {audit.get('anchoring_feasible', '—')}; "
                f"end-anchored-possible {audit.get('end_anchored_possible', '—')}; "
                f"configured-engine {audit.get('overlap_engine', '—')}"
            )
            self._set_detail_label(self.detail_overlap_lbl, "overlap", overlap_summary, overlap_full)
        else:
            statuses = [audit.get("status", "—") for audit in audit_values if audit]
            self._set_detail_label(self.detail_overlap_lbl, "overlap", self._join_values(statuses))

        for key in self._detail_paths:
            self._detail_paths[key] = []

        if self._asm_dir:
            for sample_id in sample_ids:
                sample_dir = self._asm_dir / sample_id
                contig_path = sample_dir / f"{sample_id}_paired.fasta.cap.contigs"
                singlet_path = sample_dir / f"{sample_id}_paired.fasta.cap.singlets"
                info_path = sample_dir / f"{sample_id}_paired.fasta.cap.info"
                if contig_path.exists():
                    self._detail_paths["contigs"].append(contig_path)
                if singlet_path.exists():
                    self._detail_paths["singlets"].append(singlet_path)
                if info_path.exists():
                    self._detail_paths["info"].append(info_path)

        self.detail_contigs_btn.setEnabled(bool(self._detail_paths["contigs"]))
        self.detail_singlets_btn.setEnabled(bool(self._detail_paths["singlets"]))
        self.detail_info_btn.setEnabled(bool(self._detail_paths["info"]))
        self.detail_contigs_system_btn.setEnabled(bool(self._detail_paths["contigs"]))
        self.detail_singlets_system_btn.setEnabled(bool(self._detail_paths["singlets"]))
        self.detail_info_system_btn.setEnabled(bool(self._detail_paths["info"]))
 
    # closeEvent ----------------------------
    def closeEvent(self, event):
        """
        Prevent the window from clsoing while BLAST thread is still running.
        """
        if self._thread and self._thread.isRunning():
            if self._worker:
                # ask the thread to finish, then auto-close the window
                self._worker.finished.connect(self._close_after_worker_finished, type=QtCore.Qt.ConnectionType.QueuedConnection) 
                self._thread.requestInterruption()     # politely signal
                event.ignore()                         # keep window open
                self.log_model.append("Waiting for BLAST thread to finish…")
                return
        is_max = self.isMaximized()
        normal_rect = self.normalGeometry() if is_max else self.geometry()
        self.settings.setValue(
            "window_normal_geometry",
            [normal_rect.x(), normal_rect.y(), normal_rect.width(), normal_rect.height()],
        )
        self.settings.setValue("window_start_maximized", is_max)
        self._cleanup_input_staging()
        self._flush_logs() 
        event.accept()
        root = logging.getLogger()
        if self._log_handler and self._log_handler in root.handlers:
            root.removeHandler(self._log_handler)
            self._log_handler = None

# redirects all Qt internal C++ warnings here for debugging 
def qt_message_handler(mode, _context, message):
    """Redirect Qt messages to the Python logging module."""
    level = {
        QtCore.QtMsgType.QtDebugMsg:    logging.DEBUG,
        QtCore.QtMsgType.QtInfoMsg:     logging.INFO,
        QtCore.QtMsgType.QtWarningMsg:  logging.WARNING,
        QtCore.QtMsgType.QtCriticalMsg: logging.ERROR,
        QtCore.QtMsgType.QtFatalMsg:    logging.CRITICAL,
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
    QTimer.singleShot(0, lambda: win.showMaximized() if win._start_maximized else None)

    logging.info(
        "GUI platform session: XDG_SESSION_TYPE=%s QT_QPA_PLATFORM=%s WAYLAND_DISPLAY=%s",
        os.environ.get("XDG_SESSION_TYPE", ""),
        os.environ.get("QT_QPA_PLATFORM", ""),
        os.environ.get("WAYLAND_DISPLAY", ""),
    )
    logging.info(
        "Qt backend decision: chosen_qt_backend=%s reason=%s",
        os.environ.get("QT_QPA_PLATFORM", ""),
        os.environ.get("MICROSEQ_QT_BACKEND_REASON", "unknown"),
    )
    logging.info(
        "Logging vs outputs: execution logs go to logs/microseq_<session>.log; run artifacts (qc/, asm/, hits.tsv, hits_tax.tsv) are written under each run output directory."
    )
    sys.exit(app.exec())

# allow: python -m microseq_tests.gui
if __name__ == "__main__":
    launch()
