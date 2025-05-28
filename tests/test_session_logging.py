import os, logging, time
import subprocess
from pathlib import Path
import re
import pytest

# the helper under test
import importlib.util, types, sys

ROOT = Path(__file__).resolve().parents[1]
UTIL_PATH = ROOT / "src" / "microseq_tests" / "utility" / "utils.py"
yaml_stub = types.ModuleType('yaml')
yaml_stub.safe_load = lambda fh: {}
sys.modules.setdefault('yaml', yaml_stub)

spec = importlib.util.spec_from_file_location("mutils", UTIL_PATH)
mutils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mutils)
setup_logging = mutils.setup_logging


# ──────────────────────────────────────────────────────────────
def test_auto_session_generates_unique_file(tmp_path, monkeypatch, capsys):
    """
    If no MICROSEQ_SESSION_ID is set, each invocation gets a fresh
    YYYYMMDD-HHMMSS-XXXX log file and prints a single warning.
    """
    monkeypatch.delenv("MICROSEQ_SESSION_ID", raising=False)

    f1 = setup_logging(tmp_path, force=True, console=False)
    # wait 1 second so timestamp definitely changes
    time.sleep(1)
    f2 = setup_logging(tmp_path, force=True, console=False)

    # capture the warning once (first call only)
    err_out = capsys.readouterr().err
    assert "not set – using auto session ID" in err_out
    assert f1 != f2                                # unique filenames
    pattern = re.compile(r"microseq_\d{8}-\d{6}-[0-9a-f]{4}\.log$")
    assert pattern.search(f1.name)
    assert pattern.search(f2.name)


# ──────────────────────────────────────────────────────────────
def test_explicit_session_aggregates(tmp_path, monkeypatch, capsys):
    """
    With MICROSEQ_SESSION_ID defined, multiple calls reuse the same file
    and no warning is emitted.
    """
    monkeypatch.setenv("MICROSEQ_SESSION_ID", "demo123")

    f1 = setup_logging(tmp_path, force=True, console=False)
    f2 = setup_logging(tmp_path, force=True, console=False)

    assert f1 == f2                              # same path
    err_out = capsys.readouterr().err
    assert err_out == ""                         # no warning


# ──────────────────────────────────────────────────────────────
def test_prune_only_auto_logs(tmp_path, monkeypatch):
    """
    backup_count should delete old *auto* files but never touch explicit
    session logs.
    """
    monkeypatch.delenv("MICROSEQ_SESSION_ID", raising=False)

    # create 12 fake auto-logs + 1 explicit session log
    auto_names = [f"microseq_20250101-{i:06d}-abcd.log" for i in range(12)]
    for name in auto_names:
        (tmp_path / name).touch()
    (tmp_path / "microseq_myRun.log").touch()

    # keep last 5 → 7 oldest auto files should disappear
    setup_logging(tmp_path, backup_count=5, force=True, console=False)

    remaining = {p.name for p in tmp_path.iterdir()}
    auto_remaining = sorted(n for n in remaining if n.startswith("microseq_2025"))
    assert len(auto_remaining) == 6                   # 5 old logs + new one
    assert "microseq_myRun.log" in remaining          # untouched

# ──────────────────────────────────────────────────────────────
def test_session_id_flag_respected(tmp_path, monkeypatch):
    """setup_logging should honour MICROSEQ_SESSION_ID when set later."""
    monkeypatch.delenv("MICROSEQ_SESSION_ID", raising=False)
    setup_logging(tmp_path, force=True, console=False)

    monkeypatch.setenv("MICROSEQ_SESSION_ID", "cliDemo")
    f = setup_logging(tmp_path, force=True, console=False)

    assert f.name.endswith("microseq_cliDemo.log")

# ──────────────────────────────────────────────────────────────
def test_parent_child_logging(tmp_path, monkeypatch):
    """Child process logs to same session file as parent."""
    monkeypatch.setenv("MICROSEQ_SESSION_ID", "shared123")

    log_file = setup_logging(tmp_path, force=True, console=False)
    logging.getLogger().info("parent-line")
    logging.shutdown()

    child_script = tmp_path / "child.py"
    child_script.write_text(
        """
import importlib.util, types, os, sys, logging
from pathlib import Path
ROOT = Path(os.environ['REPO_ROOT'])
UTIL_PATH = ROOT / 'src' / 'microseq_tests' / 'utility' / 'utils.py'
yaml_stub = types.ModuleType('yaml')
yaml_stub.safe_load = lambda fh: {}
sys.modules.setdefault('yaml', yaml_stub)
spec = importlib.util.spec_from_file_location('mutils', UTIL_PATH)
mutils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mutils)
setup_logging = mutils.setup_logging
setup_logging(Path(os.environ['LOG_DIR']), force=True, console=False)
logging.getLogger().info('child-line')
logging.shutdown()
"""
    )

    env = os.environ.copy()
    env['REPO_ROOT'] = str(ROOT)
    env['LOG_DIR'] = str(tmp_path)
    subprocess.run([sys.executable, str(child_script)], check=True, env=env)

    logs = log_file.read_text()
    assert 'parent-line' in logs
    assert 'child-line' in logs

