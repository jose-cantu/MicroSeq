# ── src/microseq_tests/utility/utils.py ────────────────────────────────
from __future__ import annotations

import errno
import logging
import logging.handlers
import os
import secrets
import sys
import yaml
from pathlib import Path
import datetime as dt                       

# ── locate repo root & default log dir  ────────────────────────────────
def _find_repo_root(start: Path | None = None) -> Path:
    """Walk parents until we see pyproject.toml or .git."""
    here = start or Path(__file__).resolve()
    for p in [here, *here.parents]:
        if (p / "pyproject.toml").exists() or (p / ".git").exists():
            return p
    return Path(__file__).resolve().parents[1]       # site-packages wheel

ROOT      = _find_repo_root()
LOG_ROOT  = ROOT / "logs"
CONF_PATH = ROOT / "config" / "config.yaml"

# ── tiny helpers  ──────────────────────────────────────────────────────
def load_config(path: str | Path = CONF_PATH):
    with Path(path).open() as fh:
        return yaml.safe_load(fh)

def expand_db_path(template: str) -> str:
    """
    Replace ${MICROSEQ_DB_HOME} in template with the environment variable.
    Falls back to the value in config/config.yaml, or ~/.microseq_dbs.
    """
    db_home = os.getenv("MICROSEQ_DB_HOME")
    if not db_home:
        cfg = load_config()
        db_home = cfg.get("microseq_db_home") or "~/.microseq_dbs"
    return template.replace("${MICROSEQ_DB_HOME}", os.path.expanduser(db_home))

def set_module_level(module_name: str, level: int) -> None:
    logging.getLogger(module_name).setLevel(level)

# ── main helper  ───────────────────────────────────────────────────────
def setup_logging(
    log_dir: str | Path = LOG_ROOT,
    *,
    level: int | None = None,
    console: bool = True,
    force: bool = False,
    rotate_mb: int | None = None,
    max_bytes: int | None = None,
    backup_count: int = 0,                      # keep everything by default
    session_env: str = "MICROSEQ_SESSION_ID",
    warn_if_generated: bool = True,
    log_file_prefix: str = "microseq",
) -> Path:
    """
    Create one log file called 'microseq_<SESSION_ID>.log'.

    SESSION_ID priority
    1. value of $<session_env>  (e.g. MICROSEQ_SESSION_ID)
    2. auto-generated 'YYYYMMDD-HHMMSS-<4-hex>'

    Only auto-generated (timestamp-rand) files are ever pruned,
    and only when backup_count > 0 so I kept this behavior in by default. Might just remove pruning later on. 
    """

    # ── size-rotation helper ------------------------------------------
    if max_bytes:
        rotate_bytes = max_bytes
    elif rotate_mb:
        rotate_bytes = int(rotate_mb * 1024 * 1024)
    else:
        rotate_bytes = None

    # ── decide root dir (env-var beats arg beats default) -------------
    if os.getenv("MICROSEQ_LOG_FILE"):
        logfile = Path(os.getenv("MICROSEQ_LOG_FILE")).expanduser()
        logfile.parent.mkdir(parents=True, exist_ok=True)
        root_dir = logfile.parent
    else:
        root_dir = (
            Path(log_dir).expanduser()
            if log_dir is not None
            else Path(os.getenv("MICROSEQ_LOG_DIR", LOG_ROOT)).expanduser()
        )
        root_dir.mkdir(parents=True, exist_ok=True)

        # ── choose session ID -----------------------------------------
        sess_id = os.getenv(session_env)
        if not sess_id:
            ts = dt.datetime.now().strftime("%Y%m%d-%H%M%S")
            sess_id = f"{ts}-{secrets.token_hex(2)}"
            if warn_if_generated:
                sys.stderr.write(
                    f"⚠️  {session_env} not set – using auto session ID {sess_id}\n"
                    f"   (export {session_env}=YOUR_ID to aggregate multiple commands)\n"
                )

        # ── prune old auto logs --------------------------------------
        is_auto = sess_id.count("-") == 2         # timestamp-rand pattern
        if is_auto and backup_count:
            patt  = f"{log_file_prefix}_????????-??????-*.log"
            logs  = sorted(root_dir.glob(patt))   # oldest → newest
            excess = len(logs) - backup_count
            for old in logs[:excess]:
                try:
                    old.unlink()
                except OSError:
                    pass

        logfile = root_dir / f"{log_file_prefix}_{sess_id}.log"

    # ── short-circuit if already configured ---------------------------
    root_logger = logging.getLogger()
    if root_logger.handlers and not force:
        return logfile

    root_logger.handlers.clear()
    root_logger.setLevel(level or logging.INFO)

    fmt = logging.Formatter("%(asctime)s  %(levelname)-7s  %(name)s:  %(message)s")

    # ── choose handler -------------------------------------------------
    if rotate_bytes:
        fh = logging.handlers.RotatingFileHandler(
            logfile, maxBytes=rotate_bytes, backupCount=backup_count,
            encoding="utf-8", delay=True
        )
    else:
        fh = logging.FileHandler(logfile, mode="a", encoding="utf-8", delay=True)

    fh.setFormatter(fmt)
    root_logger.addHandler(fh)

    if console:
        ch = logging.StreamHandler(sys.stderr)
        ch.setFormatter(fmt)
        root_logger.addHandler(ch)

    # ── refresh _latest symlink ---------------------------------------
    latest = logfile.parent / f"{log_file_prefix}_latest.log"
    try:
        if latest.is_symlink() or latest.exists():
            latest.unlink()
        latest.symlink_to(logfile.name)          # relative link
    except OSError as e:
        if e.errno not in (errno.EPERM, errno.EACCES, errno.EEXIST):
            raise

    root_logger.info("Logging to %s", logfile)
    return logfile

