from __future__ import annotations 
import yaml 
import os
from pathlib import Path 
import logging, datetime 
import sys
from logging.handlers import RotatingFileHandler 

ROOT = Path(__file__).resolve().parents[3]  # repo root (think going 3 steps back here ../../.. until it sees top level of repo =) )  
CONF_PATH = ROOT / "config" / "config.yaml"
def load_config(config_path: str | Path = CONF_PATH):
    config_path = Path(config_path)
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

# default log directory (can be overwritten via function arg - keep that in mind) 
LOG_ROOT = ROOT/ "logs" 


def setup_logging(log_dir: str | Path = LOG_ROOT, log_file_prefix: str ="microseq", 
                  *, level: int | None = None, max_bytes: int = 20_000_000, backup_count: int = 3, console: bool = True, force: bool = False,) -> Path:
    """
       Configure a timestamped file + console logger.

    Parameters
    log_dir : str | pathlib.Path
        Folder for log files (created if missing).
    log_file_prefix : str, default "microseq"
        The filename stem; timestamp + ".log" is appended.
    level : int | None, default logging.INFO 
        Root verbosity if not already configured.
    max_bytes : int, default 20_000_000
        Size threshold before a log file rolls over.
    backup_count : int, default 3
        How many rolled files (.log.1 …) to keep.
    console : bool, default True
        Also echo log lines to stderr.
    force : bool, default False
        Reinstall handlers even if logging was already configured
       (useful inside pytest).
    """

    Path(log_dir).mkdir(exist_ok=True)
    # allow MICROSEQ_LOG_DIR to override fefault 
    log_dir = Path(
            os.environ.get("MICROSEQ_LOG_DIR", log_dir)
    ).expanduser()
    log_dir.mkdir(exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    file_path = Path(log_dir) / f"{log_file_prefix}_{timestamp}.log" 

    root = logging.getLogger() # guarding against double-installation 
    if root.handlers and not force:
        return 
    
    root.handlers.clear()
    root.setLevel(level or logging.INFO)# this is default if caller omitted so take note 
    
    fmt = logging.Formatter("%(asctime)s  %(levelname)-7s  %(name)s:  %(message)s")

    fh = RotatingFileHandler(
            file_path,
            maxBytes=max_bytes,
            backupCount=backup_count,
            )
    fh.setFormatter(fmt)
    root.addHandler(fh) 

    if console:  # console handler 
        ch = logging.StreamHandler(sys.stderr)
        ch.setFormatter(fmt)
        root.addHandler(ch) 


    root.info("Logging to %s", file_path)
    return file_path 

def expand_db_path(template: str) -> str:
    db_home = os.environ.get("MICROSEQ_DB_HOME", os.path.expanduser("~/.microseq_dbs"))
    return template.replace("${MICROSEQ_DB_HOME}", db_home)

def set_module_level(module_name: str, level: int) -> None: 
    """
    Change verbosity of one sub-logger at runtime. 

    Example
    ______ 
    >>>> set_module_level("microseq_tests.blast", logging.ERROR)
    """
    logging.getLogger(module_name).setLevel(level) 
