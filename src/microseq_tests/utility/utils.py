# ---- src/microseq_tests/utility/utils.py --------------------- 
from __future__ import annotations 
import yaml, logging, os, sys
from datetime import datetime
from pathlib import Path 
from logging.handlers import RotatingFileHandler
import importlib.resources as pkg_resources 

# default log directory (can be overwritten via function arg - keep that in mind)
def _find_repo_root(start: Path | None = None) -> Path: 
    """
    Climb parents until we hit a file that marks the top of the project (which in this case its pyproject.toml, or .git). Fallback: 
LOG_ROOT = ROOT/ "logs" # repo-local fallback option the package root inside site-packages. 
    """
    here = start or Path(__file__).resolve() 
    for p in [here, *here.parents]:
        if (p / "pyproject.toml").exists() or (p / ".git").exists():
            return p 

    # inside a wheel / site-packages 
    return Path(__file__).resolve().parents[1] # microseq_tests 

ROOT = _find_repo_root()
LOG_ROOT = ROOT / "logs" 

def expand_db_path(template: str) -> str:
    """
    Replace ${MICROSEQ_DB_HOME} in template with the environment variable. 
    If set, otherwise with the value stored in config/config.yaml, 
    If not then it defaults to ~/.microseq_dbs 
    """
    db_home = os.getenv("MICROSEQ_DB_HOME")
    if not db_home:
        cfg = load_config() 
        db_home = cfg.get("microseq_db_home") or "~/.microseq_dbs" 
    return template.replace("${MICROSEQ_DB_HOME}", os.path.expanduser(db_home)) 

def set_module_level(module_name: str, level: int) -> None: 
    """
    Change verbosity of one sub-logger at runtime. 

    Example
    ______
    from microseq_tests.utility.utils import set_module_level 
    >>>> set_module_level("microseq_tests.blast", logging.ERROR)
    """
    logging.getLogger(module_name).setLevel(level) 


CONF_PATH = ROOT / "config" / "config.yaml"
def load_config(config_path: str | Path = CONF_PATH):
    config_path = Path(config_path)
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def setup_logging(log_dir: str | Path = LOG_ROOT, *, level: int | None = None, console: bool = True, force: bool = False, rotate_mb: int | None = None, log_file_prefix: str = "microseq", max_bytes: int | None = None, backup_count: int = 3) -> Path:
    """
    If $MICROSEQ_LOG_FILE is set -> use that exact path I recommend you do. 
    Else if $MICROSEQ_LOG_DIR is set instead -> microseq.log in that folder.
    Otherwise ./log/microseq.log inside the repo.. 
    
    If the chose file already exists, rename it to:
        microseq.log<YYYYMMDD-HHMMSS>  (rolls-on-startup). 

    If rotate_mb is provided, also adds a RotatingFileHandler that rolls every 1.5 MB (maxBytes = rotate_mb * 1 048 576,
                                                                                     backupCount = 5). 


       Configure a timestamped file + console logger.
       rotate_mb: int None = None so None was used to emulate Nextflow style logging here. 

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

    Returns
    ----------
    pathlib.Path 
    The path of the log file in use. 
    """
    # tranlates the legacy kwargs to current behavior I have it designed for rotate_mb ----
    rotate_bytes: int | None = None 
    if max_bytes is not None: # new tests use this 
        rotate_bytes = int(max_bytes) # bytes as is 
    elif rotate_mb is not None: 
        rotate_bytes = int(rotate_mb * 1024 * 1024) 
    
   # ----- pick a destination path --------------------------------------
    explicit_file = os.getenv("MICROSEQ_LOG_FILE") # this is designed to make sure explicity file always wins but can be overwritten 
    if explicit_file:
        logfile = Path(explicit_file).expanduser() 
        logfile.parent.mkdir(parents=True, exist_ok=True) 
        root_dir = logfile.parent 

    # explicit argument outranks env-var  
    else:
        root_dir = (
            Path(log_dir).expanduser() # arg provided 
            if log_dir is not None 
            else Path(os.getenv("MICROSEQ_LOG_DIR", LOG_ROOT)).expanduser()
            )
        root_dir.mkdir(parents=True, exist_ok=True)


        ts = datetime.now().strftime("%Y%m%d-%H%M%S")
        logfile = root_dir / f"{log_file_prefix}_{ts}.log"


    # ------- roll previous run a la Nextflow naming microseq.log! --------------------- 
    if (
        log_file_prefix == "microseq" 
        and rotate_bytes is None 
        and logfile.exists()
    ): 
        ts = datetime.now().strftime("%Y%m%d-%H%M%S") 
        logfile.rename(logfile.with_suffix(f".log.{ts}")) 



    # --- guard against re-init configuring root logger ----------
    root_logger = logging.getLogger() 
    if root_logger.handlers and not force:
        return logfile 

    root_logger.handlers.clear() 
    root_logger.setLevel(level or logging.INFO) 

    fmt = logging.Formatter(
            "%(asctime)s  %(levelname)-7s  %(name)s:  %(message)s"
            )

    # one or other handler (tests expect a single file handler here)  
    # size based rotation (keeps .log, .log.1, etc.) 
    if rotate_bytes:  
        fh = RotatingFileHandler(
                logfile, maxBytes=rotate_bytes, backupCount=backup_count
                )
    else:
        fh = logging.FileHandler(logfile, mode="w")
    fh.setFormatter(fmt)
    root_logger.addHandler(fh) 
    
    if console:
        ch = logging.StreamHandler(sys.stderr)
        ch.setFormatter(fmt)
        root_logger.addHandler(ch) 

    root_logger.info("Logging to %s", logfile)
    return logfile  
     



