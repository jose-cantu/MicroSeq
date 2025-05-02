# ---- src/microseq_tests/utility/utils.py --------------------- 
from __future__ import annotations 
import yaml, logging, os, sys, datetime  
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
    for parent in [here, *here.parents]:
        if (parent / "pyproject.toml").exists() \
                or (parent / "setup.cfg").exists() \
                or (parent / ".git").exists():
                    return parent 
    # installed wheel -> use the package root on site-packages 
    return pkg_resources.files("microseq_tests").resolve () 

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


def setup_logging(*, level: int | None = None, console: bool = True, force: bool = False, rotate_mb: int | None = None,) -> Path:
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
    """
    # ----- pick a destination path --------------------------------------
    explicit = os.getenv("MICROSEQ_LOG_FILE")
    if explicit:
        logfile = Path(explicit).expanduser() 

    else: 
        log_dir = Path(os.getenv("MICROSEQ_LOG_DIR", LOG_ROOT)).expanduser()
        log_dir.mkdir(parents=True, exist_ok=True)
        logfile = log_dir / "microseq.log" 


    # ------- roll previous run a la Nextflow! --------------------- 
    if logfile.exists():
        ts = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
        logfile.rename(logfile.with_suffix(f".log.{ts}"))

    logfile.parent.mkdir(parents=True, exist_ok=True) 

    # --- guard against re-init ----------
    root = logging.getLogger() 
    if root.handlers and not force:
        return logfile 

    root.handlers.clear() 
    root.setLevel(level or logging.INFO) 

    fmt = logging.Formatter(
            "%(asctime)s  %(levelname)-7s  %(name)s:  %(message)s"
            )

    # always have a primary FileHandler (overwrite mode) 
    fh = logging.FileHandler(logfile, mode="w")
    fh.setFormatter(fmt)
    root.addHandler(fh) 

    # size based rotation (keeps .log, .log.1, etc.) 
    if rotate_mb:
        rotate = RotatingFileHandler(
                logfile, maxBytes=rotate_mb * 1024 * 1024, backupCount=3
                )
        rotate.setFormatter(fmt)
        root.addHandler(rotate) 

    if console:
        ch = logging.StreamHandler(sys.stderr)
        ch.setFormatter(fmt)
        root.addHandler(ch) 

    root.info("Logging to %s", logfile)
    return logfile  
     



