from __future__ import annotations 
import yaml 
import os
from pathlib import Path 
import logging, datetime 

ROOT = Path(__file__).resolve().parents[3]  # repo root (think going 3 steps back here ../../.. until it sees top level of repo =) )  
CONF_PATH = ROOT / "config" / "config.yaml"
def load_config(config_path: str | Path = CONF_PATH):
    config_path = Path(config_path)
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

LOG_ROOT = ROOT/ "logs" 


def setup_logging(log_dir: str | Path = "logs", log_file_prefix: str ="microseq", 
    *,
                  force: bool = False,) -> None:
    """
       Configure a timestamped file + console logger.

    Parameters
    log_dir : str | pathlib.Path
        Directory to write log files (created if missing).
    log_file_prefix : str
        File name prefix, e.g. 'microseq_YYYYmmdd_HHMMSS.log'.
    force : bool, default False
        If True, re‑install handlers even when root already has some
        (useful in pytest where a capture handler is pre‑installed).
    """

    Path(log_dir).mkdir(exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    file_path = Path(log_dir) / f"{log_file_prefix}_{timestamp}.log" 

    if logging.getLogger().handlers and not force:
        return
    logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(levelname)s: %(message)s",
            handlers=[
                logging.FileHandler(file_path),
                logging.StreamHandler()
                ],
            force=force, 
            )
    logging.info(f"Logging to {file_path}") 
           
def expand_db_path(template: str) -> str:
    db_home = os.environ.get("MICROSEQ_DB_HOME", os.path.expanduser("~/.microseq_dbs"))
    return template.replace("${MICROSEQ_DB_HOME}", db_home) 
