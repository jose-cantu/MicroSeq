from __future__ import annotations 
import sys, inspect
print("DEBUG python:", sys.executable)
import microseq_tests, importlib.util
print("DEBUG pkg:", inspect.getfile(microseq_tests))
print("DEBUG spec:", importlib.util.find_spec("microseq_tests.utility"))
import os, logging, pathlib, pytest 
from microseq_tests.utility.utils import load_config, setup_logging

def test_load_config(): 
    cfg = load_config("config/config.yaml")
    assert "tools" in cfg 

def test_setup_logging(tmp_path: pathlib.Path):
    """
    tmp_path is a py test fixture that yields a fresh, auto-cleaned path. 
    """
    setup_logging(log_dir=str(tmp_path), log_file_prefix="test", force=True)
    logging.info("hello")

    logging.shutdown() 

    logs = list(tmp_path.glob("test_*.log"))
    assert logs, "no log file created" 

