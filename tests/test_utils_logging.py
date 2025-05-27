# tests/test_utils_logging.py 

import logging
import importlib.util, sys
from pathlib import Path
import pytest

# import utility module directly to avoid importing the full package
ROOT = Path(__file__).resolve().parents[1]
UTIL_PATH = ROOT / "src" / "microseq_tests" / "utility" / "utils.py"
spec = importlib.util.spec_from_file_location("mutils", UTIL_PATH)
mutils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mutils)
setup_logging = mutils.setup_logging

def test_setup_logging_handlers(tmp_path):
    log_file = setup_logging(
            log_dir=tmp_path,
            force=True,
            console=False,
            max_bytes=1_000,
            backup_count=1,
            )
    root = logging.getLogger()
    # one file handler only 
    assert len(root.handlers) == 1 
    assert log_file.exists()
    # rollover works 
    root.info("x" * 2_000) # exceed 1 kb 
    root.handlers[0].flush() 
    rotated = log_file.with_suffix(".log.1")
    assert rotated.exists() 
