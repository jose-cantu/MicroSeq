# Re-export the thin-wrapper module so callers can 
# from microseq_tests import pipeline 
from importlib import import_module as _imp 
import importlib, sys 
pipeline = _imp('.pipeline', __name__) # noqa: F401

__all__ = ['pipeline']

sys.modules["microseq"] = importlib.import_module(__name__) 
__version__ = "0.1.0-alpha2" # bump in version will update 
