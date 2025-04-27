# Re-export the thin-wrapper module so callers can 
# from microseq_tests import pipeline 
from importlib import import_module as _imp 
pipeline = _imp('.pipeline', __name__) # noqa: F401

__all__ = ['pipeline'] 
