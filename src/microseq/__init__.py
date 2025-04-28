from importlib import import_module as _imp 
_mod = _imp("microseq_tests")
globals().update(vars(_mod)) # re-export everything 
import sys 
sys.modules[__name__] = _mod # one canonical module object 
