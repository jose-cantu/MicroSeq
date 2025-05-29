from importlib.metadata import entry_points 
ALIGNER_REGISTRY = {ep.name: ep.load()
                    for ep in entry_points(group="microseq.aligners")} 

def load(name: str):
    return ALIGNER_REGISTRY[name]() # return an instance 


