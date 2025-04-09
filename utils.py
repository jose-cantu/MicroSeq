import logging 
import yaml 
import os 

def load_config(config_path="config.yaml"):
    """
    Loads a YAML config and returns the contents as a dictionary. 
    """ 
    with open(config_path, "r") as f:
        return yaml.safe_load(f) 
