import logging # built in python module for logging 
import yaml #pyYAML libary for paring YAML files 
import os # for path/OS operations 

def load_config(config_path="config.yaml"):
    """
    Loads a YAML config file and returns the contents as a dictionary.
    """
    with open(config_path, "r") as f:
        return yaml.safe_load(f)
