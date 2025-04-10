import logging 
import yaml 
import os 

def load_config(config_path="config.yaml"):
    """
    Load a YAML config file and return contents as a dictionary.
    """
    with open(config_path, "r") as f:
        return yaml.safe_load(f) 

def setup_logging(log_file="output.log"):
    """
    Sets up basic logging. Logs go to both a file and console. 
    """ 
    logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(levelname)s: %(message)s",
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler()
                ]
            ) 


