import logging 
import yaml 
import os 

def load_config(config_path="config.yaml"):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def setup_logging(log_file="output.log"):
    logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(levelname)s: %(message)s",
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler() 
                ]
            )

