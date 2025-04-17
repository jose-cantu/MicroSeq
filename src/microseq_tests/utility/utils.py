import logging 
import yaml 
import os
from pathlib import Path 
import logging, datetime 

def load_config(config_path="config.yaml"):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def setup_logging(log_dir="logs", log_file_prefix="microseq"):
    Path(log_dir).mkdir(exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    file_path = Path(log_dir) / f"{log_file_prefix}_{timestamp}.log" 

    logging.basicConfig(
            level=logging.DEBUG,
            format="%(asctime)s %(levelname)s: %(message)s",
            handlers=[
                logging.FileHandler(file_path),
                logging.StreamHandler()
                ]
            )
    logging.info(f"Logging to {file_path}") 
            
