import unittest
import os 
from utils import load_config, setup_logging 
import logging 

class TestUtils(unittest.TestCase):
    def test_load_config(self):
        config = load_config("config.yaml")
        self.assertIn("tools", config) # Just an example check to see if section or data exists in the YAML file quick and easy. If tools doesn't exit the test fails. 

    def test_setup_logging(self):
        test_log_file = "test.log"
        if os.path.exists(test_log_file):
            os.remove(test_log_file)
        setup_logging(test_log_file)
        logging.info("Test log message")
        self.assertTrue(os.path.exists(test_log_file)) # Checks that "test.log" now exists after logging message if the file doesn't exist the test fails

if __name__ == "__main__":
    unittest.main() 


