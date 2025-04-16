import unittest
import os 
from utils import load_config, setup_logging
import logging 

class TestUtils(unittest.TestCase):
    def test_load_config(self):
        config = load_config("config.yaml")
        self.assertIn("tools", config) # Just an example to check if tools is in yaml otherwise fails test

    def test_setup_logging(self):
        test_log_file = "test.log"
        if os.path.exists(test_log_file):
            os.remove(test_log_file)
        setup_logging(test_log_file)
        logging.info("Test log message")

        self.assertTrue(os.path.exists(test_log_file))

if __name__ == "__main__":
    unittest.main()
