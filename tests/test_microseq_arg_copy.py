import unittest
import subprocess

class TestHelloArgs(unittest.TestCase):
    def test_default_name(self):
        result = subprocess.check_output(["python", "hello_microseq_argparse.py"]).decode() 
        self.assertIn("Hello, MicroSeq!", result)
    def test_custom_name(self):
        result = subprocess.check_output(["python", "hello_microseq_argparse.py", "--name", "TestUser!"]).decode()
        self.assertIn("Hello, TestUser!", result)

if __name__ == "__main__":
    unittest.main() 
