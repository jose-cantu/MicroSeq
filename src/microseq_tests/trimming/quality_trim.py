import argparse
import subprocess 
from microseq_tests.utility.utils import load_config, setup_logging 

def quality_trim(input_file, output_file):
    """
    Run Trimmomatic to quality-trim reads.
    """
    config = load_config() # loads 'config.yaml'
    trimmomatic_executable = config["tools"]["trimmomatic"] 

    # Example usuage: using some deafult config-based parameters 
    cmd = [
            trimmomatic_executable, "SE", "-phred33", # quality scores are in PHRED + 33 encoding 
            input_file, output_file,
            "SLIDINGWINDOW:5:20", # sliding window algorithm used 
            "MINLEN:200"  # length bp less than 200 after QC trimmed then toss it 
            ]

    subprocess.run(cmd, check=True) # Check process "if Trimmomatic exits with a non-zero status raise process error!"  

def main():
    setup_logging() # Logging set up to capture all later messages aka Python's root logger (format, level, file handler) 
    parser = argparse.ArgumentParser(description="Quality trim using Trimmomatic")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ")
    parser.add_argument("-o", "--output", required=True, help="Output FASTQ")
    args = parser.parse_args() # validates flags being used here otherwise exits with proper flags to use 

    quality_trim(args.input, args.output)

if __name__ == "__main__":
    main() 
