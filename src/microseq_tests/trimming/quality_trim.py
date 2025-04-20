import argparse, subprocess, pathlib, logging 
from microseq_tests.utility.utils import load_config, setup_logging 

def quality_trim(input_file: str, output_file: str) -> pathlib.Path:
    """
    Run Trimmomatic single-end mode. 

    Returns here: 

    pathlib.Path
        Soo this is the absulute path to the trimmed FASTA/FASTQ file so downstream functions (assembly, GUI, tests) can chain w/o guessing. 
        """
    cfg = load_config()
    trimm = cfg["tools"]["trimmomatic"]

    cmd = [ 
        trimm, "SE", "-phred33",
        input_file, output_file,
        "SLIDINGWINDOW:5:20",
        "MINLEN:200",
        ]

    logging.info("RUN Trimmomatic: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        logging.error("Trimmomatic failed (exit %s):\n%s", e.returncode, e.stderr)
        raise 

    out_path = pathlib.Path(output_file).resolve()
    if not out_path.exists():
        raise FileNotFoundError(out_path)
    return out_path 


def main():
    setup_logging() # Logging set up to capture all later messages aka Python's root logger (format, level, file handler) 
    parser = argparse.ArgumentParser(description="Quality trim using Trimmomatic")
    parser.add_argument("-i", "--input", required=True, help="Input FASTQ")
    parser.add_argument("-o", "--output", required=True, help="Output FASTQ")
    args = parser.parse_args() # validates flags being used here otherwise exits with proper flags to use 

    quality_trim(args.input, args.output)

if __name__ == "__main__":
    main() 
