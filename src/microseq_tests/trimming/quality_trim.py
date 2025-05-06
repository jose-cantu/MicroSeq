# microseq_tests/src/microseq_tests/trimming/quality_trim.py 

from __future__ import annotations 
import logging
L = logging.getLogger(__name__)

from pathlib import Path 
import argparse, subprocess, logging 
from microseq_tests.utility.utils import load_config, setup_logging 

PathLike = str | Path 

def quality_trim(input_file: PathLike, output_file: PathLike, *, threads: int=1, **kwargs, ) -> Path:
    """
    Run Trimmomatic single-end mode and return the absolute path to trimmed FASTAQ/FASTA so down stream steps (assembly, GUI, tests) can chain without guessing. 
    threads : int 
    Forwarded to Trimmomatics's threads option. 

        """
    cfg = load_config()
    trimm = cfg["tools"]["trimmomatic"]

    cmd = [ 
        trimm, "SE", "-threads", str(threads), "-phred33",
        str(input_file), str(output_file),    # <<< cast once here 
        "SLIDINGWINDOW:5:20",
        "MINLEN:200",
        ]

    L.info("RUN Trimmomatic: %s", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True, stderr=subprocess.PIPE, text=True)
    except subprocess.CalledProcessError as e:
        L.error("Trimmomatic failed (exit %s):\n%s", e.returncode, e.stderr)
        raise 

    out_path = Path(output_file).resolve()
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
