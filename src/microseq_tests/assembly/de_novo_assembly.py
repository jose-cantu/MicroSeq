"""
de_novo_assembly.py

Wraps CAP3 for simple single-file assemblies. 

Usage from CLI module:
    de_novo_assembly("trimmed_reads.fasta", "asm/")
    
The CAP3 binary path is read from config.yaml.
"""

from pathlib import Path
import subprocess 
import logging 
from microseq_tests.utility.utils import load_config, setup_logging 

def de_novo_assembly(input_fasta: str, output_dir: str) -> None:
    """
    Run CAP3 and place results in output_dir
    """
    cfg = load_config()
    cap3_exe = cfg["tools"]["cap3"] # keeps executable path outisde code; edits once in config.yaml 
    
    in_path = Path(input_fasta).resolve()
    out_dir = Path(output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # CAP3 writes files in CWD wiht the same base name, so change CWD to out_dir and pass only basename. 
    cmd = [cap3_exe, str(in_path)]

    logging.info("RUN CAP3 -> %s", out_dir)
    subprocess.run(cmd, check=True, cwd=out_dir) #CAP3 always writes basename.cap.* in CWD; helps keep project clean 

    logging.info("CAP3 finished' contigs are in %s.cap.contigs", in_path.name) 


# optional CLI to call file directly 
if __name__ == "__main__":
    import argparse, sys
    setup_logging()
    ap = argparse.ArgumentParser(description="De-novo assembly via CAP3")
    ap.add_argument("-i", "--input", required=True, help="Trimmed FASTA here!")
    ap.add_argument("-o", "--output", required=True, help="Output directory") 
    args = ap.parse_args() 
    try:
        de_novo_assembly(args.input, args.output)
    except Exception as e:
        logging.error(e)
        sys.exit(1) 


