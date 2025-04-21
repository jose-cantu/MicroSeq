"""
de_novo_assembly.py

Wraps CAP3 for simple single-file assemblies. 

Usage from CLI module:
    de_novo_assembly("trimmed_reads.fasta", "asm/")
    
The CAP3 binary path is read from config.yaml.
"""
from __future__ import annotations 
from pathlib import Path
import subprocess 
import logging 
from microseq_tests.utility.utils import load_config, setup_logging 

PathLike = str | Path 

def de_novo_assembly(input_fasta: PathLike, output_dir: PathLike) -> Path: 
    """
    Run CAP3 on input_fasta. 

    Parameters
    input_fasta : str 
        Trimmed reads in FASTA format here. 
    output_dir: str 
        Folder where CAP3 files are written in. 

    Returns
    pathlib.Path
        Absolute path to ''*.cap.contigs''.
        """
    cfg = load_config()
    cap3_exe = cfg["tools"]["cap3"] # path comes from config file 

    in_path = Path(input_fasta).resolve()
    out_dir = Path(output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    cmd = [cap3_exe, str(in_path)] 
    logging.info("RUN CAP3: %s (cwd=%s)", " ".join(cmd), out_dir)

    try:
        subprocess.run(cmd, check=True, cwd=out_dir, stderr=subprocess.PIPE, text=True) 
    except subprocess.CalledProcessError as exc:
        logging.error("CAP3 failed (exit %s):\n%s", exc.returncode, exc.stderr)
        raise 

    contig_path = out_dir / f"{in_path.name}.cap.contigs"
    if not contig_path.exists():
        raise FileNotFoundError(contig_path)

    logging.info("CAP3 finished; contigs file: %s", contig_path) 
    return contig_path 

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


