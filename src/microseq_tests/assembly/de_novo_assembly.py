"""
microseq_tests/src/microseq_tests/assembly/de_novo_assembly.py 

de_novo_assembly.py

Wraps CAP3 for simple single-file assemblies. 

Usage from CLI module:
    de_novo_assembly("trimmed_reads.fasta", "asm/")
    
The CAP3 binary path is read from config.yaml.
"""
from __future__ import annotations 
import logging
L = logging.getLogger(__name__)

from pathlib import Path
import shutil
import subprocess
from microseq_tests.utility.utils import load_config, setup_logging 

PathLike = str | Path 

def de_novo_assembly(input_fasta: PathLike, output_dir: PathLike, *, threads: int=1, **kwargs, ) -> Path: 
    """
    Run CAP3 on ``input_fasta``.

    CAP3 itself is single-threaded; the ``threads`` argument is
    accepted for API compatibility but currently ignored.

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

    local_fasta = out_dir / in_path.name
    if in_path != local_fasta:
        shutil.copy2(in_path, local_fasta)

    cmd = [cap3_exe, local_fasta.name]
    L.info("RUN CAP3: %s (cwd=%s)", " ".join(cmd), out_dir)

    try:
        subprocess.run(cmd, check=True, cwd=out_dir, stderr=subprocess.PIPE, text=True) 
    except subprocess.CalledProcessError as exc:
        L.error("CAP3 failed (exit %s):\n%s", exc.returncode, exc.stderr)
        raise 

    contig_path = out_dir / f"{local_fasta.name}.cap.contigs"
    if not contig_path.exists():
        raise FileNotFoundError(contig_path)

    L.info("CAP3 finished; contigs file: %s", contig_path) 
    return contig_path 

# optional CLI to call file directly 
if __name__ == "__main__":
    import argparse, sys
    setup_logging()
    ap = argparse.ArgumentParser(description="De-novo assembly via CAP3")
    ap.add_argument("-i", "--input", required=True, help="Trimmed FASTA here!")
    ap.add_argument("-o", "--output", required=True, help="Output directory") 
    ap.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Ignored: CAP3 is single-threaded (flag kept for compatibility)",
    )
    args = ap.parse_args() 
    try:
        de_novo_assembly(args.input, args.output, threads=args.threads)
    except Exception as e:
        L.error(e)
        sys.exit(1) 


