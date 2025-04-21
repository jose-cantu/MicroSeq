from __future__ import annotations 
from pathlib import Path 
from microseq_tests.utility.utils import load_config, expand_db_path, setup_logging 
import subprocess, logging, argparse

PathLike = str | Path 

# db_key is the shorthand string for "gg2" or "silva" best keep it str here for future reference 
def run_blast(query_fa: PathLike, db_key: str, out_tsv: PathLike) -> None:
    cfg = load_config()
    try: 
        tmpl = cfg["databases"][db_key]["blastdb"]
    except KeyError:
        raise KeyError(
                f"[run_blast] unknown DB key '{db_key}'."
                f"Valid keys: {', '.join(cfg['databases'].keys())}"
                )

    blastdb = expand_db_path(tmpl) 

    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    outfmt = "6 qseqid sseqid pident qlen qcovhsp length evalue bitscore stitle"
    cmd=["blastn",
         "-task", "blastn", # why not? 
         "-query", str(query_fa),   # casting query_fa to str before passing to subprocess 
         "-db", blastdb, 
         "-qcov_hsp_perc", "80",     # here I am requiring ≥ 80% of query to align.... 
         "-out", str(out_tsv),
         "-outfmt", outfmt,
         ]

    logging.info("RUN BLAST:%s", " ".join(cmd))
    subprocess.run(cmd, check=True) 
