from __future__ import annotations 
from pathlib import Path 
from microseq_tests.utility.utils import load_config, expand_db_path, setup_logging 
import subprocess, logging, argparse, os

PathLike = str | Path 

# db_key is the shorthand string for "gg2" or "silva" best keep it str here for future reference 
def run_blast(query_fa: PathLike, db_key: str, out_tsv: PathLike) -> None:
    """
    Run blastn against one of the configured 16 S databases. 

    Parameters:
    query_fa : PathLike
        FASTA file to search (single sample or contigs).
    db_key : str
        Key in `config.yaml` under ``databases:`` (e.g. "gg2", "silva").
    out_tsv : PathLike
        Destination TSV file (BLAST 6 columns + extras).
    """
    cfg = load_config()
    # look up DB path in config.yaml 
    try: 
        tmpl = cfg["databases"][db_key]["blastdb"]
    except KeyError as e: # bad key 
        valid = ", ".join(cfg["databases"].keys())
        raise KeyError(
                f"[run_blast] unknown DB key '{db_key}'."
                f"Valid keys: {valid}") from e 
                 

    blastdb = expand_db_path(tmpl) # fills in $HOME etc. 

    # safety make sure query even exists before spawning BLAST .... 
    q = Path(query_fa)
    if not q.is_file():
        raise FileNotFoundError(q) 


    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)

    outfmt = "6 qseqid sseqid pident qlen qcovhsp length evalue bitscore stitle"
    cmd=["blastn",
         "-task", "blastn", # why not? 
         "-query", str(q),   # casting q to str before passing to subprocess 
         "-db", blastdb, 
         "-qcov_hsp_perc", "80",     # here I am requiring ≥ 80% of query to align.... 
         "-out", str(out_tsv),
         "-outfmt", outfmt,
         ]

    # merge env here so BLASTDB_LMDB_MAP_SIZE is kept 
    env = os.environ.copy() # start from full parent env 

    # force-set map-size 
    env.setdefault("BLASTDB_LMDB_MAP_SIZE",str(64 * 1024 ** 3))

    logging.info("RUN BLAST:%s", " ".join(cmd))
    subprocess.run(cmd, check=True, env=env) # launches blast, propogate all vars
