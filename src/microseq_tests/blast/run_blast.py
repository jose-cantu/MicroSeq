from microseq_tests.utility.utils import load_config, expand_db_path, setup_logging 
import subprocess, logging, argparse, pathlib

def run_blast(query_fa: str, db_key: str, out_tsv: str):
    cfg = load_config()
    try: 
        tmpl = cfg["databases"][db_key]["blastdb"]
    except KeyError:
        raise KeyError(f"Unknown DB key '{db_key}'")
    blastdb = expand_db_path(tmpl) 

    pathlib.Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    cmd=["blastn", "-query", query_fa, "-db", blastdb, 
         "-out", out_tsv,
         "-outfmt", "6 qseqid sseqid pident length evalue bitscore stitle"]
    logging.info("RUN BLAST:%s", "".join(cmd))
    subprocess.run(cmd, check=True) 
