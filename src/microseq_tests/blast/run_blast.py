# -- src/micrseq_tests/blast/run_blast.py ---------------
from __future__ import annotations 
import logging, os, subprocess 
from pathlib import Path 
from typing import Optional, Callable 
from Bio import SeqIO 

from microseq_tests.utility.utils import load_config, expand_db_path 

L = logging.getLogger(__name__) 
PathLike = str | Path 

# db_key is the shorthand string for "gg2" or "silva" best keep it str here for future reference 
def run_blast(query_fa: PathLike, db_key: str, out_tsv: PathLike,
              pct_id: float = 97.0, qcov: float = 80.0, max_target_seqs: int = 5, threads: int = 1, on_progress: Optional[Callable[[int], None]] = None, log_missing: PathLike | None = None,) -> None:
    """
    Run blastn against one of the configured 16 S databases.
    You will also emit a percentage progress bar I have set to make it look more a      ppealing. =) 

    Parameters:
    query_fa : PathLike
        FASTA file to search (single sample or contigs).
    db_key : str
        Key in `config.yaml` under ``databases:`` (e.g. "gg2", "silva").
    out_tsv : PathLike
        Destination TSV file (BLAST 6 columns + extras listed below).
    Thresholds applied by BLAST before writing to out_tsv. 
    log_missing 
        if given, appends isolate-ID to this file when resulting TSV contains ≤ 1 line header only -> zero hits ≥ pct_id / qcov).
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

    # compute total queries here 

    total = sum(1 for _ in SeqIO.parse(q, "fasta"))
    if total == 0:    # input is FASTQ 
        total = sum(1 for _ in SeqIO.parse(q, "fastq")) 
    if total == 0:
        raise ValueError(
            f"{q} contains no FASTA/FASTQ records - nothing to BLAST" 
            )

    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)

    outfmt = "6 qseqid sseqid pident qlen qcovhsp length evalue bitscore stitle"
    cmd=["blastn",
         "-task", "blastn", # why not? 
         "-query", str(q),   # casting q to str before passing to subprocess 
         "-db", blastdb,
         "-max_target_seqs", str(max_target_seqs), # keep only the best alignment here HSP that has the best-overall score
         "-perc_identity", str(pct_id),
         "-qcov_hsp_perc", str(qcov),  # here I am requiring ≥ 80% of query to align.... 
         "-out", str(out_tsv),
         "-outfmt", outfmt,
         "-num_threads", str(threads),
         ]

    # merge env here so BLASTDB_LMDB_MAP_SIZE is kept 
    env = os.environ.copy() # start from full parent env 

    # force-set map-size 
    env.setdefault("BLASTDB_LMDB", "0")

    # remove any pre existing map-size so LMDB can't re-activate 
    env.pop("BLASTDB_LMDB_MAP_SIZE", None) 

    L.info("RUN BLAST:%s", " ".join(cmd))

    # ------------ progress setup bar callback ------------------------ 
  
    if on_progress:
        on_progress(0) 

    proc = subprocess.Popen(
            cmd, 
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            env=env,
            )

    done, next_log = 0, 5 # log every 5% incrementally ......

    for ln in proc.stdout:
        if ln and not ln.startswith("#"): # data row, not BLAST banner 
            done += 1                         # one query finished 
            pct = int(done / total * 100) 

            if on_progress:                   # GUI / tqdm callback 
                on_progress(min(pct, 99))    # keep 100% for the end 

            if pct >= next_log: # terminal log 
                L.info("Progress %d %%", pct)
                next_log += 5

    rc = proc.wait()
    if rc:
        raise RuntimeError(f"blastn exited with code {rc}") 

    if on_progress:
        on_progress(100) # final tick in progress bar 
    L.info("BLAST finished OK -> %s", out_tsv)

    # optional no hit logging setup ------------------ 

    if log_missing:
        n_lines = 0 
        try:
            n_lines = sum(1 for _ in Path(out_tsv).open()) 
        except FileNotFoundError:
            pass 

        if n_lines <= 1:  # header-only TSV 
            iso = Path(query_fa).stem
            L.debug(
                    "[post-filter] isolate %s -> 0 hits ≥%.1f%% id / ≥%.1f%% qcov",
                    iso, pct_id, qcov, ) 
            with open(log_missing, "a") as fh: 
                fh.write(f"{iso}\n") 




