# -- src/micrseq_tests/blast/run_blast.py ---------------
from __future__ import annotations 
import logging, os, subprocess, shutil, re, functools, threading, time
try:
    from PySide6.QtCore import QThread
except Exception:  # allow running without Qt
    class QThread:
        @staticmethod
        def currentThread():
            return None
from pathlib import Path 
from typing import Optional, Callable 
from Bio import SeqIO 
from microseq_tests.utility.utils import load_config, expand_db_path
from microseq_tests.utility.progress import _tls # access to parent progress bar 
import pandas as pd
from microseq_tests.utility.id_normaliser import NORMALISERS 

L = logging.getLogger(__name__) 
PathLike = str | Path

FIELD_LIST = [ "qseqid","sseqid","pident","qlen","qcovhsp","length","evalue","bitscore","stitle"]
def header_row() -> str: 
    """Returns the TSV header line for BLAST hits including final newline. """ 
    return "\t".join(FIELD_LIST) + "\n"

# progress helper: tail the temporary TSV once per second
def _progress_tail(tmp_path: Path, total: int, callback):
    """Background reader that counts unique qseqid already written."""
    seen: set[str] = set()
    while not getattr(_progress_tail, "stop", False):
        if tmp_path.exists():
            with tmp_path.open() as fh:
                for ln in fh:
                    if ln.startswith("#") or not ln:
                        continue
                    qseqid = ln.split("\t", 1)[0]
                    seen.add(qseqid)
        if callback:
            pct = int(len(seen) / total * 100)
            callback(min(pct, 99))
        time.sleep(1)

# db_key is the shorthand string for "gg2" or "silva" best keep it str here for future reference 
def run_blast(query_fa: PathLike, db_key: str, out_tsv: PathLike, *,
              search_id: float = 97.0, search_qcov: float = 80.0, report_id: float = 97.0, report_qcov: float = 80.0, max_target_seqs: int = 5, threads: int = 1, on_progress: Optional[Callable[[int], None]] = None, log_missing: PathLike | None = None, clean_titles: bool = False, export_sweeper: bool = False, id_normaliser: str = "strip_suffix" ) -> None:
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
    thr = QThread.currentThread()
    if thr and thr.isInterruptionRequested():
        raise RuntimeError("Cancelled")

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
        total = sum(1 for _ in SeqIO.parse(q, "fastq")) # blast now accepts fastq on the offchance the user wants to use fastq instead of fasta.....  
    if total == 0:
        raise ValueError(
            f"{q} contains no FASTA/FASTQ records - nothing to BLAST" 
            )
    # ------ parent tqdm bar auto-callback -------------------- 
    parent_bar = getattr(_tls, "current", None)
    if on_progress is None and parent_bar is not None: 
        on_progress = lambda pct: parent_bar.update(
            max(0, int(pct * total / 100) - parent_bar.n) 
            ) 

    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)


    # column layout I want here 
    DEFAULT_OUTFMT = f"6 {' '.join(FIELD_LIST)}" 

    # grabbing from config or default to here if nothing in config 
    cfg_blast = cfg.get("blast", {})
    outfmt = cfg_blast.get("outfmt", DEFAULT_OUTFMT) 

    # BLAST will write to a temp file first 
    tmp_out = Path(out_tsv).with_suffix(".blasttmp")   


    cmd=["stdbuf", "-oL", "-eL", # line -buffered to see if bar updates real time 
         "blastn",
         "-task", "blastn", # why not? 
         "-query", str(q),   # casting q to str before passing to subprocess 
         "-db", blastdb,
         "-max_target_seqs", str(max_target_seqs), # keep only the best alignment here HSP that has the best-overall score
         "-perc_identity", str(search_id),
         "-qcov_hsp_perc", str(search_qcov),  # here I am requiring ≥ 80% of query to align.... 
         "-outfmt", outfmt,
         "-out", str(tmp_out), # temporary path will append blast with header 
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
    _progress_tail.stop = False # reset flag
    tailer = threading.Thread(
        target=_progress_tail,
        args=(tmp_out, total, on_progress),
        daemon=True
        )
    tailer.start()

    thr = QThread.currentThread()
    if thr and thr.isInterruptionRequested():
        raise RuntimeError("Cancelled")

    try:
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
            if thr and thr.isInterruptionRequested():
                proc.terminate()
                proc.wait()
                raise RuntimeError("Cancelled")
            if ln and not ln.startswith("#"): # data row, not BLAST banner
                done += 1                         # one query finished
                pct = int(done / total * 100)

                if on_progress:                   # GUI / tqdm callback 
                    on_progress(min(pct, 99))    # keep 100% for the end 

                if pct >= next_log: # terminal log 
                    L.info("Progress %d %%", pct)
                    next_log += 5

        rc = proc.wait()
        if thr and thr.isInterruptionRequested():
            proc.terminate()
            raise RuntimeError("Cancelled")

        if rc:
            raise RuntimeError(f"blastn exited with code {rc}") 
    
    finally:
        _progress_tail.stop = True # tell loop to exit 
        tailer.join() # wait for it to finish 


    if on_progress:
        on_progress(100) # final tick in progress bar GUI/tqdm callback 

    # ---- tic outer stage bar exactly once --------------------
    parent = getattr(_tls, "current", None)
    if parent:
        parent.update(1)


    # post-processing adding header here --------------
    header_needed  = outfmt.startswith("6 ") 
    temp_empty = tmp_out.stat().st_size == 0 

    # branch here ----- 
    if header_needed: 
        # Decide your empty-fily policy here so 
        if temp_empty:
            # zero hits -> means header-only file 
            Path(out_tsv).write_text(header_row())  
        else:
            # hits present head is appended then data is transfered over to final results file 
            with open(out_tsv, "w") as final_fh, tmp_out.open("r") as blast_fh:
                final_fh.write(header_row()) 
                shutil.copyfileobj(blast_fh, final_fh) # copies over blast results after header row is appeded first  
    else: 
        # other formats (such as 7) no custom header 
        shutil.move(tmp_out, out_tsv) 

    L.info("BLAST finished OK -> %s", out_tsv)
    if thr and thr.isInterruptionRequested():
        raise RuntimeError("Cancelled")

    # after BLAST call finishes this here will help in cleaning 
    if clean_titles:
        import re 
        # keep only Genus-Species (dropping sseqid and hitlength information from database attached to name of ID that was submitted) 
        df = pd.read_csv(out_tsv, sep="\t", names=FIELD_LIST, dtype=str)
        # strip any # banner / padding and create sample_id 
        df = df.rename(columns=lambda c: c.lstrip("# ").strip())
        # derive sample_id (example rule: drop well / driection suffix) 

        df["stitle"] = (
            df["stitle"]
            .str.split("|").str[-1] # drop gi|...|ref|.... 
            .str.lstrip(">")  # stray '>' and whitespace assuming fastq or stitle is used here.......  
            .str.strip() # this is what removes the whitespace  
            .str.replace(r"^[A-Z]{1,4}\d{5,8}(?:\.\d+){0,2}\s+", "", regex=True) # SILVA "JN193283" prefix
            )
        df.to_csv(out_tsv, sep="\t", index=False)  

    # final cleanup 
    tmp_out.unlink(missing_ok=True) 

    # optional no hit logging setup ------------------ 
    if log_missing:
        import numpy as np
        all_hits = out_tsv 
        
        # ---------- read full BLAST table, clean column whitespace --------
        hits_df = (pd.read_csv(all_hits, sep="\t", dtype=str)
                   .rename(columns=lambda c: c.lstrip("# ").strip())) 

        # add sample_id once, using the SAME normaliser here as postblast consistent       
        norm = NORMALISERS[id_normaliser]
        if "sample_id" not in hits_df.columns:
            hits_df["sample_id"] = hits_df["qseqid"].map(norm) 

        # optional sweeper file 
        if export_sweeper: 
            sweeper_path = Path(log_missing).with_name("hits_full_sweeper.tsv") 
            hits_df.to_csv(sweeper_path, sep="\t", index=False) 
            L.info("sweeper table -> %s", sweeper_path) 

        # collapse best hit per qseqid 
        hits_ok = (hits_df[["qseqid", "sample_id", 
                            "pident", "qcovhsp",
                            "length", "bitscore"]] 
                   .groupby("qseqid").first()) 


        hits_ok.columns = ["sample_id","best_pident", "best_qcov", "align_len", "bitscore"]
        hits_ok["best_pident"] = pd.to_numeric(hits_ok["best_pident"], errors="coerce")
        hits_ok["best_qcov"] = pd.to_numeric(hits_ok["best_qcov"], errors="coerce") 

        all_ids = {rec.id for rec in SeqIO.parse(query_fa, "fasta")}
        full = pd.DataFrame(index=sorted(all_ids))
        full["status"] = np.where(full.index.isin(hits_ok.index), "PASS", "FAIL")
        full = full.join(hits_ok)
        
        full["need_id"] = (
            (report_id - full["best_pident"].fillna(0)).clip(lower=0).round(2)
        )
        full["need_cov"] = (
            (report_qcov - full["best_qcov"].fillna(0)).clip(lower=0).astype("Int64")) 
        fail_id  = full["best_pident"] < report_id
        fail_cov = full["best_qcov"]   < report_qcov
        both     = fail_id & fail_cov
        no_alignment = full["best_pident"].isna() 
        full.loc[no_alignment, "reason"] = "no_alignment" 
        full.loc[ both, "reason"] = "both"
        full.loc[ fail_id & ~both, "reason"] = "low_identity"
        full.loc[ fail_cov & ~both, "reason"] = "low_qcov"
        full.loc[fail_id | fail_cov | no_alignment, "status"] = "FAIL"
        full["reason"] = full["reason"].fillna("—")

        base      = Path(log_missing)
        full_path = base.with_name("hits_full.tsv")
        hits_path = base.with_name("hits.tsv")

        cols = ["status", "reason", "need_id", "need_cov", "sample_id", "best_pident", "best_qcov", "align_len", "bitscore"]
        (full.reset_index()
             .rename(columns={"index": "qseqid"})[["qseqid"] + cols]
        ).to_csv(full_path, sep="\t", index=False)

        pass_ids = full.query("status == 'PASS'").index.astype(str)
        (pd.read_csv(all_hits, sep="\t", dtype=str).rename(columns=lambda c: c.lstrip("# ").strip()).query("qseqid in @pass_ids")).to_csv(hits_path, sep="\t", index=False)

        Path(log_missing).parent.mkdir(parents=True, exist_ok=True)
        missing_ids = full.query("status == 'FAIL'").index.tolist()
        Path(log_missing).write_text("\n".join(missing_ids) +
                                     ("\n" if missing_ids else ""))

        L.info("PASS %d | FAIL %d → %s  (full=%s , hits=%s)",
               len(pass_ids), len(missing_ids),
               log_missing, full_path, hits_path)      
    



