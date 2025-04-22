# microseq_test/src/microseq_test/post_blast_analysis.py 

"""
post_blast_analysis.py 
Take a BLAST TSV + metadata TSV => write one-column BIOM + CSV mirror. 
""" 

from pathlib import Path 
import pandas as pd 
import numpy as np 
from biom import Table
from biom.util import biom_open 
from microseq_tests.utility.utils import setup_logging # reusing logger 
from microseq_tests.utility.utils import load_config # for default DB paths here 
import logging 

setup_logging() # initialize global logging by configure as root logger  
logger = logging.getLogger(__name__) # Now this then set as the real logger by passing everything from the root logger which doesn't return anything on its own  

# -------
# constants - expose as CLI flags for later on will write them in master file CLI
IDENTITY_TH = 97.0 # % identity threshold using for species-grade hits OTU 
QCOV_TH = 80.0 # used in blast wrapper Q-coverage threshold 

# --------------

def _tax_depth(taxon: str) -> int: 
    """Count filled ranks (d__, p__, ... s__).""" 

    return sum(bool(seg.split("__")[1]) for seg in taxon.split(";"))

def _choose_best_hit(blast_df: pd.DataFrame) -> pd.DataFrame:
    """
    Best hit per sample with robustness: 
        keep only rows meeting species-level identity (≥ 97 %)
        then choose lowest e-value, deepest taxonomy, highest bitscore
        """
    # filter on indentity threshold (qcov already filtered upstream) 
    filt = blast_df["pident"] >= IDENTITY_TH 
    pruned = blast_df.loc[filt].copy() 

    # if a sample lost all hits, fall back to full table so it's not dropped 
    if pruned.empty:
        pruned = blast_df.copy() 

    pruned["tax_depth"] = pruned["taxonomy"].map(_tax_depth) 

    # deterministic ordering 
    pruned = pruned.sort_values(
            ["sample_id", "evalue", "tax_depth", "bitscore"],
            ascending=[True, True, False, False]
            )

    # Keep first row per sample 
    best = (pruned.groupby("sample_id", as_index=False)
            .first()
            .drop(columns="tax_depth") 
            )
    return best 

def biom_to_csv(biom_path: Path) -> Path: 
    """Convert one-column BIOM to Prism friendly wide CSV (samples in rows).""" 
    with biom_open(str(biom_path)) as fh: 
        table = Table.from_hdf5(fh) 

    df = table.to_dataframe(dense=True).T       # rows = samples since 
    out_csv = biom_path.with_suffix(".csv")
    df.to_csv(out_csv)
    logger.info(f"Wrote {out_csv}") 
    return out_csv 

def run(blast_tsv: Path,
        metadata_tsv: Path, 
        out_biom: Path,
        write_csv: bool = True) -> None:

        # ---- more tolerant parser: any whitespace, not just tabs ---------
        blast = pd.read_csv(blast_tsv, sep=r"\s+", engine="python")
        meta  = pd.read_csv(metadata_tsv, sep=r"\s+", engine="python")
        
        # ensure pident is numeric so the ≥97 filter works
        if "pident" in blast.columns:
            blast["pident"] = pd.to_numeric(blast["pident"], errors="coerce") 

        best = _choose_best_hit(blast)
        merged = best.merge(meta, on="sample_id", how="left")
        print(blast.head(), blast.columns) 

        # one count per sample (presence/absence) 
        counts = np.ones(len(merged), dtype=int)
        biom = Table(counts[:, None],
                     observation_ids=merged["taxonomy"],
                     sample_ids=["sample"])    

        with biom_open(str(out_biom), "w") as fh:
            biom.to_hdf5(fh, "MicroSeq")
        logger.info(f"Wrote {out_biom}")

        if write_csv:
            biom_to_csv(out_biom) 
