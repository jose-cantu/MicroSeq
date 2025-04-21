# microseq_test/src/microseq_test/post_blast_analysis.py 

"""
post_blast_analysis.py 
Take a BLAST TSV + metadata TSV => write one-column BIOM + CSV mirror. 
""" 

from pathlib import Path 
import pandas as pd 
import numpy as np 
from biom import Table 
from microseq_tests.utility.utils import setup_logging # reusing logger 
from microseq_tests.utility.utils import load_config # for default DB paths here 

logger = setup_logging(__name__)

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
    table = Table.from_hdf5(str(biom_path)) 
    df = table.to_dataframe(dense=True).T       # rows = samples since 
    out_csv = biom_path.with_suffix(".csv")
    df.to_csv(out_csv)
    logger.info(f"Wrote {out_csv}") 
    return out_csv 

def run(blast_tsv: Path,
        metadata_tsv: Path, 
        out_biom: Path,
        write_csv: bool = True) -> None: 
        blast = pd.read_csv(blast_tsv, sep="\t")
        meta  = pd.read_csv(metadata_tsv, sep="\t")

        best = _choose_best_hit(blast)
        merged = best.merge(meta, on="sample_id", how="left")

        # one count per sample (presence/absence) 
        counts = np.ones(len(merged), dtype=int)
        biom = Table(counts[:, None],
                     observation_ids=merged["taxonomy"],
                     sample_ids=merged["sample_id"])

        biom.to_hdf5(str(out_biom), "MicroSeq")
        logger.info(f"Wrote {out_biom}")

        if write_csv:
            biom_to_csv(out_biom) 
