# microseq_test/src/microseq_test/post_blast_analysis.py 

"""
post_blast_analysis.py 
Take a BLAST TSV + metadata TSV/CSV => write one-column BIOM + Prism-Friendly CSV mirror. 
""" 

from __future__ import annotations 
from pathlib import Path 
import csv, logging, pandas as pd, numpy as np 
from biom import Table
from biom.util import biom_open 
from microseq_tests.utility.utils import load_config, setup_logging # for default DB paths here
from microseq_tests.utility.io_utils import normalise_tsv 

setup_logging() # initialize global logging by configure as root logger  
logger = logging.getLogger(__name__) # Now this then set as the real logger by passing everything from the root logger which doesn't return anything on its own  

# -------
# constants - expose as CLI flags for later on will write them in master file CLI
DEFAULT_IDENTITY_TH = 97.0 # % identity threshold using for species-grade hits OTU 

# --- helper function: read table with auto delimiter detecter ---------
def _smart_read(path: Path) -> pd.DataFrame:
    """
    Robust reader that keeps the GG2 lineage intact 
    *.tsv force tab as the delimiter here so taxnomy strings may contain spaces that autosniffing for whatever reason explode them. 
    different file such as csv then fall back to older sniffer logic. 
    """
    if path.suffix.lower() == ".tsv":
        df = pd.read_csv(normalise_tsv(path), sep='\t')
        logger.info(f"Loaded {path.name} as explicit TSV rows={len(df)}")
        return df 

    sample = path.read_text(errors="replace")[:1024]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,;")
        sep = dialect.delimiter
    except csv.Error:
        sep = r"\s+"
    df = pd.read_csv(path, sep=sep, engine="python")
    logger.info(f"Loaded {path.name} with delimiter = '{sep}' rows={len(df)}")
    return df 

# ----- helper function: choose the metadata column matching BLAST sample_IDs ------- 
def _detect_sample_col(meta: pd.DataFrame, 
    blast_ids: set[str],
    preferred: str | None = None) -> str:
    """Return name of metadata column to treat as sample_id.""" 
    
    # explicit --sample-col flag to use 
    if preferred and preferred in meta.columns:
        return preferred 

    # config.yaml default 
    cfg_col = load_config().get("metadata", {}).get("sample_col")
    if cfg_col and cfg_col in meta.columns:
        return cfg_col 

    # canonical name 
    if "sample_id" in meta.columns:
        return "sample_id" 

    # fuzzy overlap > 80% 
    for col in meta.columns:
        overlap = blast_ids & set(meta[col].astype(str)) 
        if len(overlap) / len(blast_ids) >= 0.8:
            return col 

    raise ValueError("No metadata column matches BLAST sample IDs either rename to sample_id or use the --sample-col flag") 


# ---taxonomy depth function here for tie-breaking (e-values) ------------- Note this assumed you blasted against GG2 ONLY given how its parsed out 
def _tax_depth(taxon: str | float) -> int:
    """Return how many ranks are filled in the lineage string."""
    if not isinstance(taxon, str):
        return 0 # NaN or non-string means depth of 0.
    parts = [seg.split("__", 1)[-1] for seg in taxon.split(";")] # strip prefix if present 
    return sum(bool(p.strip()) for p in parts) 



# --- chosing the best hit per sample ------------

def _choose_best_hit(df: pd.DataFrame, *, identity_th: float) -> pd.DataFrame:
    """
    Best hit per sample with robustness:
        keep only rows meeting species‑level identity (≥97 %)
        then choose lowest e‑value, deepest taxonomy, highest bitscore
    """
    # filter on identity threshold (qcov already filtered upstream)
    if "pident" in df.columns:
        df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
        df = df[df["pident"] >= identity_th].copy()

    if df.empty:
        raise ValueError("No BLAST hits survive identity threshold")

    df.loc[:, "tax_depth"] = df["taxonomy"].map(_tax_depth)
    df = df.sort_values(
        ["sample_id", "evalue", "tax_depth", "bitscore"],
        ascending=[True, True, False, False]
    )
    return (
        df.groupby("sample_id", as_index=False)
          .first()
          .drop(columns="tax_depth")
    )




# ---CSV mirror for prism -----------------
def biom_to_csv(biom_path: Path) -> Path: 
    """Convert one-column BIOM to Prism friendly wide CSV (samples in rows).""" 
    with biom_open(str(biom_path)) as fh: 
        table = Table.from_hdf5(fh) 

    df = table.to_dataframe(dense=True).T       # rows = samples  
    out_csv = biom_path.with_suffix(".csv")
    df.to_csv(out_csv)
    logger.info(f"Wrote {out_csv}") 
    return out_csv 

# ----- running file API! ---------------------------------
def run(blast_tsv: Path,
        metadata_tsv: Path, 
        out_biom: Path,
        write_csv: bool = True,
        sample_col: str | None = None,
        identity_th: float = DEFAULT_IDENTITY_TH) -> None:

        # ---- more tolerant parser: any whitespace, not just tabs ---------
        blast = _smart_read(blast_tsv)
        meta = _smart_read(metadata_tsv) 

        # detect/rename sample_id column 
        col = _detect_sample_col(meta, set(blast["sample_id"].astype(str)), preferred=sample_col) 
        if col != "sample_id":
            meta = meta.rename(columns={col: "sample_id"}) 

        best = _choose_best_hit(blast, identity_th=identity_th)
        merged = best.merge(meta, on="sample_id", how="left") 
     
        print(blast.head(), blast.columns) 

        # one count per sample (presence/absence) matrix  
        mat = (
            merged.assign(count=1) # add constant 1 
                  .pivot_table(index="taxonomy", # rows = taxa 
                                  columns="sample_id", # cols = isolates 
                                  values="count",
                                  fill_value=0) # 0/1 matrix 
                  )
        biom_table = Table(
            mat.values,
            observation_ids=mat.index.tolist(),
            sample_ids = mat.columns.tolist()
            )
        with biom_open(str(out_biom), "w") as fh:
            biom_table.to_hdf5(fh, "MicroSeq")
        logger.info(f"Wrote {out_biom} (shape {biom_table.shape})") 

        if write_csv:
            biom_to_csv(out_biom) 
