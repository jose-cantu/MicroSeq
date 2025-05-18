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
from microseq_tests.utility.id_normaliser import NORMALISERS
from microseq_tests.utility.metadata_tools import resolve_duplicates
from microseq_tests.utility.taxonomy_utils import parse_lineage

try:
    from microseq_tests.utility.add_taxonomy import embed_taxonomy_from_metadata
except ImportError:
    # Fallback: leave BIOM table unchanged
    def embed_taxonomy_from_metadata(tbl, *_args, **_kw):
        return tbl

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
def _choose_best_hit(
        df: pd.DataFrame,
        *,
        identity_th: float,
        taxonomy_col: str | None = None
) -> pd.DataFrame:
    """
    Best hit per sample:
      keep only rows with pident ≥ identity_th
      pick lowest e-value
      if a taxonomy column is present, break ties by deepest lineage
      finally break ties by highest bitscore
    """
    # ── identity filter ───────────────────────────────────────────────
    if "pident" in df.columns:
        df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
        df = df[df["pident"] >= identity_th].copy()

    if df.empty:
        raise ValueError("No BLAST hits survive identity threshold")

    # ── optional taxonomy depth for tie-breaking ─────────────────────
    if taxonomy_col and taxonomy_col in df.columns:
        df["tax_depth"] = df[taxonomy_col].map(_tax_depth)
        sort_keys  = ["sample_id", "evalue", "tax_depth", "bitscore"]
        asc_flags  = [ True,       True,     False,       False     ]
    else:
        sort_keys  = ["sample_id", "evalue", "bitscore"]
        asc_flags  = [ True,       True,     False       ]

    # ── sort & collapse ───────────────────────────────────────────────
    df = df.sort_values(sort_keys, ascending=asc_flags)
    return (
        df.groupby("sample_id", as_index=False)
          .first()
          .drop(columns=[c for c in ("tax_depth",) if c in df.columns])
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
        identity_th: float = DEFAULT_IDENTITY_TH,
        *, # again force keyword args used avoids accidnetal position mistake args
        id_normaliser: str = "none",
        taxonomy_col: str = "auto",
        taxonomy_format: str = "auto",
        duplicate_policy: str = "error",
        **kw) -> None:

        # ---- more tolerant parser: any whitespace, not just tabs ---------
        blast = _smart_read(blast_tsv)
        meta = _smart_read(metadata_tsv)
        
        # ---- id-normaliser config -> pull rule from YAML 
        if id_normaliser == "config":
            cfg_norm = load_config().get("metadata", {}).get("sample_id_normaliser")
            if cfg_norm:
                logger.info("Using sample_id normaliser from config.yaml: %s", cfg_norm) 
                id_normaliser = cfg_norm
            else:
                logger.warning("--id-normaliser config requested, but "
                       "metadata.sample_id_normaliser not found in YAML; "
                       "falling back to 'none'")
        # ---------- detect/rename firest, normalise second here ........ 
        # Figure out which metadata column actually contains the sample IDs
        col = _detect_sample_col(meta, set(blast["sample_id"].astype(str)),
                                 preferred=sample_col)
        
        # stach unmodified raw id name here before normalization 
        meta["RawID"] = meta[col].astype(str) 

        if col != "sample_id":
            # if a messy sample_id column already exists, drop it first
            if "sample_id" in meta.columns:
                meta = meta.drop(columns=["sample_id"])
            meta = meta.rename(columns={col: "sample_id"})

        # --- aplying chosen ID normaliser --------------
        norm = NORMALISERS[id_normaliser]
        meta["sample_id"] = meta["sample_id"].astype(str).map(norm) 
        blast["sample_id"] = blast["sample_id"].astype(str).map(norm) 

        meta = resolve_duplicates(meta, policy=duplicate_policy, col="sample_id") 
            

        best = _choose_best_hit(blast, identity_th=identity_th, taxonomy_col=taxonomy_col)
        merged = best.merge(meta, on="sample_id", how="left")
        # after the merge ― immediately normalise the duplicate columns
        dup_suffixes = (".x", ".y")

        for base in ["taxonomy", "sseqid", "pident", "qlen",
                     "qcovhsp", "length", "evalue", "bitscore", "stitle"]:
            # find any versions of this base name, e.g. taxonomy_x / taxonomy_y
            matches = [c for c in merged.columns if c.split(".", 1)[0] == base]
            if matches:
                # prefer the *first* (left-hand) copy, drop the rest
                keep = matches[0]
                merged = merged.rename(columns={keep: base}).drop(columns=matches[1:])

        # also drop any extra SampleID copies created by the merge
        merged = merged.loc[:, ~merged.columns.str.startswith(("sample_id.", "SampleID."))]

        # created a separate blast file for providance 
        prov_out = out_biom.with_name(out_biom.stem + "_blast_provenance.csv")
        best.to_csv(prov_out, index=False)
        logger.info("BLAST provenance → %s", prov_out)

        # remove accidental duplicate-label columns created by the merge
        dup_cols = [c for c in merged.columns
                    if c.startswith("sample_id.") or c.endswith(".1")]
        if dup_cols:
            merged = merged.drop(columns=dup_cols)
        # --- resolve taxonomy column ------------------------
        if taxonomy_col == "auto":
            def looks_like_tax(col: pd.Series) -> bool: 
                return (col.astype(str)
                           .str.count(r"(d__|p__|c__|o__|f__|g__|s__)") >= 4).mean() > 0.9
            taxonomy_col = next(
                (c for c in merged.columns if looks_like_tax(merged[c])),
                None)
            if taxonomy_col is None:
                raise ValueError(
                    "Could not auto-detect taxonomy column; "
                    "use --taxonomy_col") 
            logger.info("Auto-detected taxonomy column: %s", taxonomy_col) 

            

        # Ensure the column we’ll pivot on is called exactly 'taxonomy'
        # (handles taxonomy_x / taxonomy_y after the merge)
        if taxonomy_col != "taxonomy" and taxonomy_col in merged.columns:
            merged = merged.rename(columns={taxonomy_col: "taxonomy"})
            taxonomy_col = "taxonomy"


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
            sample_ids = mat.columns.tolist(),
            type_ = "OTU table", 
            )
        # attach taxonomy list so ATIMA shows ranks 
        biom_table = embed_taxonomy_from_metadata(
            biom_table,
            merged, # merged has taxonomy_col present
            col=taxonomy_col,
            fmt=taxonomy_format,
            )

        with biom_open(str(out_biom), "w") as fh:
            biom_table.to_hdf5(fh, "MicroSeq")
        logger.info(f"Wrote {out_biom} (shape {biom_table.shape})")
        
        obs_ids = biom_table.ids(axis="observation")
        tax_df  = pd.DataFrame(
            [parse_lineage(oid, taxonomy_format) for oid in obs_ids],
            index=obs_ids,
            columns=["domain", "phylum", "class", "order", "family", "genus", "species"],
        )
        tax_df.index.name = "OTU_ID"
        tax_df.reset_index(inplace=True)

        tax_csv = out_biom.with_name(out_biom.stem + "_taxonomy_only.csv")
        tax_df.to_csv(tax_csv, index=False)
        logger.info("Wrote taxonomy-only CSV → %s", tax_csv)


        if write_csv:
            # rename first, so we’re always working with 'SampleID'
            meta_out = out_biom.with_name(out_biom.stem + "_metadata.csv") 
            meta_df = merged.rename(columns={"sample_id": "SampleID"})

            # blast columns you never want to treat as metadata
            blast_cols = ["sseqid", "pident", "qlen", "qcovhsp", "length",
                          "evalue", "bitscore", "stitle", "taxonomy"]

            # TRUE metadata = everything that is NOT a BLAST column and NOT SampleID
            meta_cols = [c for c in meta_df.columns
                         if c not in blast_cols and c != "SampleID"]

            # Build the final ordered list once, then make it unique
            ordered = ["SampleID"] + meta_cols 
                

            meta_df = meta_df.loc[:, ordered]
            meta_df.to_csv(meta_out, index=False)
            logger.info("ATIMA metadata → %s", meta_out)
            biom_to_csv(out_biom) # i'll keep the prism mirror here for now go prism users 
