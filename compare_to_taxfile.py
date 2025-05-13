#!/usr/bin/env python3
"""
compare_to_taxfile.py
Check where SILVA’s curated lineage stops at genus while the BLAST title
(stitle) appears to specify a species.

usage:
    python compare_to_taxfile.py  silva.tax.tsv  ~/.microseq_dbs/silva/taxonomy.tsv
"""
from __future__ import annotations
import sys, re
from pathlib import Path
import pandas as pd

if len(sys.argv) != 3:
    sys.exit("usage: compare_to_taxfile.py <silva.tax.tsv> <taxonomy.tsv>")

taxed_path   = Path(sys.argv[1]).expanduser().resolve()   # output of add_taxonomy
raw_tax_path = Path(sys.argv[2]).expanduser().resolve()   # SILVA taxonomy.tsv

# ── helpers ────────────────────────────────────────────────────────────
def last_rank(line):
    if not isinstance(line, str):
        return None
    parts = [p.strip() for p in line.split(";") if p.strip()]
    return parts[-1].split("__", 1)[-1] if parts else None       # strip d__/p__ etc.

species_pat = re.compile(r"\b([A-Z][a-zA-Z]+[_\s][a-z][-a-z0-9]*)\b")
def species_from_title(line):
    if not isinstance(line, str):
        return None
    m = species_pat.search(line)
    return m.group(1).replace(" ", "_") if m else None

# ── 1. load both tables ────────────────────────────────────────────────
df  = pd.read_csv(taxed_path, sep="\t", dtype=str)        # 10-column TSV from add_taxonomy
raw = pd.read_csv(raw_tax_path, sep="\t", dtype=str)      # SILVA lineage

# normalise header names
raw.columns = raw.columns.str.strip().str.lower()
id_col  = next(col for col in raw.columns if col in {"feature id", "feature_id", "sseqid"})
tax_col = next(col for col in raw.columns if col in {"taxon", "taxonomy"})

raw = raw.rename(columns={id_col: "sseqid", tax_col: "tax_raw"})

# ── 2. annotate species from both sources ─────────────────────────────
df["tax_species"] = df["taxonomy"].map(last_rank)
df["hit_species"] = df["stitle"].map(species_from_title)

df = df.merge(raw[["sseqid", "tax_raw"]], on="sseqid", how="left")
df["raw_species"] = df["tax_raw"].map(last_rank)

# ── 3. find genus-only lineages where title claims species ────────────
mask = (
    df["hit_species"].notna() &
    df["raw_species"].notna() &
    (df["raw_species"].str.count("_") < 1) &   # genus-only in taxonomy.tsv
    (df["hit_species"] != df["raw_species"])
)

subset = df.loc[mask, [
    "sample_id", "sseqid", "pident",
    "hit_species", "raw_species", "tax_raw", "stitle"
]]

# ── 4. report ──────────────────────────────────────────────────────────
print(f"{len(subset)}/{len(df)} accessions have genus-only lineage in "
      f"{raw_tax_path.name} while the BLAST title carries a species")

if not subset.empty:
    out = taxed_path.with_suffix(".title_vs_taxfile.tsv")
    subset.to_csv(out, sep="\t", index=False)
    print(f"wrote detailed list → {out}")

