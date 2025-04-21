#!/usr/bin/env python3 
"""
fetch_dbs.py 
@author José Cantú

Douwnload/ unpack reference databases for MicroSeq:
    - Greengenes2 (2022.10 backbone full-length) "gg2"
    - SILVA 138.1 SSU NR99 "silva" 
    - NCBI 16S_ribosomal_RNA "ncbi16s" 

Files are stored under $MICROSEQ_DB_HOME (default is ~/.microseq_dbs). 
Each DB direcotry ends up with BLAST-ready indices (makeblastdb -dbtype nucl). 

Run this once, then point config/config.yaml at the same MICROSEQ_DB_HOME.
"""

import os, pathlib, subprocess, urllib.request, tarfile, sys, zipfile, shutil, gzip  

DB_HOME = pathlib.Path(os.environ.get("MICROSEQ_DB_HOME", "~/.microseq_dbs")).expanduser()
DB_HOME.mkdir(parents=True, exist_ok=True) # path check 

def log(msg: str) -> None:
    print(f"[fetch_dbs] {msg}")
    
    
def dl(url: str, dest: pathlib.Path) -> None:
    """Download URL only if dest does not yet exist."""
    if dest.exists():
        log(f"good {dest.name} already present")
        return
    log(f"→ downloading {url}")
    urllib.request.urlretrieve(url, dest)
    
    
def run(cmd: list[str]) -> None:
    """Run shell command (show it, fail loud)."""
    log("+" + " ".join(cmd))
    subprocess.run(cmd, check=True)
    
    
def makeblastdb(fasta: pathlib.Path, out_prefix: pathlib.Path) -> None:
    if (out_prefix.with_suffix(".nsq")).exists():
        log(f"BLAST index for {out_prefix.name} already exists")
        return
    run(
        ["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", str(out_prefix)]
    )

def extract_member(zip_path: pathlib.Path, pattern: str, out_path: pathlib.Path):
    """
    Copy the first file in *zip_path* whose name ends‑with *pattern* to *out_path*.
    Works even when QIIME2 wraps files in random UUID folders.
    """
    with zipfile.ZipFile(zip_path) as zf:
        try:
            member = next(n for n in zf.namelist() if n.endswith(pattern))
        except StopIteration:
            raise ValueError(f"{pattern} not found inside {zip_path.name}")
        with zf.open(member) as fin, open(out_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)

    
def fetch_gg2() -> None:
    """ Download Greengenes 2 (2024‑09 backbone full‑length) and build BLAST    db. """
    gg_dir = DB_HOME / "gg2"
    gg_dir.mkdir(exist_ok=True)

    # Use HTTP to dodge the incomplete‑TLS‑chain problem
    BASE  = "http://ftp.microbio.me/greengenes_release/2024.09"
    FILES = {
        "backbone_fasta": "2024.09.backbone.full-length.fna.qza",
        "taxonomy":       "2024.09.backbone.tax.qza",
    }

    # ── download & extract dna‑sequences.fasta ────────────────────────────────
    backbone_qza = gg_dir / FILES["backbone_fasta"]
    dl(f"{BASE}/{FILES['backbone_fasta']}", backbone_qza)

    fasta = gg_dir / "dna-sequences.fasta"
    if not fasta.exists():
        log("→ extracting dna‑sequences.fasta from QZA")
        extract_member(backbone_qza, "dna-sequences.fasta", fasta)


    # ── taxonomy table ────────────────────────────────────────────────────────
    tax_qza = gg_dir / FILES["taxonomy"]
    dl(f"{BASE}/{FILES['taxonomy']}", tax_qza)
    
    taxonomy_tsv = gg_dir / "taxonomy.tsv"
    if not taxonomy_tsv.exists():
        log("→ extracting taxonomy.tsv from QZA")
        extract_member(tax_qza, "taxonomy.tsv", taxonomy_tsv)
    

    # ── build BLAST index (creates .nsq/.nin/.nhr) ────────────────────────────
    makeblastdb(fasta, gg_dir / "greengenes2_db")
# ─────────────────────────────────────────────────────────────────────────────


    
    
# SILVA 138.1 SSU NR99 (truncated fasta)
def fetch_silva() -> None:
    silva_dir = DB_HOME / "silva"
    silva_dir.mkdir(exist_ok=True)
    
    url = (
        "ftp://ftp.arb-silva.de/release_138.1/Exports/"
        "SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz"
    )
    gz = silva_dir / "SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz"
    dl(url, gz)
    
    fasta = silva_dir / "SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta"
    if not fasta.exists():
        log("→ extracting SILVA with Python gzip")
        with gzip.open(gz, "rb") as fin, open(fasta, "wb") as fout:
            shutil.copyfileobj(fin, fout)

    # Build BLAST index  
    makeblastdb(fasta, silva_dir / "silva_db")  
    
    
#  NCBI 16S_ribosomal_RNA  (pre‑built BLAST db tar.gz)
def fetch_ncbi16s() -> None:
    ncbi_dir = DB_HOME / "ncbi"
    ncbi_dir.mkdir(exist_ok=True)
    url = "ftp://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz"
    tar_file = ncbi_dir / "16S_ribosomal_RNA.tar.gz"
    dl(url, tar_file)
    
    # the tar already contains BLAST index files
    if not (ncbi_dir / "16S_ribosomal_RNA.nsq").exists():
        log("→ extracting NCBI 16S tarball")
        tarfile.open(tar_file).extractall(ncbi_dir)
        
        
# Run all fetchers
def main():
    for fn in (fetch_gg2, fetch_silva, fetch_ncbi16s):
        try:
            fn()
        except Exception as e:
            log(f"[ERROR] {fn.__name__} failed: {e}")
            sys.exit(1)
            
    log(f"[DONE] all DBs are in {DB_HOME}")
    
    
if __name__ == "__main__":
    main() 
