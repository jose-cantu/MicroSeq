#!/usr/bin/env python3
"""
microseq-setup  (fetch_dbs.py)

Download & register reference databases, write env-hook + config.yaml.

Creates
  ${DB_ROOT}/{gg2,silva,ncbi}/greengenes2_db.*
  ${LOG_DIR}/                      (empty folder – just the place for logs)

Appends to ~/.bashrc (and ­conda activate hook if inside env)

  export MICROSEQ_DB_HOME="…"
  export BLASTDB="$MICROSEQ_DB_HOME/gg2:$MICROSEQ_DB_HOME/silva:$MICROSEQ_DB_HOME/ncbi"
  export BLASTDB_LMDB=0
  export MICROSEQ_LOG_DIR="…"
"""
from __future__ import annotations

import argparse, os, sys, urllib.request, tarfile, zipfile, gzip, \
       shutil, subprocess, yaml
from pathlib import Path
from datetime import datetime

# ───────────────────────────── CLI ──────────────────────────────────────
ap = argparse.ArgumentParser(
        prog        = "microseq-setup",
        description = "Download BLAST DBs and write env-hook / config."
)
ap.add_argument("--db-root", metavar="PATH",
                help="Folder where BLAST databases will live "
                     "(default: ~/.microseq_dbs)")
ap.add_argument("--log-dir", metavar="PATH",
                help="Folder for MicroSeq run logs "
                     "(default: microseq_tests/logs)")
ap.add_argument("--quiet", action="store_true",
                help="Skip interactive prompts")
args = ap.parse_args()

# ───────────────────────── default paths ────────────────────────────────
REPO_ROOT   = Path(__file__).resolve().parents[3]          # microseq_tests/
DEFAULT_LOG = REPO_ROOT / "logs"

def ask(msg: str, default: str) -> Path:
    if args.quiet:
        return Path(default).expanduser()
    ans = input(f"{msg} [{default}] ").strip() or default
    return Path(ans).expanduser()

db_root = Path(
    args.db_root or ask("Where should databases be stored?", "~/.microseq_dbs")
).resolve()

log_dir = Path(
        args.log_dir or ask("Where should logs be stored?", str(DEFAULT_LOG)) 
).resolve()

db_root.mkdir(parents=True, exist_ok=True)
log_dir.mkdir(parents=True, exist_ok=True)

# ───────────────────────── helpers ──────────────────────────────────────
def log(msg: str) -> None:
    print(f"[setup] {msg}")

def dl(url: str, dest: Path) -> None:
    if dest.exists():
        log(f"✓ {dest.name} already present")
        return
    log(f"→ downloading {url}")
    urllib.request.urlretrieve(url, dest)

def run(cmd: list[str]) -> None:
    log("+" + " ".join(cmd))
    subprocess.run(cmd, check=True)

def makeblastdb(fasta: Path, out_prefix: Path) -> None:
    if (out_prefix.with_suffix(".nsq")).exists():
        return
    run(["makeblastdb", "-in", fasta, "-dbtype", "nucl", "-out", out_prefix])

def extract_member(zip_path: Path, pattern: str, out_path: Path) -> None:
    with zipfile.ZipFile(zip_path) as zf:
        member = next(n for n in zf.namelist() if n.endswith(pattern))
        with zf.open(member) as fin, open(out_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)

# ───────────────────── database download functions ──────────────────────
def fetch_gg2() -> None:
    gg = db_root / "gg2"; gg.mkdir(exist_ok=True)
    base = "http://ftp.microbio.me/greengenes_release/2024.09"
    qza  = gg / "2024.09.backbone.full-length.fna.qza"
    dl(f"{base}/2024.09.backbone.full-length.fna.qza", qza)

    fasta = gg / "dna-sequences.fasta"
    if not fasta.exists():
        extract_member(qza, "dna-sequences.fasta", fasta)

    tax_qza = gg / "2024.09.backbone.tax.qza"
    dl(f"{base}/2024.09.backbone.tax.qza", tax_qza)

    taxonomy = gg / "taxonomy.tsv"
    if not taxonomy.exists():
        extract_member(tax_qza, "taxonomy.tsv", taxonomy)

    makeblastdb(fasta, gg / "greengenes2_db")

def fetch_silva() -> None:
    sd = db_root / "silva"; sd.mkdir(exist_ok=True)
    gz = sd / "SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz"
    dl(
        ("ftp://ftp.arb-silva.de/release_138.1/Exports/"
         "SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta.gz"),
        gz
    )
    fasta = sd / gz.stem
    if not fasta.exists():
        log("→ extracting SILVA")
        with gzip.open(gz, "rb") as fin, open(fasta, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    makeblastdb(fasta, sd / "silva_db")

def fetch_ncbi() -> None:
    nd = db_root / "ncbi"; nd.mkdir(exist_ok=True)
    tgz = nd / "16S_ribosomal_RNA.tar.gz"
    dl("ftp://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz", tgz)
    if not (nd / "16S_ribosomal_RNA.nsq").exists():
        log("→ extracting NCBI 16S")
        tarfile.open(tgz).extractall(nd)

# ────────────────────────── main driver ─────────────────────────────────
def main() -> None:
    for fn in (fetch_gg2, fetch_silva, fetch_ncbi):
        fn()
    log(f"✔  All databases downloaded to {db_root}")

    # ---------- env-hook snippet -----------------------------------------
    snippet = f"""
# ── MicroSeq auto-export ───────────────────────────────────────────────
export MICROSEQ_DB_HOME="{db_root}"
export BLASTDB="$MICROSEQ_DB_HOME/gg2:$MICROSEQ_DB_HOME/silva:$MICROSEQ_DB_HOME/ncbi"
export BLASTDB_LMDB=0
export MICROSEQ_LOG_DIR="{log_dir}"
"""

    def append_once(target: Path, text: str) -> None:
        if target.exists() and text.strip() in target.read_text():
            return
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("a") as fh:
            fh.write(text)
        log(f"✓ exports appended to {target}")

    shell_rc = Path.home() / (os.getenv("SHELL", "/bin/bash").split("/")[-1] + "rc")
    append_once(shell_rc, snippet)

    if os.getenv("CONDA_PREFIX"):
        hook = Path(os.getenv("CONDA_PREFIX")) / "etc/conda/activate.d/microseq.sh"
        append_once(hook, snippet)

    # ---------- patch config.yaml ---------------------------------------
    cfg_path = REPO_ROOT / "config" / "config.yaml"
    cfg_path.parent.mkdir(exist_ok=True)
    cfg = yaml.safe_load(cfg_path.read_text()) if cfg_path.exists() else {}
    cfg["logging"] = {"dir": str(log_dir)}
    cfg["databases"] = {
        "gg2":   {"blastdb": "${MICROSEQ_DB_HOME}/gg2/greengenes2_db"},
        "silva": {"blastdb": "${MICROSEQ_DB_HOME}/silva/silva_db"},
        "ncbi":  {"blastdb": "${MICROSEQ_DB_HOME}/ncbi/16S_ribosomal_RNA"}
    }
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    log(f"✓ config.yaml updated at {cfg_path}")

    print("\nDone! Open a new shell (or `source ~/.bashrc`) and test it using the smoke test I have setup in github install:")

if __name__ == "__main__":
    main()

