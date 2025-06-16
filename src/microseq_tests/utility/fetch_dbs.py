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
                     "(default: MicroSeq/logs)")
ap.add_argument("--quiet", action="store_true",
                help="Skip interactive prompts")
args = ap.parse_args()

# ───────────────────────── default paths ────────────────────────────────
REPO_ROOT   = Path(__file__).resolve().parents[3]          # microseq_tests/
DEFAULT_LOG = REPO_ROOT / "logs"

def ask_path(msg: str, default: str | Path) -> Path:
    """
    Prompt until the user gives a writable folder.
    Pressing <Enter> accepts the default.
    """
    if args.quiet:                       # non-interactive install
        return Path(default).expanduser()

    while True:
        ans = input(f"{msg} [{default}]: ").strip() or default
        p   = Path(ans).expanduser()
        try:
            p.mkdir(parents=True, exist_ok=True)   # proves we can write there
            return p
        except OSError as e:
            print(f"{e}; try again.")

db_root = Path(
    args.db_root or ask_path("Where should databases be stored?", "~/.microseq_dbs")
).resolve()

log_dir = Path(
        args.log_dir or ask_path("Where should logs be stored?", str(DEFAULT_LOG)) 
).resolve()

db_root.mkdir(parents=True, exist_ok=True)
log_dir.mkdir(parents=True, exist_ok=True)

def ask_pipeline_usage() -> bool:
    """
    Returns True if the user says they will call MicroSeq from a workflow
    manager (Nextflow, Snakemake, SLURM array …).  Non-interactive installs
    default to False.
    """
    if args.quiet:
        return False 
    while True:
        ans = input(
            "Will MicroSeq ever run inside an automated pipeline on "
            "this machine or on the HPC? (y/N) "
        ).strip().lower() or "n" 
        if ans in ("y", "yes"):
            return True 
        if ans in ("n", "no"):
            return False 
        print("Please type y or n") 

pipeline_user = ask_pipeline_usage()       

# ───────────────────────── helpers ──────────────────────────────────────
def log(msg: str) -> None:
    print(f"[setup] {msg}")

def dl(url: str, dest: Path) -> None:
    if dest.exists():
        log(f"✓ {dest.name} already present")
        return
    log(f"→ downloading {url}")
    urllib.request.urlretrieve(url, dest)

def run(cmd: list) -> None: # accept Path or str 
    log("+" + " ".join(map(str, cmd))) # stringify for loggin 
    subprocess.run(list(map(str,cmd)), check=True)

def makeblastdb(fasta: Path, out_prefix: Path) -> None:
    if (out_prefix.with_suffix(".nsq")).exists():
        return
    run(["makeblastdb", "-in", fasta, "-dbtype", "nucl", "-out", out_prefix])

def extract_member(zip_path: Path, pattern: str, out_path: Path) -> None:
    with zipfile.ZipFile(zip_path) as zf:
        member = next(n for n in zf.namelist() if n.endswith(pattern))
        with zf.open(member) as fin, open(out_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)


# ------ Adding TaxonKit helper -------------------------------------
TAXONKIT_DB = Path.home() / ".taxonkit"
def ensure_taxdump() -> None: 
    """
    Download and unpack the four core NCBI tax-dump files if they are missing 
    (names.dmp, nodes.dmp, delnodes.dmp, merged.dmp). 
    """
    TAXONKIT_DB.mkdir(parents=True, exist_ok=True) # creates the folder first if it doesn't exist 

    needed = {"names.dmp", "nodes.dmp", "delnodes.dmp", "merged.dmp"}
    existing = {p.name for p in TAXONKIT_DB.iterdir() if p.is_file()}
    if needed <= existing: 
        return # already present 

    TAXONKIT_DB.mkdir(parents=True, exist_ok=True)
    dump = TAXONKIT_DB / "taxdump.tar.gz" 
    dl("ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", dump)
    
    log("-> extracting NCBI taxdump") 
    with tarfile.open(dump) as tf: 
        for m in tf.getmembers(): 
            if m.name.split("/")[-1] in needed:
                tf.extract(m, TAXONKIT_DB) 

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

    silva_tax = sd / "taxonomy.tsv"
    if not silva_tax.exists():
        gz_url = ("https://ftp.arb-silva.de/release_138.1/Exports/"
            "taxonomy/taxmap_slv_ssu_ref_138.1.txt.gz")
        tax_gz = sd / "tax_slv_ref.gz" 
        dl(gz_url, tax_gz) 

        log("-> extracting SILVA taxonomy")
        import csv 
        with gzip.open(tax_gz, "rt") as fin, silva_tax.open("w") as fout: 
            w = csv.writer(fout, delimiter="\t", lineterminator="\n")
            w.writerow(["sseqid", "taxonomy"]) 
            for i, ln in enumerate(fin):
                if i == 0: # first record = helper -> skip this 
                    continue 
                fields = ln.rstrip("\n").split("\t")
                acc = fields[0] # primary Acession column I need 
                lineage = fields[3].rstrip(";") # path semicolon delimited lineage that I need as well parsing out last ; post process 
                w.writerow([acc, lineage]) 

    # build the BLAST index only once 
    if not (sd / "silva_db.nsq").exists():  #.nsq is one makeblastdb's outputs acting as a safety guard not to redownload  
        makeblastdb(fasta, sd / "silva_db") 

def fetch_ncbi() -> None:
    nd = db_root / "ncbi"; nd.mkdir(exist_ok=True)

    # ---- downloading BLAST database ------------ 
    tgz = nd / "16S_ribosomal_RNA.tar.gz"
    dl("ftp://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz", tgz)
    if not (nd / "16S_ribosomal_RNA.nsq").exists():
        log("→ extracting NCBI 16S")
        tarfile.open(tgz).extractall(nd)


    # generating accession -> lineage taxid table created acc_taxid.tsv once 
    tax_tsv = nd / "taxonomy.tsv" 
    if tax_tsv.exists():
        log("✓ taxonomy.tsv already present, skipping TaxonKit")
        return 
    

    log("-> generating NCBI taxonomy TSV (this may take a minute please wait... =)")


    # accession -> taxid (blastdbcmd) 
    acc_tax = nd / "acc_taxid.tsv"
    blast_db = nd / "16S_ribosomal_RNA"  
    run(["blastdbcmd", "-db", blast_db, "-entry", "all", "-outfmt", "%a\t%T", "-out", acc_tax]) 

    # taxid -> lineage via TaxonKit 
    log("-> running TaxonKit lineage (TaxonKit will fetch taxonomy DB if missing)") 

    tmp_path = nd / "lineage_raw.tsv" 
    ensure_taxdump() # download nodes.dmp etc if missing
    # lineage -> reformt -> tmp_path 
    with tmp_path.open("w") as out_fh: 
        subprocess.run(
                ["taxonkit", "lineage", "-n", "-i", "2", str(acc_tax)],
                check=True, stdout=out_fh
                )


        # rewrite tmp -> two-column TSV 
    with tmp_path.open() as src, tax_tsv.open("w") as dst:
            dst.write("sseqid\ttaxonomy\n")
            for ln in src:
                cols = ln.rstrip("\n").split("\t", 3)
                if len(cols) < 3: # some rows may have been filtered out 
                    continue 
                acc, lineage = cols[0], cols[2] # take full lineage 
                lineage = lineage.split(";", 1)[-1] # removes lineage past Domain rank which I don't want here  
                dst.write(f"{acc}\t{lineage}\n") 

    tmp_path.unlink(missing_ok=True) # clean up 
           

# ────────────────────────── main driver ─────────────────────────────────
def main() -> None:
    for fn in (fetch_gg2, fetch_silva, fetch_ncbi):
        fn()
    log(f" ✓ All databases downloaded to {db_root}")

    # ---------- env-hook snippet -----------------------------------------
    snippet = f"""
# ── MicroSeq auto-export ───────────────────────────────────────────────
export MICROSEQ_DB_HOME="{db_root}"
export BLASTDB="$MICROSEQ_DB_HOME/gg2:$MICROSEQ_DB_HOME/silva:$MICROSEQ_DB_HOME/ncbi"
export BLASTDB_LMDB=0
export MICROSEQ_LOG_DIR="{log_dir}"
# export MICROSEQ_SESSION_ID=$SLURM_JOB_ID      # Slurm
# export MICROSEQ_SESSION_ID=${{workflow.runName}} # Nextflow
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
    cfg.setdefault("logging", {})
    cfg["logging"]["dir"] = str(log_dir)
    cfg["logging"]["backup_count"] = 0 # keep every log
    if pipeline_user:
        cfg["logging"]["session_env"] = "MICROSEQ_SESSION_ID"
    cfg["databases"] = {
        "gg2":   {"blastdb": "${MICROSEQ_DB_HOME}/gg2/greengenes2_db",
            "taxonomy": "${MICROSEQ_DB_HOME}/gg2/taxonomy.tsv"},
        "silva": {"blastdb": "${MICROSEQ_DB_HOME}/silva/silva_db",
            "taxonomy": "${MICROSEQ_DB_HOME}/silva/taxonomy.tsv"},
        "ncbi":  {"blastdb": "${MICROSEQ_DB_HOME}/ncbi/16S_ribosomal_RNA",
            "taxonomy": "${MICROSEQ_DB_HOME}/ncbi/taxonomy.tsv"},
    }
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    log(f"✓ config.yaml updated at {cfg_path}")

    # --------- helper template locations ---------------------------------
    template_dir = REPO_ROOT / "config" / "templates"
    sh_tpl = template_dir / "activate_microseq_session.sh.example"
    fish_tpl = template_dir / "activate_microseq_session.fish.example"

    print("\nDone! Open a new shell (or type `exec $SHELL -l` or `source ~/.bashrc` or `source ~/zshrc`) and test it using the smoke test I have setup in github install:")
    print("Session helper templates:")
    print(f"  bash/zsh: {sh_tpl}")
    print(f"  fish:     {fish_tpl}")

if __name__ == "__main__":
    main() 
