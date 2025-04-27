#!/usr/bin/env python3 
"""
gen_demo_fasta.py builds tests/data/demo.fasta with 5 short records 

Usage:
    python tests/gen_demo_fasta.py 
"""

from Bio import SeqIO
from Bio.Seq import Seq 
from pathlib import Path 
import random, itertools, sys 

SCRIPT_DIR = Path(__file__).resolve().parent # .../tests 
ROOT = SCRIPT_DIR.parent 
SRC_DIR = ROOT/ "data" / "final_fasta_output_cultivation_repo"
OUT = SCRIPT_DIR / "data" / "demo.fasta" 

N       = 5        # how many demo sequences
LENGTH  = 300      # trim to first 300 bp for speed

def main():
    if not SRC_DIR.exists():
        sys.exit(f" {SRC_DIR} does not exist – point SRC_DIR to a folder with FASTA files.")

    # gather sequences from every *.fasta in SRC_DIR
    records = list(
        itertools.chain.from_iterable(
            SeqIO.parse(fa, "fasta") for fa in SRC_DIR.glob("*.fasta")
        )
    )
    if len(records) < N:
        sys.exit(f"❌  only {len(records)} sequences found in {SRC_DIR} – need at least {N}")

    random.seed(42)
    demo = random.sample(records, N)
    for i, r in enumerate(demo, 1):
        r.id, r.description = f"demo{i}", ""
        r.seq = Seq(str(r.seq)[:LENGTH])

    OUT.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(demo, OUT, "fasta")
    print(f"✓ wrote {OUT} with {len(demo)} records")

if __name__ == "__main__":
    main()
