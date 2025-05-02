#!/usr/bin/env python3 
"""
fetch_dbs.py 
@author José Cantú

microseq-setup 

Douwnload/ unpack and register reference databases along with writing environmental hooks for MicroSeq:
    - Greengenes2 (2022.10 backbone full-length) "gg2"
    - SILVA 138.1 SSU NR99 "silva" 
    - NCBI 16S_ribosomal_RNA "ncbi16s" 

Creates:

    ${DB_ROOT}/{gg2,silva,ncbi}/greengenes2_db.*
    ${LOG_DIR}/    ---just the folder here   

Appends to ~/.bashrc and conda activate hook if inside env: 
    export MICROSEQ_DB_HOME="....."
    export BLASTDB="$MICROSEQ_DB_HOME/gg2:$MICROSEQ_DB_HOME/silva:$MICROSEQ_DB_HOME/ncbi 
    export BLASTDB_LMDB=0 
    export MICROSEQ_LOG_DIR="....." 

Also patches config/config.yaml so libary code can fall back when the environment variables are missing =) 

Files are stored under $MICROSEQ_DB_HOME (default is ~/.microseq_dbs). 
Each DB direcotry ends up with BLAST-ready indices (makeblastdb -dbtype nucl). 

This is designed to be only runned the one time..... 

"""
from __future__ import annotations 

import os, pathlib, subprocess, urllib.request, tarfile, sys, zipfile, shutil, gzip, argparse, yaml 
from pathlib import Path 
from datetime import datetime 

# ------------------------ Setting up CLI -------------------------------
ap = argparse.ArgumentParser(prog="microseq-setup")
ap.add_argument("--db-root", metavar="PATH", help="Folder where BLAST databases will live" 
                "(default: ~/.microseq_dbs)") 
ap.add_argument("--log-dir", metavar="PATH", help="Folder for MicroSeq run logs 

if __name__ == "__main__":
    main() 
