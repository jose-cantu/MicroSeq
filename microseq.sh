#!/usr/bin/env bash
##############################################################################
# blast_and_biom.sh – run MicroSeq BLAST → BIOM pipeline on trimmed FASTA
# ❶ activate conda ▸ ❷ export LMDB ▸ ❸ BLAST loop ▸ ❹ merge ▸ ❺ post-blast
##############################################################################
set -euo pipefail          # abort on error / unset var / pipe failure
set -x                       # echo each command (good for debugging)

export MICROSEQ_DB_HOME="$HOME/.microseq_dbs"   # used by config.yaml 
export BLASTDB="$MICROSEQ_DB_HOME/gg2"		# explicit for BLASt 
export BLASTDB_LMDB=0 				# disable LMDB on macOS
unset BLASTDB_LMDB_MAP_SIZE 			# safety - remove stale value 

# ── ❶ Bring the conda functions into this shell, then activate env ──────────
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate MicroSeq        # puts blastn & microseq on $PATH

python - <<'PY'
import os, subprocess, shlex, sys
print("child sees BLASTDB      =", os.environ.get("BLASTDB",  "<unset>"))
print("child sees BLASTDB_LMDB =", os.environ.get("BLASTDB_LMDB", "<unset>"))
print("child sees MAP_SIZE     =", os.environ.get("BLASTDB_LMDB_MAP_SIZE", "<unset>"))
subprocess.run(shlex.split("blastn -version"), check=True, stdout=sys.stdout)
PY

# ---- ① BLAST each FASTA (skip existing) -----------------------------------
cd ~/microseq_test
mkdir -p data/blast data/biom

for f in data/final_fasta_output_cultivation_repo/*.fasta; do
    sid=$(basename "${f%.fasta}")
    out="data/blast/${sid}.tsv"
    tmp="${out}.tmp"
    
    [[ -s "$out" ]] && { echo "[skip] $out"; continue; }
    
    microseq --workdir data blast -i "$f" -d gg2 -o "$tmp"
    
    if [[ -s "$tmp" && $(wc -l <"$tmp") -gt 1 ]]; then
        mv "$tmp" "$out"
        echo "[ok]   $out"
    else
        echo "[warn] no BLAST hit for $sid – removed" >&2
        rm -f "$tmp"
    fi
done


# ── ❹ Build BOTH files (full + slim) with fields in the right order ────────
awk '
BEGIN {
   # hdr_full = "qseqid\tsseqid\tpident\tqlen\tqcovhsp\tlength\tevalue\tbitscore\tstitle"
    hdr_slim = "sample_id\tsseqid\tevalue\tbitscore\tpident\tqcov"
   # print hdr_full >  "data/blast/all_hits_full.tsv"
    print hdr_slim  > "data/blast/all_hits.tsv"
}
FNR == 1 {                    # skip the header line of each per-sample file
    count = 0                 # reset row counter for optional 50-row cap
    next
}
$3 >= 97 {                    # keep only rows with ≥ 97 % identity
    gsub(/_trimmed$/, "", $1) # OPTIONAL: drop “_trimmed” from qseqid/sample_id
    
    ## 1 – write the full 9-column row for provenance
   # print $0 >> "data/blast/all_hits_full.tsv"
    
    ## 2 – write up to 50 rows per sample in the 6-column order postblast expects
    if (count++ < 50) {
	    # fields:     sample_id  sseqid  evalue($7)  bitscore($8)  pident($3)  qcovhsp($5)
        print $1,       $2,        $7,         $8,        $3,          $5 \
            >> "data/blast/all_hits.tsv"
    }
}
' OFS="\t" data/blast/*.tsv

# providence check 
du -h data/blast/all_hits.tsv 
wc -l data/blast/all_hits.tsv 
head data/blast/all_hits.tsv 


# 1 ⟶ addtaxonomy first
microseq add_taxonomy \
    --hits      data/blast/all_hits.tsv \
    --taxonomy  "$MICROSEQ_DB_HOME/gg2/taxonomy.tsv" \
    --out       data/blast/all_hits_with_tax.tsv

head -3 data/blast/all_hits_with_tax.tsv | column -t 


# ── ❺ Post-BLAST: build BIOM + CSV ──────────────────────────────────────────
microseq --workdir data postblast \
    -b data/blast/all_hits_with_tax.tsv \
    -m data/metadata/Final_Master_Seq_File.csv \
    --sample-col matched_isolate_id \
    -o data/biom/ocular_isolates.biom



    
# ── ❻ Quick sanity peeks ────────────────────────────────────────────────────
biom summarize-table -i data/biom/ocular_isolates.biom | head
echo "[info] first rows with taxon_names" 

