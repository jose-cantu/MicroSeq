#!/usr/bin/env bash
##############################################################################
# run_microseq.sh –  one-shot BLAST ➜ lineage ➜ BIOM pipeline
##############################################################################
set -euo pipefail
set -x                                              # echo commands for tracing

# -- 0 · env -----------------------------------------------------------------
export MICROSEQ_DB_HOME="$HOME/.microseq_dbs"
export BLASTDB="$MICROSEQ_DB_HOME/gg2"              # explicit for blastn
export BLASTDB_LMDB=0                               # macOS: disable LMDB
unset  BLASTDB_LMDB_MAP_SIZE                        # remove stale value

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate MicroSeq                             # blastn + microseq on PATH

# ---------------------------------------------------------------------------#
# ❶ BLAST each FASTA (skip if the .tsv already exists and is non-empty)
# ---------------------------------------------------------------------------#
cd  ~/microseq_tests
mkdir -p data/blast data/biom

for f in data/final_fasta_output_cultivation_repo/*.fasta; do
    sid=$(basename "${f%.fasta}")
    out="data/blast/${sid}.tsv"
    tmp="${out}.tmp"

    [[ -s "$out" ]] && { echo "[skip] $out"; continue; }

    microseq --workdir data blast  -i "$f" -d gg2 -o "$tmp"

    if [[ -s "$tmp" && $(wc -l <"$tmp") -gt 1 ]]; then
        mv "$tmp" "$out"
        echo "[ok]   $out"
    else
        echo "[warn] no BLAST hit for $sid – removed" >&2
        rm -f "$tmp"
    fi
done

# ---------------------------------------------------------------------------#
# ❷ Merge per-sample TSVs → slim & full tables (cap slim at 50 rows/sample)
# ---------------------------------------------------------------------------#
awk '
BEGIN {
    hdr_full = "qseqid\tsseqid\tpident\tqlen\tqcovhsp\tlength\tevalue\tbitscore\tstitle"
    hdr_slim = "sample_id\tsseqid\tevalue\tbitscore\tpident\tqcov"
    print hdr_full >  "data/blast/all_hits_full.tsv"
    print hdr_slim  > "data/blast/all_hits.tsv"
}
FNR == 1 { count = 0; next }               # skip header of each file
$3 >= 97 {                                 # keep pident ≥ 97 %
    gsub(/_trimmed$/, "", $1)              # strip “_trimmed” if present

    # -- full file ----------------------------------------------------------
    print $0 >> "data/blast/all_hits_full.tsv"

    # -- slim file (≤50 rows / sample) --------------------------------------
    if (count++ < 50) {
        # sample_id  sseqid  evalue($7)  bitscore($8)  pident($3)  qcovhsp($5)
        print $1,     $2,     $7,         $8,          $3,         $5 \
             >> "data/blast/all_hits.tsv"
    }
}
' OFS="\t"  data/blast/*.tsv

# optional: compress the full table to save space (~3-5 MB)
gzip -f data/blast/all_hits_full.tsv

du -h  data/blast/all_hits.tsv
wc -l  data/blast/all_hits.tsv | xargs echo "[rows]"
head   data/blast/all_hits.tsv

# ---------------------------------------------------------------------------#
# ❸ Append Greengenes-2 lineage (taxon_name) – creates *_with_tax.tsv
# ---------------------------------------------------------------------------#
microseq add_taxonomy \
        --hits      data/blast/all_hits.tsv \
        --taxonomy  "$MICROSEQ_DB_HOME/gg2/taxonomy.tsv" \
        --out       data/blast/all_hits_with_tax.tsv

head -3 data/blast/all_hits_with_tax.tsv | column -t

# ---------------------------------------------------------------------------#
# ❹ Post-BLAST  ➜ BIOM  (+ CSV mirror)
# ---------------------------------------------------------------------------#
microseq --workdir data postblast \
        -b data/blast/all_hits_with_tax.tsv \
        -m data/metadata/Final_Master_Seq_File.csv \
        --sample-col matched_isolate_id \
        -o ocular_isolates.biom

# ---------------------------------------------------------------------------#
# ❺ Quick sanity peeks
# ---------------------------------------------------------------------------#
biom summarize-table -i data/biom/ocular_isolates.biom | head
echo "[tax preview]"
head -5 data/blast/all_hits_with_tax.tsv | cut -f1-3
##############################################################################

