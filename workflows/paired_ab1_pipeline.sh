#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  paired_ab1_pipeline.sh <input_dir> <db_key> <blast_threads> [out_dir]

Arguments:
  <input_dir>       Directory containing raw paired AB1 files
  <db_key>          gg2 | silva | ncbi
  <blast_threads>   Positive integer; used only for BLAST
  [out_dir]         Optional output directory
                    Default: <input_dir>_microseq

Examples:
  paired_ab1_pipeline.sh /home/jason/MicroSeq/tests/paired_single_pair_ab1_demo_run gg2 4
  paired_ab1_pipeline.sh /home/jason/MicroSeq/tests/paired_ab1_demo_run/10292025_1080497 gg2 4
EOF
}

INPUT_DIR="${1:-}"
DB_KEY="${2:-}"
BLAST_THREADS="${3:-}"
OUT_DIR="${4:-${INPUT_DIR%/}_microseq}"
DUP_POLICY="${DUP_POLICY:-error}"
FWD_PATTERN="${FWD_PATTERN:-27F}"
REV_PATTERN="${REV_PATTERN:-1492R}"

# Create a one run scoped session id for this wrapper I made 
# If the user already exported MICROSEQ_SESSION_ID outside the script, keep it then.
RUN_SESSION_ID="${MICROSEQ_SESSION_ID:-$(date +%Y%m%d-%H%M%S)-paired-wrapper}"

# Export it once so every later MicroSeq command inherits the same session id 
export MICROSEQ_SESSION_ID="$RUN_SESSION_ID"

# Show session id being used for this whole wrapper run for my user 
echo "MICROSEQ_SESSION_ID=$MICROSEQ_SESSION_ID"

[[ -n "$INPUT_DIR" && -n "$DB_KEY" && -n "$BLAST_THREADS" ]] || { usage; exit 1; }
[[ -d "$INPUT_DIR" ]] || { echo "Error: input_dir not found: $INPUT_DIR" >&2; exit 1; }
[[ "$DB_KEY" =~ ^(gg2|silva|ncbi)$ ]] || { echo "Error: db_key must be one of: gg2 | silva | ncbi" >&2; exit 1; }
[[ "$BLAST_THREADS" =~ ^[0-9]+$ && "$BLAST_THREADS" -ge 1 ]] || { echo "Error: blast_threads must be a positive integer" >&2; exit 1; }
command -v microseq >/dev/null 2>&1 || { echo "Error: microseq not found in PATH" >&2; exit 1; }

echo "[1/6] trim AB1 -> QC-passed reads"
microseq  --workdir  "$OUT_DIR" trim -i "$INPUT_DIR" --sanger 

FASTQ_STAGING_DIR="$OUT_DIR/passed_qc_fastq_primer_trim"
if [[ ! -d "$FASTQ_STAGING_DIR" ]] || [[ -z "$(find "$FASTQ_STAGING_DIR" -type f -name '*.fastq' -print -quit)" ]]; then
  FASTQ_STAGING_DIR="$OUT_DIR/passed_qc_fastq"
fi

[[ -d "$FASTQ_STAGING_DIR" ]] || { echo "Error: expected FASTQ staging dir not found after trim: $FASTQ_STAGING_DIR" >&2; exit 1; }
[[ -n "$(find "$FASTQ_STAGING_DIR" -type f -name '*.fastq' -print -quit)" ]] || { echo "Error: no FASTQ files found in $FASTQ_STAGING_DIR" >&2; exit 1; }

echo "[2/6] stage paired FASTA"
microseq stage-paired-fasta \
  -i "$FASTQ_STAGING_DIR" \
  -o "$OUT_DIR/qc/paired_fasta"

PAIRED_FASTA_DIR="$OUT_DIR/qc/paired_fasta"
[[ -d "$PAIRED_FASTA_DIR" ]] || { echo "Error: expected paired FASTA staging dir not found: $PAIRED_FASTA_DIR" >&2; exit 1; }

# write canonical pairing report using MIcroSeq 
microseq pairing-report \
  -i "$PAIRED_FASTA_DIR" \
  -o "$OUT_DIR/qc/pairing_report.tsv" \
  --dup-policy "$DUP_POLICY" \
  --fwd-pattern "$FWD_PATTERN" \
  --rev-pattern "$REV_PATTERN" 


echo "[3/6] paired assembly"
microseq assembly --mode paired \
  -i "$PAIRED_FASTA_DIR" \
  -o "$OUT_DIR/asm" \
  --dup-policy error \
  --fwd-pattern "27F" \
  --rev-pattern "1492R"

# Write canonical assembly summary 
microseq assembly-summary \
  --asm-dir "$OUT_DIR/asm" \
  --pairing-input-dir "$PAIRED_FASTA_DIR" \
  -o "$OUT_DIR/asm/assembly_summary.tsv" \
  --dup-policy "$DUP_POLICY" \
  --fwd-pattern "$FWD_PATTERN" \
  --rev-pattern "$REV_PATTERN"

# Write canonical overlap audit using MicroSeq 
microseq overlap-audit \
  -i "$PAIRED_FASTA_DIR" \
  -o "$OUT_DIR/qc/overlap_audit.tsv" \
  --dup-policy "$DUP_POLICY" \
  --fwd-pattern "$FWD_PATTERN" \
  --rev-pattern "$REV_PATTERN" \
  --pretrim-input-dir "$OUT_DIR/qc/paired_fasta_pretrim" \
  --primer-trim-report "$OUT_DIR/qc/primer_trim_report.tsv" 

echo "[4/6] build BLAST input FASTA"
BLAST_INPUT="$OUT_DIR/asm/blast_inputs.fasta"
# Write canonical BLAST input FASTA + manifest using MicroSeq's 
microseq blast-inputs \
  --asm-dir "$OUT_DIR/asm" \
  --pairing-input-dir "$PAIRED_FASTA_DIR" \
  --output-fasta "$OUT_DIR/asm/blast_inputs.fasta" \
  --output-tsv "$OUT_DIR/asm/blast_inputs.tsv" \
  --dup-policy "$DUP_POLICY" \
  --fwd-pattern "$FWD_PATTERN" \
  --rev-pattern "$REV_PATTERN"


echo "[5/6] BLAST"
microseq blast \
  -i "$BLAST_INPUT" \
  -d "$DB_KEY" \
  -o "$OUT_DIR/hits.tsv" \
  --threads "$BLAST_THREADS"

echo "[6/6] taxonomy join"
microseq add_taxonomy \
  -i "$OUT_DIR/hits.tsv" \
  -d "$DB_KEY" \
  -o "$OUT_DIR/hits_tax.tsv"

cat <<EOF

Done.

Key outputs:
  trim summary      : $OUT_DIR/qc/trim_summary.tsv
  pairing report    : $OUT_DIR/qc/pairing_report.tsv
  assembly summary  : $OUT_DIR/asm/assembly_summary.tsv
  blast payload map : $OUT_DIR/asm/blast_inputs.tsv
  blast hits        : $OUT_DIR/hits.tsv
  hits + taxonomy   : $OUT_DIR/hits_tax.tsv
EOF
