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

[[ -n "$INPUT_DIR" && -n "$DB_KEY" && -n "$BLAST_THREADS" ]] || { usage; exit 1; }
[[ -d "$INPUT_DIR" ]] || { echo "Error: input_dir not found: $INPUT_DIR" >&2; exit 1; }
[[ "$DB_KEY" =~ ^(gg2|silva|ncbi)$ ]] || { echo "Error: db_key must be one of: gg2 | silva | ncbi" >&2; exit 1; }
[[ "$BLAST_THREADS" =~ ^[0-9]+$ && "$BLAST_THREADS" -ge 1 ]] || { echo "Error: blast_threads must be a positive integer" >&2; exit 1; }
command -v microseq >/dev/null 2>&1 || { echo "Error: microseq not found in PATH" >&2; exit 1; }

echo "[1/4] trim AB1 -> QC-passed reads"
microseq trim -i "$INPUT_DIR" --sanger --workdir "$OUT_DIR"

PAIRED_DIR="$OUT_DIR/qc/paired_fasta"
[[ -d "$PAIRED_DIR" ]] || { echo "Error: expected paired FASTA staging dir not found: $PAIRED_DIR" >&2; exit 1; }

echo "[2/4] paired assembly"
microseq assembly --mode paired \
  -i "$PAIRED_DIR" \
  -o "$OUT_DIR/asm" \
  --dup-policy error

BLAST_INPUT="$OUT_DIR/asm/blast_inputs.fasta"
[[ -f "$BLAST_INPUT" ]] || { echo "Error: expected BLAST input FASTA not found: $BLAST_INPUT" >&2; exit 1; }

echo "[3/4] BLAST"
microseq blast \
  -i "$BLAST_INPUT" \
  -d "$DB_KEY" \
  -o "$OUT_DIR/hits.tsv" \
  --threads "$BLAST_THREADS"

echo "[4/4] taxonomy join"
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
