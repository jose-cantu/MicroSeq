#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  paired_ab1_pipeline.sh <input_dir> <db_key> <blast_threads> [out_dir] \
    [--identity N] \
    [--qcov N] \
    [--max-hits N]

Arguments:
  <input_dir>       Directory containing raw paired AB1 files
  <db_key>          gg2 | silva | ncbi
  <blast_threads>   Positive integer; used for BLAST
  [out_dir]         Optional output directory
                    Default: <input_dir>_microseq

Optional BLAST overrides:
  --identity N      Override CLI default BLAST identity
  --qcov N          Override CLI default BLAST query coverage
  --max-hits N      Override CLI default BLAST max_target_seqs

Examples:
  paired_ab1_pipeline.sh /home/jason/MicroSeq/tests/paired_single_pair_ab1_demo_run gg2 4
  paired_ab1_pipeline.sh /home/jason/MicroSeq/tests/paired_ab1_demo_run/10292025_1080497 gg2 4 run_strict
  paired_ab1_pipeline.sh /home/jason/MicroSeq/tests/paired_ab1_demo_run/10292025_1080497 gg2 4 \
    --identity 99 --qcov 90 --max-hits 5
  paired_ab1_pipeline.sh /home/jason/MicroSeq/tests/paired_ab1_demo_run/10292025_1080497 gg2 4 run_strict \
    --identity 99 --qcov 90 --max-hits 10
EOF
}

# ---- required positionals -------------------------------------------------
INPUT_DIR="${1:-}"                                  # req: raw paired AB1 directory
DB_KEY="${2:-}"                                     # req: database key
BLAST_THREADS="${3:-}"                              # req: blast threads
shift $(( $# >= 3 ? 3 : $# ))                       # consume required positionals if present

# ---- optional positional out_dir ------------------------------------------
OUT_DIR="${INPUT_DIR%/}_microseq"                   # default output dir
if [[ $# -gt 0 && "${1:-}" != --* ]]; then          # only treat 4th positional as out_dir if it is not a flag
  OUT_DIR="$1"
  shift
fi

# ---- optional named BLAST overrides ---------------------------------------
IDENTITY=""                                         # empty => inherit CLI canonical default
QCOV=""                                             # empty => inherit CLI canonical default
MAX_HITS=""                                         # empty => inherit CLI canonical default

# ---- other wrapper policy knobs -------------------------------------------
DUP_POLICY="${DUP_POLICY:-error}"
FWD_PATTERN="${FWD_PATTERN:-27F}"
REV_PATTERN="${REV_PATTERN:-1492R}"

# ---- parse remaining named flags ------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --identity)
      [[ $# -ge 2 ]] || { echo "Error: --identity requires a value" >&2; usage; exit 1; }
      IDENTITY="$2"
      shift 2
      ;;
    --qcov)
      [[ $# -ge 2 ]] || { echo "Error: --qcov requires a value" >&2; usage; exit 1; }
      QCOV="$2"
      shift 2
      ;;
    --max-hits)
      [[ $# -ge 2 ]] || { echo "Error: --max-hits requires a value" >&2; usage; exit 1; }
      MAX_HITS="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Error: unknown argument: $1" >&2
      usage
      exit 1
      ;;
  esac
done

# ---- one-run scoped session id --------------------------------------------
RUN_SESSION_ID="${MICROSEQ_SESSION_ID:-$(date +%Y%m%d-%H%M%S)-paired-wrapper}"
export MICROSEQ_SESSION_ID="$RUN_SESSION_ID"
echo "MICROSEQ_SESSION_ID=$MICROSEQ_SESSION_ID"

# ---- validation -----------------------------------------------------------
[[ -n "$INPUT_DIR" && -n "$DB_KEY" && -n "$BLAST_THREADS" ]] || { usage; exit 1; }
[[ -d "$INPUT_DIR" ]] || { echo "Error: input_dir not found: $INPUT_DIR" >&2; exit 1; }
[[ "$DB_KEY" =~ ^(gg2|silva|ncbi)$ ]] || { echo "Error: db_key must be one of: gg2 | silva | ncbi" >&2; exit 1; }
[[ "$BLAST_THREADS" =~ ^[0-9]+$ && "$BLAST_THREADS" -ge 1 ]] || {
  echo "Error: blast_threads must be a positive integer" >&2
  exit 1
}
command -v microseq >/dev/null 2>&1 || { echo "Error: microseq not found in PATH" >&2; exit 1; }

# Optional light numeric validation here.
# Keep canonical range enforcement in the CLI/backend layer.
if [[ -n "$IDENTITY" && ! "$IDENTITY" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "Error: --identity must be numeric" >&2
  exit 1
fi
if [[ -n "$QCOV" && ! "$QCOV" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
  echo "Error: --qcov must be numeric" >&2
  exit 1
fi
if [[ -n "$MAX_HITS" && ! "$MAX_HITS" =~ ^[0-9]+$ ]]; then
  echo "Error: --max-hits must be an integer" >&2
  exit 1
fi

echo "[1/6] trim AB1 -> QC-passed reads"
microseq --workdir "$OUT_DIR" trim -i "$INPUT_DIR" --sanger

FASTQ_STAGING_DIR="$OUT_DIR/passed_qc_fastq_primer_trim"
if [[ ! -d "$FASTQ_STAGING_DIR" ]] || [[ -z "$(find "$FASTQ_STAGING_DIR" -type f -name '*.fastq' -print -quit)" ]]; then
  FASTQ_STAGING_DIR="$OUT_DIR/passed_qc_fastq"
fi

[[ -d "$FASTQ_STAGING_DIR" ]] || { echo "Error: expected FASTQ staging dir not found after trim: $FASTQ_STAGING_DIR" >&2; exit 1; }
[[ -n "$(find "$FASTQ_STAGING_DIR" -type f -name '*.fastq' -print -quit)" ]] || {
  echo "Error: no FASTQ files found in $FASTQ_STAGING_DIR" >&2
  exit 1
}

echo "[2/6] stage paired FASTA"
microseq stage-paired-fasta \
  -i "$FASTQ_STAGING_DIR" \
  -o "$OUT_DIR/qc/paired_fasta"

PAIRED_FASTA_DIR="$OUT_DIR/qc/paired_fasta"
[[ -d "$PAIRED_FASTA_DIR" ]] || { echo "Error: expected paired FASTA staging dir not found: $PAIRED_FASTA_DIR" >&2; exit 1; }

echo "[2.5/6] pairing report"
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
  --dup-policy "$DUP_POLICY" \
  --fwd-pattern "$FWD_PATTERN" \
  --rev-pattern "$REV_PATTERN"

echo "[3.3/6] assembly summary"
microseq assembly-summary \
  --asm-dir "$OUT_DIR/asm" \
  --pairing-input-dir "$PAIRED_FASTA_DIR" \
  -o "$OUT_DIR/asm/assembly_summary.tsv" \
  --dup-policy "$DUP_POLICY" \
  --fwd-pattern "$FWD_PATTERN" \
  --rev-pattern "$REV_PATTERN"

echo "[3.6/6] overlap audit"
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
microseq blast-inputs \
  --asm-dir "$OUT_DIR/asm" \
  --pairing-input-dir "$PAIRED_FASTA_DIR" \
  --output-fasta "$OUT_DIR/asm/blast_inputs.fasta" \
  --output-tsv "$OUT_DIR/asm/blast_inputs.tsv" \
  --dup-policy "$DUP_POLICY" \
  --fwd-pattern "$FWD_PATTERN" \
  --rev-pattern "$REV_PATTERN"

echo "[5/6] BLAST"
blast_cmd=(
  microseq blast
  -i "$BLAST_INPUT"
  -d "$DB_KEY"
  -o "$OUT_DIR/hits.tsv"
  --threads "$BLAST_THREADS"
)

[[ -n "$IDENTITY" ]] && blast_cmd+=( --identity "$IDENTITY" )           # opt: override canonical default
[[ -n "$QCOV" ]] && blast_cmd+=( --qcov "$QCOV" )                       # opt: override canonical default
[[ -n "$MAX_HITS" ]] && blast_cmd+=( --max_target_seqs "$MAX_HITS" )    # opt: wrapper name -> CLI name

printf 'Running:'                               # show exact effective BLAST command
printf ' %q' "${blast_cmd[@]}"
printf '\n'

"${blast_cmd[@]}"

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
  blast input sequences : $OUT_DIR/asm/blast_inputs.tsv
  blast hits        : $OUT_DIR/hits.tsv
  hits + taxonomy   : $OUT_DIR/hits_tax.tsv
EOF
