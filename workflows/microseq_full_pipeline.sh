#!/usr/bin/env bash
# MicroSeq full workflow (no vsearch, no BIOM/postblast).
#
# Usage:
#   microseq_full.sh <input> <db_key> <threads> [out_dir]
#
# <input>  : FASTA file OR FASTQ file/dir OR AB1 dir
# <db_key> : gg2 | silva | ncbi  (must exist in your config.yaml databases section)
# <threads>: threads for BLAST
# [out_dir]: optional output directory (default: microseq_run_YYYYmmdd_HHMMSS)
#
# Notes:
# - If input is FASTA, trimming/conversion are skipped.
# - If input is FASTQ/AB1, outputs are written into out_dir/ and out_dir/qc/.
#

# @AUTHOR JC 
# Date: 2025/01/19

input=$1
db_key=$2
threads=$3
out_dir=${4:-microseq_run_$(date +%Y%m%d_%H%M%S)}

# ---------- validation ----------
if [[ -z "${input}" || -z "${db_key}" || -z "${threads}" ]]; then
  echo "Usage: $0 <input> <db_key> <threads> [out_dir]" 1>&2
  exit 1
fi

if [[ ! -e "${input}" ]]; then
  echo "Error: input path does not exist: ${input}" 1>&2
  exit 1
fi

if ! [[ "${threads}" =~ ^[0-9]+$ ]]; then
  echo "Error: threads must be an integer: ${threads}" 1>&2
  exit 1
fi

command -v microseq >/dev/null 2>&1 || { echo "Error: microseq not found in PATH" 1>&2; exit 1; }

mkdir -p "${out_dir}"
mkdir -p "${out_dir}/qc"

# ---------- helpers ----------
is_fasta() {
  [[ "$1" =~ \.(fa|fna|fasta|fas)$ ]]
}

pick_fastq_dir() {
  # MicroSeq sometimes writes to qc/ and sometimes to passed_qc_fastq/ depending on your trim logic.
  local qc_dir="${out_dir}/qc"
  local alt_dir="${out_dir}/passed_qc_fastq"

  shopt -s nullglob
  local qc_fastq=("${qc_dir}"/*.fastq "${qc_dir}"/*.fq)
  local alt_fastq=("${alt_dir}"/*.fastq "${alt_dir}"/*.fq)
  shopt -u nullglob

  if (( ${#qc_fastq[@]} > 0 )); then
    echo "${qc_dir}"
    return
  fi
  if (( ${#alt_fastq[@]} > 0 )); then
    echo "${alt_dir}"
    return
  fi

  echo "Error: no FASTQ files found after trimming in ${qc_dir} or ${alt_dir}" 1>&2
  exit 1
}

# ---------- stage 0: set stable output paths ----------
reads_fasta="${out_dir}/reads.fasta"
hits_tsv="${out_dir}/hits.tsv"
hits_tax_tsv="${out_dir}/hits_tax.tsv"

# ---------- stage 1: trim / convert ----------
if is_fasta "${input}"; then
  echo "[1/3] Input is FASTA: skipping trim + fastq-to-fasta"
  reads_fasta="$(cd "$(dirname "${input}")" && pwd)/$(basename "${input}")"
else
  echo "[1/3] Trim/QC"
  # NOTE: --workdir is a GLOBAL arg in your argparse (if you kept it that way), so it must come before the subcommand:
  microseq --workdir "${out_dir}" trim -i "${input}"

  echo "[2/3] FASTQ -> FASTA"
  fastq_dir="$(pick_fastq_dir)"
  microseq fastq-to-fasta -i "${fastq_dir}" -o "${reads_fasta}"
fi

# ---------- stage 2: BLAST ----------
echo "[3/3] BLAST"
microseq blast -i "${reads_fasta}" -d "${db_key}" -o "${hits_tsv}" --threads "${threads}"

# ---------- stage 3: taxonomy join ----------
echo "[4/3] Add taxonomy"
microseq add_taxonomy -i "${hits_tsv}" -d "${db_key}" -o "${hits_tax_tsv}"

echo ""
echo "Done."
echo "Outputs:"
echo "  FASTA:      ${reads_fasta}"
echo "  BLAST hits: ${hits_tsv}"
echo "  Hits+tax:   ${hits_tax_tsv}"

