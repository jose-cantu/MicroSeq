#!/usr/bin/env bash
# MicroSeq full workflow with vsearch collapse + reference chimera filtering (no BIOM/postblast).
#
# Usage:
#   microseq_full_with_vsearch.sh <input> <db_key> <threads> [out_dir]
#

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
command -v vsearch >/dev/null 2>&1 || { echo "Error: vsearch not found in PATH (required for this workflow)" 1>&2; exit 1; }

mkdir -p "${out_dir}"
mkdir -p "${out_dir}/qc"

# ---------- helpers ----------
is_fasta() { [[ "$1" =~ \.(fa|fna|fasta|fas)$ ]]; }

pick_fastq_dir() {
  local qc_dir="${out_dir}/qc"
  local alt_dir="${out_dir}/passed_qc_fastq"
  shopt -s nullglob
  local qc_fastq=("${qc_dir}"/*.fastq "${qc_dir}"/*.fq)
  local alt_fastq=("${alt_dir}"/*.fastq "${alt_dir}"/*.fq)
  if (( ${#qc_fastq[@]} > 0 )); then echo "${qc_dir}"; return; fi
  if (( ${#alt_fastq[@]} > 0 )); then echo "${alt_dir}"; return; fi
  echo "Error: no FASTQ files found after trimming in ${qc_dir} or ${alt_dir}" 1>&2
  exit 1
}

strip_size_headers() {
  # Removes ';size=123;' tokens from FASTA headers.
  # Works on GNU/Linux and macOS.
  local in_fa=$1
  local out_fa=$2
  perl -pe 'if(/^>(\S+)/){$h=$1;$h=~s/;size=\d+;?//g;$_=">$h\n"}' "${in_fa}" > "${out_fa}"
}

# ---------- stage 0: output paths ----------
reads_fasta="${out_dir}/reads.fasta"
collapsed_fasta="${out_dir}/qc/replicates_collapsed.fasta"
nonchimera_fasta="${out_dir}/qc/nonchimeras.fasta"
blast_fasta="${out_dir}/qc/nonchimeras_clean.fasta"
hits_tsv="${out_dir}/hits.tsv"
hits_tax_tsv="${out_dir}/hits_tax.tsv"

# ---------- stage 1: trim / convert ----------
if is_fasta "${input}"; then
  echo "[1/5] Input is FASTA: skipping trim + fastq-to-fasta"
  reads_fasta="$(cd "$(dirname "${input}")" && pwd)/$(basename "${input}")"
else
  echo "[1/5] Trim/QC"
  microseq --workdir "${out_dir}" trim -i "${input}"

  echo "[2/5] FASTQ -> FASTA"
  fastq_dir="$(pick_fastq_dir)"
  microseq fastq-to-fasta -i "${fastq_dir}" -o "${reads_fasta}"
fi

# ---------- stage 2: per-sample replicate collapse ----------
echo "[3/5] vsearch collapse (per-sample)"
microseq vsearch-collapse -i "${reads_fasta}" -o "${collapsed_fasta}" --threads "${threads}"

# ---------- stage 3: reference chimera filtering ----------
echo "[4/5] vsearch chimera (reference)"
microseq vsearch-chimera -i "${collapsed_fasta}" -o "${nonchimera_fasta}" -d "${db_key}" --threads "${threads}"

# ---------- stage 3.5: clean IDs for BLAST ----------
echo "[4.5/5] Clean FASTA IDs for BLAST (strip ;size=)"
strip_size_headers "${nonchimera_fasta}" "${blast_fasta}"

# ---------- stage 4: BLAST ----------
echo "[5/5] BLAST"
microseq blast -i "${blast_fasta}" -d "${db_key}" -o "${hits_tsv}" --threads "${threads}"

# ---------- stage 5: taxonomy join ----------
echo "[6/5] Add taxonomy"
microseq add_taxonomy -i "${hits_tsv}" -d "${db_key}" -o "${hits_tax_tsv}"

echo ""
echo "Done."
echo "Outputs:"
echo "  FASTA:      ${reads_fasta}"
echo "  Collapsed:  ${collapsed_fasta}"
echo "  Nonchim:    ${nonchimera_fasta}"
echo "  BLAST FASTA:${blast_fasta}"
echo "  BLAST hits: ${hits_tsv}"
echo "  Hits+tax:   ${hits_tax_tsv}"

