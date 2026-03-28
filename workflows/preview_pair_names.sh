#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  preview_pair_names.sh <input_dir>

What it does:
  - recursively scans sequence-like files under <input_dir>
  - checks whether forward/reverse primer labels look canonically delimited
  - reports CURRENT / STATUS / PRIMER_LABEL / ISSUE / SUGGESTED
  - optionally renames files in place
  - recommends the actual MicroSeq paired preview step afterward

Current default primer labels:
  Forward: 27F, 8F, 515F
  Reverse: 1492R, 806R, 926R

Notes:
  - this script ignores well matching
  - this script is a filename hygiene helper, not the authoritative pairing step
  - after preview/rename, run the actual MicroSeq preview:
      microseq assembly --mode paired -i <input_dir> -o <preview_out> --dup-policy error --preview-pairs
EOF
}

INPUT_DIR="${1:-}"
[[ -n "$INPUT_DIR" ]] || { usage; exit 1; }
[[ -d "$INPUT_DIR" ]] || { echo "Error: input_dir not found: $INPUT_DIR" >&2; exit 1; }

# GUI-like broad discovery for common sequence-bearing file types.
mapfile -t FILES < <(
  find "$INPUT_DIR" -type f \
    \( -iname '*.ab1' -o -iname '*.seq' -o -iname '*.fastq' -o -iname '*.fq' \
       -o -iname '*.fasta' -o -iname '*.fa' -o -iname '*.fas' \) \
    | LC_ALL=C sort
)

(( ${#FILES[@]} > 0 )) || {
  echo "Error: no supported files found under $INPUT_DIR" >&2
  exit 1
}

FWD_LABELS=(27F 8F 515F)
REV_LABELS=(1492R 806R 926R)

TMP_MAP="$(mktemp)"
trap 'rm -f "$TMP_MAP"' EXIT

# Print nothing at all when there are no matches.
emit_found_lines() {
  local -n arr_ref=$1
  if (( ${#arr_ref[@]} > 0 )); then
    printf '%s\n' "${arr_ref[@]}"
  fi
}

# Canonical delimiter policy for this helper:
# primer labels should be separated by "_" or "-" or appear at string boundaries.
# "+" is intentionally NOT treated as canonical, so A+1492R... becomes FIXABLE.
detect_safe_primer_label() {
  local stem="$1"
  local found=()
  local label

  for label in "${FWD_LABELS[@]}"; do
    if [[ "$stem" =~ (^|[_-])${label}([_-]|$) ]]; then
      found+=("forward:${label}")
    fi
  done

  for label in "${REV_LABELS[@]}"; do
    if [[ "$stem" =~ (^|[_-])${label}([_-]|$) ]]; then
      found+=("reverse:${label}")
    fi
  done

  emit_found_lines found
}

# Relaxed detection: substring-only. Used to identify likely fixable cases.
detect_relaxed_primer_label() {
  local stem="$1"
  local found=()
  local label

  for label in "${FWD_LABELS[@]}"; do
    if [[ "$stem" == *"$label"* ]]; then
      found+=("forward:${label}")
    fi
  done

  for label in "${REV_LABELS[@]}"; do
    if [[ "$stem" == *"$label"* ]]; then
      found+=("reverse:${label}")
    fi
  done

  emit_found_lines found
}

canonicalize_not_delimited() {
  local stem="$1"
  local primer_label="$2"

  local prefix="${stem%%${primer_label}*}"
  local suffix="${stem#*${primer_label}}"

  prefix="${prefix%_}"
  prefix="${prefix%-}"
  suffix="${suffix#_}"
  suffix="${suffix#-}"

  local newstem
  if [[ -n "$suffix" ]]; then
    newstem="${prefix}_${primer_label}_${suffix}"
  else
    newstem="${prefix}_${primer_label}"
  fi

  # light cleanup
  while [[ "$newstem" == *__* ]]; do newstem="${newstem//__/_}"; done
  while [[ "$newstem" == *--* ]]; do newstem="${newstem//--/-}"; done
  newstem="${newstem#_}"
  newstem="${newstem#-}"
  newstem="${newstem%_}"
  newstem="${newstem%-}"

  printf '%s\n' "$newstem"
}

printf '%-50s %-10s %-18s %-38s %s\n' "CURRENT" "STATUS" "PRIMER_LABEL" "ISSUE" "SUGGESTED"
printf '%-50s %-10s %-18s %-38s %s\n' "-------" "------" "------------" "-----" "---------"

NEEDS_APPLY=0

for path in "${FILES[@]}"; do
  base="$(basename "$path")"
  ext=".${base##*.}"
  stem="${base%.*}"

  mapfile -t SAFE < <(detect_safe_primer_label "$stem")
  mapfile -t RELAXED < <(detect_relaxed_primer_label "$stem")

  STATUS="OK"
  PRIMER_LABEL=""
  ISSUE=""
  SUGGESTED=""

  if (( ${#SAFE[@]} == 1 )); then
    PRIMER_LABEL="${SAFE[0]}"
    ISSUE="primer label canonically delimited"
  elif (( ${#SAFE[@]} > 1 )); then
    STATUS="REVIEW"
    PRIMER_LABEL="$(IFS=,; echo "${SAFE[*]}")"
    ISSUE="ambiguous forward/reverse primer labels"
  elif (( ${#RELAXED[@]} == 1 )); then
    STATUS="FIXABLE"
    PRIMER_LABEL="${RELAXED[0]}"
    ISSUE="primer label not delimited"

    PRIMER_LABEL_VALUE="${RELAXED[0]#*:}"
    NEWSTEM="$(canonicalize_not_delimited "$stem" "$PRIMER_LABEL_VALUE")"
    SUGGESTED="${NEWSTEM}${ext}"
    NEEDS_APPLY=1
  elif (( ${#RELAXED[@]} > 1 )); then
    STATUS="REVIEW"
    PRIMER_LABEL="$(IFS=,; echo "${RELAXED[*]}")"
    ISSUE="multiple primer label matches"
  else
    STATUS="REVIEW"
    PRIMER_LABEL="NA"
    ISSUE="no forward/reverse primer label detected"
  fi

  printf '%-50s %-10s %-18s %-38s %s\n' "$base" "$STATUS" "$PRIMER_LABEL" "$ISSUE" "$SUGGESTED"
  printf '%s\t%s\t%s\n' "$path" "$base" "$SUGGESTED" >> "$TMP_MAP"
done

echo
echo "Interpretation:"
echo "  OK      = primer label found and clearly delimited"
echo "  FIXABLE = primer label found but jammed into the sample name; suggested rename generated"
echo "  REVIEW  = missing or ambiguous primer label; inspect manually"

echo
echo "Recommended authoritative MicroSeq preview:"
echo "  microseq assembly --mode paired -i \"$INPUT_DIR\" -o \"${INPUT_DIR%/}_pair_preview\" --dup-policy error --preview-pairs"
echo
echo "Then run:"
echo "  bash workflows/paired_ab1_pipeline.sh \"$INPUT_DIR\" gg2 4"

if (( NEEDS_APPLY == 0 )); then
  echo
  echo "No auto-fixable names detected."
  exit 0
fi

echo
read -r -p "Apply suggested FIXABLE renames in place? [y/N] " RESP
if [[ ! "$RESP" =~ ^[Yy]$ ]]; then
  echo "Rename step skipped."
  exit 0
fi

DUPES="$(awk -F'\t' '$3 != "" {print $3}' "$TMP_MAP" | sort | uniq -d || true)"
if [[ -n "$DUPES" ]]; then
  echo "Error: rename collisions detected. Refusing to rename in place." >&2
  echo "$DUPES" >&2
  exit 1
fi

RENAME_LOG="$INPUT_DIR/rename_map.tsv"
: > "$RENAME_LOG"
printf "old_path\told_name\tnew_name\n" >> "$RENAME_LOG"

while IFS=$'\t' read -r SRC OLD NEW; do
  [[ -n "$NEW" ]] || continue
  DST="$(dirname "$SRC")/$NEW"
  [[ "$SRC" == "$DST" ]] && continue

  if [[ -e "$DST" ]]; then
    echo "Error: destination already exists: $DST" >&2
    exit 1
  fi

  mv -- "$SRC" "$DST"
  printf "%s\t%s\t%s\n" "$SRC" "$OLD" "$NEW" >> "$RENAME_LOG"
done < "$TMP_MAP"

echo
echo "Renames applied."
echo "Audit log: $RENAME_LOG"
echo
echo "Now run:"
echo "  microseq assembly --mode paired -i \"$INPUT_DIR\" -o \"${INPUT_DIR%/}_pair_preview\" --dup-policy error --preview-pairs"
echo
echo "Then run:"
echo "  bash workflows/paired_ab1_pipeline.sh \"$INPUT_DIR\" gg2 4"
