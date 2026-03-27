#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  preview_pair_names.sh <input_dir>

What it does:
  - scans *.ab1 files in <input_dir>
  - checks whether forward/reverse primer labels look canonically delimited
  - reports CURRENT / STATUS / PRIMER_LABEL / ISSUE / SUGGESTED
  - optionally renames files in place

Current default primer labels:
  Forward: 27F, 8F, 515F
  Reverse: 1492R, 806R, 926R

Notes:
  - this script ignores well matching
  - it is a filename hygiene helper, not the authoritative pairing detector
  - after renaming, confirm with:
      microseq assembly --mode paired -i <input_dir> -o <preview_out> --dup-policy error --dry-run
EOF
}

INPUT_DIR="${1:-}"
[[ -n "$INPUT_DIR" ]] || { usage; exit 1; }
[[ -d "$INPUT_DIR" ]] || { echo "Error: input_dir not found: $INPUT_DIR" >&2; exit 1; }

shopt -s nullglob
FILES=( "$INPUT_DIR"/*.ab1 "$INPUT_DIR"/*.AB1 )
shopt -u nullglob

(( ${#FILES[@]} > 0 )) || { echo "Error: no .ab1 files found in $INPUT_DIR" >&2; exit 1; }

FWD_LABELS=(27F 8F 515F)
REV_LABELS=(1492R 806R 926R)

TMP_MAP="$(mktemp)"
trap 'rm -f "$TMP_MAP"' EXIT

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

  printf '%s\n' "${found[@]}"
}

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

  printf '%s\n' "${found[@]}"
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

  newstem="${newstem//__/_}"
  newstem="${newstem//--/-}"
  newstem="${newstem#_}"
  newstem="${newstem#-}"
  newstem="${newstem%_}"
  newstem="${newstem%-}"

  printf '%s\n' "$newstem"
}

printf '%-42s %-10s %-18s %-38s %s\n' "CURRENT" "STATUS" "PRIMER_LABEL" "ISSUE" "SUGGESTED"
printf '%-42s %-10s %-18s %-38s %s\n' "-------" "------" "------------" "-----" "---------"

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

  printf '%-42s %-10s %-18s %-38s %s\n' "$base" "$STATUS" "$PRIMER_LABEL" "$ISSUE" "$SUGGESTED"
  printf '%s\t%s\t%s\n' "$path" "$base" "$SUGGESTED" >> "$TMP_MAP"
done

echo
echo "Interpretation:"
echo "  OK      = primer label found and clearly delimited"
echo "  FIXABLE = primer label found but jammed into the sample name; suggested rename generated"
echo "  REVIEW  = missing or ambiguous primer label; inspect manually"

echo
echo "Recommended authoritative check after preview:"
echo "  microseq assembly --mode paired -i \"$INPUT_DIR\" -o \"${INPUT_DIR%/}_pair_preview\" --dup-policy error --dry-run"

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
echo "Now run the pipeline bash script I have paired_ab1_pipeline.sh"
echo "you can also do a dry run before hand as a sanity check for assembly like so but I'll leave that up to you as the user..."
echo "  microseq assembly --mode paired -i \"$INPUT_DIR\" -o \"${INPUT_DIR%/}_pair_preview\" --dup-policy error --dry-run"
