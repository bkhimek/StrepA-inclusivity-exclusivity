#!/usr/bin/env bash
set -euo pipefail

LIST="${LIST:-demo_accessions.txt}"
OUTDIR="${OUTDIR:-pipelines/panaroo/2_panaroo/demo}"
THREADS="${THREADS:-8}"

# Build array of GFFs
mapfile -t GFFS < <(awk '{printf "pipelines/panaroo/1_annotate_prokka/%s/%s.gff\n",$1,$1}' "$LIST")

# Sanity checks
if [[ ${#GFFS[@]} -lt 2 ]]; then
  echo "Need at least 2 GFFs. Check $LIST and Step 1 outputs." >&2
  exit 1
fi

missing=0
for g in "${GFFS[@]}"; do
  if [[ ! -s "$g" ]]; then
    echo "[MISSING] $g"
    ((missing++)) || true
  fi
done
if [[ $missing -gt 0 ]]; then
  echo "Abort: $missing GFF(s) missing. Run Prokka first for all accessions in $LIST." >&2
  exit 1
fi

echo "[INFO] Running Panaroo on ${#GFFS[@]} GFFs"
mkdir -p "$OUTDIR"

# IMPORTANT: pass files space-separated (array expansion), not comma-joined
conda run -n panaroo_env panaroo \
  -i "${GFFS[@]}" \
  -o "$OUTDIR" \
  --clean-mode strict \
  --core_threshold 0.95 \
  --threads "$THREADS"

echo "[OK] Panaroo demo finished -> $OUTDIR"
