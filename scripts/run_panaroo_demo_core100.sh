#!/usr/bin/env bash
set -euo pipefail

LIST="${LIST:-demo_accessions.txt}"
OUTDIR="${OUTDIR:-pipelines/panaroo/2_panaroo/demo_core100}"
THREADS="${THREADS:-8}"

# Build array of GFFs from the list
mapfile -t GFFS < <(awk '{printf "pipelines/panaroo/1_annotate_prokka/%s/%s.gff\n",$1,$1}' "$LIST")

# Sanity checks
if [[ ${#GFFS[@]} -lt 2 ]]; then
  echo "Need at least 2 GFFs. Check $LIST and Step 1 outputs." >&2
  exit 1
fi
missing=0
for g in "${GFFS[@]}"; do
  [[ -s "$g" ]] || { echo "[MISSING] $g"; ((missing++)) || true; }
done
(( missing == 0 )) || { echo "Abort: $missing GFF(s) missing."; exit 1; }

echo "[INFO] Running Panaroo with core_threshold=1.0 on ${#GFFS[@]} GFFs"
mkdir -p "$OUTDIR"

# This enforces 100% presence
conda run -n panaroo_env panaroo \
  -i "${GFFS[@]}" \
  -o "$OUTDIR" \
  --clean-mode strict \
  --core_threshold 1.0 \
  --threads "$THREADS"

echo "[OK] Panaroo core=100% done -> $OUTDIR"
