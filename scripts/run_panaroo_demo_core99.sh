#!/usr/bin/env bash
set -euo pipefail
LIST="${LIST:-demo_accessions.txt}"
OUTDIR="${OUTDIR:-pipelines/panaroo/2_panaroo/demo_core99}"
THREADS="${THREADS:-8}"

# Build list of GFFs from Step 1
mapfile -t GFFS < <(awk '{printf "pipelines/panaroo/1_annotate_prokka/%s/%s.gff\n",$1,$1}' "$LIST")

# Sanity checks
missing=0; for g in "${GFFS[@]}"; do [[ -s "$g" ]] || { echo "[MISSING] $g"; ((missing++))||true; }; done
(( missing == 0 )) || { echo "Abort: $missing GFF(s) missing."; exit 1; }

echo "[INFO] Panaroo core=0.99 with alignment=core, aligner=mafft on ${#GFFS[@]} GFFs"
mkdir -p "$OUTDIR"

conda run -n panaroo_env panaroo \
  -i "${GFFS[@]}" \
  -o "$OUTDIR" \
  --clean-mode strict \
  --core_threshold 0.99 \
  --alignment core \
  --aligner mafft \
  --threads "$THREADS"

echo "[OK] Panaroo done -> $OUTDIR"
