#!/usr/bin/env bash
set -euo pipefail

ENV="roary_env"
IN="pipelines/panaroo/3_inclusivity/demo_core100_identity.tsv"
OUT="pipelines/panaroo/3_inclusivity/demo_core100_candidates.txt"
THR="${THR:-98.0}"   # percent identity cutoff

echo "[INFO] Filtering at ≥${THR}% identity"
echo "[INFO] Input  -> $IN"
echo "[INFO] Output -> $OUT"
mkdir -p "$(dirname "$OUT")"

if [[ ! -s "$IN" ]]; then
  echo "[ERR] Input file missing or empty: $IN" >&2
  echo "      Run scripts/run_inclusivity_demo_core100.sh first." >&2
  exit 1
fi

if [[ "${CONDA_DEFAULT_ENV:-}" == "$ENV" ]]; then
  python -u scripts/filter_high_identity_genes.py --input "$IN" --output "$OUT" --thr "$THR"
else
  conda run -n "$ENV" python -u scripts/filter_high_identity_genes.py --input "$IN" --output "$OUT" --thr "$THR"
fi

echo "[OK] Wrote $OUT"
