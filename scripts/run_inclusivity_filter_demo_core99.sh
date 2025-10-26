#!/usr/bin/env bash
set -euo pipefail
ENV="roary_env"
IN="pipelines/panaroo/3_inclusivity/demo_core99_identity.tsv"
OUT="pipelines/panaroo/3_inclusivity/demo_core99_candidates.tsv"
THR="${THR:-98.0}"   # percent identity cutoff

echo "[INFO] Filtering core=99% at ≥${THR}% identity"
echo "[INFO] Input  -> $IN"
echo "[INFO] Output -> $OUT"
mkdir -p "$(dirname "$OUT")"

[[ -s "$IN" ]] || { echo "[ERR] Missing input: $IN"; exit 1; }

if [[ "${CONDA_DEFAULT_ENV:-}" == "$ENV" ]]; then
  python -u scripts/filter_high_identity_genes.py --input "$IN" --output "$OUT" --thr "$THR"
else
  conda run -n "$ENV" python -u scripts/filter_high_identity_genes.py --input "$IN" --output "$OUT" --thr "$THR"
fi

echo "[OK] Wrote $OUT"
