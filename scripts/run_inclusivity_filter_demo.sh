# scripts/run_inclusivity_filter_demo.sh
#!/usr/bin/env bash
set -euo pipefail

ENV="roary_env"
IN="pipelines/panaroo/3_inclusivity/demo_core_identity.tsv"
OUT="pipelines/panaroo/3_inclusivity/demo_inclusivity_candidates.tsv"
THR="${THR:-0.98}"

echo "[INFO] Filtering inclusivity candidates at threshold ${THR}"
echo "[INFO] Input  -> $IN"
echo "[INFO] Output -> $OUT"

if [[ ! -s "$IN" ]]; then
  echo "[ERR] Input file not found or empty: $IN" >&2
  exit 1
fi

if [[ "${CONDA_DEFAULT_ENV:-}" == "$ENV" ]]; then
  python -u scripts/filter_high_identity_genes.py --input "$IN" --thr "$THR" --output "$OUT"
else
  conda run -n "$ENV" python -u scripts/filter_high_identity_genes.py --input "$IN" --thr "$THR" --output "$OUT"
fi

echo "[OK] Wrote $OUT"
