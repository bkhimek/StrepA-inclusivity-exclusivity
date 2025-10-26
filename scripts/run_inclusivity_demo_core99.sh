#!/usr/bin/env bash
set -euo pipefail
ENV="roary_env"
ALIGN="pipelines/panaroo/2_panaroo/demo_core99/core_gene_alignment.aln"
OUT="pipelines/panaroo/3_inclusivity/demo_core99_identity.tsv"

echo "[INFO] Inclusivity (identity) on: $ALIGN"
echo "[INFO] Output -> $OUT"
mkdir -p "$(dirname "$OUT")"

if [[ "${CONDA_DEFAULT_ENV:-}" == "$ENV" ]]; then
  python -u scripts/calculate_identity_with_names.py --alignment "$ALIGN" --output "$OUT"
else
  conda run -n "$ENV" python -u scripts/calculate_identity_with_names.py --alignment "$ALIGN" --output "$OUT"
fi

echo "[OK] Wrote $OUT"
