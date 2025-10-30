# scripts/run_inclusivity.sh
#!/usr/bin/env bash
set -euo pipefail

ENV="roary_env"
ALIGN="pipelines/panaroo/2_panaroo/demo/core_gene_alignment.aln"
OUT="pipelines/panaroo/3_inclusivity/demo_core_identity.tsv"

echo "[INFO] Inclusivity (identity) on: $ALIGN"
echo "[INFO] Output -> $OUT"

# If the requested env is already active, run python directly.
if [[ "${CONDA_DEFAULT_ENV:-}" == "$ENV" ]]; then
  echo "[INFO] Using already-active conda env: $ENV"
  python scripts/calculate_identity_with_names.py \
    --alignment "$ALIGN" \
    --output "$OUT"
else
  # Otherwise, do not 'conda activate' in a non-interactive script; use conda run
  echo "[INFO] Using conda run with env: $ENV"
  conda run -n "$ENV" python scripts/calculate_identity_with_names.py \
    --alignment "$ALIGN" \
    --output "$OUT"
fi

echo "[OK] Wrote $OUT"
