#!/usr/bin/env bash
set -euo pipefail
ENV="roary_env"

IN="pipelines/panaroo/5_exclusivity/core100_vs_nonpyogenes.tsv"
OUT="pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv"
MIN_ID="${MIN_ID:-85}"
MIN_QCOV="${MIN_QCOV:-80}"

[[ -s "$IN" ]] || { echo "[ERR] Missing BLAST hits: $IN"; exit 1; }
mkdir -p "$(dirname "$OUT")"

if [[ "${CONDA_DEFAULT_ENV:-}" == "$ENV" ]]; then
  python -u scripts/summarize_blast_exclusivity.py \
    --blast_tsv "$IN" --out_tsv "$OUT" \
    --min_pident "$MIN_ID" --min_qcovs "$MIN_QCOV"
else
  conda run -n "$ENV" python -u scripts/summarize_blast_exclusivity.py \
    --blast_tsv "$IN" --out_tsv "$OUT" \
    --min_pident "$MIN_ID" --min_qcovs "$MIN_QCOV"
fi

echo "[OK] Summary -> $OUT"
