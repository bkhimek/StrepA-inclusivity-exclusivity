#!/usr/bin/env bash
set -euo pipefail

ENV="roary_env"

CONS="pipelines/panaroo/4_consensus/demo_core99_consensus.fasta"
EXCL="pipelines/panaroo/5_exclusivity/demo_core99_exclusivity.tsv"
OUTDIR="pipelines/panaroo/6_reports/pass_consensus_split"

mkdir -p "$OUTDIR"

echo "[INFO] Extracting PASS consensus genes..."
if [[ "${CONDA_DEFAULT_ENV:-}" == "$ENV" ]]; then
  python -u scripts/extract_pass_consensus.py \
    --consensus_fasta "$CONS" \
    --exclusivity_tsv "$EXCL" \
    --outdir "$OUTDIR"
else
  conda run -n "$ENV" python -u scripts/extract_pass_consensus.py \
    --consensus_fasta "$CONS" \
    --exclusivity_tsv "$EXCL" \
    --outdir "$OUTDIR"
fi

echo "[OK] PASS consensus FASTAs are in: $OUTDIR"
ls -lh "$OUTDIR" | head
