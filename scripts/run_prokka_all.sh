#!/usr/bin/env bash
set -euo pipefail

CPUS=${CPUS:-8}
IN_DIR="data/genomes_genomic"
OUT_ROOT="pipelines/panaroo/1_annotate_prokka"

shopt -s nullglob
# pick up both .fna and .fna.gz if present in the future
for f in "$IN_DIR"/*genomic.fna "$IN_DIR"/*genomic.fna.gz; do
  [[ -e "$f" ]] || continue
  base="$(basename "$f")"
  base="${base%.gz}"
  base="${base%.fna}"          # e.g., GCF_000006785.2.genomic
  acc="${base%.genomic}"       # e.g., GCF_000006785.2

  OUTDIR="${OUT_ROOT}/${acc}"
  GFF="${OUTDIR}/${acc}.gff"

  if [[ -s "$GFF" ]]; then
    echo "[SKIP] ${acc} already annotated -> ${GFF}"
    continue
  fi

  echo "[RUN] Prokka for ${acc}"
  conda run -n prokka_env prokka \
    --outdir "$OUTDIR" \
    --prefix "$acc" \
    --cpus "$CPUS" --kingdom Bacteria \
    "$f"
done

echo "[OK] Prokka all done"
