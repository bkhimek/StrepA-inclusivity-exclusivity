#!/usr/bin/env bash
set -euo pipefail

LIST="demo_accessions.txt"
IN_DIR="data/genomes_genomic"
OUT_ROOT="pipelines/panaroo/1_annotate_prokka"
CPUS=${CPUS:-8}

if [[ ! -s "$LIST" ]]; then
  echo "Missing or empty $LIST" >&2
  exit 1
fi

while read -r acc; do
  [[ -z "$acc" ]] && continue
  fasta="${IN_DIR}/${acc}.genomic.fna"
  [[ -f "$fasta" ]] || { echo "[MISS] $fasta not found"; continue; }

  outdir="${OUT_ROOT}/${acc}"
  gff="${outdir}/${acc}.gff"

  if [[ -s "$gff" ]]; then
    echo "[SKIP] ${acc} already annotated -> ${gff}"
    continue
  fi

  echo "[RUN] Prokka for ${acc}"
  conda run -n prokka_env prokka \
    --outdir "$outdir" \
    --prefix "$acc" \
    --cpus "$CPUS" --kingdom Bacteria \
    "$fasta"
done < "$LIST"

echo "[OK] Prokka subset done"
