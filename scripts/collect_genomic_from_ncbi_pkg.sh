#!/usr/bin/env bash
set -euo pipefail

# ----- CONFIG -----
# Windows/OneDrive source (NCBI Datasets unzipped "data" folder)
PKG_ROOT="/mnt/c/Users/krist/OneDrive/Documents/Bioinformatics/Bioinfo_instructions/StrepA_inclusivity_exclusivity/S_pyogenes_genomes/ncbi_dataset/data"

# Destination inside your repo (Linux side)
DEST="data/genomes_genomic"
LOG="${DEST}/_missing_genomic.log"
# -------------------

mkdir -p "$DEST"
: > "$LOG"

shopt -s nullglob

for dir in "$PKG_ROOT"/GCF_*; do
  [[ -d "$dir" ]] || continue
  acc="$(basename "$dir")"
  fna=""
  # look for genomic FASTA (compressed or not)
  for cand in "$dir"/*_genomic.fna.gz "$dir"/*_genomic.fna; do
    [[ -f "$cand" ]] && { fna="$cand"; break; }
  done

  if [[ -z "$fna" ]]; then
    echo "[MISSING] ${acc} has no *_genomic.fna[.gz]" | tee -a "$LOG"
    continue
  fi

  if [[ "$fna" == *.gz ]]; then
    cp -f "$fna" "$DEST/${acc}.genomic.fna.gz"
  else
    cp -f "$fna" "$DEST/${acc}.genomic.fna"
  fi
  echo "[OK] ${acc} -> ${DEST}/${acc}.genomic.fna*"
done

echo "Done. Destination: $DEST"
echo "Missing log: $LOG"
