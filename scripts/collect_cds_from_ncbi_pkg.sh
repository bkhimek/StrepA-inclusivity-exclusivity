#!/usr/bin/env bash
set -euo pipefail

# ----- CONFIG -----
# Windows/OneDrive source (NCBI Datasets unzipped "data" folder)
PKG_ROOT="/mnt/c/Users/krist/OneDrive/Documents/Bioinformatics/StrepA_inclusivity_exclusivity/S_pyogenes_genomes/ncbi_dataset/data"

# Destination inside your repo (Linux side)
DEST="data/genomes_cds"
LOG="${DEST}/_missing_cds.log"

# Dry run: set to 1 to print actions only; set to 0 to actually copy
DRYRUN=${DRYRUN:-0}
# -------------------

mkdir -p "$DEST"
: > "$LOG"

shopt -s nullglob

found_any=0
for dir in "$PKG_ROOT"/GCF_*; do
  [[ -d "$dir" ]] || continue
  found_any=1
  acc="$(basename "$dir")"
  cds=""
  # common NCBI Datasets filenames
  for cand in "$dir"/cds_from_genomic.fna.gz "$dir"/cds_from_genomic.fna; do
    [[ -f "$cand" ]] && { cds="$cand"; break; }
  done

  if [[ -z "$cds" ]]; then
    echo "[MISSING] ${acc} has no cds_from_genomic.fna[.gz]" | tee -a "$LOG"
    continue
  fi

  if [[ "$cds" == *.gz ]]; then
    out="${DEST}/${acc}.cds.fna.gz"
  else
    out="${DEST}/${acc}.cds.fna"
  fi

  echo "[OK] ${acc} -> ${out}"
  if [[ "$DRYRUN" -eq 0 ]]; then
    cp -f "$cds" "$out"
  fi
done

if [[ "$found_any" -eq 0 ]]; then
  echo "No GCF_* folders found under: $PKG_ROOT" >&2
  exit 1
fi

echo "Done. Destination: $DEST"
echo "Missing log: $LOG"
