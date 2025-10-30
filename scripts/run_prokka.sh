#!/usr/bin/env bash
set -euo pipefail

# Prokka via conda run (no need to activate in the script)
# Adjust --cpus if you want fewer/more
CPUS=${CPUS:-8}

# Genome 1
conda run -n prokka_env prokka \
  --outdir pipelines/panaroo/1_annotate_prokka/GCF_000006785.2 \
  --prefix GCF_000006785.2 \
  --cpus $CPUS --kingdom Bacteria \
  data/genomes_genomic/GCF_000006785.2.genomic.fna

# Genome 2
conda run -n prokka_env prokka \
  --outdir pipelines/panaroo/1_annotate_prokka/GCF_000007285.1 \
  --prefix GCF_000007285.1 \
  --cpus $CPUS --kingdom Bacteria \
  data/genomes_genomic/GCF_000007285.1.genomic.fna

echo "[OK] Prokka demo finished"
