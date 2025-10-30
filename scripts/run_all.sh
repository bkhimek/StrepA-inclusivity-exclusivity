#!/usr/bin/env bash
set -euo pipefail

DATE_TAG="$(date +%F)"
RUN_DIR="results/runs/run-${DATE_TAG}"
IN_GFF_DIR="pipelines/panaroo/1_annotate_prokka"         # your Prokka output root
PANAROO_OUT="pipelines/panaroo/2_panaroo/all"
SPLIT_DIR="${PANAROO_OUT}/core_gene_alignment.aln.split"

mkdir -p "$RUN_DIR" logs

echo "[1/6] Panaroo (all genomes)…"
conda run -n panaroo_env panaroo \
  -i ${IN_GFF_DIR}/*/*.gff \
  -o "${PANAROO_OUT}" \
  --clean-mode strict \
  --core_threshold 1.00 \
  --threads 8 2>&1 | tee "logs/panaroo_${DATE_TAG}.log"

echo "[2/6] Split per-gene & align (MAFFT)…"
conda run -n roary_env python -u scripts/split_core100_from_panaroo.py \
  --gpa   "${PANAROO_OUT}/gene_presence_absence.csv" \
  --cds   "${PANAROO_OUT}/combined_DNA_CDS.fasta" \
  --outdir "${SPLIT_DIR}" \
  --min_presence 1.00 \
  --threads 8 2>&1 | tee "logs/split_${DATE_TAG}.log"

echo "[3/6] Per-gene identity…"
ID_TSV="pipelines/panaroo/3_inclusivity/all_core100_identity.tsv"
conda run -n roary_env python -u scripts/calculate_identity_with_names.py \
  --split_dir "${SPLIT_DIR}" \
  --output    "${ID_TSV}" 2>&1 | tee "logs/identity_${DATE_TAG}.log"

echo "[4/6] Filter ≥98% identity…"
CAND_TSV="pipelines/panaroo/3_inclusivity/all_core100_candidates.tsv"
conda run -n roary_env python -u scripts/filter_high_identity_genes.py \
  --input "${ID_TSV}" --output "${CAND_TSV}" --thr 98.0 2>&1 | tee "logs/filter_${DATE_TAG}.log"

echo "[5/6] Consensus & SNPs…"
conda run -n roary_env bash scripts/run_consensus.sh 2>&1 | tee "logs/consensus_${DATE_TAG}.log"

echo "[6/6] Exclusivity BLAST + summary…"
conda run -n roary_env bash scripts/run_exclusivity.sh 2>&1 | tee "logs/exclusivity_${DATE_TAG}.log"
conda run -n roary_env bash scripts/run_exclusivity_summarize.sh 2>&1 | tee -a "logs/exclusivity_${DATE_TAG}.log"

# Collect compact outputs into RUN_DIR
cp -f pipelines/panaroo/3_inclusivity/all_core100_identity.tsv           "${RUN_DIR}/"
cp -f pipelines/panaroo/3_inclusivity/all_core100_candidates.tsv        "${RUN_DIR}/"
cp -f pipelines/panaroo/4_consensus/core100_pergene.tsv            "${RUN_DIR}/pergene_snps.tsv"
cp -f pipelines/panaroo/4_consensus/core100_snps.tsv               "${RUN_DIR}/perstrain_snps.tsv"
cp -f pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv      "${RUN_DIR}/exclusivity_summary.tsv"

# Update 'latest' symlink
rm -f results/latest && ln -s "runs/run-${DATE_TAG}" results/latest

echo "[OK] Full run complete -> ${RUN_DIR}"
