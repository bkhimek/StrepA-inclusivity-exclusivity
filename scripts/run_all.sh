#!/usr/bin/env bash
set -euo pipefail

# === Config ===
DATE_TAG="$(date +%F)"
RUN_DIR="results/runs/run-${DATE_TAG}"
LOG_DIR="logs"
THREADS="${THREADS:-8}"

# Envs
PANAROO_ENV="panaroo_env"
ANALYSIS_ENV="roary_env"   # Biopython, pandas, mafft, blastn in PATH

# Paths
IN_GFF_DIR="pipelines/panaroo/1_annotate_prokka"
PANAROO_OUT="pipelines/panaroo/2_panaroo/core100"
SPLIT_DIR="${PANAROO_OUT}/core_gene_alignment.aln.split"

INCL_DIR="pipelines/panaroo/3_inclusivity"
CONS_DIR="pipelines/panaroo/4_consensus"
EXCL_DIR="pipelines/panaroo/5_exclusivity"
PASS_DIR="pipelines/panaroo/6_reports/pass_consensus_fasta"  # optional

mkdir -p "$RUN_DIR" "$LOG_DIR" "$INCL_DIR" "$CONS_DIR" "$EXCL_DIR"

# Thresholds
CORE_THRESHOLD="1.00"   # 100% presence
IDENTITY_THR="98.0"     # ≥98% pairwise identity for inclusivity
MIN_PIDENT="85"         # exclusivity reject rule (pident)
MIN_QCOVS="80"          # exclusivity reject rule (qcovs)

# BLAST DB to non-pyogenes (adjust if yours differs)
BLAST_DB="pipelines/panaroo/5_exclusivity/blastdb/streptococcus_non_pyogenes_db"

echo "[1/6] Panaroo (core=${CORE_THRESHOLD})…"
conda run -n "$PANAROO_ENV" panaroo \
  -i ${IN_GFF_DIR}/*/*.gff \
  -o "$PANAROO_OUT" \
  --clean-mode strict \
  --core_threshold "$CORE_THRESHOLD" \
  --threads "$THREADS" 2>&1 | tee "${LOG_DIR}/panaroo_${DATE_TAG}.log"

echo "[2/6] Split per-gene & align (MAFFT)…"
conda run -n "$ANALYSIS_ENV" python -u scripts/split_core100_from_panaroo.py \
  --gpa    "${PANAROO_OUT}/gene_presence_absence.csv" \
  --cds    "${PANAROO_OUT}/combined_DNA_CDS.fasta" \
  --outdir "${SPLIT_DIR}" \
  --min_presence "$CORE_THRESHOLD" \
  --threads "$THREADS" 2>&1 | tee "${LOG_DIR}/split_${DATE_TAG}.log"

echo "[3/6] Per-gene identity (from split alignments)…"
conda run -n "$ANALYSIS_ENV" python -u scripts/calculate_identity_with_names.py \
  --split_dir "${SPLIT_DIR}" \
  --output    "${INCL_DIR}/core100_identity.tsv"

echo "[4/6] Filter ≥${IDENTITY_THR}% identity…"
conda run -n "$ANALYSIS_ENV" python -u scripts/filter_high_identity_genes.py \
  --input  "${INCL_DIR}/core100_identity.tsv" \
  --output "${INCL_DIR}/core100_candidates.tsv" \
  --thr    "$IDENTITY_THR"

echo "[5/6] Consensus & SNPs (per-strain + per-gene)…"
conda run -n "$ANALYSIS_ENV" python -u scripts/build_consensus_from_split.py \
  --candidates "${INCL_DIR}/core100_candidates.tsv" \
  --split_dir  "${SPLIT_DIR}" \
  --out_fasta  "${CONS_DIR}/core100_consensus.fasta"

conda run -n "$ANALYSIS_ENV" python -u scripts/summarize_snps_vs_consensus.py \
  --consensus  "${CONS_DIR}/core100_consensus.fasta" \
  --split_dir  "${SPLIT_DIR}" \
  --out_perstrain "${CONS_DIR}/core100_snps.tsv" \
  --out_pergene   "${CONS_DIR}/core100_pergene.tsv"

echo "[6/6] Exclusivity: BLAST vs non-pyogenes + summarize…"
# Raw BLAST hits
conda run -n "$ANALYSIS_ENV" bash -c "
  mkdir -p '${EXCL_DIR}'
  blastn -query '${CONS_DIR}/core100_consensus.fasta' \
         -db '${BLAST_DB}' \
         -out '${EXCL_DIR}/core100_vs_nonpyogenes.tsv' \
         -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs sscinames staxids stitle' \
         -max_target_seqs 5 -num_threads '${THREADS}'
"

# Summary with neighbor names + PASS/REJECT
conda run -n "$ANALYSIS_ENV" python -u scripts/summarize_blast_exclusivity.py \
  --blast_tsv "${EXCL_DIR}/core100_vs_nonpyogenes.tsv" \
  --out_tsv   "${EXCL_DIR}/core100_exclusivity.tsv" \
  --min_pident "${MIN_PIDENT}" \
  --min_qcovs  "${MIN_QCOVS}"

# (Optional) Extract individual FASTAs for PASS genes
conda run -n "$ANALYSIS_ENV" python -u scripts/extract_pass_consensus.py \
  --consensus_fasta "${CONS_DIR}/core100_consensus.fasta" \
  --exclusivity_tsv "${EXCL_DIR}/core100_exclusivity.tsv" \
  --out_dir         "${PASS_DIR}"

# Stash key outputs for this run
cp -f "${INCL_DIR}/core100_identity.tsv"       "${RUN_DIR}/"
cp -f "${INCL_DIR}/core100_candidates.tsv"     "${RUN_DIR}/"
cp -f "${CONS_DIR}/core100_consensus.fasta"    "${RUN_DIR}/"
cp -f "${CONS_DIR}/core100_pergene.tsv"        "${RUN_DIR}/pergene_snps.tsv"
cp -f "${CONS_DIR}/core100_snps.tsv"           "${RUN_DIR}/perstrain_snps.tsv"
cp -f "${EXCL_DIR}/core100_exclusivity.tsv"    "${RUN_DIR}/exclusivity_summary.tsv"

echo "[DONE] Run folder: ${RUN_DIR}"
