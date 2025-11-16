#!/usr/bin/env bash
set -euo pipefail

# -----------------------------
# Defaults (override via env)
# -----------------------------
DATE_TAG="${DATE_TAG:-$(date +%F_%H%M%S)}"

# Conda envs
PANAROO_ENV="${PANAROO_ENV:-panaroo_env}"
ROARY_ENV="${ROARY_ENV:-roary_env}"

# Threads
THREADS="${THREADS:-8}"

# Paths
IN_GFF_DIR="${IN_GFF_DIR:-pipelines/panaroo/1_annotate_prokka}"     # Prokka outputs (one subdir per genome)
PANAROO_OUT="${PANAROO_OUT:-pipelines/panaroo/2_panaroo/core100}"   # Panaroo output folder
SPLIT_DIR="${SPLIT_DIR:-${PANAROO_OUT}/core_gene_alignment.aln.split}"

INCL_DIR="${INCL_DIR:-pipelines/panaroo/3_inclusivity}"
CONS_DIR="${CONS_DIR:-pipelines/panaroo/4_consensus}"
EXCL_DIR="${EXCL_DIR:-pipelines/panaroo/5_exclusivity}"

LOG_DIR="${LOG_DIR:-logs}"
RUN_DIR="${RUN_DIR:-results/runs/run-${DATE_TAG}}"

# Thresholds
CORE_THRESHOLD="${CORE_THRESHOLD:-1.00}"  # fraction (1.00 = 100% presence)
IDENTITY_THR="${IDENTITY_THR:-98.0}"      # ≥98% avg pairwise identity
MIN_PIDENT="${MIN_PIDENT:-85}"            # exclusivity REJECT rule: % identity
MIN_QCOVS="${MIN_QCOVS:-80}"              # exclusivity REJECT rule: % query coverage

# BLAST DB (non-pyogenes panel)
BLAST_DB="${BLAST_DB:-pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/streptococcus_non_pyogenes_db}"

# -----------------------------
# Create dirs & print config
# -----------------------------
mkdir -p "$LOG_DIR" "$RUN_DIR" "$INCL_DIR" "$CONS_DIR" "$EXCL_DIR"

echo "[VARS] DATE_TAG=$DATE_TAG"
echo "[VARS] PANAROO_ENV=$PANAROO_ENV  ROARY_ENV=$ROARY_ENV  THREADS=$THREADS"
echo "[VARS] IN_GFF_DIR=$IN_GFF_DIR"
echo "[VARS] PANAROO_OUT=$PANAROO_OUT"
echo "[VARS] SPLIT_DIR=$SPLIT_DIR"
echo "[VARS] INCL_DIR=$INCL_DIR  CONS_DIR=$CONS_DIR  EXCL_DIR=$EXCL_DIR"
echo "[VARS] LOG_DIR=$LOG_DIR  RUN_DIR=$RUN_DIR"
echo

# Guard critical vars are non-empty
for v in DATE_TAG PANAROO_ENV ROARY_ENV THREADS IN_GFF_DIR PANAROO_OUT SPLIT_DIR INCL_DIR CONS_DIR EXCL_DIR LOG_DIR RUN_DIR; do
  eval "val=\${$v}"
  [[ -n "$val" ]] || { echo "[ERR] Variable $v is empty. Aborting."; exit 1; }
done

# [1/6] Run Panaroo (core=100%)
# -----------------------------
echo "[1/6] Panaroo (core=${CORE_THRESHOLD})…"

# Build a robust list of input GFFs whether flat or nested
shopt -s nullglob
GFFS=( "${IN_GFF_DIR}"/*/*.gff "${IN_GFF_DIR}"/*.gff )
shopt -u nullglob

if ((${#GFFS[@]}==0)); then
  echo "[ERR] No .gff files found in ${IN_GFF_DIR} (flat or subfolders)."
  exit 1
fi
echo "[INFO] Panaroo will process ${#GFFS[@]} GFFs"

# Run Panaroo with explicit env
conda run -n "$PANAROO_ENV" panaroo \
  -i "${GFFS[@]}" \
  -o "$PANAROO_OUT" \
  --clean-mode strict \
  --core_threshold "$CORE_THRESHOLD" \
  --threads "$THREADS" 2>&1 | tee "${LOG_DIR}/panaroo_${DATE_TAG}.log"

# quick existence check
for f in gene_presence_absence.csv gene_data.csv combined_DNA_CDS.fasta; do
  [[ -s "${PANAROO_OUT}/${f}" ]] || { echo "[ERR] Missing ${PANAROO_OUT}/${f}"; exit 1; }
done

# -------------------------------------------------
# [2/6] Split per-gene & align (MAFFT from CDS map)
# -------------------------------------------------
echo "[2/6] Split per-gene & align (MAFFT)…"
conda run -n "$ROARY_ENV" python -u scripts/split_core100_from_panaroo.py \
  --gpa       "${PANAROO_OUT}/gene_presence_absence.csv" \
  --gene_data "${PANAROO_OUT}/gene_data.csv" \
  --cds       "${PANAROO_OUT}/combined_DNA_CDS.fasta" \
  --outdir    "${SPLIT_DIR}" \
  --min_presence "${CORE_THRESHOLD}" \
  --threads   "$THREADS" 2>&1 | tee "${LOG_DIR}/split_${DATE_TAG}.log"

# Count aligned per-gene FASTAs
N_FASTA=$(find "${SPLIT_DIR}" -maxdepth 1 -type f -name "*.fasta" | wc -l || echo 0)
echo "[INFO] Per-gene aligned FASTAs: ${N_FASTA}"
if [[ "${N_FASTA}" -eq 0 ]]; then
  echo "[ERR] No per-gene FASTAs produced in ${SPLIT_DIR}. Aborting."
  exit 1
fi

# --------------------------------------------
# [3/6] Inclusivity: per-gene identity & filter
# --------------------------------------------
echo "[3/6] Inclusivity identities & filter (≥${IDENTITY_THR}%)…"
# 3A identities
conda run -n "$ROARY_ENV" python -u scripts/calculate_identity_with_names.py \
  --split_dir "${SPLIT_DIR}" \
  --output    "${INCL_DIR}/core100_identity.tsv"

# 3B filter
conda run -n "$ROARY_ENV" python -u scripts/filter_high_identity_genes.py \
  --input  "${INCL_DIR}/core100_identity.tsv" \
  --output "${INCL_DIR}/core100_candidates.tsv" \
  --thr    "${IDENTITY_THR}"

# Make a gene list for consensus step
tail -n +2 "${INCL_DIR}/core100_candidates.tsv" > "${INCL_DIR}/core100_candidates.genes.txt" || true

# ----------------------------------------------------
# [4/6] Consensus building & SNP summaries (per gene)
# ----------------------------------------------------
echo "[4/6] Consensus & SNP summaries…"
# Build consensus from split dir, restricted to candidates
conda run -n "$ROARY_ENV" python -u scripts/build_consensus_from_split.py \
  --split_dir "${SPLIT_DIR}" \
  --out_fasta "${CONS_DIR}/core100_consensus.fasta" \
  --genes_txt "${INCL_DIR}/core100_candidates.genes.txt"

# SNP summaries
conda run -n "$ROARY_ENV" python -u scripts/summarize_snps_vs_consensus.py \
  --split_dir       "${SPLIT_DIR}" \
  --consensus_fasta "${CONS_DIR}/core100_consensus.fasta" \
  --out_tsv         "${CONS_DIR}/core100_snps.tsv" \
  --pergene_tsv     "${CONS_DIR}/core100_pergene.tsv"

# -------------------------------------------------------
# [5/6] Exclusivity: BLAST consensus vs exclusion panel
# -------------------------------------------------------
echo "[5/6] Exclusivity BLAST…"
if [[ -e "${BLAST_DB}.nin" || -e "${BLAST_DB}.nhr" ]]; then
  conda run -n "$ROARY_ENV" blastn -task megablast -db "$BLAST_DB" \
    -query "${CONS_DIR}/core100_consensus.fasta" \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs sscinames staxids stitle" \
    -max_target_seqs 5 -evalue 1e-10 -num_threads "$THREADS" \
    > "${EXCL_DIR}/core100_vs_nonpyogenes.tsv"

  conda run -n "$ROARY_ENV" python -u scripts/summarize_blast_exclusivity.py \
    --blast_tsv "${EXCL_DIR}/core100_vs_nonpyogenes.tsv" \
    --out_tsv   "${EXCL_DIR}/core100_exclusivity.tsv" \
    --min_pident "${MIN_PIDENT}" \
    --min_qcovs  "${MIN_QCOVS}"
else
  echo "[WARN] BLAST DB not found at ${BLAST_DB}.* — skipping BLAST & summary."
  : > "${EXCL_DIR}/core100_exclusivity.tsv"
fi

# -----------------------------------------------------------
# [6/6] Extract PASS consensus sequences into per-gene FASTAs
# -----------------------------------------------------------
echo "[6/6] Extract PASS consensus FASTAs…"
conda run -n "$ROARY_ENV" python -u scripts/extract_pass_consensus.py \
  --consensus       "${CONS_DIR}/core100_consensus.fasta" \
  --exclusivity_tsv "${EXCL_DIR}/core100_exclusivity.tsv" \
  --outdir          "${EXCL_DIR}/PASS_FASTAs"

# -----------------------------
# Copy key deliverables to RUN_DIR
# -----------------------------
cp -f "${INCL_DIR}/core100_identity.tsv"            "${RUN_DIR}/"
cp -f "${INCL_DIR}/core100_candidates.tsv"          "${RUN_DIR}/"
cp -f "${CONS_DIR}/core100_consensus.fasta"         "${RUN_DIR}/"
cp -f "${CONS_DIR}/core100_snps.tsv"                "${RUN_DIR}/perstrain_snps.tsv"
cp -f "${CONS_DIR}/core100_pergene.tsv"             "${RUN_DIR}/pergene_snps.tsv" || true
cp -f "${EXCL_DIR}/core100_exclusivity.tsv"         "${RUN_DIR}/exclusivity_summary.tsv" || true

echo
echo "[DONE] Run completed -> ${RUN_DIR}"
echo "       Candidates:   ${INCL_DIR}/core100_candidates.tsv"
echo "       Consensus:    ${CONS_DIR}/core100_consensus.fasta"
echo "       BLAST sum:    ${EXCL_DIR}/core100_exclusivity.tsv"
echo "       PASS FASTAs:  ${EXCL_DIR}/PASS_FASTAs/"
