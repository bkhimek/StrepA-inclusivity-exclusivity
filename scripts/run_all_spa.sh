#!/usr/bin/env bash
set -euo pipefail

# ========= Config (override by exporting before running) =========
DATE_TAG="${DATE_TAG:-$(date +%F_%H%M%S)}"

PANAROO_ENV="${PANAROO_ENV:-panaroo_env}"
ROARY_ENV="${ROARY_ENV:-roary_env}"
THREADS="${THREADS:-8}"

# Prokka GFFs for *S. pyogenes* (391 genomes you just annotated)
IN_GFF_DIR="${IN_GFF_DIR:-pipelines/panaroo/1_annotate_prokka_spa}"

# Panaroo outputs (core 100%)
PANAROO_OUT="${PANAROO_OUT:-pipelines/panaroo/2_panaroo/core100_spa}"
SPLIT_DIR="${SPLIT_DIR:-${PANAROO_OUT}/core_gene_alignment.aln.split}"

# Downstream
INCL_DIR="${INCL_DIR:-pipelines/panaroo/3_inclusivity}"
CONS_DIR="${CONS_DIR:-pipelines/panaroo/4_consensus}"
EXCL_DIR="${EXCL_DIR:-pipelines/panaroo/5_exclusivity}"
LOG_DIR="${LOG_DIR:-logs}"
RUN_DIR="${RUN_DIR:-results/runs/run-${DATE_TAG}}"

# Thresholds
CORE_THRESHOLD="1.00"     # 100% presence
IDENTITY_THR="98.0"       # ≥98% avg pairwise identity
MIN_PIDENT="85"           # exclusivity: reject if (pident >= MIN_PIDENT AND qcovs >= MIN_QCOVS)
MIN_QCOVS="80"

# Non-pyogenes BLAST DB (auto-rebuilds if missing)
EXCL_DB_DIR="${EXCL_DB_DIR:-${EXCL_DIR}/blastdb/non_pyogenes}"
EXCL_DB="${EXCL_DB_DIR}/streptococcus_non_pyogenes_db"
EXCL_COMBINED="${EXCL_DB_DIR}/streptococcus_non_pyogenes_combined.fna"

# Optional: copy outputs to OneDrive (set ONEDRIVE_OUT before running)
ONEDRIVE_OUT="${ONEDRIVE_OUT:-}"

# ========= Prep =========
mkdir -p "$PANAROO_OUT" "$SPLIT_DIR" "$INCL_DIR" "$CONS_DIR" "$EXCL_DIR" "$EXCL_DB_DIR" "$LOG_DIR" "$RUN_DIR"

echo "[VARS] DATE_TAG=$DATE_TAG"
echo "[VARS] PANAROO_ENV=$PANAROO_ENV  ROARY_ENV=$ROARY_ENV  THREADS=$THREADS"
echo "[VARS] IN_GFF_DIR=$IN_GFF_DIR"
echo "[VARS] PANAROO_OUT=$PANAROO_OUT"
echo "[VARS] SPLIT_DIR=$SPLIT_DIR"
echo "[VARS] INCL_DIR=$INCL_DIR  CONS_DIR=$CONS_DIR  EXCL_DIR=$EXCL_DIR"
echo "[VARS] LOG_DIR=$LOG_DIR  RUN_DIR=$RUN_DIR"
echo

# ========= 1) Panaroo (core=100%) =========
echo "[1/8] Panaroo (core=${CORE_THRESHOLD})…"

# Build a glob that matches your layout (subfolders with GFFs).
GLOB_A="${IN_GFF_DIR}/*.gff"
GLOB_B="${IN_GFF_DIR}/*/*.gff"

if compgen -G "$GLOB_A" > /dev/null; then
  GFF_GLOB="$GLOB_A"
elif compgen -G "$GLOB_B" > /dev/null; then
  GFF_GLOB="$GLOB_B"
else
  echo "[ERR] No .gff files found under ${IN_GFF_DIR}"
  exit 1
fi

conda run -n "$PANAROO_ENV" panaroo \
  -i $GFF_GLOB \
  -o "$PANAROO_OUT" \
  --clean-mode strict \
  --core_threshold "$CORE_THRESHOLD" \
  --threads "$THREADS" 2>&1 | tee "${LOG_DIR}/panaroo_${DATE_TAG}.log"

for f in gene_presence_absence.csv gene_data.csv combined_DNA_CDS.fasta; do
  [[ -s "${PANAROO_OUT}/${f}" ]] || { echo "[ERR] Missing ${PANAROO_OUT}/${f}"; exit 1; }
done

# ========= 2) Split per-gene + MAFFT align =========
echo "[2/8] Split per-gene & align (MAFFT)…"

conda run -n "$ROARY_ENV" python -u scripts/split_core100_from_panaroo.py \
  --gpa       "${PANAROO_OUT}/gene_presence_absence.csv" \
  --gene_data "${PANAROO_OUT}/gene_data.csv" \
  --cds       "${PANAROO_OUT}/combined_DNA_CDS.fasta" \
  --outdir    "$SPLIT_DIR" \
  --min_presence "$CORE_THRESHOLD" \
  --threads   "$THREADS" 2>&1 | tee "${LOG_DIR}/split_${DATE_TAG}.log"

# ========= 3) Per-gene identity =========
echo "[3/8] Compute per-gene identity…"
conda run -n "$ROARY_ENV" python -u scripts/calculate_identity_with_names.py \
  --split_dir "$SPLIT_DIR" \
  --output    "${INCL_DIR}/core100_identity.tsv"

[[ -s "${INCL_DIR}/core100_identity.tsv" ]] || { echo "[ERR] identity table empty"; exit 1; }

# ========= 4) Filter ≥98% & dedupe .raw =========
echo "[4/8] Filter ≥${IDENTITY_THR}% and dedupe *.raw…"

conda run -n "$ROARY_ENV" python -u scripts/filter_inclusivity_candidates.py \
  --input  "${INCL_DIR}/core100_identity.tsv" \
  --output "${INCL_DIR}/core100_candidates.tsv" \
  --threshold "${IDENTITY_THR}"

# Drop .raw duplicates; keep first column only when we just need gene list later
awk 'NR==1 || $1 !~ /\.raw$/' "${INCL_DIR}/core100_candidates.tsv" > "${INCL_DIR}/core100_candidates.nodup.tsv"

# ========= 5) Build consensus FASTA =========
echo "[5/8] Build consensus FASTA from split alignments…"
conda run -n "$ROARY_ENV" python -u scripts/build_consensus_from_split.py \
  --split_dir "$SPLIT_DIR" \
  --genes_txt "${INCL_DIR}/core100_candidates.nodup.tsv" \
  --out_fasta "${CONS_DIR}/core100_consensus.fasta"

[[ -s "${CONS_DIR}/core100_consensus.fasta" ]] || { echo "[ERR] consensus FASTA empty"; exit 1; }

# ========= 6) BLAST vs non-pyogenes DB =========
echo "[6/8] BLAST exclusivity…"
# Build DB if missing
if [[ ! -s "${EXCL_DB}.nsq" && ! -s "${EXCL_DB}.00.nsq" ]]; then
  echo "[INFO] BLAST DB missing; (re)building…"
  if [[ ! -s "$EXCL_COMBINED" ]]; then
    # Try to resurrect from known backup location if available
    BK=$(find backups -type f -name "streptococcus_non_pyogenes_combined.fna" | head -n 1 || true)
    if [[ -n "${BK:-}" ]]; then
      cp -f "$BK" "$EXCL_COMBINED"
    else
      echo "[ERR] Could not find ${EXCL_COMBINED} to make BLAST DB."
      exit 1
    fi
  fi
  conda run -n "$ROARY_ENV" makeblastdb \
    -in "$EXCL_COMBINED" \
    -dbtype nucl \
    -out "$EXCL_DB" \
    -parse_seqids
fi

conda run -n "$ROARY_ENV" blastn -task megablast -db "$EXCL_DB" \
  -query "${CONS_DIR}/core100_consensus.fasta" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs sscinames staxids stitle" \
  -max_target_seqs 5 -evalue 1e-10 -num_threads "$THREADS" \
  > "${EXCL_DIR}/core100_vs_nonpyogenes.tsv"

# ========= 7) Summarize exclusivity & extract PASS FASTAs =========
echo "[7/8] Summarize exclusivity & extract PASS FASTAs…"
conda run -n "$ROARY_ENV" python -u scripts/summarize_blast_exclusivity.py \
  --blast_tsv "${EXCL_DIR}/core100_vs_nonpyogenes.tsv" \
  --out_tsv   "${EXCL_DIR}/core100_exclusivity.tsv" \
  --min_pident "${MIN_PIDENT}" \
  --min_qcovs  "${MIN_QCOVS}"

conda run -n "$ROARY_ENV" python -u scripts/extract_pass_consensus.py \
  --consensus      "${CONS_DIR}/core100_consensus.fasta" \
  --exclusivity_tsv "${EXCL_DIR}/core100_exclusivity.tsv" \
  --outdir         "${EXCL_DIR}/PASS_FASTAs"

# ========= 8) Collect & (optional) copy to OneDrive =========
echo "[8/8] Collecting artifacts…"
cp -f "${INCL_DIR}/core100_identity.tsv"                  "${RUN_DIR}/"
cp -f "${INCL_DIR}/core100_candidates.tsv"               "${RUN_DIR}/"
cp -f "${INCL_DIR}/core100_candidates.nodup.tsv"         "${RUN_DIR}/"
cp -f "${CONS_DIR}/core100_consensus.fasta"              "${RUN_DIR}/"
cp -f "${EXCL_DIR}/core100_vs_nonpyogenes.tsv"           "${RUN_DIR}/"
cp -f "${EXCL_DIR}/core100_exclusivity.tsv"              "${RUN_DIR}/"
mkdir -p "${RUN_DIR}/PASS_FASTAs"
cp -f ${EXCL_DIR}/PASS_FASTAs/*.fasta "${RUN_DIR}/PASS_FASTAs/" 2>/dev/null || true

if [[ -n "${ONEDRIVE_OUT}" ]]; then
  mkdir -p "${ONEDRIVE_OUT}"
  cp -f "${RUN_DIR}/core100_exclusivity.tsv" "${ONEDRIVE_OUT}/"
  mkdir -p "${ONEDRIVE_OUT}/PASS_FASTAs"
  cp -f "${RUN_DIR}/PASS_FASTAs/"*.fasta "${ONEDRIVE_OUT}/PASS_FASTAs/" 2>/dev/null || true
fi

echo
echo "[DONE] Run completed -> ${RUN_DIR}"
echo "       Candidates:   ${INCL_DIR}/core100_candidates.tsv"
echo "       Consensus:    ${CONS_DIR}/core100_consensus.fasta"
echo "       BLAST sum:    ${EXCL_DIR}/core100_exclusivity.tsv"
echo "       PASS FASTAs:  ${EXCL_DIR}/PASS_FASTAs/"
