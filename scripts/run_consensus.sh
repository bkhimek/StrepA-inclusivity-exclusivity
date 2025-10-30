#!/usr/bin/env bash
set -euo pipefail

# ---------- Config ----------
ENV="roary_env"

# Use the core=99% run + split per-gene alignments
SPLIT_DIR="pipelines/panaroo/2_panaroo/core100/core_gene_alignment.aln.split"

# Use the inclusivity candidates from Step 3 (core100 + ≥98% identity)
CAND_TSV="pipelines/panaroo/3_inclusivity/core100_candidates.tsv"

# Outputs
OUT_DIR="pipelines/panaroo/4_consensus"
CONS_OUT="${OUT_DIR}/core100_consensus.fasta"
SNP_OUT="${OUT_DIR}/core100_snps.tsv"
PERGENE_OUT="${OUT_DIR}/core100_pergene.tsv"

# ---------- Checks ----------
[[ -d "$SPLIT_DIR" ]] || { echo "[ERR] Missing split dir: $SPLIT_DIR"; exit 1; }
[[ -s "$CAND_TSV" ]]  || { echo "[ERR] Missing candidates TSV: $CAND_TSV"; exit 1; }
mkdir -p "$OUT_DIR"

echo "[INFO] Consensus -> $CONS_OUT"
echo "[INFO] SNPs (per-strain) -> $SNP_OUT"
echo "[INFO] SNPs (per-gene)   -> $PERGENE_OUT"

# ---------- Step 4A: build consensus for candidate genes ----------
if [[ "${CONDA_DEFAULT_ENV:-}" == "$ENV" ]]; then
  python -u scripts/build_consensus_from_split.py \
    --split_dir "$SPLIT_DIR" \
    --out_fasta "$CONS_OUT" \
    --genes_txt "$CAND_TSV"
else
  conda run -n "$ENV" python -u scripts/build_consensus_from_split.py \
    --split_dir "$SPLIT_DIR" \
    --out_fasta "$CONS_OUT" \
    --genes_txt "$CAND_TSV"
fi

# ---------- Step 4B: SNP summary vs consensus ----------
if [[ "${CONDA_DEFAULT_ENV:-}" == "$ENV" ]]; then
  python -u scripts/summarize_snps_vs_consensus.py \
    --split_dir "$SPLIT_DIR" \
    --consensus_fasta "$CONS_OUT" \
    --out_tsv "$SNP_OUT" \
    --pergene_tsv "$PERGENE_OUT"
else
  conda run -n "$ENV" python -u scripts/summarize_snps_vs_consensus.py \
    --split_dir "$SPLIT_DIR" \
    --consensus_fasta "$CONS_OUT" \
    --out_tsv "$SNP_OUT" \
    --pergene_tsv "$PERGENE_OUT"
fi

echo "[OK] Wrote:"
ls -lh "$CONS_OUT" "$SNP_OUT" "$PERGENE_OUT"
