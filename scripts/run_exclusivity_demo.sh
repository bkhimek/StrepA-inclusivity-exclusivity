#!/usr/bin/env bash
set -euo pipefail

ENV="roary_env"

CAND_FASTA="pipelines/panaroo/4_consensus/demo_core99_consensus.fasta"
DB="databases/streptococcus_non_pyogenes_db"

OUT_DIR="pipelines/panaroo/5_exclusivity"
mkdir -p "$OUT_DIR"

BLAST_OUT="${OUT_DIR}/demo_core99_vs_nonpyogenes.tsv"

[[ -s "$CAND_FASTA" ]] || { echo "[ERR] Missing candidate FASTA: $CAND_FASTA"; exit 1; }
for ext in nin nsq nhr; do
  [[ -f "${DB}.${ext}" ]] || { echo "[ERR] Missing BLAST DB component: ${DB}.${ext}"; exit 1; }
done

# Include qlen, qcovs, and stitle (subject title)
OUTFMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs sscinames staxids stitle"

if [[ "${CONDA_DEFAULT_ENV:-}" == "$ENV" ]]; then
  blastn -task blastn -query "$CAND_FASTA" -db "$DB" \
         -evalue 1e-10 -word_size 11 -dust yes \
         -max_target_seqs 25 -outfmt "$OUTFMT" > "$BLAST_OUT"
else
  conda run -n "$ENV" blastn -task blastn -query "$CAND_FASTA" -db "$DB" \
         -evalue 1e-10 -word_size 11 -dust yes \
         -max_target_seqs 25 -outfmt "$OUTFMT" > "$BLAST_OUT"
fi

echo "[OK] BLAST exclusivity search complete -> $BLAST_OUT"
