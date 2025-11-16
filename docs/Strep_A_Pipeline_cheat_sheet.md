Strep A pipeline cheat-sheet
0) Environments you’ll use
# check envs
conda env list

# activate the main one (has Biopython, pandas, BLAST+ etc.)
conda activate roary_env

# panaroo-only pieces are run through conda run -n panaroo_env inside the wrapper

1) One-liner run (Strep A, all genomes)
# GFFs produced by Prokka (one subfolder per genome):
ls pipelines/panaroo/1_annotate_prokka_spa/*/*.gff | wc -l  # expect ~391

# Run end-to-end (core=100%, identity≥98%, BLAST vs non-pyogenes DB)
bash scripts/run_all_spa.sh

What the wrapper does (step-by-step)

Panaroo (core=100%): builds pangenome

input: pipelines/panaroo/1_annotate_prokka_spa/*/*.gff

output dir: pipelines/panaroo/2_panaroo/core100_spa/

key files:

gene_presence_absence.csv

combined_DNA_CDS.fasta

gene_data.csv

Split per-gene & MAFFT align (only genes present in 100%):

helper: scripts/split_core100_from_panaroo.py

output folder: pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split/*.fasta

Per-gene identity (avg pairwise %ID per alignment):

helper: scripts/calculate_identity_with_names.py

output: pipelines/panaroo/3_inclusivity/core100_identity.tsv

Filter inclusivity (≥98% avg identity):

helper: scripts/filter_inclusivity_candidates.py

output: pipelines/panaroo/3_inclusivity/core100_candidates.tsv

wrapper also writes a de-duplicated list: core100_candidates.nodup.tsv (drops “*.raw” dupes)

Consensus FASTA for PASS (inclusive) genes & SNP summaries:

helpers:

scripts/build_consensus_from_split.py → pipelines/panaroo/4_consensus/core100_consensus.fasta

scripts/summarize_snps_vs_consensus.py → pipelines/panaroo/4_consensus/{core100_snps.tsv, core100_pergene.tsv}

Exclusivity BLAST vs near-neighbors (non-pyogenes) and summary:

DB path (editable in wrapper):
pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/streptococcus_non_pyogenes_db

run: blastn → core100_vs_nonpyogenes.tsv

summarize: scripts/summarize_blast_exclusivity.py → core100_exclusivity.tsv (PASS/REJECT + best neighbor)

Extract PASS consensus FASTAs (only markers that passed exclusivity gate):

helper: scripts/extract_pass_consensus.py

output folder: pipelines/panaroo/5_exclusivity/PASS_FASTAs/ (one FASTA per gene)

Run folder snapshot:

copies a few key artifacts into results/runs/run-YYYY-MM-DD_HHMMSS/

2) Key inputs & outputs

Inputs

pipelines/panaroo/1_annotate_prokka_spa/*/*.gff — Prokka GFFs for all S. pyogenes genomes.

BLAST DB at pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/streptococcus_non_pyogenes_db.*
(If missing: rebuild with makeblastdb—see Troubleshooting.)

Primary outputs

pipelines/panaroo/3_inclusivity/core100_identity.tsv

pipelines/panaroo/3_inclusivity/core100_candidates.tsv (+ .nodup.tsv)

pipelines/panaroo/4_consensus/core100_consensus.fasta

pipelines/panaroo/4_consensus/core100_snps.tsv, core100_pergene.tsv

pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv (PASS/REJECT)

pipelines/panaroo/5_exclusivity/PASS_FASTAs/ (one FASTA per PASS gene)

results/runs/run-*/ snapshot

3) Sanity checks (quick)
# GFF count
find pipelines/panaroo/1_annotate_prokka_spa -name '*.gff' | wc -l

# split FASTAs present?
find pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split -name '*.fasta' | wc -l

# identity table should be > header-only
wc -l pipelines/panaroo/3_inclusivity/core100_identity.tsv

# candidates should be non-empty (and nodup file should not contain *.raw rows)
head pipelines/panaroo/3_inclusivity/core100_candidates.tsv
grep -c '\.raw$' pipelines/panaroo/3_inclusivity/core100_candidates.nodup.tsv

# consensus contains 1 header per PASS gene
grep -c '^>' pipelines/panaroo/4_consensus/core100_consensus.fasta

# exclusivity summary counts
awk -F'\t' 'NR>1{c[$5]++} END{ for(k in c)print k":",c[k] }' pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv

# PASS FASTAs count
find pipelines/panaroo/5_exclusivity/PASS_FASTAs -name '*.fasta' | wc -l

4) Common fixes

A. “Database memory map file error” (BLAST)

# Rebuild DB from combined FASTA (adjust path if your combined .fna lives elsewhere)
DBDIR="pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes"
mkdir -p "$DBDIR"
makeblastdb -in "$DBDIR/streptococcus_non_pyogenes_combined.fna" \
  -dbtype nucl -out "$DBDIR/streptococcus_non_pyogenes_db" -parse_seqids

# check:
blastdbcmd -db "$DBDIR/streptococcus_non_pyogenes_db" -info | head


B. Empty identity table (only header)
Usually means zero split FASTAs (earlier step failed). Re-run split step or the full wrapper.

C. Duplicate “*.raw” genes
The wrapper creates core100_candidates.nodup.tsv and uses that downstream. If running manually, use the nodup file.

D. Panaroo “Error reading prokka input!”
Your -i glob doesn’t match. For a flat folder of GFFs use:

panaroo -i pipelines/panaroo/1_annotate_prokka_spa/*.gff ...


For a nested layout (one subdir per genome) use:

panaroo -i pipelines/panaroo/1_annotate_prokka_spa/*/*.gff ...


Make sure the wrapper’s IN_GFF_DIR matches your structure.

5) Handy partial re-runs
# Only split + identity + candidates (reuse existing panaroo outputs)
python scripts/split_core100_from_panaroo.py \
  --gpa pipelines/panaroo/2_panaroo/core100_spa/gene_presence_absence.csv \
  --gene_data pipelines/panaroo/2_panaroo/core100_spa/gene_data.csv \
  --cds pipelines/panaroo/2_panaroo/core100_spa/combined_DNA_CDS.fasta \
  --outdir pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split \
  --min_presence 1.00 --threads 8

python scripts/calculate_identity_with_names.py \
  --split_dir pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split \
  --output pipelines/panaroo/3_inclusivity/core100_identity.tsv

python scripts/filter_inclusivity_candidates.py \
  --input pipelines/panaroo/3_inclusivity/core100_identity.tsv \
  --output pipelines/panaroo/3_inclusivity/core100_candidates.tsv \
  --threshold 98.0

python scripts/drop_raw_duplicates.py \
  --input  pipelines/panaroo/3_inclusivity/core100_candidates.tsv \
  --output pipelines/panaroo/3_inclusivity/core100_candidates.nodup.tsv

6) Copy outputs to Windows (WSL)
# PASS FASTAs
cp -v pipelines/panaroo/5_exclusivity/PASS_FASTAs/*.fasta \
   /mnt/c/Users/krist/OneDrive/Documents/Bioinformatics/Bioinfo_instructions/StrepA_inclusivity_exclusivity/outputs/PASS_FASTAs/

# exclusivity summary
cp -v pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv \
   /mnt/c/Users/krist/OneDrive/Documents/Bioinformatics/Bioinfo_instructions/StrepA_inclusivity_exclusivity/outputs/

7) Adapting to another species (e.g., S. dysgalactiae)

Folders

Put GFFs under a separate root, e.g.
pipelines/panaroo/1_annotate_prokka_sdy/*/*.gff

Run (override input & outputs cleanly):

IN_GFF_DIR="pipelines/panaroo/1_annotate_prokka_sdy" \
PANAROO_OUT="pipelines/panaroo/2_panaroo/core100_sdy" \
bash scripts/run_all_spa.sh


This reuses the same logic but writes to core100_sdy trees.
If you want a different exclusivity DB (e.g., “non-dysgalactiae”), point BLAST_DB in the wrapper to a new DB path and build it with makeblastdb.

Tip: for cross-species studies, keep species-specific wrappers (run_all_spa.sh, run_all_sdy.sh) that do nothing but set the three variables and then source a common driver.

8) Git & housekeeping
# commit code + docs; large outputs are ignored by .gitignore
git add scripts/*.py scripts/*.sh docs/*.md
git commit -m "SPA: full core100 pipeline wrapper + docs; stable exclusivity"
git push origin main

# back up blast DBs or huge inputs (ignored) under backups/
mkdir -p backups/$(date +%F_%H%M%S)
cp -r pipelines/panaroo/5_exclusivity/blastdb backups/$(date +%F_%H%M%S)/

Quick reference of helper scripts

split_core100_from_panaroo.py — picks core genes by presence (GPA), extracts per-gene sequences from gene_data.csv/combined_DNA_CDS.fasta, aligns with MAFFT.

calculate_identity_with_names.py — avg pairwise %ID per alignment.

filter_inclusivity_candidates.py — keeps genes ≥ threshold (e.g., 98%).

drop_raw_duplicates.py — strips duplicate “*.raw” entries from tables.

build_consensus_from_split.py — makes one consensus per gene from split alignments.

summarize_snps_vs_consensus.py — per-strain and per-gene SNP summaries vs consensus.

summarize_blast_exclusivity.py — reads BLAST tabular output, annotates nearest neighbor & PASS/REJECT (pident/qcovs gate).

extract_pass_consensus.py — exports one FASTA per PASS gene to a folder.
