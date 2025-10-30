Runbook: Strep A – All Genomes Analysis

This runbook describes the workflow for full-scale analysis of Streptococcus pyogenes genomes using Prokka + Panaroo + custom scripts.
It extends the demo workflow and standardizes outputs for inclusivity and exclusivity analysis.

Runbook: Strep A – All Genomes Analysis

This runbook describes the workflow for full-scale analysis of Streptococcus pyogenes genomes using Prokka + Panaroo + custom scripts.
It extends the demo workflow and standardizes outputs for inclusivity and exclusivity analysis.

## Step 0. Genome collection

Purpose: Download and organize complete genomes from NCBI.
Inputs: List of accession IDs (e.g., data/accessions/strepA_all.txt).
Outputs: Genome FASTA files in data/genomes_genomic/.
Example:
bash scripts/collect_genomic_from_ncbi_pkg.sh data/accessions/strepA_all.txt

## Step 1. Genome annotation (Prokka)

Inputs: Genome FASTAs from data/genomes_genomic/.
Outputs: Annotated .gff files in pipelines/panaroo/1_annotate_prokka/.
Example:
bash scripts/run_prokka_all.sh

## Step 2. Pan-genome construction (Panaroo)

Purpose: Build pan-genome and core alignment with core threshold = 100%.
Inputs: Annotated .gff files.
Outputs:
pipelines/panaroo/2_panaroo/all/core_gene_alignment.aln
pipelines/panaroo/2_panaroo/all/gene_presence_absence.csv
pipelines/panaroo/2_panaroo/all/combined_DNA_CDS.fasta

Example (inside run_all.sh):
panaroo \
  -i pipelines/panaroo/1_annotate_prokka/*/*.gff \
  -o pipelines/panaroo/2_panaroo/all \
  --clean-mode strict \
  --core_threshold 1.0 \
  --threads 8

Step 3. Inclusivity analysis
## 3A. Split and align per-gene
python scripts/split_core100_from_panaroo.py \
  --gpa pipelines/panaroo/2_panaroo/all/gene_presence_absence.csv \
  --cds pipelines/panaroo/2_panaroo/all/combined_DNA_CDS.fasta \
  --outdir pipelines/panaroo/2_panaroo/all/core_gene_alignment.aln.split \
  --min_presence 1.0 \
  --threads 8

## 3B. Compute pairwise identity
python scripts/calculate_identity_with_names.py \
  --split_dir pipelines/panaroo/2_panaroo/all/core_gene_alignment.aln.split \
  --output    pipelines/panaroo/3_inclusivity/all_core100_identity.tsv

## 3C. Filter ≥98% identity
bash scripts/run_inclusivity_filter_core100.sh

Outputs:
pipelines/panaroo/3_inclusivity/all_core100_identity.tsv
pipelines/panaroo/3_inclusivity/all_core100_candidates.tsv

## Step 4. Consensus & SNPs
bash scripts/run_consensus_core100.sh
Outputs:

pipelines/panaroo/4_consensus/core100_consensus.fasta
pipelines/panaroo/4_consensus/core100_snps.tsv
pipelines/panaroo/4_consensus/core100_pergene.tsv

Optional: extract PASS-only consensus FASTAs into per-gene files:
bash scripts/run_extract_pass_consensus.sh

## Step 5. Exclusivity analysis

5A. BLAST against non-pyogenes Streptococci
bash scripts/run_exclusivity_core100.sh

Input DB: streptococcus_non_pyogenes_db (prebuilt with makeblastdb).
Output: pipelines/panaroo/5_exclusivity/core100_vs_nonpyogenes.tsv.

5B. Summarize best hits

bash scripts/run_exclusivity_summarize.sh

Output: pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv
with columns:
Gene   neighbor_species   best_pident   best_qcovs   decision

Step 6. Results & reporting

Inclusivity: all_core100_candidates.tsv → list of highly conserved genes.
Consensus FASTAs: core100_consensus.fasta → reference for SNP mapping.
Exclusivity: core100_exclusivity.tsv → best matches in non-pyogenes.

PASS-only FASTAs: folder with candidate-specific sequences for downstream primer/probe design.

Example: Run all steps at once
The pipeline is wrapped in:

bash scripts/run_all.sh

Logs are written into logs/.

🔄 Extending to Other Species

If you want to analyze another species (say Streptococcus agalactiae):

1. Genome collection

Prepare a new accession list, e.g. data/accessions/strep_agalactiae.txt.

Download genomes with:

bash scripts/collect_genomic_from_ncbi_pkg.sh data/accessions/strep_agalactiae.txt


2. Folder organization

Keep outputs in species-specific subfolders, e.g.:

pipelines/panaroo/2_panaroo/strep_agalactiae/
pipelines/panaroo/3_inclusivity/strep_agalactiae/
pipelines/panaroo/4_consensus/strep_agalactiae/
pipelines/panaroo/5_exclusivity/strep_agalactiae/


3. Scripts

Most scripts take --gpa, --cds, --output arguments → you just point them to the correct input/output.

The only strict thing: adjust core threshold (e.g., 1.0 for 100% presence or 0.95 for relaxed core).

4. Exclusivity

Build a non-target database excluding the species of interest.
Example: for S. agalactiae analysis, build DB from all Streptococci except S. agalactiae.

makeblastdb -in streptococcus_non_agalactiae.fna -dbtype nucl -out streptococcus_non_agalactiae_db


Update exclusivity scripts to use the correct DB.

🗂️ Tips for Project Organization

data/ → raw input (accession lists, FASTAs, BLAST DBs).

pipelines/ → heavy intermediate outputs.

results/runs/run-YYYY-MM-DD/ → symlink/copy of summary TSVs + FASTAs.

docs/ → runbooks, reports, figures.

scripts/ → wrappers and Python helpers.

✅ With this structure you can:

Run Strep A full analysis via run_all.sh.

Add new species by creating a new accession list + DB, without touching core scripts.

Compare inclusivity/exclusivity results across species in results/runs/.













