Runbook: Strep A – All Genomes Analysis

This runbook describes the workflow for full-scale analysis of Streptococcus pyogenes genomes using Prokka + Panaroo + custom scripts.
It extends the demo workflow and standardizes outputs for inclusivity and exclusivity analysis.

Runbook: All-Genomes Analysis (Wrapper-driven)

This runbook documents the end-to-end pipeline for running inclusivity → consensus/SNPs → exclusivity over all downloaded genomes of Streptococcus pyogenes, driven by a single wrapper script.

It also explains how to adapt the same workflow to another species (e.g., Streptococcus dysgalactiae).

0) Folder layout (expected)
StrepA-inclusivity-exclusivity/
├─ data/
│  └─ all_genomes_raw/                     # raw .fna (genomic FASTA) for the target species
├─ pipelines/
│  └─ panaroo/
│     ├─ 1_annotate_prokka_spa/            # Prokka outputs (one subfolder per genome)
│     ├─ 2_panaroo/
│     │  └─ core100_spa/                    # Panaroo outputs for the target species
│     ├─ 3_inclusivity/
│     ├─ 4_consensus/
│     ├─ 5_exclusivity/
│     │  └─ blastdb/
│     │     └─ non_pyogenes/                # BLAST DB of exclusion panel (near-neighbors)
│     └─ 6_reports/
├─ results/
│  └─ runs/                                 # timestamped run folders with final artifacts
├─ scripts/                                 # wrapper + helper scripts
└─ docs/


You already have this structure; just keep new analyses consistent (one Panaroo output folder per species/threshold).

1) Environments (conda)

prokka_env: Prokka

panaroo_env: Panaroo + cd-hit

roary_env: Biopython, pandas, BLAST+ (makeblastdb, blastn), mafft

Activate as needed; the wrapper will call the right env per step using conda run -n <env> ….

2) One-command run (the wrapper)

From repo root:

# Option A: just run with defaults (core=1.00 on spa folders)
bash scripts/run_all.sh

# Option B: override inputs/outputs (e.g., for a different species run)
IN_GFF_DIR="pipelines/panaroo/1_annotate_prokka_spa" \
PANAROO_OUT="pipelines/panaroo/2_panaroo/core100_spa" \
bash scripts/run_all.sh


On success you’ll see something like:

[DONE] Run completed -> results/runs/run-YYYY-MM-DD_HHMMSS
       Candidates:   pipelines/panaroo/3_inclusivity/core100_candidates.tsv
       Consensus:    pipelines/panaroo/4_consensus/core100_consensus.fasta
       BLAST sum:    pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv
       PASS FASTAs:  pipelines/panaroo/5_exclusivity/PASS_FASTAs/

3) What the wrapper does (step-by-step)
Step 1 — Panaroo (core threshold = 100%)

Input: Prokka GFFs under pipelines/panaroo/1_annotate_prokka_spa/*/*.gff

Runs: panaroo with --clean-mode strict --core_threshold 1.00

Output:
pipelines/panaroo/2_panaroo/core100_spa/{gene_presence_absence.csv, combined_DNA_CDS.fasta, gene_data.csv}

Helper: none (external tool).

Step 2 — Split per-gene & MAFFT align

Runs: scripts/split_core100_from_panaroo.py

Uses gene_presence_absence.csv + gene_data.csv to map cluster → CDS IDs, pulls sequences from combined_DNA_CDS.fasta, groups by cluster, and MAFFT-aligns each gene.

Output: pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split/*.fasta

Helper: split_core100_from_panaroo.py

Step 3 — Per-gene average identity

Runs: scripts/calculate_identity_with_names.py

Computes average pairwise identity per aligned gene file.

Output: pipelines/panaroo/3_inclusivity/core100_identity.tsv
(columns: Gene_File, Gene_Name, Average_Identity(%))

Helper: calculate_identity_with_names.py

Step 4 — Filter ≥98% identity (inclusivity)

Runs: scripts/filter_inclusivity_candidates.py
(Your wrapper may call the pre-existing run_inclusivity_filter_core100.sh—both are equivalent.)

Output: pipelines/panaroo/3_inclusivity/core100_candidates.tsv
(genes that meet the identity threshold, potentially with .raw duplicates)

Helper: filter_inclusivity_candidates.py (or run_inclusivity_filter_core100.sh)

We then deduplicate any *.raw companion rows into core100_candidates.nodup.tsv.

Step 5 — Build consensus sequences

Runs: scripts/build_consensus_from_split.py with --genes_txt core100_candidates.nodup.tsv

Input: split alignments (Step 2) + the filtered gene list

Output: pipelines/panaroo/4_consensus/core100_consensus.fasta
(one consensus per kept gene)

Helper: build_consensus_from_split.py

Step 6 — Exclusivity BLAST & summary

DB: pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/streptococcus_non_pyogenes_db

If missing: wrapper rebuilds it with makeblastdb from streptococcus_non_pyogenes_combined.fna

Runs: blastn on the consensus FASTA; then:

scripts/summarize_blast_exclusivity.py to produce:
pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv
(columns: Gene, neighbor, best_pident, best_qcovs, decision)

scripts/extract_pass_consensus.py to export PASS genes into
pipelines/panaroo/5_exclusivity/PASS_FASTAs/ as individual FASTA files (one per gene; no .raw duplicates).

Helpers: summarize_blast_exclusivity.py, extract_pass_consensus.py

4) Quick checks
# Panaroo key files
ls -lh pipelines/panaroo/2_panaroo/core100_spa/{gene_presence_absence.csv,combined_DNA_CDS.fasta,gene_data.csv}

# Split alignments present?
find pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split -name "*.fasta" | wc -l

# Inclusivity tables not empty?
wc -l pipelines/panaroo/3_inclusivity/core100_identity.tsv
wc -l pipelines/panaroo/3_inclusivity/core100_candidates.tsv
wc -l pipelines/panaroo/3_inclusivity/core100_candidates.nodup.tsv

# Consensus FASTA has many entries?
grep -c '^>' pipelines/panaroo/4_consensus/core100_consensus.fasta

# Exclusivity summary & PASS FASTAs
head -5 pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv
find pipelines/panaroo/5_exclusivity/PASS_FASTAs -name "*.fasta" | wc -l

5) Helper scripts (one-line purpose)

split_core100_from_panaroo.py — map Panaroo clusters to CDSs using gene_data.csv, create per-gene MAFFT alignments.

calculate_identity_with_names.py — compute average pairwise identity per gene alignment.

filter_inclusivity_candidates.py — keep genes with identity ≥ threshold (default 98%).

build_consensus_from_split.py — build consensus sequences for a list of genes from the split alignments.

summarize_blast_exclusivity.py — parse BLAST tabular output with species columns (includes sscinames, stitle), pick best off-target hit per gene, and mark PASS/REJECT using min_pident and min_qcovs.

extract_pass_consensus.py — export only PASS genes to one-FASTA-per-gene in PASS_FASTAs/ and avoid .raw duplicates.

6) Adapting to a different species (e.g., S. dysgalactiae)

Download genomes (e.g., NCBI datasets) to data/all_genomes_raw/ for the new species.

Annotate with Prokka (new output root to keep things separate). Example:

# in prokka_env
mkdir -p pipelines/panaroo/1_annotate_prokka_sdy
parallel -j 8 '
  BASENAME=$(basename {} .fna)
  prokka \
    --outdir pipelines/panaroo/1_annotate_prokka_sdy/${BASENAME} \
    --prefix ${BASENAME} \
    --force --cpus 1 {}
' ::: data/all_genomes_raw/*.fna


Choose a new Panaroo output folder, e.g. pipelines/panaroo/2_panaroo/core100_sdy.

Run the wrapper overriding input GFF folder + Panaroo out:

IN_GFF_DIR="pipelines/panaroo/1_annotate_prokka_sdy" \
PANAROO_OUT="pipelines/panaroo/2_panaroo/core100_sdy" \
bash scripts/run_all.sh


Exclusion BLAST DB

Either reuse the existing non-pyogenes panel if appropriate, or build a new exclusion DB tailored to your new target (e.g., other streptococci, oropharyngeal flora).

Place combined FASTA at:
pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/streptococcus_non_pyogenes_combined.fna

The wrapper will build/refresh the DB automatically (or you can run makeblastdb yourself).

Thresholds

Inclusivity identity can stay at ≥98% to stay strict across the species.

Exclusivity decision defaults: min_pident=85, min_qcovs=80. Raise pident (e.g., 90–95) to be more conservative, or lower qcovs to allow near hits on shorter segments.

Outputs
The wrapper will emit the same artifacts under pipelines/panaroo/… and a timestamped results/runs/run-YYYY-MM-DD_HHMMSS folder with summarized copies.

7) Common pitfalls (and fixes we already implemented)

CDS/cluster naming mismatches between Panaroo tables and FASTA headers → fixed by using gene_data.csv mapping in split_core100_from_panaroo.py.

.raw duplicates in candidates/consensus → we deduplicate lists and avoid exporting .raw in final PASS FASTAs.

BLAST DB missing/corrupt → the wrapper rebuilds from …/non_pyogenes/streptococcus_non_pyogenes_combined.fna when needed.

Globbing GFF paths → wrapper supports flat or nested GFF organization by reading ${IN_GFF_DIR}/*/*.gff. Put Prokka outputs in …/strain_id/strain_id.gff and you’re good.

8) Reproducibility & versioning

After a successful run:

git add docs/runbook_panaroo_Strep_A_All_Genomes_Analysis.md \
        scripts/run_all.sh scripts/*.py \
        results/runs/run-YYYY-MM-DD_HHMMSS/*
git commit -m "All-genomes run: core=100%; wrapper-driven; results + docs"
git push origin main


Large intermediates are ignored by .gitignore; summarized TSV/FASTA/plots in results/runs/… are kept.

TL;DR

Run everything with: bash scripts/run_all.sh

Check outputs in results/runs/run-*/

Adapt to new species by swapping Prokka GFF folder and Panaroo output root, and (optionally) adjusting the exclusion BLAST DB.








