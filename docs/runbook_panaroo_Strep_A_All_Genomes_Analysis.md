# Panaroo Runbook — *Streptococcus pyogenes* (Strep A)  
## All Genomes Inclusivity–Exclusivity Analysis

---

## Table of Contents
- [Overview](#overview)
- [Scope](#scope)
- [Expected Folder Layout](#expected-folder-layout)
- [Software & Environments](#software--environments)
- [Workflow Summary](#workflow-summary)
- [Step-by-Step Pipeline](#step-by-step-pipeline)
  - [Step 0 — Input Preparation](#step-0--input-preparation)
  - [Step 1 — Prokka Annotation](#step-1--prokka-annotation)
  - [Step 2 — Panaroo Pangenome Analysis](#step-2--panaroo-pangenome-analysis)
  - [Step 3 — Inclusivity Analysis](#step-3--inclusivity-analysis)
  - [Step 4 — Consensus / SNP Analysis](#step-4--consensus--snp-analysis)
  - [Step 5 — Exclusivity Analysis](#step-5--exclusivity-analysis)
  - [Step 6 — Reports & Outputs](#step-6--reports--outputs)
- [Common Pitfalls & Troubleshooting](#common-pitfalls--troubleshooting)
- [Notes & Best Practices](#notes--best-practices)

---

## Overview

This runbook documents the **end-to-end pipeline** for large-scale analysis of  
*Streptococcus pyogenes* (Strep A) genomes using:

- **Prokka** for genome annotation  
- **Panaroo** for pangenome construction  
- **Custom scripts** for inclusivity, consensus/SNP, and exclusivity analysis  

The workflow is **wrapper-driven**, standardized, and intended for reproducible
analysis across **all downloaded Strep A genomes**.

---

## Scope

This runbook covers:

- Full-genome inclusivity analysis
- Core-gene consensus and SNP detection
- Exclusivity screening against non-*S. pyogenes* genomes
- Generation of standardized outputs for downstream reporting

This document **extends the demo workflow** and formalizes the structure and
expected outputs for production-scale runs.

---

## Expected Folder Layout

The pipeline assumes the following repository structure:

```text
StrepA-inclusivity-exclusivity/
├─ data/
│  ├─ all_genomes_raw/          # Raw input FASTA (.fna)
│  └─ non_pyogenes/             # Non–S. pyogenes genomes for exclusivity
│
├─ pipelines/
│  └─ panaroo/
│     ├─ 1_annotate_prokka_spa/
│     │  └─ gff/
│     ├─ 2_panaroo/
│     │  └─ core100_spa/
│     ├─ 3_inclusivity/
│     ├─ 4_consensus_snps/
│     ├─ 5_exclusivity/
│     │  └─ blastdb/
│     └─ 6_reports/
│
├─ scripts/
│  ├─ run_all_genomes_wrapper.sh
│  └─ analysis_helpers/
│
└─ docs/
   └─ runbook_panaroo_Strep_A_All_Genomes_Analysis.md
⚠️ Important:
The wrapper script expects this structure. Renaming folders will break paths.

Software & Environments
Recommended environment setup (Conda):

prokka

panaroo

blast

mafft

python ≥ 3.8

Example:

conda activate panaroo_env
Workflow Summary
Prepare and validate genome FASTA inputs

Annotate genomes with Prokka

Run Panaroo to generate a core pangenome

Identify inclusivity candidates

Generate consensus sequences and SNPs

Perform exclusivity analysis via BLAST

Collect standardized reports

Step-by-Step Pipeline
Step 0 — Input Preparation
Goal: Ensure all genome FASTA files are valid and consistent.

Input: .fna files

Location: data/all_genomes_raw/

ls data/all_genomes_raw/*.fna | wc -l
Tips

Use consistent naming (no spaces)

Avoid mixed assemblies (chromosome + contigs in one file)

Step 1 — Prokka Annotation
Goal: Annotate all genomes uniformly before Panaroo.

prokka \
  --outdir pipelines/panaroo/1_annotate_prokka_spa \
  --prefix spa \
  data/all_genomes_raw/*.fna
Outputs

.gff files used as Panaroo input

Stored in 1_annotate_prokka_spa/gff/

Step 2 — Panaroo Pangenome Analysis
Goal: Build the pangenome and extract 100% core genes.

panaroo \
  -i pipelines/panaroo/1_annotate_prokka_spa/gff/*.gff \
  -o pipelines/panaroo/2_panaroo/core100_spa \
  --clean-mode strict \
  --core_threshold 1.00 \
  --threads 8
Key Parameters

--clean-mode strict → removes annotation artifacts

--core_threshold 1.00 → genes present in 100% of genomes

Step 3 — Inclusivity Analysis
Goal: Identify genes consistently present across all Strep A genomes.

bash scripts/run_inclusivity.sh \
  pipelines/panaroo/2_panaroo/core100_spa \
  pipelines/panaroo/3_inclusivity
Outputs

Core gene presence tables

Candidate inclusivity markers

Step 4 — Consensus / SNP Analysis
Goal: Build consensus sequences and detect SNP variation.

bash scripts/run_consensus_snps.sh \
  pipelines/panaroo/3_inclusivity \
  pipelines/panaroo/4_consensus_snps
Notes

MAFFT used for alignments

SNP density can guide marker selection

Step 5 — Exclusivity Analysis
Goal: Confirm candidates are absent from non-S. pyogenes genomes.

Build BLAST database:

makeblastdb \
  -in data/non_pyogenes/non_pyogenes.fna \
  -dbtype nucl \
  -out pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes
Run BLAST screening:

bash scripts/run_exclusivity_blast.sh \
  pipelines/panaroo/4_consensus_snps \
  pipelines/panaroo/5_exclusivity
Step 6 — Reports & Outputs
Goal: Generate standardized outputs for review.

Outputs saved in:

text
pipelines/panaroo/6_reports/
Includes:

Final inclusivity–exclusivity tables

Candidate gene lists

Summary statistics

## Common Pitfalls & Troubleshooting

### ❌ Panaroo fails or produces empty core genes
- Check GFF consistency
- Ensure Prokka versions are uniform
- Avoid mixed genome formats (chromosome + contigs)

### ❌ Too few core genes
- Confirm `--core_threshold 1.00` is appropriate
- Inspect input genome quality
- Check for failed annotations

### ❌ BLAST returns many hits
- Check database composition
- Ensure non-*S. pyogenes* genomes only
- Tighten identity and coverage thresholds

### ❌ Wrapper script fails
- Verify folder structure matches the runbook
- Ensure executable permissions:

```bash
chmod +x scripts/*.sh
