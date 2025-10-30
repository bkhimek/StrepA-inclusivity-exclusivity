# Runbook: Panaroo-based Strep A Pipeline

This runbook documents the workflow for finding **Streptococcus pyogenes–specific marker genes** using a pan-genome approach.  

The pipeline combines **Prokka** (annotation) + **Panaroo** (pangenome construction) with custom downstream scripts for:

- **Inclusivity** → identify genes present in almost all *S. pyogenes* strains with high sequence identity  
- **Consensus + SNP mapping** → generate representative consensus sequences and assess intra-species variation  
- **Exclusivity** → filter out genes with close matches in non-*S. pyogenes* species using BLAST  

The final output is a set of **candidate marker genes** that are:  
- (i) present in nearly all *S. pyogenes* strains (inclusivity), and  
- (ii) absent from or divergent in non-*S. pyogenes* bacteria (exclusivity).  

🔧 Scripts overview

All helper scripts (Bash + Python) are kept in scripts/.
This table shows which script belongs to which pipeline step:

| **Step** | **Purpose**                                 | **Script(s)**                                                                                                                                                            | **Type**      |
| -------- | ------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------- |
| Step 0   | Genome collection & cleanup                 | `collect_genomic_from_ncbi_pkg.sh`                                                                                                                                       | Bash          |
| Step 1   | Genome annotation (Prokka)                  | `run_prokka.sh`, `run_prokka_subset.sh`, `run_prokka_all.sh`                                                                                                        | Bash          |
| Step 2   | Pan-genome (Panaroo)                        | `run_panaroo.sh`, `run_panaroo_demo_core99.sh`, `run_panaroo_demo_core100.sh`                                                                                       | Bash          |
| Step 3   | Inclusivity analysis (identity + filtering) | `split_core99_from_panaroo.py`, `calculate_identity_with_names.py`, `filter_high_identity_genes.py`, wrappers (`run_inclusivity_demo*.sh`, `run_inclusivity_filter*.sh`) | Python + Bash |
| Step 4   | Consensus & SNP mapping                     | `build_consensus_from_split.py`, `summarize_snps_vs_consensus.py`, wrapper `run_consensus.sh`                                                                       | Python + Bash |
| Step 5   | Exclusivity BLAST & summary                 | `blast_consensus_exclusivity.py`, `summarize_blast_exclusivity.py`, wrapper `run_exclusivity.sh` + `run_exclusivity_summarize.sh`                                   | Python + Bash |
| Step 6   | Extract final PASS consensus FASTAs         | `extract_pass_consensus.py`, wrapper `run_extract_pass_consensus.sh`                                                                                                     | Python + Bash |

flowchart TD
  A0[Step 0: Download & clean genomes\n(data/genomes_genomic/*.fna)] --> A1
  A1[Step 1: Prokka annotate\nscripts/run_prokka.sh\n→ pipelines/panaroo/1_annotate_prokka/<ACC>/*.gff] --> A2
  A2[Step 2: Panaroo (core threshold)\n95%: scripts/run_panaroo_demo95.sh\n99%: scripts/run_panaroo_demo_core99.sh\n→ core_gene_alignment.aln + gene_presence_absence.csv] --> A2b
  A2b[Split concatenated alignment\nscripts/split_core_alignment.py\n→ core_gene_alignment.aln.split/*.fasta] --> A3a
  A3a[Step 3A: Identity per gene\nscripts/run_inclusivity_demo_core99.sh\n→ 3_inclusivity/*_identity.tsv] --> A3b
  A3b[Step 3B: Filter ≥98% identity\nscripts/run_inclusivity_filter_demo_core99.sh\n→ 3_inclusivity/*_candidates.tsv] --> A4
  A4[Step 4: Consensus + SNP mapping\nscripts/run_consensus.sh\n→ 4_consensus/*.fasta + SNP tables] --> A5
  A5[Step 5: Exclusivity BLAST vs non-pyogenes\nscripts/run_exclusivity.sh + summarize_blast_exclusivity.py\n→ 5_exclusivity/*_exclusivity.tsv] --> A6
  A6[Step 6: Extract PASS consensus FASTAs\nscripts/run_extract_pass_consensus.sh\n→ 6_reports/pass_consensus_split/*.fasta]

Notes:
- Bash wrappers run heavy external tools (Prokka, Panaroo, BLAST).
- Python scripts perform downstream analysis (identity filtering, consensus, SNP mapping, BLAST summarization).
- All scripts live in scripts/ for consistency.
- Steps 0–5 now tested end-to-end with real outputs.

📊 Example outputs per step
| **Step** | **Main outputs**                                                                    | **Folder(s)**                                                                             |
| -------- | ----------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------- |
| Step 0   | Raw genome FASTAs                                                                   | `data/genomes/`                                                                           |
| Step 1   | Prokka GFF annotations                                                              | `pipelines/panaroo/1_annotate_prokka/`                                                    |
| Step 2   | Panaroo pan-genome outputs (incl. core alignment, gene presence/absence, CDS fasta) | `pipelines/panaroo/2_panaroo/demo_core99/`                                                |
| Step 3   | Inclusivity: per-gene identity + filtered candidates                                | `pipelines/panaroo/3_inclusivity/`<br> (`*_identity.tsv`, `*_candidates.tsv`)             |
| Step 4   | Consensus FASTA (all candidates) + SNP tables (per-strain, per-gene)                | `pipelines/panaroo/4_consensus/`<br> (`*_consensus.fasta`, `*_snps.tsv`, `*_pergene.tsv`) |
| Step 5   | BLAST results against non-*S. pyogenes* db + exclusivity decision table             | `pipelines/panaroo/5_exclusivity/`<br> (`*_vs_nonpyogenes.tsv`, `*_exclusivity.tsv`)      |
| Step 6   | Final marker FASTAs (only PASS genes)                                               | `pipelines/panaroo/6_reports/pass_consensus_split/*.fasta`                                |


## Step 0. Genome download & preparation

- **Purpose**: Obtain and prepare RefSeq (GCF_) genome FASTA files for analysis.  
- Two common options are supported: **CDS FASTA** (coding sequences only) and **Genomic FASTA** (whole assemblies).  
- For Panaroo and Prokka, **Genomic FASTA is recommended**.

---

### 0A) Manual download on NCBI (browser)
1. Go to the NCBI Datasets page for *Streptococcus pyogenes* (taxon 1314) and filter **Assembly level = Complete genomes**.  
2. Click **Download package**.  
3. In the **Download Package** dialog:  
   - **Select file source**: All  
   - **Select file types**:  
     - If doing CDS workflow: **Genomic coding sequences (FASTA)**  
     - If doing genomic workflow: **Genomic sequences (FASTA)**  
4. Save the zip file, e.g. `S_pyogenes_ncbi_dataset.zip`.  
5. Unzip into a stable folder on Windows (example below).

---

### 0B) Windows → WSL path translation (example)
If the package unzips into Windows here:
C:\Users\krist\OneDrive\Documents\Bioinformatics\Bioinfo_instructions\StrepA_inclusivity_exclusivity\S_pyogenes_genomes\ncbi_dataset\data
… then the equivalent path in WSL is:
/mnt/c/Users/krist/OneDrive/Documents/Bioinformatics/Bioinfo_instructions/StrepA_inclusivity_exclusivity/S_pyogenes_genomes/ncbi_dataset/data

Check it in WSL:
```bash
ls -l "/mnt/c/Users/krist/OneDrive/Documents/Bioinformatics/Bioinfo_instructions/StrepA_inclusivity_exclusivity/S_pyogenes_genomes/ncbi_dataset/data" | head

### 0C) Collecting CDS FASTA (optional workflow)

Use the script:
scripts/collect_cds_from_ncbi_pkg.sh
It looks for cds_from_genomic.fna[.gz] inside GCF_* folders and copies them to:
data/genomes_cds/
with filenames like GCF_XXXXXX.cds.fna[.gz].
Run:bash scripts/collect_cds_from_ncbi_pkg.sh

ls -l data/genomes_cds | head
echo "Total CDS genomes:" $(ls data/genomes_cds/*.fna* 2>/dev/null | wc -l)
sed -n '1,15p' data/genomes_cds/_missing_cds.log

### 0D) Collecting Genomic FASTA (recommended workflow)
Use the script:
scripts/collect_genomic_from_ncbi_pkg.sh
It looks for *_genomic.fna[.gz] inside GCF_* folders and copies them to:
data/genomes_genomic/
with filenames like GCF_XXXXXX.genomic.fna[.gz].

Run:bash scripts/collect_genomic_from_ncbi_pkg.sh
ls -l data/genomes_genomic | head
echo "Total genomic genomes:" $(ls data/genomes_genomic/*.fna* 2>/dev/null | wc -l)
sed -n '1,15p' data/genomes_genomic/_missing_genomic.log

### 0E) Keep FASTA out of GitHub

Large FASTA files are ignored via .gitignore.
Only the scripts and logs are committed to GitHub.

### Outputs of Step 0: 
data/genomes_genomic/ (preferred) or data/genomes_cds/ with clean, renamed FASTA files.

A log file (_missing_genomic.log or _missing_cds.log) listing any missing sequences.

---

## Step 1. Genome annotation (Prokka)

- **Purpose**: Annotate raw genome assemblies with gene features.

- **Process**:
  1. Run Prokka on each genome FASTA file.  
  2. Store all resulting `.gff` annotation files in one folder.

- **Inputs**: Genome FASTA files from `data/genomes/`.  
- **Outputs**: Annotated `.gff` files in `pipelines/panaroo/1_annotate_prokka/`.

- **Example command**:
```bash
prokka --outdir annotated_genomes/STRAIN_NAME \
       --prefix STRAIN_NAME \
       data/genomes/STRAIN_NAME.fasta

## Step 2. Pan-genome construction (Panaroo)

Purpose: Build the pan-genome and core-genome alignment for S. pyogenes strains.
Inputs: Prokka GFFs from Step 1.
Outputs (demo core=95% default):
- pipelines/panaroo/2_panaroo/demo/core_gene_alignment.aln
- pipelines/panaroo/2_panaroo/demo/gene_presence_absence.csv
- other QC files in that folder

Example (demo):
conda activate panaroo_env
bash scripts/run_panaroo.sh
conda deactivate

## Step 2 (alt): “Near-universal” core (≥99% presence)
Use this when you want genes present in ~all strains.
Command:
conda activate panaroo_env
bash scripts/run_panaroo_demo_core99.sh
conda deactivate

Expected outputs:
pipelines/panaroo/2_panaroo/demo_core99/core_gene_alignment.aln
pipelines/panaroo/2_panaroo/demo_core99/gene_presence_absence.csv
pipelines/panaroo/2_panaroo/demo_core99/combined_DNA_CDS.fasta

(Summary showed ~1430 core genes at ≥99% presence.)


## Step 3. Inclusivity analysis

Goal: Find genes present in ≈all S. pyogenes strains that also show high sequence identity (e.g., ≥98%).

Inputs: Panaroo outputs from Step 2 (core_gene_alignment.aln, gene_presence_absence.csv, combined_DNA_CDS.fasta).

Outputs:

- pipelines/panaroo/3_inclusivity/demo_core99_identity.tsv — per-gene average pairwise identity
- pipelines/panaroo/3_inclusivity/demo_core99_candidates.tsv — genes passing the identity threshold

Scripts used:
- scripts/split_core99_from_panaroo.py — creates per-gene MAFFT alignments using GPA + gene_data mapping
- scripts/calculate_identity_with_names.py — computes per-gene % identity from split FASTAs
- scripts/filter_high_identity_genes.py — filters by identity cutoff

---

### Step 3A — Prepare per-gene alignments (core=99%) (split with MAFFT)

conda activate roary_env   # env with Biopython + mafft available
python scripts/split_core99_from_panaroo.py \
  --gpa       pipelines/panaroo/2_panaroo/demo_core99/gene_presence_absence.csv \
  --cds       pipelines/panaroo/2_panaroo/demo_core99/combined_DNA_CDS.fasta \
  --outdir    pipelines/panaroo/2_panaroo/demo_core99/core_gene_alignment.aln.split \
  --min_presence 0.99 \
  --threads 8
conda deactivate

Checks:

find pipelines/panaroo/2_panaroo/demo_core99/core_gene_alignment.aln.split -type f -name "*.fasta" | wc -l
ls -lh pipelines/panaroo/2_panaroo/demo_core99/core_gene_alignment.aln.split | head

Observed: 1430 per-gene aligned FASTAs written.

# Step 3B — Compute per-gene identities

conda activate roary_env
python scripts/calculate_identity_with_names.py \
  --split_dir pipelines/panaroo/2_panaroo/demo_core99/core_gene_alignment.aln.split \
  --output    pipelines/panaroo/3_inclusivity/demo_core99_identity.tsv
conda deactivate

Checks:
ls -lh pipelines/panaroo/3_inclusivity/demo_core99_identity.tsv
head -5 pipelines/panaroo/3_inclusivity/demo_core99_identity.tsv
- Expect header: Gene_File  Gene_Name  Average_Identity(%)

Observed: file present (~31 KB); identities e.g., COQ5_1 99.27, IMPDH 99.75, …

# Step 3C — Filter genes by identity threshold (≥98%)

bash scripts/run_inclusivity_filter_demo_core99.sh

Checks:
wc -l  pipelines/panaroo/3_inclusivity/demo_core99_candidates.tsv
head -10 pipelines/panaroo/3_inclusivity/demo_core99_candidates.tsv

Observed: 1142 genes (one name per line) passed the ≥98% identity filter from the core≥99% set.


If you want to version this in GitHub now:

```bash
git add docs/runbook_panaroo.md pipelines/panaroo/3_inclusivity/*.tsv scripts/*.py
git commit -m "Step 3 complete: split+identity (core99) and ≥98% candidates (n=1142)"
git push origin main

**Provenance**
- Panaroo core threshold used: **0.99**
- Identity filter: **≥98%**
- Result counts:
  - Core≥99% genes: ~1430 (from `summary_statistics.txt`)
  - Passing identity (≥98%): **1142** (from `demo_core99_candidates.tsv`)
- Environments:
  - `panaroo_env` → panaroo v1.5.x, mafft v7.52x
  - `roary_env`  → Biopython + mafft present
- Files tracked in Git:
  - `pipelines/panaroo/3_inclusivity/demo_core99_identity.tsv`
  - `pipelines/panaroo/3_inclusivity/demo_core99_candidates.tsv`

## ## Step 4 — Consensus + SNP summary

**Goal
From the high-identity inclusivity candidates (Step 3B), build consensus sequences for each gene and quantify per-strain and per-gene SNP variation.

**Inputs
From Step 2 (Panaroo, core=99%):
pipelines/panaroo/2_panaroo/demo_core99/core_gene_alignment.aln
From Step 3B (filtered inclusivity candidates):
pipelines/panaroo/3_inclusivity/demo_core99_candidates.tsv

**Scripts
scripts/build_consensus_from_split.py
Builds consensus FASTA per candidate gene from the split core alignments.
scripts/map_SNPs_vs_consensus.py
Counts SNPs per strain vs. consensus and aggregates to per-gene stats.

Wrapper: scripts/run_consensus.sh
Orchestrates both steps above and handles inputs/outputs.

Run:
conda activate roary_env
bash scripts/run_consensus.sh
conda deactivate


**Outputs
pipelines/panaroo/4_consensus/demo_core99_consensus.fasta — consensus sequences (one per gene)
pipelines/panaroo/4_consensus/demo_core99_snps.tsv — SNPs per (gene × strain)
pipelines/panaroo/4_consensus/demo_core99_pergene.tsv — per-gene SNP summary (mean, max SNPs)

**Checks
ls -lh pipelines/panaroo/4_consensus/
head -5 pipelines/panaroo/4_consensus/demo_core99_pergene.tsv
head -5 pipelines/panaroo/4_consensus/demo_core99_snps.tsv
grep -m 5 "^>" pipelines/panaroo/4_consensus/demo_core99_consensus.fasta


**QC / interpretation
# Top 15 most conserved genes (lowest mean SNPs)
(head -1 && tail -n +2 | sort -k3,3n -k4,4n | head -15) < pipelines/panaroo/4_consensus/demo_core99_pergene.tsv

# Genes with worrying divergence (max SNPs > 10)
awk 'NR==1 || $4>10' pipelines/panaroo/4_consensus/demo_core99_pergene.tsv | head

# Count of very stable genes (mean ≤1 SNP, max ≤3 SNPs)
awk 'NR>1 && $3<=1 && $4<=3 {c++} END{print c+0}' pipelines/panaroo/4_consensus/demo_core99_pergene.tsv

## Step 5. Exclusivity analysis (BLAST)

Goal: Ensure that candidate markers are specific to Streptococcus pyogenes by removing genes that also appear (with high similarity) in non-S. pyogenes genomes.

 - Inputs:
Inclusivity consensus genes (from Step 4)
Non-S. pyogenes BLAST database (pre-built, e.g. streptococcus_non_pyogenes_db.*)
BLAST search results vs that database

- Outputs:
pipelines/panaroo/5_exclusivity/demo_core99_vs_nonpyogenes.tsv — raw BLAST hits
pipelines/panaroo/5_exclusivity/demo_core99_exclusivity.tsv — summarized PASS/REJECT per gene with nearest neighbor species

- Scripts used:
scripts/run_exclusivity.sh — runs BLAST of inclusivity candidates vs non-S. pyogenes db
scripts/summarize_blast_exclusivity.py — parses BLAST TSV, extracts best hits per gene, assigns decision
scripts/run_exclusivity_summarize.sh — wrapper to run the summarizer

## Step 5A — Run BLAST search

conda activate roary_env
bash scripts/run_exclusivity.sh
conda deactivate

Produces raw BLAST output:
pipelines/panaroo/5_exclusivity/demo_core99_vs_nonpyogenes.tsv

## Step 5B — Summarize exclusivity results

conda activate roary_env
bash scripts/run_exclusivity_summarize.sh
conda deactivate

Produces summary table:
pipelines/panaroo/5_exclusivity/demo_core99_exclusivity.tsv

Columns:
Gene — candidate gene
neighbor — closest non-S. pyogenes species
best_pident — best % identity
best_qcovs — query coverage %
decision — PASS (exclusive to S. pyogenes) or REJECT (too similar to others

Gene    neighbor                                    best_pident  best_qcovs  decision
COQ5_1  Streptococcus dysgalactiae subsp. equisimilis   97.548      99       REJECT
COQ5_2  Streptococcus agalactiae                        65.867      98       PASS
IMPDH   Streptococcus canis                             82.555      100      PASS
accB    Streptococcus dysgalactiae                      84.449      99       PASS

- Count how many passed vs rejected
awk -F'\t' 'NR>1 {count[$5]++} END{for (k in count) print k,count[k]}' \
  pipelines/panaroo/5_exclusivity/demo_core99_exclusivity.tsv
- Preview some PASS genes
awk -F'\t' '$5=="PASS"{print $0}' pipelines/panaroo/5_exclusivity/demo_core99_exclusivity.tsv | head -10


## Step 6. Extract consensus FASTAs for PASS genes

# Goal: From the consensus sequences (Step 4) and exclusivity summary (Step 5), extract individual FASTA files only for genes that passed exclusivity (decision = PASS).

Inputs:
Consensus FASTA: pipelines/panaroo/4_consensus/demo_core99_consensus.fasta
Exclusivity summary: pipelines/panaroo/5_exclusivity/demo_core99_exclusivity.tsv

# Outputs:
One FASTA file per PASS gene under pipelines/panaroo/6_reports/pass_consensus_split/

# Script:
scripts/extract_pass_consensus.py (Python)
Wrapper: scripts/run_extract_pass_consensus.sh (Bash)

Command:
bash scripts/run_extract_pass_consensus.sh

ls -lh pipelines/panaroo/6_reports/pass_consensus_split | head
- Example output:
 -rw-r--r-- 1 user group  800 Oct 25 14:10 accB.fasta
 -rw-r--r-- 1 user group  745 Oct 25 14:10 COQ5_2.fasta
 -rw-r--r-- 1 user group  912 Oct 25 14:10 IMPDH.fasta


Final Deliverables

At the end of the full workflow (Steps 0–6), the pipeline produces a shortlist of candidate diagnostic markers for Streptococcus pyogenes.

✅ Main final outputs

FASTA files (per gene):
Located in: pipelines/panaroo/6_reports/pass_consensus_split/*.fasta

Each file contains the consensus sequence for a single gene that:

is present in ≥99% of S. pyogenes strains
has ≥98% pairwise identity
passes exclusivity check (no high-identity matches in non-S. pyogenes species)
Master TSV with decisions

pipelines/panaroo/5_exclusivity/demo_core99_exclusivity.tsv

Contains the per-gene summary with:
 - closest non-S. pyogenes match (species name)
 - best % identity and query coverage
 - final decision (PASS or REJECT)

🧾 Interpretation

- Genes labeled PASS are candidate inclusivity/exclusivity markers.
- Genes labeled REJECT are excluded because of close similarity to non-S. pyogenes species.
- The FASTA files in pass_consensus_split/ represent the final marker panel and can be used directly for downstream assay design (primer/probe design, etc.).

