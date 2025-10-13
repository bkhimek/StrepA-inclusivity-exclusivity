# Runbook: Panaroo-based Strep A Pipeline

This runbook captures the main workflow for inclusivity/exclusivity analysis of *Streptococcus pyogenes* using Prokka + Panaroo.

It is meant as a **living document**:  
- Every step shows *what is done*, *which scripts are used*, and *where outputs go*.  
- As you refine/correct scripts, update the notes here so this always reflects the working path.

---

## Step 0. Genome download & preparation

- **Purpose**: Obtain and prepare RefSeq (GCF_) CDS FASTA files for analysis, with clear, repeatable steps.

### 0A) Manual download on NCBI (browser)
1. Go to the NCBI Datasets page for *Streptococcus pyogenes* (taxon 1314) and filter **Assembly level = Complete genomes**.
2. Click **Download package**.
3. In the **Download Package** dialog:
   - **Select file source**: All
   - **Select file types**: **Genomic coding sequences (FASTA)** (CDS)
4. Save the zip file, e.g. `S_pyogenes_ncbi_dataset.zip`.

> This package typically contains both **GCF_** (RefSeq) and **GCA_** (GenBank) folders.

### 0B) Unzip to Windows folder (kept as your raw archive)
- Unzip the package to a stable location on Windows, e.g.:
C:\Users\krist\OneDrive\Documents\Bioinformatics\StrepA_inclusivity_exclusivity\S_pyogenes_genomes\ncbi_dataset\

After unzip, the path that contains all the per-accession folders is:

C:\Users\krist\OneDrive\Documents\Bioinformatics\StrepA_inclusivity_exclusivity\S_pyogenes_genomes\ncbi_dataset\data


### 0C) Windows → WSL path translation (example)
- Any `C:\...` Windows path is seen in WSL as `/mnt/c/...`.

**Your exact example:**
- Windows: C:\Users\krist\OneDrive\Documents\Bioinformatics\StrepA_inclusivity_exclusivity\S_pyogenes_genomes\ncbi_dataset\data
- WSL: /mnt/c/Users/krist/OneDrive/Documents/Bioinformatics/StrepA_inclusivity_exclusivity/S_pyogenes_genomes/ncbi_dataset/data


To sanity-check in WSL:
```bash
ls -l "/mnt/c/Users/krist/OneDrive/Documents/Bioinformatics/StrepA_inclusivity_exclusivity/S_pyogenes_genomes/ncbi_dataset/data" | head
You should see many GCF_* and GCA_* folders.

### 0D) Collect only RefSeq CDS into the repo (scripted cleanup)
We keep the heavy raw package on Windows, and copy only what we need (RefSeq GCF_ CDS FASTA) into the repo on the Linux side.
1. Script location in repo: scripts/collect_cds_from_ncbi_pkg.sh
2. Script config inside the file (already set for you):
PKG_ROOT="/mnt/c/Users/krist/OneDrive/Documents/Bioinformatics/StrepA_inclusivity_exclusivity/S_pyogenes_genomes/ncbi_dataset/data"
DEST="data/genomes_cds"
3. Dry-run (prints what would be copied): DRYRUN=1 bash scripts/collect_cds_from_ncbi_pkg.sh
4. Real run (copies files): bash scripts/collect_cds_from_ncbi_pkg.sh
The script:
- iterates only over GCF_* folders (RefSeq),
- looks for cds_from_genomic.fna[.gz],
- copies to data/genomes_cds/ as GCF_xxxxxxx.cds.fna[.gz],
- logs any missing CDS to data/genomes_cds/_missing_cds.log.

### 0E) Verify results

ls -l data/genomes_cds | head
echo "Total CDS files:" $(ls data/genomes_cds/*.fna* 2>/dev/null | wc -l)
sed -n '1,15p' data/genomes_cds/_missing_cds.log

### Keep FASTA out of GitHub (.gitignore)
We track scripts/configs/docs in Git, not large FASTA files.
.gitignore (root of repo) should include:
# Ignore raw genome data
data/genomes_cds/
data/genomes_genomic/
# Commit the .gitignore change if you added those lines:
git add .gitignore
git commit -m "Update .gitignore: ignore genome FASTA folders"
git push origin main

### 0G) Track the collector script in Git (safe, small) 
git add scripts/collect_cds_from_ncbi_pkg.sh
git commit -m "Add script: collect RefSeq (GCF) CDS FASTA from NCBI package"
git push origin main

### Outputs of Step 0: Clean, renamed RefSeq CDS FASTA files in data/genomes_cds/ on the Linux side; missing entries listed in data/genomes_cds/_missing_cds.log.

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

- **Purpose**: Build the pan-genome and core-genome alignment for *S. pyogenes* strains.

- **Process**:
  1. Collect the annotated `.gff` files from Step 1 (Prokka).  
  2. Run Panaroo on all `.gff` files with strict cleaning mode.  
  3. Outputs include the core genome alignment and a presence/absence matrix.

- **Inputs**: Annotated genome `.gff` files from `pipelines/panaroo/1_annotate_prokka/`.  
- **Outputs**:  
  - `panaroo_out/core_gene_alignment.aln` (core genome alignment)  
  - `panaroo_out/gene_presence_absence.csv` (gene presence/absence matrix)  
  - Other QC outputs in `panaroo_out/`

- **Example command**:
```bash
panaroo -i annotated_genomes/*.gff \
        -o panaroo_out \
        --clean-mode strict

## Step 3. Inclusivity analysis

- **Purpose**: Identify genes present in nearly all *S. pyogenes* strains (≥98% identity).

- **Process**:
  1. Calculate pairwise identity across core genes.  
  2. Filter genes that meet high-identity thresholds.  
  3. Collect candidate inclusive genes.

- **Inputs**: Core alignment from Step 2.  
- **Outputs**: Inclusivity candidate table.  

- **Scripts**:  
  - `calculate_identity_with_names.py`  
  - `filter_high_identity_genes.py`

- **Example command**:
```bash
python pipelines/panaroo/3_inclusivity/calculate_identity_with_names.py
python pipelines/panaroo/3_inclusivity/filter_high_identity_genes.py

## Step 4. Consensus building & SNP mapping

- **Purpose**: Build consensus sequences of candidate inclusive genes and evaluate SNP variation.

- **Process**:
  1. Generate consensus sequences from alignments of high-identity genes (from Step 3).  
  2. Map SNPs against each consensus sequence.  
  3. Summarize SNP positions and variability.  

- **Inputs**: Candidate inclusive gene alignments (Step 3).  
- **Outputs**: Consensus FASTA files, SNP summary tables.  

- **Scripts**:  
  - `build_consensus_sequences.py`  
  - `map_SNPs_vs_consensus.py`  
  - `split_SNP_summary.py`

- **Example command**:
```bash
python pipelines/panaroo/4_consensus/build_consensus_sequences.py \
       --input inclusivity_candidates.fasta \
       --output consensus_sequences.fasta

python pipelines/panaroo/4_consensus/map_SNPs_vs_consensus.py \
       --consensus consensus_sequences.fasta \
       --alignment core_gene_alignment.aln

python pipelines/panaroo/4_consensus/split_SN

## Step 5. Exclusivity analysis (BLAST)

- **Purpose**: Verify that candidate genes are unique to *S. pyogenes* and not present in near-neighbor species.

- **Process**:
  1. Run BLAST searches of consensus sequences (from Step 4) against a panel of near-neighbor genomes.  
  2. Parse BLAST results to identify unique or specific hits.  
  3. Summarize exclusivity scores for each candidate gene.  

- **Inputs**: Consensus sequences from Step 4, near-neighbor genomes (downloaded in Step 0).  
- **Outputs**: BLAST reports, exclusivity summary tables.  

- **Scripts**:  
  - `blast_consensus_exclusivity.py`  
  - `blast_output_analysis.py`  
  - `summarize_blast_results.py`

- **Example command**:
```bash
python pipelines/panaroo/5_exclusivity/blast_consensus_exclusivity.py \
       --input consensus_sequences.fasta \
       --db data/near_neighbors_db

python pipelines/panaroo/5_exclusivity/blast_output_analysis.py \
       --input blast_results.tsv \
       --output exclusivity_analysis.tsv

python pipelines/panaroo/5_exclusivity/summarize_blast_results.py \
       --input exclusivity_analysis.tsv \
       --outdir exclusivity_summaries/

## Step 6. Reporting

- **Purpose**: Generate human-readable reports, tables, and plots that summarize inclusivity and exclusivity results.

- **Process**:
  1. Collect outputs from inclusivity (Step 3), consensus/SNP analysis (Step 4), and exclusivity (Step 5).  
  2. Combine into a unified summary.  
  3. Format results into tables, plots, and a final report.  

- **Inputs**: Results from Steps 3–5.  
- **Outputs**: Final reports stored in `pipelines/panaroo/6_reports/`.  

- **Scripts** (examples, depending on implementation):  
  - `generate_summary_tables.py`  
  - `plot_results.py`  
  - `compile_final_report.py`

- **Example command**:
```bash
python pipelines/panaroo/6_reports/generate_summary_tables.py \
       --inclusivity inclusivity_candidates.tsv \
       --exclusivity exclusivity_analysis.tsv \
       --snps snp_summaries/ \
       --output summary_tables/

python pipelines/panaroo/6_reports/plot_results.py \
       --input summary_tables/ \
       --outdir plots/

python pipelines/panaroo/6_reports/compile_final_report.py \
       --tables summary_tables/ \
       --plots plots/ \
       --output final_report.pdf

## Status markers

Use this checklist to track progress as each step is verified.  
Mark items with `[x]` when complete.

- [ ] Step 0: Genome download & preparation verified  
- [ ] Step 1: Genome annotation (Prokka) verified  
- [ ] Step 2: Pan-genome construction (Panaroo) verified  
- [ ] Step 3: Inclusivity analysis scripts working  
- [ ] Step 4: Consensus building & SNP mapping working  
- [ ] Step 5: Exclusivity analysis (BLAST) working  
- [ ] Step 6: Reporting pipeline working
