# Runbook: Panaroo-based Strep A Pipeline

This runbook captures the main workflow for inclusivity/exclusivity analysis of *Streptococcus pyogenes* using Prokka + Panaroo.

It is meant as a **living document**:  
- Every step shows *what is done*, *which scripts are used*, and *where outputs go*.  
- As you refine/correct scripts, update the notes here so this always reflects the working path.

---

## Step 0. Genome download & preparation

- **Purpose**: Obtain and prepare genome FASTA files for analysis.

- **Process**:
  1. **Download genomes**  
     - Assemblies of *Streptococcus pyogenes* and near-neighbor species are downloaded manually from NCBI Assembly.  
     - Manual download is currently more reliable than scripted batch download.
  2. **Folder structure from NCBI**  
     - Each strain is provided as a folder.  
     - Inside each folder is a file named `GCF_XXXXXX...` (generic name).
  3. **Renaming**  
     - The FASTA file is renamed to carry the strain name (taken from the folder name).  
     - This ensures filenames remain unique and informative.
  4. **Filtering**  
     - Folders without genomic sequence are removed.  
     - Only assemblies with complete sequence files are kept.
  5. **Final organization**  
     - All renamed FASTA files are moved into a single folder: `data/genomes/`.

- **Inputs**: Genome folders as downloaded from NCBI.  
- **Outputs**: Clean, renamed FASTA files in `data/genomes/`.

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
