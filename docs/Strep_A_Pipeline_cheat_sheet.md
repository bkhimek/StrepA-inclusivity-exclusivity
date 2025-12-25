# Strep A Pipeline Cheat Sheet

## 0) Environments You’ll Use

### Activate environments
```bash
# list available conda environments
conda env list

# activate main environment (with Biopython, pandas, BLAST+ etc.)
conda activate roary_env
For Panaroo-specific steps inside the wrapper:

conda run -n panaroo_env
1) One-Liner Run (Strep A, All Genomes)
Quick sanity check
# Count Prokka GFFs (one subfolder per genome)
ls pipelines/panaroo/1_annotate_prokka_spa/*/*.gff | wc -l
# expect ~391
Run full pipeline
bash scripts/run_all_spa.sh
2) What the Wrapper Does (Step-by-Step)
Panaroo (core = 100%)
Input: pipelines/panaroo/1_annotate_prokka_spa/*/*.gff

Output: pipelines/panaroo/2_panaroo/core100_spa/

Key files:

gene_presence_absence.csv

combined_DNA_CDS.fasta

gene_data.csv

Split per-gene & MAFFT alignment
# per-gene FASTA alignments
pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split/*.fasta
Per-gene identity (avg pairwise %ID)
Helper:

scripts/calculate_identity_with_names.py
Output:

pipelines/panaroo/3_inclusivity/core100_identity.tsv
Filter for ≥98% identity
scripts/filter_inclusivity_candidates.py
Produces:

core100_candidates.tsv

core100_candidates.nodup.tsv (de-duplicated)

Consensus & SNP summaries
Helpers:

scripts/build_consensus_from_split.py
scripts/summarize_snps_vs_consensus.py
Outputs:

core100_consensus.fasta

core100_snps.tsv

core100_pergene.tsv

Exclusivity BLAST vs near neighbors
BLAST DB:
pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/streptococcus_non_pyogenes_db.*

Run BLASTN → core100_vs_nonpyogenes.tsv

Summarize:

scripts/summarize_blast_exclusivity.py
→ core100_exclusivity.tsv

Extract PASS consensus FASTAs
scripts/extract_pass_consensus.py
→ pipelines/panaroo/5_exclusivity/PASS_FASTAs/

3) Key Inputs & Outputs
Inputs
pipelines/panaroo/1_annotate_prokka_spa/*/*.gff — Prokka GFFs for all S. pyogenes

BLAST DB at:
pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/streptococcus_non_pyogenes_db.*

If the BLAST DB is missing, rebuild it using makeblastdb.

4) Sanity Checks (Quick)
# GFF count
find pipelines/panaroo/1_annotate_prokka_spa -name '*.gff' | wc -l

# Split FASTAs present?
find pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split -name '*.fasta' | wc -l

# Identity table should have more than header
wc -l pipelines/panaroo/3_inclusivity/core100_identity.tsv

# Candidates should be non-empty
head pipelines/panaroo/3_inclusivity/core100_candidates.tsv

# Check duplicates removed
grep -c '\.raw$' pipelines/panaroo/3_inclusivity/core100_candidates.nodup.tsv

# Consensus headers (one per PASS gene)
grep -c '^>' pipelines/panaroo/4_consensus/core100_consensus.fasta

# Exclusivity summary counts
awk -F'\t' 'NR>1{c[$5]++} END{for(k in c)print k":",c[k]}' pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv

# PASS FASTAs count
find pipelines/panaroo/5_exclusivity/PASS_FASTAs -name '*.fasta' | wc -l
5) Common Fixes
A) BLAST “Database memory map file error”
Rebuild DB from combined FASTA (adjust path if needed):

DBDIR="pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes"
mkdir -p "$DBDIR"
makeblastdb \
  -in "$DBDIR/streptococcus_non_pyogenes_combined.fna" \
  -dbtype nucl \
  -out "$DBDIR/streptococcus_non_pyogenes_db" \
  -parse_seqids

# check with:
blastdbcmd -db "$DBDIR/streptococcus_non_pyogenes_db" -info | head
B) Empty identity table
Usually means zero split FASTAs → rerun split step or full wrapper.

C) Duplicate “*.raw” genes
Use core100_candidates.nodup.tsv downstream.

D) Panaroo “Error reading prokka input!”
If your GFF path layout differs:

# flat folder
panaroo -i pipelines/panaroo/1_annotate_prokka_spa/*.gff ...

# nested layout
panaroo -i pipelines/panaroo/1_annotate_prokka_spa/*/*.gff ...
Make sure the wrapper’s IN_GFF_DIR matches your structure.

6) Handy Partial Re-runs
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
  --input pipelines/panaroo/3_inclusivity/core100_candidates.tsv \
  --output pipelines/panaroo/3_inclusivity/core100_candidates.nodup.tsv
7) Output Copy Tips (WSL → Windows)
# PASS FASTAs
cp -v pipelines/panaroo/5_exclusivity/PASS_FASTAs/*.fasta \
  /mnt/c/Users/krist/OneDrive/Documents/.../outputs/PASS_FASTAs/

# exclusivity summary
cp -v pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv \
  /mnt/c/Users/krist/OneDrive/Documents/.../outputs/
8) Adapting to Another Species
For example, S. dysgalactiae:

Put GFFs under separate root:
pipelines/panaroo/1_annotate_prokka_sdy/*/*.gff

Set:

IN_GFF_DIR="pipelines/panaroo/1_annotate_prokka_sdy"
PANAROO_OUT="pipelines/panaroo/2_panaroo/core100_sdy"
bash scripts/run_all_spa.sh
Build a new non-species BLAST DB as required.

9) Git & Housekeeping
git add scripts/*.py scripts/*.sh docs/*.md
git commit -m "SPA: pipeline cheat sheet + docs"
git push origin main
Back up large inputs (ignored via .gitignore):

mkdir -p backups/$(date +%F_%H%M%S)
cp -r pipelines/panaroo/5_exclusivity/blastdb backups/$(date +%F_%H%M%S)/
Quick Reference — Helper Scripts
Script	Purpose
split_core100_from_panaroo.py	Extract & align core genes
calculate_identity_with_names.py	Compute avg %ID of alignments
filter_inclusivity_candidates.py	Filter genes by identity threshold
drop_raw_duplicates.py	Remove .raw duplicates
build_consensus_from_split.py	Build consensus FASTA per gene
summarize_snps_vs_consensus.py	SNP summary vs consensus
summarize_blast_exclusivity.py	Summarize BLAST exclusivity results
extract_pass_consensus.py	Export PASS gene FASTAs
