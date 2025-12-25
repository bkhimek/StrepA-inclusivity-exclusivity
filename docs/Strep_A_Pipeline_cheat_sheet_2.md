# Strep A Panaroo Pipeline — Cheat Sheet (Part 2)

---

## 0) Quick Glossary — Key Paths & Scripts

**Folders & Outputs**
- **Raw genomes:** `data/all_genomes_raw/*.fna`
- **Prokka GFFs:** `pipelines/panaroo/1_annotate_prokka_spa/*/*.gff`
- **Panaroo outputs:** `pipelines/panaroo/2_panaroo/core100_spa/`
- **Per-gene alignments:** `.../core_gene_alignment.aln.split/*.fasta`
- **Inclusivity tables:**  
  - `pipelines/panaroo/3_inclusivity/core100_identity.tsv`  
  - `pipelines/panaroo/3_inclusivity/core100_candidates.tsv`
- **Consensus & SNPs:**  
  - `pipelines/panaroo/4_consensus/core100_consensus.fasta`  
  - `pipelines/panaroo/4_consensus/core100_pergene.tsv`  
  - `pipelines/panaroo/4_consensus/core100_snps.tsv`
- **Exclusivity tables:**  
  - Raw BLAST: `pipelines/panaroo/5_exclusivity/core100_vs_nonpyogenes.tsv`  
  - Summary: `pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv`
- **PASS FASTAs:** `pipelines/panaroo/5_exclusivity/PASS_FASTAs/`

**Main Wrapper**
```bash
bash scripts/run_all_spa.sh
Helper Scripts

split_core100_from_panaroo.py

calculate_identity_with_names.py

filter_inclusivity_candidates.py

build_consensus_from_split.py

summarize_snps_vs_consensus.py

summarize_blast_exclusivity.py

extract_pass_consensus.py

1) Minimal “From Zero to Results” Run
Run the full pipeline
# Starts from raw FASTA, runs Prokka → Panaroo → identity → consensus → exclusivity
bash scripts/run_all_spa.sh
2) Copy Outputs (WSL → Windows OneDrive)
PASS FASTAs
mkdir -p "/mnt/c/Users/krist/OneDrive/Documents/Bioinformatics/Bioinfo_instructions/StrepA_inclusivity_exclusivity/outputs/PASS_FASTAs"
cp -f pipelines/panaroo/5_exclusivity/PASS_FASTAs/*.fasta \
  "/mnt/c/Users/krist/OneDrive/Documents/Bioinformatics/Bioinfo_instructions/StrepA_inclusivity_exclusivity/outputs/PASS_FASTAs/"
Exclusivity Summary Table
cp -f pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv \
  "/mnt/c/Users/krist/OneDrive/Documents/Bioinformatics/Bioinfo_instructions/StrepA_inclusivity_exclusivity/outputs/"
3) One-Line Sanity Checks
# Prokka GFF count
find pipelines/panaroo/1_annotate_prokka_spa -name '*.gff' | wc -l

# Panaroo expected outputs
ls -lh pipelines/panaroo/2_panaroo/core100_spa/{gene_presence_absence.csv,combined_DNA_CDS.fasta,gene_data.csv}

# Check split FASTA existence & strain count
f=$(find pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split -name '*.fasta' ! -name '*.raw.fasta' | head -1)
echo "[FILE] $f"; grep -c '^>' "$f"

# Identity table should be larger than just the header
wc -l pipelines/panaroo/3_inclusivity/core100_identity.tsv

# First few candidate lines
head -5 pipelines/panaroo/3_inclusivity/core100_candidates.tsv

# Consensus sequences count
grep -c '^>' pipelines/panaroo/4_consensus/core100_consensus.fasta

# BLAST summary counts
awk -F'\t' 'NR>1 && $5=="PASS"{c++} END{print "PASS", c+0}' pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv
awk -F'\t' 'NR>1 && $5=="REJECT"{c++} END{print "REJECT", c+0}' pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv
4) Common Fixes & Gotchas
A) BLAST DB Errors (“memory map file error”)
Rebuild non-pyogenes BLAST DB:

mkdir -p pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes
cp -f backups/*/blastdb/non_pyogenes/streptococcus_non_pyogenes_combined.fna \
  pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/

makeblastdb \
  -in pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/streptococcus_non_pyogenes_combined.fna \
  -dbtype nucl \
  -out pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/streptococcus_non_pyogenes_db \
  -parse_seqids

blastdbcmd -db pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/streptococcus_non_pyogenes_db -info | head
5) Handling Duplicates (.raw sequences)
If a gene and gene.raw both appear in candidates:

python - <<'PY'
import pandas as pd
p="pipelines/panaroo/3_inclusivity/core100_candidates.tsv"
df=pd.read_csv(p,sep="\t")
df=df[~df['Gene_File'].str.endswith('.raw', na=False)]
df.to_csv(p.replace(".tsv",".nodup.tsv"),sep="\t",index=False)
print("[OK] wrote:", p.replace(".tsv",".nodup.tsv"), "rows:", len(df))
PY
Then rebuild consensus using the *.nodup.tsv list.

6) “Error reading prokka input!” Troubleshooting
If Panaroo complains that your GFF glob pattern is wrong:

For nested subfolders:

-i pipelines/panaroo/1_annotate_prokka_spa/*/*.gff
For flat structure:

-i pipelines/panaroo/1_annotate_prokka_spa/*.gff
Match the glob to how your GFFs are organized.

7) Threshold Tuning
Panaroo core presence:
--core_threshold 1.00 is strict. For more permissive cores try 0.99 or 0.98.

Inclusivity identity threshold:
≥98% (default).
Change with --threshold in filter_inclusivity_candidates.py.

Exclusivity “reject” filters:
Default: percent identity ≥85 AND query coverage ≥80.
Tweak in wrapper or in summarize_blast_exclusivity.py.

8) Run Individual Steps (Handy for Resuming)
# Split per-gene (MAFFT)
python scripts/split_core100_from_panaroo.py \
  --gpa pipelines/panaroo/2_panaroo/core100_spa/gene_presence_absence.csv \
  --gene_data pipelines/panaroo/2_panaroo/core100_spa/gene_data.csv \
  --cds pipelines/panaroo/2_panaroo/core100_spa/combined_DNA_CDS.fasta \
  --outdir pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split \
  --min_presence 1.00 --threads 8

# Identity table
python scripts/calculate_identity_with_names.py \
  --split_dir pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split \
  --output pipelines/panaroo/3_inclusivity/core100_identity.tsv

# Filter ≥98%
python scripts/filter_inclusivity_candidates.py \
  --input pipelines/panaroo/3_inclusivity/core100_identity.tsv \
  --output pipelines/panaroo/3_inclusivity/core100_candidates.tsv \
  --threshold 98.0

# Consensus
python scripts/build_consensus_from_split.py \
  --split_dir pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split \
  --genes_txt pipelines/panaroo/3_inclusivity/core100_candidates.nodup.tsv \
  --out_fasta pipelines/panaroo/4_consensus/core100_consensus.fasta

# SNP summaries
python scripts/summarize_snps_vs_consensus.py \
  --split_dir pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split \
  --consensus_fasta pipelines/panaroo/4_consensus/core100_consensus.fasta \
  --out_tsv pipelines/panaroo/4_consensus/core100_snps.tsv \
  --pergene_tsv pipelines/panaroo/4_consensus/core100_pergene.tsv

# BLAST exclusivity
DB="pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/streptococcus_non_pyogenes_db"
blastn -task megablast -db "$DB" \
  -query pipelines/panaroo/4_consensus/core100_consensus.fasta \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs sscinames staxids stitle" \
  -max_target_seqs 5 -evalue 1e-10 -num_threads 8 \
  > pipelines/panaroo/5_exclusivity/core100_vs_nonpyogenes.tsv

python scripts/summarize_blast_exclusivity.py \
  --blast_tsv pipelines/panaroo/5_exclusivity/core100_vs_nonpyogenes.tsv \
  --out_tsv pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv \
  --min_pident 85 --min_qcovs 80

# Extract PASS FASTAs
python scripts/extract_pass_consensus.py \
  --consensus pipelines/panaroo/4_consensus/core100_consensus.fasta \
  --exclusivity_tsv pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv \
  --outdir pipelines/panaroo/5_exclusivity/PASS_FASTAs
9) Switching Species (e.g., S. dysgalactiae)
Place raw .fna in a separate folder (e.g., data/sdysg/)

Run Prokka into a new root:

conda activate prokka_env
mkdir -p pipelines/panaroo/1_annotate_prokka_sdysg
parallel -j 8 \
  'B=$(basename {} .fna) prokka --outdir pipelines/panaroo/1_annotate_prokka_sdysg/${B} \
   --prefix ${B} --force --cpus 1 {}' ::: data/sdysg/*.fna
conda deactivate
In wrapper:

IN_GFF_DIR="pipelines/panaroo/1_annotate_prokka_sdysg"
PANAROO_OUT="pipelines/panaroo/2_panaroo/core100_sdysg"
bash scripts/run_all_spa.sh
Build a “non-target” exclusion database accordingly.

10) Git & Docs Quickies
git add docs/*.md scripts/*.sh scripts/*.py
git commit -m "Docs: update pipeline cheat sheet 2"
git push origin main
