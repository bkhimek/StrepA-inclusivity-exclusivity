🧪 Strep A Panaroo Pipeline — Cheat-Sheet
0) Quick glossary (folders & scripts)

Genomes (raw): data/all_genomes_raw/*.fna

Prokka GFFs: pipelines/panaroo/1_annotate_prokka_spa/*/*.gff

Panaroo outputs: pipelines/panaroo/2_panaroo/core100_spa/

Per-gene alignments: .../core_gene_alignment.aln.split/*.fasta

Inclusivity tables: pipelines/panaroo/3_inclusivity/core100_identity.tsv, core100_candidates.tsv

Consensus & SNPs: pipelines/panaroo/4_consensus/core100_consensus.fasta, core100_pergene.tsv, core100_snps.tsv

Exclusivity:
raw hits → pipelines/panaroo/5_exclusivity/core100_vs_nonpyogenes.tsv
summary → pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv
PASS FASTAs → pipelines/panaroo/5_exclusivity/PASS_FASTAs/

Wrapper (all steps): scripts/run_all_spa.sh

Key helpers (used by wrapper):
split_core100_from_panaroo.py, calculate_identity_with_names.py, filter_inclusivity_candidates.py, build_consensus_from_split.py, summarize_snps_vs_consensus.py, summarize_blast_exclusivity.py, extract_pass_consensus.py

1) Minimal “from zero to results” run

Runs Prokka ➜ Panaroo (core=100%) ➜ per-gene MAFFT ➜ identity ➜ filter (≥98%) ➜ consensus & SNPs ➜ BLAST to non-pyogenes ➜ PASS FASTAs.

# 1) Prokka annotate (from raw .fna)
conda activate prokka_env
mkdir -p pipelines/panaroo/1_annotate_prokka_spa
parallel -j 8 '
  B=$(basename {} .fna)
  prokka --outdir pipelines/panaroo/1_annotate_prokka_spa/${B} \
         --prefix ${B} --force --cpus 1 {}
' ::: data/all_genomes_raw/*.fna
conda deactivate

# 2) Run the full pipeline
bash scripts/run_all_spa.sh


Outputs are summarized at the end and mirrored into results/runs/run-YYYY-MM-DD_HHMMSS/.

2) Windows copy (WSL ⇄ OneDrive)

Update the Windows path if needed.

# PASS FASTAs
mkdir -p "/mnt/c/Users/krist/OneDrive/Documents/Bioinformatics/Bioinfo_instructions/StrepA_inclusivity_exclusivity/outputs/PASS_FASTAs"
cp -f pipelines/panaroo/5_exclusivity/PASS_FASTAs/*.fasta \
      "/mnt/c/Users/krist/OneDrive/Documents/Bioinformatics/Bioinfo_instructions/StrepA_inclusivity_exclusivity/outputs/PASS_FASTAs/"

# Exclusivity table
cp -f pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv \
      "/mnt/c/Users/krist/OneDrive/Documents/Bioinformatics/Bioinfo_instructions/StrepA_inclusivity_exclusivity/outputs/"

3) One-line sanity checks
# Prokka count
find pipelines/panaroo/1_annotate_prokka_spa -name '*.gff' | wc -l

# Panaroo must produce these 3:
ls -lh pipelines/panaroo/2_panaroo/core100_spa/{gene_presence_absence.csv,combined_DNA_CDS.fasta,gene_data.csv}

# Split FASTAs exist & have right strain count?
f=$(find pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split -name '*.fasta' ! -name '*.raw.fasta' | head -1)
echo "[FILE] $f"; grep -c '^>' "$f"

# Inclusivity table has many rows (not just header)
wc -l pipelines/panaroo/3_inclusivity/core100_identity.tsv
head -5 pipelines/panaroo/3_inclusivity/core100_candidates.tsv

# Consensus produced?
grep -c '^>' pipelines/panaroo/4_consensus/core100_consensus.fasta

# BLAST summary counts
awk -F'\t' 'NR>1 && $5=="PASS"{c++} END{print "PASS", c+0}' pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv
awk -F'\t' 'NR>1 && $5=="REJECT"{c++} END{print "REJECT", c+0}' pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv

4) Common fixes & gotchas
A) BLAST DB errors (“memory map file error”)

Rebuild the database from your combined FASTA:

mkdir -p pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes
cp -f backups/*/blastdb/non_pyogenes/streptococcus_non_pyogenes_combined.fna \
      pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/

makeblastdb \
  -in pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/streptococcus_non_pyogenes_combined.fna \
  -dbtype nucl \
  -out pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/streptococcus_non_pyogenes_db \
  -parse_seqids

blastdbcmd -db pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/streptococcus_non_pyogenes_db -info | head

B) Duplicate “.raw” sequences in candidates

If you see gene and gene.raw both in candidates:

python - <<'PY'
import pandas as pd
p="pipelines/panaroo/3_inclusivity/core100_candidates.tsv"
df=pd.read_csv(p,sep="\t")
df=df[~df['Gene_File'].str.endswith('.raw', na=False)]
df.to_csv(p.replace(".tsv",".nodup.tsv"),sep="\t",index=False)
print("[OK] wrote:", p.replace(".tsv",".nodup.tsv"), "rows:", len(df))
PY


Then rebuild consensus using the .nodup.tsv list (see §6).

C) “Error reading prokka input!”

This usually means the -i glob is wrong. If your GFFs are in a flat folder:

# In run_all_spa.sh, use:
-i pipelines/panaroo/1_annotate_prokka_spa/*/*.gff   # nested
# …or…
-i pipelines/panaroo/1_annotate_prokka_spa/*.gff     # flat


Match the pattern to your actual layout.

5) Tuning thresholds (quick)

Core presence (Panaroo): --core_threshold 1.00 → strict core. For permissive core, try 0.99 or 0.98.

Inclusivity identity: default ≥98%. Change with --threshold in filter_inclusivity_candidates.py.

Exclusivity “reject” rule: default pident ≥85 AND qcovs ≥80. Tweak in scripts/run_all_spa.sh (MIN_PIDENT, MIN_QCOVS) or pass to summarize_blast_exclusivity.py.

6) Run individual steps (handy if resuming)
# Split per-gene (MAFFT)
python scripts/split_core100_from_panaroo.py \
  --gpa       pipelines/panaroo/2_panaroo/core100_spa/gene_presence_absence.csv \
  --gene_data pipelines/panaroo/2_panaroo/core100_spa/gene_data.csv \
  --cds       pipelines/panaroo/2_panaroo/core100_spa/combined_DNA_CDS.fasta \
  --outdir    pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split \
  --min_presence 1.00 --threads 8

# Identities
python scripts/calculate_identity_with_names.py \
  --split_dir pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split \
  --output    pipelines/panaroo/3_inclusivity/core100_identity.tsv

# Filter ≥98%
python scripts/filter_inclusivity_candidates.py \
  --input  pipelines/panaroo/3_inclusivity/core100_identity.tsv \
  --output pipelines/panaroo/3_inclusivity/core100_candidates.tsv \
  --threshold 98.0

# (optional) drop .raw
# -> produces core100_candidates.nodup.tsv

# Consensus
python scripts/build_consensus_from_split.py \
  --split_dir pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split \
  --genes_txt pipelines/panaroo/3_inclusivity/core100_candidates.nodup.tsv \
  --out_fasta pipelines/panaroo/4_consensus/core100_consensus.fasta

# SNP summaries
python scripts/summarize_snps_vs_consensus.py \
  --split_dir      pipelines/panaroo/2_panaroo/core100_spa/core_gene_alignment.aln.split \
  --consensus_fasta pipelines/panaroo/4_consensus/core100_consensus.fasta \
  --out_tsv        pipelines/panaroo/4_consensus/core100_snps.tsv \
  --pergene_tsv    pipelines/panaroo/4_consensus/core100_pergene.tsv

# BLAST to non-pyogenes
DB="pipelines/panaroo/5_exclusivity/blastdb/non_pyogenes/streptococcus_non_pyogenes_db"
blastn -task megablast -db "$DB" \
  -query pipelines/panaroo/4_consensus/core100_consensus.fasta \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs sscinames staxids stitle" \
  -max_target_seqs 5 -evalue 1e-10 -num_threads 8 \
  > pipelines/panaroo/5_exclusivity/core100_vs_nonpyogenes.tsv

python scripts/summarize_blast_exclusivity.py \
  --blast_tsv pipelines/panaroo/5_exclusivity/core100_vs_nonpyogenes.tsv \
  --out_tsv   pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv \
  --min_pident 85 --min_qcovs 80

# Extract PASS FASTAs
python scripts/extract_pass_consensus.py \
  --consensus       pipelines/panaroo/4_consensus/core100_consensus.fasta \
  --exclusivity_tsv pipelines/panaroo/5_exclusivity/core100_exclusivity.tsv \
  --outdir          pipelines/panaroo/5_exclusivity/PASS_FASTAs

7) Switching species (e.g., S. dysgalactiae)

Put raw .fna in a new folder, e.g. data/sdysg/.

Prokka to a new root:

conda activate prokka_env
mkdir -p pipelines/panaroo/1_annotate_prokka_sdysg
parallel -j 8 '
  B=$(basename {} .fna)
  prokka --outdir pipelines/panaroo/1_annotate_prokka_sdysg/${B} \
         --prefix ${B} --force --cpus 1 {}
' ::: data/sdysg/*.fna
conda deactivate


Run pipeline with species-specific OUT folder:

IN_GFF_DIR="pipelines/panaroo/1_annotate_prokka_sdysg" \
PANAROO_OUT="pipelines/panaroo/2_panaroo/core100_sdysg" \
bash scripts/run_all_spa.sh


Exclusion DB: build a “non-sdysg” database (or reuse your “non-target streptococci” DB). Update BLAST_DB in the wrapper or pass it to the BLAST step.

8) Typical “why did X happen?” answers

Candidates fewer than expected: strict core (1.00) with many genomes can drop variable genes; try --core_threshold 0.99 if biologically acceptable.

Identity table tiny (41 bytes): no split FASTAs were found; re-run the split step; ensure *.fasta actually exist in the .split folder.

No species names in exclusivity: ensure -outfmt includes sscinames and/or stitle; our wrapper already uses this; otherwise rebuild BLAST step.

Duplicates (.raw): use the .nodup filter step before consensus.

“Error reading prokka input!”: your -i glob doesn’t match layout; switch between /*/*.gff (nested) and /*.gff (flat).

9) Git & docs quickies
# Add updated runbooks & scripts
git add docs/*.md scripts/*.sh scripts/*.py
git commit -m "Docs: all-genomes runbook + wrapper; pipeline stabilized"
git push origin main
