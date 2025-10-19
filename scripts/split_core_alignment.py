#!/usr/bin/env python3

from Bio import AlignIO
import os

# Set input and output paths
input_alignment = "/home/himek/s_pyogenes_project/panaroo_output_old/core_gene_alignment.aln"
output_folder = "/home/himek/s_pyogenes_project/panaroo_output_old/core_gene_alignment.aln.split/"

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Read the full core alignment
alignment = AlignIO.read(input_alignment, "fasta")

# How many strains (sequences) and how long the alignment is
n_sequences = len(alignment)
alignment_length = alignment.get_alignment_length()

# Read gene names from core_genes.csv
gene_list = []
with open("/home/himek/s_pyogenes_project/panaroo_output_old/core_genes.csv", "r") as f:
    next(f)  # skip header
    for line in f:
        gene_name = line.strip().split(",")[0]
        gene_list.append(gene_name)

# Check if gene number matches the number of windows we expect
expected_genes = alignment_length // 300  # Rough guess: ~300 bp per gene
print(f"Total alignment length: {alignment_length}")
print(f"Number of core genes expected: ~{len(gene_list)}")

# Now split the alignment into per-gene alignments
window_size = alignment_length // len(gene_list)

for idx, gene in enumerate(gene_list):
    start = idx * window_size
    end = (idx + 1) * window_size
    sliced_alignment = alignment[:, start:end]
    output_file = os.path.join(output_folder, f"{gene}.fasta")
    AlignIO.write(sliced_alignment, output_file, "fasta")

print(f"✅ Done! Core gene alignments saved to {output_folder}")
