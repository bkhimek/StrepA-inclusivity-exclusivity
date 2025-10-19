#!/usr/bin/env python3

from Bio import AlignIO
import os
import csv

# Correct Paths (now using panaroo_output)
gene_presence_absence_file = "/home/himek/s_pyogenes_project/panaroo_output/gene_presence_absence.csv"
split_folder = "/home/himek/s_pyogenes_project/panaroo_output/core_gene_alignment.aln.split/"
output_file = "/home/himek/s_pyogenes_project/core_gene_identity_with_names.csv"

# Read gene names and their corresponding file names
gene_name_mapping = {}
with open(gene_presence_absence_file, "r") as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        gene_name = row["Gene"]
        gene_name_mapping[gene_name] = gene_name

# Output list
output_lines = []

# Walk through split alignments
for file in os.listdir(split_folder):
    if file.endswith(".fasta") or file.endswith(".fa") or file.endswith(".aln"):
        gene = os.path.splitext(file)[0]
        path = os.path.join(split_folder, file)
        alignment = AlignIO.read(path, "fasta")
        
        # Calculate percent identity
        sequences = [str(rec.seq) for rec in alignment]
        length = len(sequences[0])
        identity_count = 0
        compare_count = 0
        
        for i in range(len(sequences)):
            for j in range(i+1, len(sequences)):
                seq1 = sequences[i]
                seq2 = sequences[j]
                matches = sum((c1 == c2) for c1, c2 in zip(seq1, seq2))
                identity = (matches / length) * 100
                identity_count += identity
                compare_count += 1
        
        avg_identity = identity_count / compare_count if compare_count > 0 else 0
        
        output_lines.append([gene, gene_name_mapping.get(gene, ""), f"{avg_identity:.2f}"])

# Save output
with open(output_file, "w", newline="") as outcsv:
    writer = csv.writer(outcsv)
    writer.writerow(["Gene_File", "Gene_Name", "Average_Identity(%)"])
    writer.writerows(output_lines)

print(f"✅ Core gene identity analysis finished! Results saved to: {output_file}")
