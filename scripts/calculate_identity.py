#!/usr/bin/env python3

from Bio import AlignIO
import os

# Set paths
split_folder = "/home/himek/s_pyogenes_project/panaroo_output/core_gene_alignment.aln.split/"
output_file = "/home/himek/s_pyogenes_project/core_gene_identity_summary.csv"

# Prepare output
out = open(output_file, "w")
out.write("Gene,Percent_Identity\n")

# Loop over all split files
for filename in sorted(os.listdir(split_folder)):
    if filename.endswith(".fasta") or filename.endswith(".fa"):
        filepath = os.path.join(split_folder, filename)
        alignment = AlignIO.read(filepath, "fasta")
        aln_len = alignment.get_alignment_length()
        total_sites = 0
        identical_sites = 0

        # Compare column by column
        for i in range(aln_len):
            column = alignment[:, i]
            if "-" in column:
                continue  # Ignore columns with gaps
            total_sites += 1
            if len(set(column)) == 1:
                identical_sites += 1

        # Calculate percent identity
        if total_sites > 0:
            percent_identity = (identical_sites / total_sites) * 100
        else:
            percent_identity = 0.0

        # Save result
        gene_name = filename.replace(".fasta", "").replace(".fa", "")
        out.write(f"{gene_name},{percent_identity:.2f}\n")

out.close()

print("Finished! Output saved to:", output_file)
