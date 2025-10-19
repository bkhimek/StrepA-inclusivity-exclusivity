#!/usr/bin/env python3

import csv

# Input file (your full table)
input_csv = "/home/himek/s_pyogenes_project/core_gene_identity_with_names.csv"

# Output file (just selected gene names)
output_txt = "/home/himek/s_pyogenes_project/selected_core_genes.txt"

# Filtering threshold
identity_cutoff = 98.0

# Process
selected_genes = []
with open(input_csv, "r") as infile:
    reader = csv.DictReader(infile)
    for row in reader:
        try:
            identity = float(row["Average_Identity(%)"])
            if identity >= identity_cutoff:
                selected_genes.append(row["Gene_File"])
        except ValueError:
            continue  # skip bad rows

# Save output
with open(output_txt, "w") as outfile:
    for gene in selected_genes:
        outfile.write(gene + "\n")

print(f"✅ Selected {len(selected_genes)} core genes with ≥{identity_cutoff}% identity.")
print(f"✅ Saved to: {output_txt}")
