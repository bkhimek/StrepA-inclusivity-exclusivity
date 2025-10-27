#!/usr/bin/env python3

import csv

input_csv = "/home/himek/s_pyogenes_project/core_gene_SNP_summary.csv"

# How many parts to split into
n_splits = 4

# Read all data
with open(input_csv, "r") as infile:
    reader = list(csv.reader(infile))
    header = reader[0]
    rows = reader[1:]

# Determine size of each part
chunk_size = len(rows) // n_splits + 1

# Write each part
for i in range(n_splits):
    part_rows = rows[i * chunk_size : (i + 1) * chunk_size]
    output_file = f"/home/himek/s_pyogenes_project/core_gene_SNP_summary_part_{i+1}.csv"
    with open(output_file, "w", newline="") as outfile:
        writer = csv.writer(outfile)
        writer.writerow(header)
        writer.writerows(part_rows)
    print(f"✅ Created: {output_file} with {len(part_rows)} SNP records.")

