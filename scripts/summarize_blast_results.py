#!/usr/bin/env python

import pandas as pd

# Define BLAST column names (outfmt 6 custom)
columns = [
    "qseqid",    # Query sequence ID
    "sseqid",    # Subject sequence ID
    "pident",    # % identity
    "length",    # Alignment length
    "mismatch",  # Number of mismatches
    "gapopen",   # Number of gap openings
    "qstart",    # Start of alignment in query
    "qend",      # End of alignment in query
    "sstart",    # Start of alignment in subject
    "send",      # End of alignment in subject
    "evalue",    # Expect value
    "bitscore",  # Bit score
    "staxids",   # Subject Taxonomy ID(s)
    "sacc",      # Subject accession(s)
    "stitle"     # Subject title (contains organism name etc.)
]

blast_file = "exclusivity_blast_results.tsv"

print("Reading BLAST results...")
df = pd.read_csv(blast_file, sep="\t", header=None, names=columns)

print(f"Total hits found: {len(df)}")

# Optional: exclude exact matches (100% identity)
df = df[df["pident"] < 100]

# Sort and pick best hit for each query (lowest evalue, highest identity)
df_sorted = df.sort_values(["qseqid", "evalue", "pident"], ascending=[True, True, False])
best_hits = df_sorted.groupby("qseqid").first().reset_index()

# Save the summary CSV with clear headers
output_file = "exclusivity_blast_summary.csv"
best_hits.to_csv(output_file, index=False)

print(f"✅ Summary saved to {output_file} with headers.")
