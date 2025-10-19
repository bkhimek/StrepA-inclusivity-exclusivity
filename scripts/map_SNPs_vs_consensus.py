#!/usr/bin/env python3

from Bio import AlignIO, SeqIO
import csv
import os

# Paths
split_folder = "/home/himek/s_pyogenes_project/panaroo_output/core_gene_alignment.aln.split/"
consensus_folder = "/home/himek/s_pyogenes_project/consensus_sequences/"
output_csv = "/home/himek/s_pyogenes_project/core_gene_SNP_summary.csv"

# Prepare output list
output_lines = []
output_lines.append(["Gene", "Strain", "Position", "Consensus_Base", "Strain_Base"])

# Process each consensus file
for file in os.listdir(consensus_folder):
    if not file.endswith("_consensus.fasta"):
        continue

    gene = file.replace("_consensus.fasta", "")
    consensus_path = os.path.join(consensus_folder, file)
    alignment_path = os.path.join(split_folder, gene + ".fasta")

    if not os.path.exists(alignment_path):
        print(f"⚠️ Alignment for {gene} not found, skipping.")
        continue

    # Read consensus
    consensus_record = SeqIO.read(consensus_path, "fasta")
    consensus_seq = str(consensus_record.seq)

    # Read alignment
    alignment = AlignIO.read(alignment_path, "fasta")

    # Compare each strain to consensus
    for record in alignment:
        strain_seq = str(record.seq)
        strain_id = record.id

        for pos, (cons_base, strain_base) in enumerate(zip(consensus_seq, strain_seq), start=1):
            if cons_base != strain_base:
                output_lines.append([gene, strain_id, pos, cons_base, strain_base])

print(f"✅ Comparison complete. Found {len(output_lines)-1} SNPs.")

# Write output CSV
with open(output_csv, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(output_lines)

print(f"✅ SNP summary saved to: {output_csv}")
