#!/usr/bin/env python3

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

# Paths
split_folder = "/home/himek/s_pyogenes_project/panaroo_output/core_gene_alignment.aln.split/"
selected_genes_file = "/home/himek/s_pyogenes_project/selected_core_genes.txt"
output_folder = "/home/himek/s_pyogenes_project/consensus_sequences/"

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Read selected gene names
with open(selected_genes_file, "r") as f:
    selected_genes = [line.strip() for line in f if line.strip()]

# Process each selected gene
for gene in selected_genes:
    fasta_path = os.path.join(split_folder, gene + ".fasta")
    if not os.path.exists(fasta_path):
        print(f"⚠️ Warning: Alignment for {gene} not found, skipping.")
        continue

    alignment = AlignIO.read(fasta_path, "fasta")

    # Build consensus sequence
    consensus = ""
    alignment_length = alignment.get_alignment_length()
    for i in range(alignment_length):
        column = [rec.seq[i] for rec in alignment]
        base_counts = {}
        for base in column:
            base_counts[base] = base_counts.get(base, 0) + 1
        consensus += max(base_counts, key=base_counts.get)

    # Create SeqRecord for the consensus
    consensus_record = SeqRecord(Seq(consensus),
                                 id=gene,
                                 description="consensus sequence")

    # Save consensus FASTA for this gene (INSIDE THE LOOP)
    output_path = os.path.join(output_folder, gene + "_consensus.fasta")
    with open(output_path, "w") as out_fasta:
        SeqIO.write(consensus_record, out_fasta, "fasta")

print(f"✅ Consensus sequences saved to: {output_folder}")
