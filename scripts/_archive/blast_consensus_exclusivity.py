#!/usr/bin/env python3

import os
import subprocess
from Bio import SeqIO

# Parameters
consensus_folder = "/home/himek/s_pyogenes_project/consensus_sequences/"
blast_db = "/home/himek/blastdb/nt_prok"
output_folder = "/home/himek/s_pyogenes_project/exclusivity_blast_results/"

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Loop over all consensus FASTA files
fasta_files = [f for f in os.listdir(consensus_folder) if f.endswith(".fasta")]

for fasta_file in fasta_files:
    gene_name = fasta_file.replace("_consensus.fasta", "")
    fasta_path = os.path.join(consensus_folder, fasta_file)
    out_file = os.path.join(output_folder, f"{gene_name}_blast_results.csv")

    # BLAST command
    blast_cmd = [
        "blastn",
        "-query", fasta_path,
        "-db", blast_db,
        "-outfmt", "10 qseqid sseqid pident length evalue stitle",  # CSV output
        "-max_target_seqs", "10",
        "-evalue", "1e-5"
    ]

    print(f"🔎 Running BLAST for {gene_name}...")
    with open(out_file, "w") as out:
        subprocess.run(blast_cmd, stdout=out)

    print(f"✅ Results saved to {out_file}")

print("\n🎯 All BLAST searches completed. Check the exclusivity_blast_results folder.")
