#!/usr/bin/env python3
import argparse
from pathlib import Path
from Bio import SeqIO
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--consensus_fasta", required=True, help="Consensus FASTA with all genes")
    ap.add_argument("--exclusivity_tsv", required=True, help="Exclusivity TSV with PASS/REJECT calls")
    ap.add_argument("--outdir", required=True, help="Output directory for split FASTA files")
    args = ap.parse_args()

    # Read PASS genes from exclusivity TSV
    excl = pd.read_csv(args.exclusivity_tsv, sep="\t")
    pass_genes = set(excl.loc[excl["decision"]=="PASS", "Gene"])
    print(f"[INFO] Found {len(pass_genes)} PASS genes")

    # Ensure output folder exists
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Read consensus fasta and write only PASS genes
    count = 0
    for rec in SeqIO.parse(args.consensus_fasta, "fasta"):
        if rec.id in pass_genes:
            out_path = outdir / f"{rec.id}.fasta"
            SeqIO.write(rec, out_path.open("w"), "fasta")
            count += 1

    print(f"[OK] Wrote {count} PASS consensus FASTA files -> {outdir}")

if __name__ == "__main__":
    main()
