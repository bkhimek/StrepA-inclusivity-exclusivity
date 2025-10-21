#!/usr/bin/env python3
import argparse
from pathlib import Path
from Bio import AlignIO
import csv

def pct_identity_of_alignment(aln):
    seqs = [str(rec.seq) for rec in aln]
    if len(seqs) < 2:
        return 100.0 if len(seqs) == 1 else 0.0
    L = len(seqs[0])
    total = 0.0
    pairs = 0
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            matches = sum(1 for a, b in zip(seqs[i], seqs[j]) if a == b)
            total += (matches / L) * 100.0
            pairs += 1
    return total / pairs if pairs else 0.0

def main():
    ap = argparse.ArgumentParser(
        description="Compute average pairwise identity (%) per gene from split core alignments."
    )
    ap.add_argument("--split_dir", required=True,
                    help="Directory of per-gene FASTA/ALN files (core_gene_alignment.aln.split/)")
    ap.add_argument("--output", required=True,
                    help="Output TSV path")
    args = ap.parse_args()

    split_dir = Path(args.split_dir)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    files = sorted([p for p in split_dir.iterdir()
                    if p.suffix.lower() in (".fasta", ".fa", ".aln")])

    rows = []
    for fp in files:
        gene = fp.stem
        try:
            aln = AlignIO.read(str(fp), "fasta")
        except Exception:
            continue
        pid = pct_identity_of_alignment(aln)
        rows.append((gene, gene, f"{pid:.2f}"))

    with out_path.open("w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["Gene_File", "Gene_Name", "Average_Identity(%)"])
        w.writerows(rows)

    print(f"✅ Core gene identity analysis finished! Results saved to: {out_path}")
    print(f"✅ Genes processed: {len(rows)}")

if __name__ == "__main__":
    main()
