#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
from Bio import SeqIO
import subprocess

def detect_sample_cols(df):
    cols = list(df.columns)
    if "Non-unique Gene name" in cols:
        start = cols.index("Non-unique Gene name") + 1
    elif "No. isolates" in cols:
        start = cols.index("No. isolates") + 1
    else:
        # Fallback: assume first 3 are metadata
        start = 3
    sample_cols = cols[start:]
    # keep only columns that look like presence/absence (strings / non-empty)
    return [c for c in sample_cols]

def run_mafft(in_fa: Path, out_fa: Path, threads: int):
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with out_fa.open("w") as out:
        subprocess.run(
            ["mafft", "--thread", str(threads), "--auto", str(in_fa)],
            stdout=out,
            check=True
        )

def main():
    ap = argparse.ArgumentParser(
        description="Split Panaroo core (>=min_presence) clusters into per-gene aligned FASTAs using MAFFT."
    )
    ap.add_argument("--gpa", required=True, help="gene_presence_absence.csv from Panaroo")
    ap.add_argument("--cds", required=True, help="combined_DNA_CDS.fasta from Panaroo")
    ap.add_argument("--outdir", required=True, help="Output directory for per-gene alignments (*.fasta)")
    ap.add_argument("--min_presence", type=float, default=0.99, help="Fraction of genomes required (e.g., 0.99)")
    ap.add_argument("--threads", type=int, default=8, help="MAFFT threads")
    args = ap.parse_args()

    gpa = pd.read_csv(args.gpa)
    sample_cols = detect_sample_cols(gpa)
    n_samples = len(sample_cols)
    if n_samples == 0:
        raise SystemExit("Could not detect sample columns in gene_presence_absence.csv")

    # presence fraction: count non-empty cells across sample columns
    def row_presence_frac(row):
        vals = row[sample_cols].astype(str)
        non_empty = (vals != "") & (vals != "nan")
        return non_empty.sum() / n_samples

    gpa["presence_frac"] = gpa.apply(row_presence_frac, axis=1)

    # choose clusters meeting threshold
    if "Gene" in gpa.columns:
        cluster_col = "Gene"
    else:
        cluster_col = gpa.columns[0]

    selected = gpa.loc[gpa["presence_frac"] >= args.min_presence, cluster_col].astype(str).tolist()
    print(f"[INFO] core≥{args.min_presence*100:.0f}% clusters: {len(selected)}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # index CDS fasta by first token of header (Panaroo uses cluster names like group_1234|sample|…)
    cds_index = {}
    for rec in SeqIO.parse(args.cds, "fasta"):
        first_token = rec.id.split("|")[0]
        cds_index.setdefault(first_token, []).append(rec)

    # for each selected cluster, write raw fasta then align with MAFFT
    written = 0
    missing = 0
    for gene in selected:
        recs = cds_index.get(gene, [])
        if len(recs) < 2:
            # need at least 2 sequences to align meaningfully
            missing += 1
            continue
        raw_path = outdir / f"{gene}.raw.fa"
        aln_path = outdir / f"{gene}.fasta"
        SeqIO.write(recs, raw_path, "fasta")
        run_mafft(raw_path, aln_path, args.threads)
        raw_path.unlink(missing_ok=True)
        written += 1

    print(f"[OK] wrote {written} aligned per-gene FASTAs to: {outdir}")
    if missing:
        print(f"[NOTE] skipped {missing} clusters with <2 sequences in CDS (likely parsing/naming differences)")

if __name__ == "__main__":
    main()
