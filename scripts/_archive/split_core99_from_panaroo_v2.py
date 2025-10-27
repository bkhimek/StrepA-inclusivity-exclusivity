#!/usr/bin/env python3
import argparse
from pathlib import Path
import subprocess
import pandas as pd
from Bio import SeqIO

def run_mafft(in_fa: Path, out_fa: Path, threads: int):
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    with out_fa.open("w") as out:
        subprocess.run(
            ["mafft", "--thread", str(threads), "--auto", str(in_fa)],
            stdout=out,
            check=True
        )

def sanitize(name: str) -> str:
    # safe filename from Gene name
    return "".join(c if c.isalnum() or c in ("_", "-", ".") else "_" for c in name)

def detect_sample_cols(gpa: pd.DataFrame):
    ann_cols = ["Gene", "Non-unique Gene name", "Annotation"]
    return [c for c in gpa.columns if c not in ann_cols]

def main():
    ap = argparse.ArgumentParser(
        description="Split Panaroo core (>=min_presence) into per-gene aligned FASTAs using GPA + gene_data mapping."
    )
    ap.add_argument("--gpa", required=True, help="gene_presence_absence.csv")
    ap.add_argument("--gene_data", required=True, help="gene_data.csv (has annotation_id <-> clustering_id)")
    ap.add_argument("--cds", required=True, help="combined_DNA_CDS.fasta")
    ap.add_argument("--outdir", required=True, help="Output dir for per-gene alignments (*.fasta)")
    ap.add_argument("--min_presence", type=float, default=0.99, help="Fraction of genomes required (e.g., 0.99)")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    # Load GPA and detect sample columns
    gpa = pd.read_csv(args.gpa)
    sample_cols = detect_sample_cols(gpa)
    n_samples = len(sample_cols)
    if n_samples == 0:
        raise SystemExit("No sample columns detected in GPA.")

    # presence fraction over sample columns (non-empty cell = present)
    def presence_frac(row):
        vals = row[sample_cols].astype(str)
        non_empty = (vals != "") & (vals != "nan") & vals.notna()
        return non_empty.sum() / n_samples

    gpa["presence_frac"] = gpa.apply(presence_frac, axis=1)

    # Select clusters meeting presence threshold
    gene_col = "Gene" if "Gene" in gpa.columns else gpa.columns[0]
    core_rows = gpa.loc[gpa["presence_frac"] >= args.min_presence, [gene_col] + sample_cols]

    # Load gene_data: annotation_id (locus tag) -> clustering_id (0_0_###)
    gd = pd.read_csv(args.gene_data, usecols=["annotation_id", "clustering_id"])
    ann2cid = dict(zip(gd["annotation_id"].astype(str), gd["clustering_id"].astype(str)))

    # Index CDS by clustering_id
    cid2rec = {}
    for rec in SeqIO.parse(args.cds, "fasta"):
        cid2rec[str(rec.id)] = rec

    written = 0
    skipped = 0
    for _, row in core_rows.iterrows():
        gene_name = str(row[gene_col])
        # Collect locus tags from all non-empty sample cells
        locus_tags = []
        for c in sample_cols:
            v = str(row[c])
            if v and v != "nan" and v.strip():
                # Cells usually contain one locus tag (e.g., LJEDHPMG_00803)
                # If multiple (rare), split on delimiters
                parts = [p.strip() for p in v.replace(";", ",").split(",") if p.strip()]
                locus_tags.extend(parts)

        # Map locus tags -> clustering_ids present in CDS
        cids = []
        for lt in locus_tags:
            cid = ann2cid.get(lt)
            if cid and cid in cid2rec:
                cids.append(cid)

        # Need >=2 sequences to align
        if len(set(cids)) < 2:
            skipped += 1
            continue

        raw_path = outdir / f"{sanitize(gene_name)}.raw.fa"
        aln_path = outdir / f"{sanitize(gene_name)}.fasta"

        SeqIO.write([cid2rec[c] for c in sorted(set(cids))], raw_path, "fasta")
        run_mafft(raw_path, aln_path, args.threads)
        raw_path.unlink(missing_ok=True)
        written += 1

    print(f"[OK] wrote {written} aligned per-gene FASTAs to: {outdir}")
    print(f"[NOTE] skipped {skipped} clusters (insufficient mapped sequences)")

if __name__ == "__main__":
    main()
