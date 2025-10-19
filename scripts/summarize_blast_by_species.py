import pandas as pd

def process_blast(blast_file, output_file):
    df = pd.read_csv(blast_file)

    if "stitle" not in df.columns:
        raise ValueError(f"'stitle' column not found in {blast_file}!")

    # Ensure stitle is treated as string, even if missing values exist
    df["stitle"] = df["stitle"].astype(str)

    # Extract species name (first two words)
    df["species"] = df["stitle"].str.extract(r'^(\w+ \w+)', expand=False)

    summary = df.groupby(["qseqid", "species"]).agg(
        num_hits=("pident", "count"),
        best_identity=("pident", "max"),
        mean_identity=("pident", "mean")
    ).reset_index()

    summary.to_csv(output_file, index=False)
    print(f"Summary saved to: {output_file}")

# Process both files
print("Processing exclusivity_blast_summary.csv...")
process_blast("exclusivity_blast_summary.csv", "exclusivity_blast_summary_species.csv")

print("Processing missing_genes_blast_results.tsv...")
process_blast("missing_genes_blast_results.tsv", "missing_genes_blast_summary_species.csv")
