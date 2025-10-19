import pandas as pd

def summarize_blast(input_file, output_file):
    cols = ["qseqid","sseqid","pident","length","mismatch","gapopen",
            "qstart","qend","sstart","send","evalue","bitscore","staxids","sacc","stitle"]
    df = pd.read_csv(input_file, sep="\t", names=cols)

    # Extract species name from stitle (first two words)
    df["species"] = df["stitle"].str.extract(r'^(\w+ \w+)')

    # Get best hit per query per species (highest % identity)
    best_hits = df.sort_values(["qseqid","species","pident"], ascending=[True,True,False]) \
                  .drop_duplicates(subset=["qseqid","species"])

    # Save summary
    best_hits.to_csv(output_file, index=False)

# Summarize both BLAST results
summarize_blast("marker_vs_s_pyogenes.tsv", "marker_vs_s_pyogenes_summary.csv")
summarize_blast("marker_vs_non_pyogenes.tsv", "marker_vs_non_pyogenes_summary.csv")

print("✅ Summaries created: marker_vs_s_pyogenes_summary.csv and marker_vs_non_pyogenes_summary.csv")
