import pandas as pd
import glob
import os

input_folder = "exclusivity_blast_results/"
output_folder = "exclusivity_blast_results_filtered/"
os.makedirs(output_folder, exist_ok=True)

exclude_terms = ["Streptococcus pyogenes", "S. pyogenes", "Streptococcus_pyogenes"]

summary = []

for file in glob.glob(os.path.join(input_folder, "*.csv")):
    df = pd.read_csv(file)

    # Filter out rows containing excluded species names
    filtered_df = df[~df.apply(
        lambda row: row.astype(str).str.contains('|'.join(exclude_terms), case=False).any(),
        axis=1
    )]

    output_file = os.path.join(output_folder, os.path.basename(file))
    filtered_df.to_csv(output_file, index=False)

    # Add to summary
    gene_name = os.path.basename(file).replace(".csv", "")
    summary.append({"Gene": gene_name, "Remaining_hits": len(filtered_df)})

    print(f"Filtered {file} → {output_file}, remaining hits: {len(filtered_df)}")

# Save summary CSV
summary_df = pd.DataFrame(summary)
summary_csv = os.path.join(output_folder, "exclusivity_summary.csv")
summary_df.to_csv(summary_csv, index=False)

print(f"\n✅ Summary saved to: {summary_csv}")
