import pandas as pd

# Load list of all genes
with open('selected_core_genes.txt') as f:
    all_genes = set(line.strip() for line in f if line.strip())

# Load the BLAST summary CSV
df = pd.read_csv('exclusivity_blast_summary.csv')

# Get the genes that had BLAST hits
genes_with_hits = set(df['qseqid'])

# Find genes with no hits
missing_genes = sorted(all_genes - genes_with_hits)

# Save the missing genes to a CSV
missing_df = pd.DataFrame({'Gene': missing_genes})
missing_df.to_csv('missing_genes.csv', index=False)

print(f"Found {len(missing_genes)} genes with no BLAST hits. Saved to missing_genes.csv.")
