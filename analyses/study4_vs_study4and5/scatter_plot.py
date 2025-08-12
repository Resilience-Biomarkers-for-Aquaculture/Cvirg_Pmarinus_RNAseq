import pandas as pd
import matplotlib.pyplot as plt

# Load the merged TSV with common genes
df = pd.read_csv("common_study5_study4and5_DN_treatment_control_inj_block.deseq2.results_filtered.tsv", sep="\t")

# Pivot to align log2FoldChange values by gene_id from both files
pivot = df.pivot(index="gene_id", columns="source_file", values="log2FoldChange")

# Drop any genes missing in either file (shouldn't happen if merge was complete)
pivot.dropna(inplace=True)

# Create scatter plot
plt.figure(figsize=(6, 6))
plt.scatter(pivot.iloc[:, 0], pivot.iloc[:, 1], alpha=0.7)

# Add reference lines
plt.axhline(0, color='gray', linestyle='--')
plt.axvline(0, color='gray', linestyle='--')

# Label axes based on filenames
plt.xlabel(f"log2FC in {pivot.columns[0]}")
plt.ylabel(f"log2FC in {pivot.columns[1]}")
plt.title("Common Genes: log2FoldChange Comparison")
plt.grid(True)
plt.tight_layout()

# Save to file
plt.savefig("log2fc_comparison.png", dpi=300)
print("Plot saved to log2fc_comparison.png")
