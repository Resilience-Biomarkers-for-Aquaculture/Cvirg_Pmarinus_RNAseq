"""
Direct from https://chatgpt.com/share/67ad0a07-f728-800d-9133-7a7f275a345aimport
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from matplotlib.lines import Line2D

### 1. Load Data ###
# Update paths as necessary
metadata_path = "path/to/further_augmented_metadata.csv"
expression_path = "path/to/merged_gene_counts.tsv"

# Load metadata
metadata = pd.read_csv(metadata_path)

# Load gene expression data
expression = pd.read_csv(expression_path, sep="\t")

# --------------------------
# 2️⃣ Prepare Data for Analysis
# --------------------------

# Select relevant metadata columns
metadata_filtered = metadata[['Study', 'Experiment', 'BREED', 'control', 'Trait']]
metadata_filtered.rename(columns={'Experiment': 'sample_id'}, inplace=True)

# Prepare gene expression data
expression.set_index('gene_id', inplace=True)  # Set gene_id as index
expression_transposed = expression.T  # Transpose so samples are rows
expression_transposed.index.name = 'sample_id'
expression_transposed.reset_index(inplace=True)

# Merge metadata with expression data
merged_data = metadata_filtered.merge(expression_transposed, on='sample_id', how='inner')


# --------------------------
# 3️⃣ Compute Gene-Level Differences (Averaged for Top Gene Selection)
# --------------------------

def compute_group_means(group):
    control_samples = group[group['control'] == 1].iloc[:, 4:]
    treatment_samples = group[group['control'] == 0].iloc[:, 4:]

    if control_samples.empty or treatment_samples.empty:
        return pd.Series(np.nan, index=control_samples.columns)

    return np.abs(control_samples.mean(numeric_only=True) - treatment_samples.mean(numeric_only=True))


# Compute group-level mean differences
gene_differences_grouped = merged_data.groupby(['Study', 'BREED']).apply(compute_group_means).reset_index()

# Drop NaN rows
gene_differences_grouped_clean = gene_differences_grouped.dropna(axis=0)

# Compute mean and std for selecting top 500 genes
numeric_cols = gene_differences_grouped_clean.select_dtypes(include=[float, int])
gene_stats = pd.DataFrame({
    'mean': numeric_cols.mean(),
    'std': numeric_cols.std()
})
gene_stats['consistency'] = np.abs(gene_stats['mean']) / (gene_stats['std'] + 1e-6)

# Select top 500 consistently varying genes
top_500_genes = gene_stats.sort_values(by='consistency', ascending=False).head(500)
top_gene_names = top_500_genes.index.tolist()

# --------------------------
# 4️⃣ Compute Per-Sample Differences for PCA
# --------------------------

per_sample_differences_list = []

for (study, breed), group in merged_data.groupby(['Study', 'BREED']):
    try:
        control_samples = group[group['control'] == 1].iloc[:, 4:]  # Extract control samples
        treatment_samples = group[group['control'] == 0]  # Extract treatment samples with metadata

        if control_samples.empty or treatment_samples.empty:
            continue  # Skip if no control or treatment samples exist

        # Compute mean of control samples for this group
        control_mean = control_samples.mean(numeric_only=True)

        # Process treatment samples in smaller batches to avoid memory overload
        chunk_size = 10
        for i in range(0, len(treatment_samples), chunk_size):
            treatment_chunk = treatment_samples.iloc[i:i + chunk_size, :].copy()

            # Compute per-sample differences (treatment - control mean)
            gene_diff = treatment_chunk.iloc[:, 4:] - control_mean

            # Combine metadata with computed differences
            result = pd.concat([treatment_chunk[['Study', 'BREED', 'sample_id', 'Trait']], gene_diff], axis=1)
            per_sample_differences_list.append(result)
    
    except MemoryError:
        print(f"Memory error processing {study} - {breed}. Skipping this group.")
        continue

# Concatenate processed data
per_sample_differences = pd.concat(per_sample_differences_list, ignore_index=True)

# Keep only the top 500 selected genes
pca_data = per_sample_differences[['Study', 'BREED', 'sample_id', 'Trait'] + top_gene_names].dropna()

# --------------------------
# 5️⃣ Perform PCA
# --------------------------

# Extract numerical data for PCA
X_pca = pca_data[top_gene_names]

# Perform PCA
pca = PCA(n_components=2)
pca_transformed = pca.fit_transform(X_pca)

# Add PCA results to dataframe
pca_data['PC1'] = pca_transformed[:, 0]
pca_data['PC2'] = pca_transformed[:, 1]

# --------------------------
# 6️⃣ Generate Improved PCA Plot with Colors & Markers
# --------------------------

# Define distinct colors for studies
study_colors = {
    study: color for study, color in zip(
        pca_data['Study'].unique(), ['red', 'blue', 'green', 'purple', 'orange', 'cyan']
    )
}

# Define marker styles for Trait
trait_markers = {
    "sensitive": "o",
    "tolerant": "x"
}

# Plot PCA results with study colors and trait markers
plt.figure(figsize=(10, 6))

for _, row in pca_data.iterrows():
    plt.scatter(
        row['PC1'], row['PC2'],
        color=study_colors[row['Study']],  # Color by study
        marker=trait_markers.get(str(row['Trait']).lower(), 'o'),  # Default to round if Trait is missing
        alpha=0.7, s=50  # Increase point size for visibility
    )

plt.xlabel("Principal Component 1")
plt.ylabel("Principal Component 2")
plt.title("PCA of Top 500 Genes (Per-Sample Differences)")

# Create a custom legend
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Sensitive', markerfacecolor='black', markersize=8),
    Line2D([0], [0], marker='x', color='w', label='Tolerant', markerfacecolor='black', markersize=8)
] + [
    Line2D([0], [0], marker='o', color=color, label=study, markersize=8) for study, color in study_colors.items()
]

plt.legend(handles=legend_elements, loc='upper right', fontsize='small', frameon=True)
plt.show()
