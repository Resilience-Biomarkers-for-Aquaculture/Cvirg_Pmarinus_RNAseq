"""
Directly from ChatGPT
https://chatgpt.com/share/67acfb36-4a34-800d-abc6-5d93e89cffba
"""

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import linkage, dendrogram
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns

# Load Data
gene_counts_path = "merged_gene_counts.tsv"  # Replace with your actual path
metadata_path = "further_augmented_metadata.csv"

gene_counts_df = pd.read_csv(gene_counts_path, sep='\t', index_col=0)
metadata_df = pd.read_csv(metadata_path)

# Column names
study_col = "Study"
family_col = "BREED"  # Assuming "BREED" represents family/species
treatment_col = "treatment"
experiment_col = "Experiment"

# Normalize Data (Counts Per Million - CPM)
total_counts_per_sample = gene_counts_df.sum(axis=0)
cpm_df = (gene_counts_df / total_counts_per_sample) * 1_000_000

# Filter out low-expression genes (CPM < 10 across all samples)
filtered_cpm_df = cpm_df[(cpm_df.sum(axis=1) >= 10)]

# Step 1: PCA within Each Study & Family and Compute Treatment-Control Distances
treatment_control_distances = {}

for study in metadata_df[study_col].unique():
    study_metadata = metadata_df[metadata_df[study_col] == study]
    study_samples = study_metadata[experiment_col].values

    # Subset the gene data for this study
    study_gene_data = filtered_cpm_df[study_samples].dropna(axis=0, how='all').T

    # Perform PCA
    pca = PCA(n_components=2)
    study_pca_result = pca.fit_transform(study_gene_data)

    # Create PCA DataFrame
    study_pca_df = pd.DataFrame(study_pca_result, columns=['PC1', 'PC2'], index=study_gene_data.index)
    study_pca_df = study_pca_df.merge(study_metadata[[experiment_col, family_col, treatment_col]], 
                                      left_index=True, right_on=experiment_col)

    for family in study_pca_df[family_col].unique():
        family_data = study_pca_df[study_pca_df[family_col] == family]
        
        control_samples = family_data[family_data[treatment_col] == "Control"]
        treatment_samples = family_data[family_data[treatment_col] != "Control"]
        
        if not control_samples.empty and not treatment_samples.empty:
            distances = pairwise_distances(treatment_samples[['PC1', 'PC2']], control_samples[['PC1', 'PC2']])
            avg_distances = np.mean(np.min(distances, axis=1))
            treatment_control_distances[(study, family)] = avg_distances

# Step 2: Compute Treatment-Control Distances in Full Gene Space
full_gene_space_distances = {}

for study in metadata_df[study_col].unique():
    study_metadata = metadata_df[metadata_df[study_col] == study]
    study_samples = study_metadata[experiment_col].values

    study_gene_data = filtered_cpm_df[study_samples].dropna(axis=0, how='all').T

    for family in study_metadata[family_col].unique():
        family_samples = study_metadata[study_metadata[family_col] == family][experiment_col].values
        family_data = study_gene_data.loc[family_samples]

        control_samples = study_metadata[(study_metadata[family_col] == family) & (study_metadata[treatment_col] == "Control")][experiment_col].values
        treatment_samples = study_metadata[(study_metadata[family_col] == family) & (study_metadata[treatment_col] != "Control")][experiment_col].values

        if len(control_samples) > 0 and len(treatment_samples) > 0:
            distances = pairwise_distances(family_data.loc[treatment_samples], family_data.loc[control_samples])
            avg_distances = np.mean(np.min(distances, axis=1))
            full_gene_space_distances[(study, family)] = avg_distances

# Step 3: Compute Gene-Level Treatment-Control Differences
gene_expression_differences = []

for study in metadata_df[study_col].unique():
    study_metadata = metadata_df[metadata_df[study_col] == study]
    study_samples = study_metadata[experiment_col].values

    study_gene_data = filtered_cpm_df[study_samples].dropna(axis=0, how='all')

    for family in study_metadata[family_col].unique():
        family_samples = study_metadata[study_metadata[family_col] == family][experiment_col].values
        family_data = study_gene_data[family_samples]

        control_samples = study_metadata[(study_metadata[family_col] == family) & (study_metadata[treatment_col] == "Control")][experiment_col].values
        treatment_samples = study_metadata[(study_metadata[family_col] == family) & (study_metadata[treatment_col] != "Control")][experiment_col].values

        if len(control_samples) > 0 and len(treatment_samples) > 0:
            avg_control_expression = family_data[control_samples].mean(axis=1)
            avg_treatment_expression = family_data[treatment_samples].mean(axis=1)
            expression_difference = avg_treatment_expression - avg_control_expression

            for gene, diff in expression_difference.items():
                gene_expression_differences.append({
                    "Study": study,
                    "Family": family,
                    "Gene": gene,
                    "Expression_Difference": diff
                })

gene_diff_df = pd.DataFrame(gene_expression_differences)

# Step 4: Identify Consistently Regulated Genes
gene_consistency_df = gene_diff_df.groupby("Gene")["Expression_Difference"].agg(["mean", "std"])
gene_consistency_df["abs_mean"] = gene_consistency_df["mean"].abs()
gene_consistency_df["consistency_score"] = gene_consistency_df["abs_mean"] / (gene_consistency_df["std"] + 1e-6)
top_consistent_genes = gene_consistency_df.nlargest(20, "consistency_score")

# Step 5: Cluster Gene Responses
gene_pivot_df = gene_diff_df.pivot(index="Gene", columns=["Study", "Family"], values="Expression_Difference").fillna(0)
linkage_matrix = linkage(gene_pivot_df, method="ward")

plt.figure(figsize=(10, 5))
dendrogram(linkage_matrix, labels=gene_pivot_df.index, leaf_rotation=90, leaf_font_size=8)
plt.title("Hierarchical Clustering of Genes Based on Expression Differences")
plt.ylabel("Distance")
plt.show()

# Step 6: PCA on Gene Differences
pca_genes = PCA(n_components=2)
pca_genes_result = pca_genes.fit_transform(gene_pivot_df)

pca_genes_df = pd.DataFrame(pca_genes_result, columns=["PC1", "PC2"], index=gene_pivot_df.index)

plt.figure(figsize=(8, 6))
sns.scatterplot(x="PC1", y="PC2", data=pca_genes_df)
plt.title("PCA of Gene-Level Treatment-Control Differences")
plt.xlabel("PC1")
plt.ylabel("PC2")
plt.show()
