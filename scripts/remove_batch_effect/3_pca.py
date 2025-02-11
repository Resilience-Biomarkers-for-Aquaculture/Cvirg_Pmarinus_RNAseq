# Reload necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import seaborn as sns
from sklearn.decomposition import PCA

matplotlib.use("Agg")  # Force non-GUI mode (prevents Qt errors)

# Load batch-corrected expression data
working_dir = "/mnt/c/Users/inter/Downloads/GMGI/Perkinsus batch correction"
output_dir = os.path.join(working_dir, "pca_plots")
os.makedirs(output_dir, exist_ok=True)  # Create output directory if not exists

batch_corrected_file = os.path.join(working_dir, "batch_corrected_expression.csv")
batch_corrected_data = pd.read_csv(batch_corrected_file, index_col=0)

# Load metadata
metadata_file = os.path.join(working_dir, "augmented_metadata.csv")
metadata_new = pd.read_csv(metadata_file)

# Transpose the expression matrix (samples as rows, genes as columns)
batch_corrected_transposed = batch_corrected_data.T

# Perform PCA on batch-corrected data
pca = PCA(n_components=2)
pca_result = pca.fit_transform(batch_corrected_transposed)

# Merge PCA results with metadata
pca_df_corrected = pd.DataFrame(pca_result, columns=["PC1", "PC2"], index=batch_corrected_transposed.index)
pca_df_corrected = pca_df_corrected.merge(metadata_new[['Experiment', 'Study', 'control', 'Collection_Interval_Days', 'Trait' ]], left_index=True, right_on='Experiment')

# Create a colormap for the unique Study values
unique_studies = pca_df_corrected["Study"].astype('category').cat.categories
study_colors = plt.cm.viridis(np.linspace(0, 1, len(unique_studies)))

# Create a mapping of Study to colors
study_color_map = dict(zip(unique_studies, study_colors))

def pca_plot(pca, pca_df, palette, style_selector, unique_studies):
    # Plot PCA after batch correction (colored by Study) with legend
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Study",
                    palette=palette,
                    style=style_selector, alpha=0.7)

    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.2f}% Variance)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.2f}% Variance)')
    plt.title(f"PCA (colored by study, styled by {style_selector}")

    # Create legend
    legend_patches = [plt.Line2D([0], [0], marker='o', color='w',
                                 markerfacecolor=study_color_map[study],
                                 markersize=8, label=study) for study in unique_studies]
    plt.legend(handles=legend_patches, title="Study", loc="best")

    plt.savefig(os.path.join(output_dir, f"PCA_Batch_Corrected_{style_selector}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()  # Close plot to avoid GUI issues

pca_plot(pca, pca_df_corrected, study_color_map, 'control', unique_studies)
pca_plot(pca, pca_df_corrected, study_color_map, 'Collection_Interval_Days', unique_studies)
pca_plot(pca, pca_df_corrected, study_color_map, 'Trait', unique_studies)

