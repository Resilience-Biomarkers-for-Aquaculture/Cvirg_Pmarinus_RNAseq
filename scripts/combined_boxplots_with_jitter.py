#!/usr/bin/env python3
"""
Combined Boxplots with Jitter Points for Specific Genes Analysis

Creates comprehensive boxplot visualizations with overlaid jitter points showing 
individual sample expression values for genes of interest across multiple datasets.
Excludes Dataset 2 from analysis.

Features:
- Combined boxplots with dataset-specific coloring
- Individual jitter points for all samples (black)
- Statistical summaries and fold change calculations
- Study-grouped visualizations
- Individual dataset breakdowns

Author: Analysis Pipeline
Date: 2025-10-16
"""

import pandas as pd
import math
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

def load_specific_datasets_excluding_study2(metadata_file, genes_of_interest):
    """
    Load and combine data from datasets 1, 3, 4, 5 (excluding study 2).
    
    Args:
        metadata_file (str): Path to metadata file
        genes_of_interest (list): List of genes to extract
        
    Returns:
        tuple: (combined_gene_data, combined_metadata)
    """
    print("Loading data from datasets 1, 3, 4, 5 (excluding Study 2)...")
    
    # Load metadata
    metadata = pd.read_csv(metadata_file)
    
    # Dataset file mapping (excluding dataset 2)
    dataset_files = {
        1: 'data/rnaseq_gene_counts/salmon.merged.gene_counts_length_scaled_dataset1.tsv',
        3: 'data/rnaseq_gene_counts/salmon.merged.gene_counts_length_scaled_dataset3.tsv',
        4: 'data/rnaseq_gene_counts/salmon.merged.gene_counts_length_scaled_dataset4.tsv',
        5: 'data/rnaseq_gene_counts/salmon.merged.gene_counts_length_scaled_dataset5.tsv'
    }
    
    # BioProject mapping
    dataset_bioproject_map = {
        1: 'PRJNA894694',
        3: 'PRJNA778545',
        4: 'PRJNA691949',
        5: 'PRJNA590205'
    }
    
    combined_data = []
    combined_metadata = []
    
    for dataset_num, file_path in dataset_files.items():
        try:
            print(f"Loading Dataset {dataset_num}...")
            data = pd.read_csv(file_path, sep='\t', index_col=0)
            
            # Filter for genes of interest
            available_genes = set(genes_of_interest).intersection(set(data.index))
            if not available_genes:
                print(f"  Warning: No genes of interest found in Dataset {dataset_num}")
                continue
                
            gene_data = data.loc[list(available_genes)]
            
            # Get metadata for this dataset
            bioproject = dataset_bioproject_map[dataset_num]
            dataset_metadata = metadata[metadata['BioProject'] == bioproject].copy()
            dataset_metadata['Dataset'] = dataset_num
            
            # Find overlapping samples
            data_samples = set(gene_data.columns)
            meta_samples = set(dataset_metadata['Experiment'])
            overlap = data_samples.intersection(meta_samples)
            
            if overlap:
                # Filter to overlapping samples
                gene_data = gene_data[list(overlap)]
                dataset_metadata = dataset_metadata[dataset_metadata['Experiment'].isin(overlap)]
                
                print(f"  Dataset {dataset_num}: {len(overlap)} samples, {len(available_genes)} genes")
                print(f"  BioProject: {bioproject}")
                
                combined_data.append(gene_data)
                combined_metadata.append(dataset_metadata)
            else:
                print(f"  Warning: No sample overlap for Dataset {dataset_num}")
                
        except FileNotFoundError:
            print(f"  Warning: File not found for Dataset {dataset_num}: {file_path}")
        except Exception as e:
            print(f"  Error loading Dataset {dataset_num}: {e}")
    
    if not combined_data:
        raise ValueError("No datasets loaded successfully")
    
    # Combine all data
    combined_gene_data = pd.concat(combined_data, axis=1)
    combined_metadata_df = pd.concat(combined_metadata, ignore_index=True)
    
    return combined_gene_data, combined_metadata_df

def prepare_plot_data(gene_data, metadata, genes_of_interest):
    """
    Prepare data for plotting by melting into long format.
    
    Args:
        gene_data (pd.DataFrame): Gene expression data
        metadata (pd.DataFrame): Sample metadata
        genes_of_interest (list): List of genes
        
    Returns:
        pd.DataFrame: Long format data ready for plotting
    """
    print("Preparing plot data...")
    
    # Create sample info lookup
    sample_info = metadata.set_index('Experiment')[['Trait', 'BioProject', 'Dataset']]
    
    # Dataset labels for cleaner display
    dataset_labels = {
        1: 'Dataset 1',
        3: 'Dataset 3', 
        4: 'Dataset 4',
        5: 'Dataset 5'
    }
    
    plot_data = []
    
    for gene in genes_of_interest:
        if gene not in gene_data.index:
            print(f"Warning: {gene} not found in data, skipping...")
            continue
            
        for sample in gene_data.columns:
            if sample in sample_info.index:
                sample_meta = sample_info.loc[sample]
                plot_data.append({
                    'Gene': gene,
                    'Sample': sample,
                    'Expression': gene_data.loc[gene, sample],
                    'Trait': sample_meta['Trait'],
                    'BioProject': sample_meta['BioProject'],
                    'Dataset': sample_meta['Dataset'],
                    'Dataset_Label': dataset_labels.get(sample_meta['Dataset'], f"Dataset {sample_meta['Dataset']}")
                })
    
    return pd.DataFrame(plot_data)

def create_combined_boxplots_with_jitter(plot_df, output_dir):
    """
    Create combined boxplots with jitter points for all genes.
    
    Args:
        plot_df (pd.DataFrame): Plot data in long format with columns:
                               ['Gene', 'Sample', 'Expression', 'Trait', 'Dataset_Label']
        output_dir (str): Output directory path for saving plots
        
    Returns:
        matplotlib.figure.Figure: The created figure object
    """
    print("Creating combined boxplots with jitter points...")
    
    # Set up the plot style
    plt.style.use('default')
    sns.set_palette("Set2")
    
    # Create figure with subplots for each gene
    n_genes = len(plot_df['Gene'].unique())
    n_cols = 3
    n_rows = (n_genes + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 6*n_rows))
    if n_rows == 1:
        axes = [axes] if n_cols == 1 else axes
    else:
        axes = axes.flatten()
    
    # Color palette for datasets (using seaborn's Set2)
    
    genes = sorted(plot_df['Gene'].unique())
    
    for i, gene in enumerate(genes):
        ax = axes[i]
        gene_data = plot_df[plot_df['Gene'] == gene]
        
        # Create boxplot with dataset as hue
        box_plot = sns.boxplot(
            data=gene_data, 
            x='Trait', 
            y='Expression',
            hue='Dataset_Label',
            ax=ax,
            palette='Set2',
            showfliers=False  # Don't show outliers since we'll show all points as jitter
        )
        
        # Add jitter points (strip plot) colored black
        strip_plot = sns.stripplot(
            data=gene_data,
            x='Trait', 
            y='Expression',
            hue='Dataset_Label',
            ax=ax,
            color='black',
            size=4,
            alpha=0.8,
            dodge=True,  # Separate by hue to match boxplot positions
            jitter=0.3
        )
        
        # Customize subplot
        ax.set_title(f'{gene}', fontsize=14, fontweight='bold', pad=20)
        ax.set_xlabel('Trait', fontsize=12)
        ax.set_ylabel('Length-Scaled Counts', fontsize=12)
        
        # Add sample counts (overall totals)
        total_sensitive = gene_data[gene_data['Trait'] == 'sensitive'].shape[0]
        total_tolerant = gene_data[gene_data['Trait'] == 'tolerant'].shape[0]

        count_text = f"n_sensitive={total_sensitive}, n_tolerant={total_tolerant}"
        ax.text(0.02, 0.98, count_text, transform=ax.transAxes, 
               fontsize=10, verticalalignment='top', 
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        # Calculate and display fold change
        sensitive_mean = gene_data[gene_data['Trait'] == 'sensitive']['Expression'].mean()
        tolerant_mean = gene_data[gene_data['Trait'] == 'tolerant']['Expression'].mean()

        if sensitive_mean > 0 and tolerant_mean > 0:
         fold_change = tolerant_mean / sensitive_mean
         log2_fc = math.log2(fold_change)
         fc_text = f"FC: {fold_change:.2f} (Log2: {log2_fc:.2f})"
         ax.text(0.98, 0.98, fc_text, transform=ax.transAxes, 
             fontsize=9, verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        
        # Keep legend for boxplot colors but remove stripplot legend
        handles, labels = ax.get_legend_handles_labels()
        # Only keep the first 4 handles (boxplot colors), remove stripplot handles
        if len(handles) > 4:
            # Position legend just below the fold change text in upper right
            ax.legend(handles[:4], labels[:4], title='Dataset', 
                     bbox_to_anchor=(1.0, 0.85), loc='upper right', fontsize=8)
        else:
            ax.legend(title='Dataset', bbox_to_anchor=(1.0, 0.85), 
                     loc='upper right', fontsize=8)
    
    # Hide empty subplots
    for i in range(len(genes), len(axes)):
        axes[i].set_visible(False)
    
    # Main title with more space from facet titles
    fig.suptitle('Gene Expression by Trait and Dataset', 
                 fontsize=16, fontweight='bold', y=0.95)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.88)
    
    # Save single PNG output (user requested only this file)
    out_png = f"{output_dir}/combined_boxplots_with_jitter_specific.png"
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    print(f"Saved combined plot to {out_png}")
    
    plt.close()
    
    return fig

# Note: study-grouped and individual study plotting helpers removed to keep this
# script minimal. Re-enable or restore from version control if you need them later.

# Note: Summary statistics helper removed to keep script minimal. Restore from
# version control if you need CSV outputs or fold-change tables.
    
def main():
    """Main function."""
    # Your specific genes of interest
    genes_of_interest = ['LOC111103066', 'LOC111133512', 'LOC111109496', 
                        'LOC111116766', 'LOC111125640', 'LOC111102778']
    
    metadata_file = 'data/merged_metadata.csv'
    output_dir = Path('analyses/20251030_gene_visualization_output')
    output_dir.mkdir(exist_ok=True)
    
    print("=== Combined Boxplots with Jitter Points (Specific Genes, Study 2 Excluded) ===")
    print(f"Genes: {genes_of_interest}")
    print(f"Output directory: {output_dir}")
    
    try:
        # Load data
        gene_data, metadata = load_specific_datasets_excluding_study2(metadata_file, genes_of_interest)
        
        print(f"\nCombined dataset summary:")
        print(f"Total samples: {gene_data.shape[1]}")
        print(f"Total genes: {gene_data.shape[0]}")
        print(f"Trait distribution: {metadata['Trait'].value_counts().to_dict()}")
        print(f"Dataset distribution: {metadata['Dataset'].value_counts().to_dict()}")
        
        # Prepare plot data
        plot_df = prepare_plot_data(gene_data, metadata, genes_of_interest)

        # Create only the combined violin/boxplot PNG per user request
        create_combined_boxplots_with_jitter(plot_df, output_dir)

        print(f"\n=== Analysis Complete ===")
        print(f"Created: combined_boxplots_with_jitter_specific.png in {output_dir}")

    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main())