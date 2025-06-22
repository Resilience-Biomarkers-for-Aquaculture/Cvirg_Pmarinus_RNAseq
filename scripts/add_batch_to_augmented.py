import pandas as pd
import numpy as np

# File paths
augmented_metadata_path = "../data/augmented_metadata_with_batch.csv"
sra_run_table_path = "../data/SraRunTable (4).csv"  # adjust if needed
output_path = "../data/augmented_metadata_with_batch2.csv"

# Load data
augmented_metadata = pd.read_csv(augmented_metadata_path)
sra_run_table = pd.read_csv(sra_run_table_path)

# Extract relevant columns and rename to avoid clash
sra_batch = sra_run_table[['Experiment', 'batch']].rename(columns={'batch': 'batch_from_sra'})

# Merge SRA batch into metadata
merged = pd.merge(augmented_metadata, sra_batch, on='Experiment', how='left')

# Initialize 'batch' if it doesn't exist
if 'batch' not in merged.columns:
    merged['batch'] = np.nan

# Fill in missing values in 'batch' with 'batch_from_sra'
merged['batch'] = merged['batch'].combine_first(merged['batch_from_sra'])

# Drop the temporary column
merged = merged.drop(columns=['batch_from_sra'])

# Move 'batch' to right of 'Trait' if both exist
if 'Trait' in merged.columns and 'batch' in merged.columns:
    cols = list(merged.columns)
    cols.remove('batch')
    trait_index = cols.index('Trait')
    cols.insert(trait_index + 1, 'batch')
    merged = merged[cols]
else:
    raise KeyError("Required columns 'Trait' or 'batch' not found for column reordering.")

# Save to output CSV with headers
merged.to_csv(output_path, index=False)
print(f"Saved output to {output_path} with headers.")
