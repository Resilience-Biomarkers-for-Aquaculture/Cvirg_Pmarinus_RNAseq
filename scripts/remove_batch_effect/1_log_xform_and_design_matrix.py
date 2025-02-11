import pandas as pd
import numpy as np
import os

working_dir = "/mnt/c/Users/inter/Downloads/GMGI/Perkinsus batch correction"
# === 1️⃣ Load Expression Data (TSV) ===
expression_file = os.path.join(working_dir, "merged_gene_counts.tsv")
expression_data = pd.read_csv(expression_file, sep="\t", index_col=0)  # gene_id as index

# === 2️⃣ Load Metadata (CSV) ===
metadata_file = os.path.join(working_dir, "augmented_metadata.csv")
metadata = pd.read_csv(metadata_file)

# === 3️⃣ Create design Matrix ===
design = metadata[['control', 'Trait', 'Collection_Interval_Days']].copy()

# One-hot encode categorical variables
design_encoded = pd.get_dummies(design, columns=['Trait'], drop_first=True)

# Ensure design have the same row order as expression_data
design_encoded.index = metadata['Experiment']  # Assuming 'Experiment' matches sample names
design_encoded = design_encoded.astype(int)

# Save the design matrix
design_output_tsv = os.path.join(working_dir, "design_matrix.tsv")
design_encoded.to_csv(design_output_tsv, sep="\t")

print(f"✅ design matrix saved as: {design_output_tsv}")

# === 4️⃣ Apply Log Transformation (Exclude Gene ID) ===
log_expression = np.log1p(expression_data.astype(float))  # Convert numeric values before applying log

# === 6️⃣ Save the Processed Log-Transformed Expression Data (CSV & TSV) ===
log_expression_output_tsv = os.path.join(working_dir, "log_transformed_expression.tsv")

log_expression.to_csv(log_expression_output_tsv, sep="\t")

print(f"✅ Log-transformed expression data saved as: {log_expression_output_tsv}")

print("🎯 Data is now fully prepared for batch correction in R!")
