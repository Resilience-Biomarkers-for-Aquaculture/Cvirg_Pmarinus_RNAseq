""" This is a one-off script to add the treatment column from SraRunTable (4).csv
which was missing because the treatment heading was misspelled
"""
import os
import pandas as pd

output_dir = "/home/syost/git/Cvirg_Pmarinus_RNAseq/data/"  # Update this path as needed
# Load the metadata CSV file
first_metadata_file = os.path.join(output_dir, "updated_metadata.csv")
updated_metadata_path = os.path.join(output_dir, "updated_metadata_2.csv")

# Load the first metadata file
df1 = pd.read_csv(first_metadata_file)

# Load the second metadata file
second_metadata_file = os.path.join(output_dir, "SraRunTable (4).csv")
df2 = pd.read_csv(second_metadata_file)

# Store the original treatment column to compare later
original_treatment = df1['treatment'].copy()

# Merge the treatment column based on matching Experiment in both files
df1 = df1.merge(df2[['Experiment', 'treatment']], on='Experiment', how='left', suffixes=('', '_new'))

# Replace the old treatment column with the new one where matches are found
df1['treatment'] = df1['treatment_new'].combine_first(df1['treatment'])

# Drop the auxiliary column
df1.drop(columns=['treatment_new'], inplace=True)

# Count how many rows had their treatment value changed
changed_rows = (df1['treatment'] != original_treatment) & df1['treatment'].notna()
num_changed = changed_rows.sum()

# Save the updated DataFrame
df1.to_csv(os.path.join(output_dir, "updated_first_metadata.csv"), index=False)

# Print the report
print(f"Updated first metadata file saved as 'updated_first_metadata.csv'.")
print(f"Number of rows where the treatment value was changed: {num_changed}")
