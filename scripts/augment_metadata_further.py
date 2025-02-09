""" Note that before running this, I took the output of augment_metadata.py and did the following:
ran temp_merge_treatment.py
Added the species as the Family (it was blank before) for PRJNA778545.
Changed 'Injection' values to 'Injected' for uniformity.
In the former 'collection_time' column, I fixed the values.
They were 0.0, 3.0, and 24.0.
I corrected the values to 0, 30, and 1 (i.e. 24 hours) respectively, according to the paper.
Renamed the column to 'Collection_Interval_Days'.
I hope that captures all the ad-hoc changes I made.
"""
import os
import pandas as pd

output_dir = "/home/syost/git/Cvirg_Pmarinus_RNAseq/data/"  # Update this path as needed
# Load the metadata CSV file
updated_metadata_path = os.path.join(output_dir, "augmented_metadata.csv")
further_updated_metadata_path = os.path.join(output_dir, "further_augmented_metadata.csv")
merged_df = pd.read_csv(updated_metadata_path)

# Convert Collection_Date to datetime format for calculations
merged_df['Collection_Date'] = pd.to_datetime(merged_df['Collection_Date'], errors='coerce')

# Fill in Collection_Interval_Days where missing
for bioproject in merged_df['BioProject'].unique():
    # Get earliest date for the BioProject
    min_date = merged_df.loc[merged_df['BioProject'] == bioproject, 'Collection_Date'].min()
    # Calculate difference in days and fill missing values
    merged_df.loc[(merged_df['BioProject'] == bioproject) & (merged_df['Collection_Interval_Days'].isna()), 'Collection_Interval_Days'] = \
        (merged_df['Collection_Date'] - min_date).dt.days

# Update Collection_Interval_Days for rows where BioProject is PRJNA894694
merged_df.loc[merged_df['BioProject'] == 'PRJNA894694', 'Collection_Interval_Days'] = 7.0

# Calculate the number of days in 14 months (assuming an average month length of 30.44 days)
days_in_14_months = round(14 * 30.44)

# Update Collection_Interval_Days for rows where BioProject is PRJNA604121
merged_df.loc[merged_df['BioProject'] == 'PRJNA604121', 'Collection_Interval_Days'] = days_in_14_months


# Define BioProject to Study mapping
study_mapping = {
    "PRJNA590205": "P&S 2020",
    "PRJNA691949": "P&S 2021",
    "PRJNA894694": "P&S 2023",
    "PRJNA778545": "Chan et al. 2021",
    "PRJNA604121": "Johnson et al. 2020"
}

# Create the Study column using mapping
merged_df['Study'] = merged_df['BioProject'].map(study_mapping)

# Reorder columns: Move Study to the beginning and BioProject next to it
column_order = ['Study', 'BioProject'] + [col for col in merged_df.columns if col not in ['Study', 'BioProject']]
merged_df = merged_df[column_order]

# Save updated metadata file
merged_df.to_csv(further_updated_metadata_path, index=False)

print(f"Processing complete! Updated metadata saved to: {further_updated_metadata_path}")
