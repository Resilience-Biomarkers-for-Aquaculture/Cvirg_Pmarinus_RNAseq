import os
import pandas as pd

output_dir = "/home/syost/git/Cvirg_Pmarinus_RNAseq/data/"  # Update this path as needed
# Load the metadata CSV file
metadata_file_path = os.path.join(output_dir, "merged_metadata.csv")
updated_metadata_path = os.path.join(output_dir, "updated_metadata.csv")
metadata_df = pd.read_csv(metadata_file_path)

# Manually defined second table (replace this with a CSV read if available)
trait_data = pd.DataFrame([
    ["P&S 2020", "PRJNA590205", "242", "susceptible"],
    ["P&S 2020", "PRJNA590205", "266", "resistant"],
    ["P&S 2021", "PRJNA691949", "286", "unknown"],
    ["P&S 2023", "PRJNA894694", "084", "tolerant"],
    ["P&S 2023", "PRJNA894694", "089", "tolerant"],
    ["P&S 2023", "PRJNA894694", "090", "sensitive"],
    ["P&S 2023", "PRJNA894694", "120", "sensitive"],
    ["Johnson et al. 2020", "PRJNA604121", "A", "resistant"],
    ["Johnson et al. 2020", "PRJNA604121", "B", "resistant"],
    ["Johnson et al. 2020", "PRJNA604121", "C", "susceptible"],
    ["Johnson et al. 2020", "PRJNA604121", "D", "susceptible"]
], columns=["Study", "BioProject", "Family", "Trait"])

# Ensure Family is treated as a string
trait_data["Family"] = trait_data["Family"].astype(str)

# Create a dictionary mapping BioProject to Family lengths
family_lengths = trait_data.set_index("BioProject")["Family"].astype(str).apply(len).to_dict()

# Extract the corresponding number of digits from BREED dynamically
metadata_df["BREED_suffix"] = metadata_df.apply(
    lambda row: str(row["BREED"])[-family_lengths.get(row["BioProject"], 0):], axis=1
)

# Merge using BioProject and dynamically extracted Family ID from BREED
merged_df = metadata_df.merge(
    trait_data[["BioProject", "Family", "Trait"]],
    how="left",
    left_on=["BioProject", "BREED_suffix"],
    right_on=["BioProject", "Family"],
    suffixes=("", "_trait1")
)

# Perform the second merge: Match BioProject and Infection level
merged_df = merged_df.merge(trait_data[["BioProject", "Family", "Trait"]], 
                            how='left', left_on=["BioProject", "infection"],
                            right_on=["BioProject", "Family"], suffixes=('', '_trait2'))

# Assign TRSS correctly, ensuring we prioritize the second table's Trait column
merged_df["TRSS"] = merged_df["Trait_trait1"].combine_first(merged_df["Trait_trait2"])

# Ensure TRSS values are locked before proceeding further
trss_values_after_fix = merged_df["TRSS"].notna().sum()

# Update the TRS column based on BioProject PRJNA778545 and Trait values
condition = merged_df["BioProject"] == "PRJNA778545"
merged_df.loc[condition, "TRSS"] = merged_df.loc[condition, "Trait"].map({"tolerant": "resistant", "sensitive": "susceptible"}).fillna("TBD")
merged_df.loc[condition, "BREED"] = merged_df.loc[condition, "Trait"].map({"tolerant": "Crassostrea gigas", "sensitive": "Crassostrea virginica"}).fillna("TBD")

# Append the 'infection' value to 'BREED', prefixed with an underscore, where 'infection' is not NaN
merged_df.loc[merged_df['infection'].notna(), 'BREED'] = merged_df['BREED'] + '_' + merged_df['infection']


# Ensure the TRSS column is not accidentally removed
columns_to_keep = list(metadata_df.columns) + ["TRSS"]
columns_to_keep.remove("BREED_suffix")
columns_to_keep.remove("infection")
merged_df = merged_df[columns_to_keep]

# Save updated metadata file
merged_df.to_csv(updated_metadata_path, index=False)

print(f"Processing complete! Updated metadata saved to: {updated_metadata_path}")
