import pandas as pd

# Load the metadata CSV file
metadata_file_path = "/home/syost/git/Cvirg_Pmarinus_RNAseq/data/merged_metadata.csv"  # Change this if needed
metadata_df = pd.read_csv(metadata_file_path)

# Manually defined second table (replace this with a CSV read if available)
trait_data = pd.DataFrame([
    ["P&S 2020", "PRJNA590205", "242", "susceptible"],
    ["P&S 2020", "PRJNA590205", "266", "resistant"],
    ["P&S 2021", "PRJNA691949", "286", "unknown"],
    ["P&S 2023", "PRJNA894694", "84", "tolerant"],
    ["P&S 2023", "PRJNA894694", "89", "tolerant"],
    ["P&S 2023", "PRJNA894694", "90", "sensitive"],
    ["P&S 2023", "PRJNA894694", "120", "sensitive"],
    ["Chan et al. 2021", "PRJNA778545", "N/A (C. virginica)", "susceptible"],
    ["Chan et al. 2021", "PRJNA778545", "N/A (C. gigas)", "resistant"],
    ["Johnson et al. 2020", "PRJNA604121", "A", "resistant"],
    ["Johnson et al. 2020", "PRJNA604121", "B", "resistant"],
    ["Johnson et al. 2020", "PRJNA604121", "C", "susceptible"],
    ["Johnson et al. 2020", "PRJNA604121", "D", "susceptible"]
], columns=["Study", "BioProject", "Family", "Trait"])

# Ensure Family is treated as a string
trait_data["Family"] = trait_data["Family"].astype(str)

# Determine the number of digits in Family and extract corresponding final digits from BREED
trait_data["Family_digits"] = trait_data["Family"].str.len()

def extract_last_n_digits(breed, max_digits):
    """Extracts the last n digits from BREED based on Family length."""
    breed_str = str(breed)
    for digits in sorted(trait_data["Family_digits"].unique(), reverse=True):
        if len(breed_str) >= digits:
            return breed_str[-digits:]
    return None

# Apply extraction function
metadata_df["BREED_extracted"] = metadata_df["BREED"].apply(lambda x: extract_last_n_digits(x, trait_data["Family_digits"].max()))

# Perform the first merge: Match BioProject and dynamically extracted digits from BREED
merged_df = metadata_df.merge(trait_data[["BioProject", "Family", "Trait", "Family_digits"]], 
                              how='left', left_on=["BioProject", "BREED_extracted"],
                              right_on=["BioProject", "Family"], suffixes=('', '_trait1'))

# Perform the second merge: Match BioProject and Infection level
merged_df = merged_df.merge(trait_data[["BioProject", "Family", "Trait"]], 
                            how='left', left_on=["BioProject", "infection"],
                            right_on=["BioProject", "Family"], suffixes=('', '_trait2'))

# Assign TRSS values from the second table only
merged_df["TRSS"] = merged_df["Trait_trait1"].combine_first(merged_df["Trait_trait2"])

# For rows with BioProject 'PRJNA778545', set TRSS to 'TBD'
merged_df.loc[merged_df["BioProject"] == "PRJNA778545", "TRSS"] = "TBD"

# Drop unnecessary columns while keeping the original Trait column
merged_df.drop(columns=["BREED_extracted", "Family", "Family_trait2", "Trait_trait1", "Trait_trait2", "Family_digits"], inplace=True)

# Save updated metadata file
updated_metadata_path = "updated_metadata_final.csv"
merged_df.to_csv(updated_metadata_path, index=False)

print(f"Processing complete! Updated metadata saved to: {updated_metadata_path}")
