import csv

# Input file paths
samplesheet_path = "rnaseq_difabundance_study4_s3_samplesheet.csv"
metadata_path = "augmented_metadata.csv"
output_path = "joined_output.csv"

# Read metadata file into a dictionary keyed by 'Experiment'
metadata_dict = {}
with open(metadata_path, newline='') as metafile:
    reader = csv.DictReader(metafile)
    for row in reader:
        metadata_dict[row["Experiment"]] = row

# Read sample file and join with metadata
with open(samplesheet_path, newline='') as samplefile, open(output_path, "w", newline='') as outfile:
    sample_reader = csv.DictReader(samplefile)
    
    # Determine output fieldnames
    sample_fields = sample_reader.fieldnames
    # meta_fields = [f for f in next(iter(metadata_dict.values())).keys() if f != "Experiment"]
    meta_fields = ['treatment']
    output_fields = sample_fields + meta_fields

    writer = csv.DictWriter(outfile, fieldnames=output_fields)
    writer.writeheader()
    
    for row in sample_reader:
        experiment = row["sample"]
        if experiment in metadata_dict:
            combined = row.copy()
            for key in meta_fields:
                combined[key] = metadata_dict[experiment][key]
            writer.writerow(combined)
