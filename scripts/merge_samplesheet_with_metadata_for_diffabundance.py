import csv

# Input file paths
samplesheet_path = "rnaseq_difabundance_fullset_s3_samplesheet.csv"
metadata_path = "../augmented_metadata.csv"
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
    
    # Remove unwanted fields from sample_fields
    sample_fields = [f for f in sample_reader.fieldnames if f not in ("fastq_1", "fastq_2")]
    
    # Only selected metadata fields
    meta_fields = ['treatment', 'batch', 'Collection_Interval_Days', 'TRSS', 'Study']
    
    # Final output field order
    output_fields = sample_fields + meta_fields

    writer = csv.DictWriter(outfile, fieldnames=output_fields)
    writer.writeheader()
    
    for row in sample_reader:
        experiment = row["sample"]
        if experiment in metadata_dict:
            # Copy row but remove the unwanted fields
            combined = {k: v for k, v in row.items() if k in sample_fields}
            
            # Add metadata fields
            for key in meta_fields:
                combined[key] = metadata_dict[experiment][key]
            
            writer.writerow(combined)
