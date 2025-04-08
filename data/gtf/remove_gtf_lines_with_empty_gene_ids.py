import gzip

input_path = "GCF_002022765.2_C_virginica-3.0_genomic.gtf.gz"
output_path = "cleaned_GCF_002022765.2_C_virginica-3.0_genomic.gtf.gz"
removed_lines = []

with gzip.open(input_path, 'rt') as infile, gzip.open(output_path, 'wt') as outfile:
    for line_number, line in enumerate(infile, 1):
        if line.startswith("#"):
            outfile.write(line)
            continue

        columns = line.strip().split('\t')
        if len(columns) == 9:
            attributes = columns[8]
            if 'gene_id' not in attributes or 'gene_id ""' in attributes:
                removed_lines.append(line_number)
                continue

        outfile.write(line)

# Print the line numbers of removed lines
print("Removed lines with empty or missing gene_id:")
for ln in removed_lines:
    print(ln)

# Optional: save removed_lines to a text file if needed
with open("removed_line_numbers.txt", "w") as f:
    for ln in removed_lines:
        f.write(f"{ln}\n")
