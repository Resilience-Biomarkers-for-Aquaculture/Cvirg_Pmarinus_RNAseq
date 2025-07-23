import sys
import csv
from collections import defaultdict

def merge_tsv_files(file_list, output_file):
    gene_info = []
    sample_data = defaultdict(dict)
    all_samples = []

    for file_index, file_path in enumerate(file_list):
        with open(file_path, 'r', newline='') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            sample_ids = header[2:]

            # Check for duplicate sample IDs
            for sid in sample_ids:
                if sid in all_samples:
                    raise ValueError(f"Duplicate sample ID '{sid}' found in {file_path}")
                all_samples.append(sid)

            for row_index, row in enumerate(reader):
                if len(row) < 2:
                    raise ValueError(f"Row {row_index+2} in {file_path} has fewer than 2 columns")

                gene_id, gene_name = row[:2]
                values = row[2:]

                # On first file, store gene list
                if file_index == 0:
                    gene_info.append((gene_id, gene_name))
                else:
                    # For other files, validate that gene list matches
                    if row_index >= len(gene_info):
                        raise ValueError(f"Extra gene row in {file_path} at line {row_index+2}")
                    ref_gene_id, ref_gene_name = gene_info[row_index]
                    if (gene_id, gene_name) != (ref_gene_id, ref_gene_name):
                        raise ValueError(
                            f"Gene mismatch at row {row_index+2} in {file_path}: "
                            f"expected ({ref_gene_id}, {ref_gene_name}), got ({gene_id}, {gene_name})"
                        )

                for sid, val in zip(sample_ids, values):
                    sample_data[row_index][sid] = val

    # Check that all files had the same number of gene rows
    for file_path in file_list:
        with open(file_path, 'r') as f:
            num_data_rows = sum(1 for _ in f) - 1  # subtract header
            if num_data_rows != len(gene_info):
                raise ValueError(f"{file_path} has {num_data_rows} gene rows; expected {len(gene_info)}")

    # Write merged output
    with open(output_file, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['gene_id', 'gene_name'] + all_samples)

        for i, (gene_id, gene_name) in enumerate(gene_info):
            row = [gene_id, gene_name]
            for sid in all_samples:
                row.append(sample_data[i].get(sid, 'NA'))
            writer.writerow(row)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python merge_gene_tsv.py <output_file> <input1.tsv> <input2.tsv> ...")
        sys.exit(1)

    output_file = sys.argv[1]
    input_files = sys.argv[2:]
    merge_tsv_files(input_files, output_file)
