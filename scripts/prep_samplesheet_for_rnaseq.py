import csv

input_file = 'samplesheet.csv'
output_file = 'samplesheet_rnaseq.csv'

with open(input_file, newline='') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.DictReader(infile)
    fieldnames = ['sample', 'fastq_1', 'fastq_2', 'strandedness']
    writer = csv.DictWriter(outfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in reader:
        writer.writerow({
            'sample': row['sample'],
            'fastq_1': row['fastq_1'],
            'fastq_2': row['fastq_2'],
            'strandedness': 'auto'
        })
