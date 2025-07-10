from goatools.obo_parser import GODag
import csv

def filter_gmt_by_go_depth(gmt_in, obo_file, gmt_out, min_depth=3):
    # Load GO ontology
    go_dag = GODag(obo_file)

    total, kept = 0, 0

    with open(gmt_in, 'r') as fin, open(gmt_out, 'w', newline='') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')

        for row in reader:
            if not row:
                continue
            go_id = row[0]
            total += 1

            # Check GO term depth
            if go_id in go_dag and go_dag[go_id].depth >= min_depth:
                writer.writerow(row)
                kept += 1

    print(f"Total GO terms in input: {total}")
    print(f"GO terms with depth ≥ {min_depth}: {kept}")
    print(f"Filtered GMT written to: {gmt_out}")

# Example usage
if __name__ == "__main__":
    gmt_in = "4706222.C_virginica_named_filtered.gmt"
    obo_file = "go-basic.obo"
    gmt_out = "4706222.C_virginica_named_filtered_depth3.gmt"
    filter_gmt_by_go_depth(gmt_in, obo_file, gmt_out, min_depth=3)
