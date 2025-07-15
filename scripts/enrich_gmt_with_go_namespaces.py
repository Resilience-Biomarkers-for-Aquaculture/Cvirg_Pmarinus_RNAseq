from goatools.obo_parser import GODag
import csv

def enrich_gmt_with_go_namespaces(gmt_in, obo_file, gmt_out):
    # Load GO DAG
    go_dag = GODag(obo_file)

    with open(gmt_in, 'r') as fin, open(gmt_out, 'w', newline='') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')

        for row in reader:
            if not row or not row[0].startswith("GO:"):
                continue
            go_id = row[0]
            if go_id in go_dag:
                term_name = go_dag[go_id].name
                namespace = go_dag[go_id].namespace  # biological_process, molecular_function, cellular_component
                label = f"{term_name} [{namespace}]"
            else:
                label = row[1]  # fallback to existing label if GO ID missing
            writer.writerow([go_id, label] + row[2:])

    print(f"Updated GMT written to {gmt_out}")

# Example usage
if __name__ == "__main__":
    gmt_in = "../data/4706222.C_virginica_named_filtered_depth3.gmt"
    obo_file = "../data/go-basic.obo"
    gmt_out = "../data/4706222.C_virginica_named_filtered_depth3_enrichednames.gmt"
    enrich_gmt_with_go_namespaces(gmt_in, obo_file, gmt_out)
