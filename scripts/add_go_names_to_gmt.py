from goatools.obo_parser import GODag
import csv

def add_go_names_to_gmt(input_gmt, obo_file, output_gmt):
    go_dag = GODag(obo_file)

    with open(input_gmt, 'r') as fin, open(output_gmt, 'w', newline='') as fout:
        reader = csv.reader(fin, delimiter='\t')
        writer = csv.writer(fout, delimiter='\t')

        for row in reader:
            go_id = row[0]
            term_name = go_dag[go_id].name if go_id in go_dag else row[1]
            writer.writerow([go_id, term_name] + row[2:])

# Example usage
add_go_names_to_gmt("4706222.C_virginica.gmt", "go-basic.obo", "4706222.C_virginica_named.gmt")
