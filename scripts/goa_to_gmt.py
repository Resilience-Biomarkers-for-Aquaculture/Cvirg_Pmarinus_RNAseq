from collections import defaultdict
import csv

def parse_goa_to_gmt(goa_file, gmt_file, use_symbol=True, min_genes=10, exclude_codes={"ND", "IEA"}):
    go_terms = defaultdict(set)

    with open(goa_file, 'r') as f:
        for line in f:
            if line.startswith("!"):
                continue
            cols = line.strip().split('\t')
            if len(cols) < 7:
                continue

            gene_id = cols[2] if use_symbol else cols[1]  # symbol or UniProt
            go_id = cols[4]
            evidence = cols[6]

            if evidence in exclude_codes:
                continue

            go_terms[go_id].add(gene_id)

    with open(gmt_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        for go_id, genes in go_terms.items():
            if len(genes) >= min_genes:
                writer.writerow([go_id, go_id] + sorted(genes))

# Example usage
parse_goa_to_gmt("../data/4706222.C_virginica.goa", "../data/4706222.C_virginica.gmt", use_symbol=True)
