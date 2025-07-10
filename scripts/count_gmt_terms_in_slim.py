from goatools.obo_parser import GODag
import csv

def count_gmt_go_slim_matches(gmt_file, go_basic_obo, go_slim_obo):
    # Load the full and slim ontologies
    go_full = GODag(go_basic_obo)
    go_slim = GODag(go_slim_obo)

    # Set of GO slim terms (explicitly defined)
    slim_ids = {term.id for term in go_slim.values()}

    # Read GO terms from your GMT file (first column of each row)
    gmt_go_ids = set()
    with open(gmt_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if row:  # non-empty
                gmt_go_ids.add(row[0])

    # Count direct matches
    direct_matches = gmt_go_ids & slim_ids

    # Optionally: also count indirect matches (e.g., descendants of slim terms)
    # For that, you'd need to check if any GMT GO term is a child of a slim term
    # That logic is more involved and not needed for this direct test

    # Report
    print(f"Total GO terms in GMT: {len(gmt_go_ids)}")
    print(f"GO terms matching GO slim (direct): {len(direct_matches)}")
    print(f"Percentage: {100 * len(direct_matches) / len(gmt_go_ids):.2f}%")

    # Optional: list matching terms
    print("\nExample matches:")
    for go_id in list(direct_matches)[:10]:
        term_name = go_full[go_id].name if go_id in go_full else "(unknown)"
        print(f"{go_id}\t{term_name}")

# Example usage
if __name__ == "__main__":
    gmt_file = "4706222.C_virginica_named_filtered.gmt"
    go_basic_obo = "go-basic.obo"
    go_slim_obo = "goslim_generic.obo"
    count_gmt_go_slim_matches(gmt_file, go_basic_obo, go_slim_obo)
