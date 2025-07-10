def filter_gmt_by_size(gmt_in, gmt_out, min_size=10, max_size=500):
    with open(gmt_in) as fin, open(gmt_out, 'w') as fout:
        for line in fin:
            fields = line.strip().split('\t')
            genes = fields[2:]
            if min_size <= len(genes) <= max_size:
                fout.write(line)

if __name__ == "__main__":
    filter_gmt_by_size("4706222.C_virginica_named.gmt", "4706222.C_virginica_named_filtered.gmt", 0, 500);
