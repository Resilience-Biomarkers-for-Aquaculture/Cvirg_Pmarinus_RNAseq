import pandas as pd
import sys
import os

def load_and_annotate(filepath):
    df = pd.read_csv(filepath, sep="\t")
    df["source_file"] = os.path.basename(filepath)
    return df

def main(file1, file2, output):
    df1 = load_and_annotate(file1)
    df2 = load_and_annotate(file2)

    # Find common gene_ids
    common_ids = set(df1["gene_id"]) & set(df2["gene_id"])

    # Filter both DataFrames
    df1_common = df1[df1["gene_id"].isin(common_ids)]
    df2_common = df2[df2["gene_id"].isin(common_ids)]

    # Combine and sort
    combined = pd.concat([df1_common, df2_common], ignore_index=True)
    combined.sort_values(by="gene_id", inplace=True)

    # Write output
    combined.to_csv(output, sep="\t", index=False)
    print(f"Written {len(combined)} rows to {output}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python merge_common_genes.py file1.tsv file2.tsv output.tsv")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
