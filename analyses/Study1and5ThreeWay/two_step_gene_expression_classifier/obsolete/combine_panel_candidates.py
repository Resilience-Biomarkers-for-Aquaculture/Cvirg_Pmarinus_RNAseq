# This is from
# https://chatgpt.com/share/68bc4d13-e418-800d-9960-f42e37d9f98b
import glob, pandas as pd, pathlib

paths = sorted(glob.glob("results/loso_*/panel_candidates.txt"))
print("Found candidate lists:", paths)

gene_counts = {}
for path in paths:
    with open(path) as f:
        for g in map(str.strip, f):
            if g:
                gene_counts[g] = gene_counts.get(g, 0) + 1

# Make a table
df = pd.DataFrame.from_dict(gene_counts, orient="index", columns=["n_folds_found"])
df = df.sort_values("n_folds_found", ascending=False)
df.to_csv("panel_candidates_combined.csv")
print(df.head())
print("Wrote panel_candidates_combined.csv")
keep = df[df["n_folds_found"] >= 2].index
pd.Series(keep).to_csv("panel_candidates_tier12_combined.txt",
                       index=False, header=False)
print(f"Wrote {len(keep)} combined candidates")
print("Wrote panel_candidates_tier12_combined.txt")