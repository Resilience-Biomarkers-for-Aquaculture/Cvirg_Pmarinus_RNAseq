import pandas as pd, json, glob
records = []
for f in glob.glob("results/loso_*/test_metrics_overall.json"):
    d = json.load(open(f))
    d["fold"] = f.split("/")[-2].replace("loso_", "")
    records.append(d)
df = pd.DataFrame(records)
print(df)
weighted = (df["auroc"] * df["n_test"]).sum() / df["n_test"].sum()
print("Weighted mean AUROC:", round(weighted, 3))
