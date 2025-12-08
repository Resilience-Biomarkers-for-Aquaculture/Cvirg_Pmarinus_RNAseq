import pandas as pd, numpy as np
df = pd.read_csv("results/loso_summary_metrics.csv")
print(df[["batch","n_test","auroc","aupr","brier"]].round(3))

w = df["n_test"].astype(float)
print("Weighted AUROC:", np.average(df["auroc"], weights=w))
print("Weighted AUPR :", np.average(df["aupr"],  weights=w))
