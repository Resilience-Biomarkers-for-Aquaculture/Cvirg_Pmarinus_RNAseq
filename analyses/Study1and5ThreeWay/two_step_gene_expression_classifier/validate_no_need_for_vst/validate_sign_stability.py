import pandas as pd
from pathlib import Path
coef = pd.read_csv("results/loso_coefficients_summary.csv")
# Inspect only genes in your chosen final panel if you’ve fixed m:
panel = set(pd.read_csv(next(Path("results").glob("loso_*/selected_genes_*_m*.txt")), header=None)[0])
stab = coef[coef["gene"].isin(panel)].copy()
print(stab[["gene","mean","sd","n_folds","sign_pos_frac","sign_neg_frac"]].sort_values("gene"))
