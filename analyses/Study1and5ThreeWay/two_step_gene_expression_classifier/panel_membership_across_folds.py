# Step 1.B (2): Panel membership stability across LOSO folds
# Measures overlap (Jaccard index) of selected genes per fold.

from pathlib import Path
import pandas as pd
import itertools
from scipy.stats import kendalltau

# --- CONFIG ---
root = Path("results_08_Dec_2025")
panel_size = 24  # leave None to auto-detect from filenames

# --- Load per-fold selected gene lists ---
gene_sets = {}
if panel_size is None:
    glob_str = "*_m*"
else:
    glob_str = f"*_m{panel_size}"
for sel_file in root.glob(f"loso_*/selected_genes_{glob_str}.txt"):
    fold = sel_file.parent.name.replace("loso_", "")
    m_str = sel_file.stem.split("_m")[-1]
    if panel_size is None:
        try:
            panel_size = int(m_str)
        except ValueError:
            panel_size = None
    genes = [ln.strip() for ln in sel_file.open() if ln.strip()]
    print(f"using {sel_file} for fold {fold:15s}")
    gene_sets[fold] = set(genes)
    print(f"{fold:15s}: {len(genes)} genes loaded")

folds = sorted(gene_sets)
print(f"\nDetected {len(folds)} folds; panel_size={panel_size}")

# --- Compute pairwise Jaccard overlaps ---
pairs = list(itertools.combinations(folds, 2))
jac_rows = []
for f1, f2 in pairs:
    g1, g2 = gene_sets[f1], gene_sets[f2]
    inter = len(g1 & g2)
    union = len(g1 | g2)
    jaccard = inter / union if union > 0 else float("nan")
    jac_rows.append({"fold1": f1, "fold2": f2, "jaccard": jaccard})
jac = pd.DataFrame(jac_rows).sort_values(["fold1","fold2"]).reset_index(drop=True)

print("\nPairwise Jaccard overlaps:")
print(jac)

mean_jacc = jac["jaccard"].mean()
min_jacc = jac["jaccard"].min()
print(f"\nMean Jaccard: {mean_jacc:.2f},  Min Jaccard: {min_jacc:.2f}")

# --- Optional: rank-correlation stability using final_consensus_genes.csv ---
try:
    cons = pd.read_csv(root / "final_consensus_genes.csv")
    ranks = cons.set_index("gene")["rank"]
    # Compute τ only on genes present in both folds
    rank_rows = []
    for f1, f2 in pairs:
        shared = list(gene_sets[f1] & gene_sets[f2])
        if len(shared) < 3:
            continue
        tau, _ = kendalltau(ranks.loc[shared], ranks.loc[shared])
        rank_rows.append({"fold1": f1, "fold2": f2, "tau": tau, "n_shared": len(shared)})
    if rank_rows:
        rank_df = pd.DataFrame(rank_rows)
        print("\nApproximate rank correlations (Kendall τ on shared genes):")
        print(rank_df)
    else:
        print("\n[Info] Too few shared genes for rank correlation.")
except Exception as e:
    print(f"\n[Info] Rank correlation skipped: {e}")

# --- Simple interpretive thresholds ---
if mean_jacc >= 0.7:
    interp = "High overlap → panel structure likely stable (fine-tuning expected)."
elif mean_jacc >= 0.5:
    interp = "Moderate overlap → partial stability; new data may reorder genes."
else:
    interp = "Low overlap → panel membership unstable; new data likely to alter gene set."
print(f"\nInterpretation: {interp}")

from collections import Counter
all_genes = sum((list(s) for s in gene_sets.values()), [])
freq = Counter(all_genes)
freq_df = pd.DataFrame.from_dict(freq, orient="index", columns=["fold_count"]).sort_values("fold_count", ascending=False)
print(freq_df.head(20))

