import pandas as pd
import numpy as np

coef = pd.read_csv("results_08_Dec_2025/loso_coefficients_summary.csv")

# CV only makes sense if we have variance across folds
coef["cv"] = coef["sd"] / (coef["mean"].abs() + 1e-9)

# Eligibility: appeared in at least 2 folds and has finite sd
eligible = (coef["n_folds"] >= 2) & np.isfinite(coef["sd"])
cv_elig = coef.loc[eligible, "cv"].replace([np.inf, -np.inf], np.nan)

n_total = len(coef)
n_elig  = cv_elig.notna().sum()
n_missing = n_total - n_elig

print(f"CV eligibility: {n_elig} / {n_total} genes (appeared in ≥2 folds with finite SD);"
      f" {n_missing} lack cross-fold SD (likely selected in only 1 fold).")

if n_elig > 0:
    med_cv = float(np.nanmedian(cv_elig))
    pct_le_0_5 = float((cv_elig <= 0.5).mean() * 100.0)
    pct_gt_1   = float((cv_elig > 1.0).mean() * 100.0)
    print(f"Among CV-eligible genes: median CV = {med_cv:.2f}, "
          f"%CV ≤ 0.5 = {pct_le_0_5:.1f}%, %CV > 1 = {pct_gt_1:.1f}%")
else:
    print("No genes are CV-eligible; cannot assess magnitude stability.")

# Directional stability can still be assessed on all genes:
coef["sign_stable"] = ((coef["sign_pos_frac"] >= 2/3) | (coef["sign_neg_frac"] >= 2/3))
print(f"Genes sign-stable (≥2/3 same direction): {coef['sign_stable'].sum()} / {n_total}")
