#!/usr/bin/env python3
# OBSOLETE:
# This is a precursor to lasso_prune_onefold.py
# that didn't do stability selection and pruning.
# train_eval.py
# Train/evaluate classifier on panel for a single LOSO fold

import argparse, json
from pathlib import Path
import pandas as pd
import joblib
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss, log_loss
import matplotlib.pyplot as plt
from sklearn.calibration import calibration_curve

# ---------------------------
# CLI args
# ---------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--counts", required=True, help="Counts matrix (genes x samples)")
parser.add_argument("--meta", required=True, help="Metadata CSV")
parser.add_argument("--panel", required=True, help="Panel candidate gene list")
parser.add_argument("--outdir", required=True, help="Output directory for fold results")
parser.add_argument("--seed", type=int, default=12345, help="Random seed")
parser.add_argument("--holdout_batch", required=True, help="Which batch is test in this fold")
args = parser.parse_args()

# ---------------------------
# Load data
counts = pd.read_csv(args.counts, sep="\t", index_col=0).T  # samples x genes
meta   = pd.read_csv(args.meta)
meta.columns = [c.lower() for c in meta.columns]
meta = meta.set_index("sample")

panel_genes = pd.read_csv(args.panel, header=None)[0].tolist()
outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

# Align
common = counts.index.intersection(meta.index)
counts = counts.loc[common]
meta   = meta.loc[common]

X = counts[panel_genes]
y = meta["condition"].map({"tolerant":1, "sensitive":0}).values
batches = meta["batch"]

print("[DEBUG] Unique batches in meta:", meta["batch"].unique())
print("[DEBUG] Holdout batch name:", args.holdout_batch)
print("[DEBUG] n_train:", (meta["batch"] != args.holdout_batch).sum(),
      "n_test:", (meta["batch"] == args.holdout_batch).sum())

# ---------------------------
# Train/test split
train_mask = batches != args.holdout_batch
test_mask  = batches == args.holdout_batch

X_train, y_train = X[train_mask], y[train_mask]
X_test,  y_test  = X[test_mask], y[test_mask]

# ---------------------------
# Preprocess + model
scaler = StandardScaler().fit(X_train)
Xz_train = scaler.transform(X_train)
Xz_test  = scaler.transform(X_test)

clf = LogisticRegression(
    penalty=None, solver="lbfgs", max_iter=5000,
    class_weight="balanced", random_state=args.seed
).fit(Xz_train, y_train)

# Save artifacts
joblib.dump(scaler, outdir / "scaler.joblib")
joblib.dump(clf,    outdir / "model.joblib")

# Save coefficients
coef = pd.Series(clf.coef_.ravel(), index=panel_genes)
# This DataFrame has:
# One row per gene (feature)
# One extra row labeled "INTERCEPT" with the model's bias term
coef_df = pd.DataFrame({
    "feature": list(coef.index) + ["INTERCEPT"],
    "coefficient": list(coef.values) + [float(clf.intercept_[0])]
})
coef_df.to_csv(outdir / "model_coefficients.csv", index=False)

# ---------------------------
# Evaluate
prob_test = clf.predict_proba(Xz_test)[:,1]
metrics = {
    "n_test": int(len(y_test)),
    "auroc": float(roc_auc_score(y_test, prob_test)),
    "aupr":  float(average_precision_score(y_test, prob_test)),
    "brier": float(brier_score_loss(y_test, prob_test)),
    "log_loss": float(log_loss(y_test, prob_test))
}
with open(outdir / "test_metrics_overall.json", "w") as f:
    json.dump(metrics, f, indent=2)

# Calibration curve plot
prob_true, prob_pred = calibration_curve(y_test, prob_test, n_bins=5)
plt.plot(prob_pred, prob_true, marker="o"); plt.plot([0,1],[0,1],"--")
plt.xlabel("Predicted"); plt.ylabel("Observed")
plt.title(f"Calibration ({args.holdout_batch} test)")
plt.savefig(outdir / "calibration_plot.png")
plt.close()
