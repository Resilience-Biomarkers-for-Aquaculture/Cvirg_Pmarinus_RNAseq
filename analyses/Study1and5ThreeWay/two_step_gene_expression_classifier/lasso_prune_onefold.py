#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lasso_prune_onefold.py
Performs the same per-fold workflow as lasso_prune_loso_size_nested.py,
but for a single holdout batch, using the same CLI arguments as train_eval.py.

Steps:
  1. Load expression matrix, metadata, and candidate gene list
  2. Split train/test by holdout batch
  3. Perform stability selection on training data
  4. Redundancy prune the top-ranked genes
  5. Evaluate performance for various panel sizes on that fold
  6. Save per-size metrics, selected gene lists, and a summary JSON
"""

import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.calibration import CalibratedClassifierCV

# ---------------------------
# CLI args (same as train_eval.py)
# ---------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--counts", required=True, help="Counts matrix (genes x samples)")
parser.add_argument("--meta", required=True, help="Metadata CSV")
parser.add_argument("--panel", required=True, help="Candidate gene list (Tier1+Tier2 etc.)")
parser.add_argument("--outdir", required=True, help="Output directory for fold results")
parser.add_argument("--seed", type=int, default=12345)
parser.add_argument("--holdout_batch", required=True, help="Which batch is test in this fold")
args = parser.parse_args()

outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)

# ---------------------------
# Load data
# ---------------------------
counts = pd.read_csv(args.counts, sep="\t", index_col=0).T
meta = pd.read_csv(args.meta)
meta.columns = [c.lower() for c in meta.columns]
meta = meta.set_index("sample")

panel_genes = pd.read_csv(args.panel, header=None)[0].astype(str).tolist()

common = counts.index.intersection(meta.index)
counts = counts.loc[common]
meta = meta.loc[common]

# Filter candidate genes that exist in count matrix
panel_genes = [g for g in panel_genes if g in counts.columns]
if not panel_genes:
    raise ValueError("No candidate genes found in counts matrix.")

X_df = counts[panel_genes]
y = meta["condition"].map({"tolerant": 1, "sensitive": 0}).values
batch = meta["batch"].astype(str)

if False:
    # -------------------------------------------------------
    # Optional: mean-center each batch to remove additive batch offsets
    # -------------------------------------------------------
    # Diagnostic: batch-wise means before centering
    # Note this would only need to be done once, so really shouldn't be in
    # per-fold script, but leaving here for convenience.
    # E.g. could just run standalone:
    # python lasso_prune_onefold.py \
    # --counts=../../../data/rnaseq_gene_counts/merged_gene_counts.tsv \
    # --meta=../../../data/differential_abundance_sheets/rnaseq_diffabundance_study1and5_samplesheet_filled_with_study5_alldays.csv \
    # --panel=results/loso_P&S_2023/panel_candidates.txt \
    # --outdir=results/loso_P&S_2023 --seed=12346 --holdout_batch=P&S 2023
    orig_means = X_df.groupby(batch).mean()

    # Apply per-batch centering
    X_batch_centered = (
        X_df.groupby(batch, group_keys=False)
            .apply(lambda g: g - g.mean())
    )

    # Diagnostic: verify centering worked
    centered_means = X_batch_centered.groupby(batch).mean()
    diff = centered_means.abs().mean().mean()
    print(f"[DEBUG] Avg abs post-centering batch mean: {diff:.4e}")
    # Got:  Avg abs post-centering batch mean: 1.8866e-14

    # Save diagnostics to do a PCA plot later if desired
    orig_means.to_csv(outdir / f"batch_means_before_{args.holdout_batch}.csv")
    centered_means.to_csv(outdir / f"batch_means_after_{args.holdout_batch}.csv")
 
    # Optionally replace X_df with centered version for modeling
    X_df = X_batch_centered.copy()
    # -------------------------------------------------------

train_mask = meta["batch"] != args.holdout_batch
test_mask = ~train_mask
X_tr, X_te = X_df.loc[train_mask], X_df.loc[test_mask]
y_tr, y_te = y[train_mask], y[test_mask]
batch_tr, batch_te = batch.loc[train_mask], batch.loc[test_mask]

print(f"[INFO] Holdout batch: {args.holdout_batch}")
print(f"[INFO] Train samples: {len(X_tr)}, Test samples: {len(X_te)}")

# ---------------------------
# Helper functions
# ---------------------------
def composite_strata(y, batch):
    return np.array([f"{b}__{c}" for b, c in zip(batch, y)])

def fit_l1_logreg(X, y, C=0.1, l1_ratio=None, seed=0):
    penalty = "l1" if l1_ratio is None else "elasticnet"
    clf = LogisticRegression(
        penalty=penalty, l1_ratio=l1_ratio, C=C,
        solver="saga", max_iter=5000, n_jobs=-1,
        class_weight="balanced", random_state=seed
    )
    clf.fit(X, y)
    return clf

def stability_selection(X_df, y, batch, n_resamples=300, subsample_frac=0.7,
                        C_grid=(0.05, 0.1, 0.2), l1_ratio=None, random_state=123):
    rng = np.random.RandomState(random_state)
    genes = X_df.columns.to_numpy()
    sel_counts = pd.Series(0.0, index=genes)
    coef_sums = pd.Series(0.0, index=genes)
    strata = composite_strata(y, batch)
    sss = StratifiedShuffleSplit(n_splits=n_resamples, train_size=subsample_frac,
                                 random_state=random_state)
    X = X_df.to_numpy()
    for i, (tr, _) in enumerate(sss.split(np.zeros_like(strata), strata), 1):
        scaler = StandardScaler().fit(X[tr])
        Xt = scaler.transform(X[tr])
        best_coef, best_nz = None, 10**9
        """ This method selects the C that gives the sparsest model per resample
        It produced too few selected features in practice.
        for C in C_grid:
            coef = fit_l1_logreg(Xt, y[tr], C=C, l1_ratio=l1_ratio,
                                 seed=rng.randint(1_000_000_000)).coef_.ravel()
            nz = np.flatnonzero(coef != 0)
            if len(nz) < best_nz:
                best_nz = len(nz)
                best_coef = coef
        """
        # Choose middle C instead (no searching for sparsest)
        C = C_grid[len(C_grid)//2]
        model = fit_l1_logreg(Xt, y[tr], C=C, l1_ratio=l1_ratio,
                              seed=rng.randint(1_000_000_000))
        best_coef = model.coef_.ravel()

        nz = np.flatnonzero(best_coef != 0)
        if len(nz):
            sel_counts.iloc[nz] += 1.0
            coef_sums.iloc[nz] += np.abs(best_coef[nz])
        if i % 50 == 0:
            print(f"[DEBUG] Stability resample {i}/{n_resamples}")
    stab = pd.DataFrame({
        "gene": genes,
        "selection_prob": (sel_counts / n_resamples).values,
        "mean_abs_coef_when_selected": (coef_sums / np.maximum(sel_counts, 1.0)).values
    }).sort_values(["selection_prob", "mean_abs_coef_when_selected"], ascending=False)
    return stab

def redundancy_prune(X_df, genes, corr_threshold=0.9):
    kept = []
    corr = X_df[genes].corr().abs()
    for g in genes:
        if kept:
            if (corr.loc[g, kept] > corr_threshold).any():
                continue
        kept.append(g)
    return kept

def loso_eval(X_tr, y_tr, X_te, y_te, calibrate=True):
    scaler = StandardScaler().fit(X_tr)
    X_tr = scaler.transform(X_tr)
    X_te = scaler.transform(X_te)
    base = LogisticRegression(penalty=None, max_iter=5000, n_jobs=-1)
    if calibrate:
        n_cv = min(5, max(2, int(min(np.bincount(y_tr)))))
        clf = CalibratedClassifierCV(base, method="sigmoid", cv=n_cv)
    else:
        clf = base
    clf.fit(X_tr, y_tr)
    p = clf.predict_proba(X_te)[:, 1]
    return {
        "n_test": len(y_te),
        "auroc": roc_auc_score(y_te, p),
        "aupr": average_precision_score(y_te, p),
        "brier": brier_score_loss(y_te, p)
    }

# ---------------------------
# Nested feature selection on training fold
# ---------------------------
print(f"[INFO] Running nested feature selection for held-out batch: {args.holdout_batch}")
stab = stability_selection(X_tr, y_tr, batch_tr, random_state=args.seed)
stab.to_csv(outdir / f"stability_selection_{args.holdout_batch}.csv", index=False)

ordered_genes = stab["gene"].tolist()
TARGET_SIZES = [6, 8, 10, 12, 16, 24]

records = []
for m in TARGET_SIZES:
    topm = ordered_genes[:max(m*2, m)]
    pruned = redundancy_prune(X_tr, topm, corr_threshold=0.9)
    # genes dropped = those in topm but not in pruned_keep
    dropped = [g for g in topm if g not in pruned]
  
    if len(pruned) > m:
        pruned = pruned[:m]

    metrics = loso_eval(X_tr[pruned].to_numpy(), y_tr, X_te[pruned].to_numpy(), y_te)
    metrics.update({"heldout": args.holdout_batch, "m": m, "genes_kept": len(pruned)})
    records.append(metrics)

    pd.Series(pruned).to_csv(outdir / f"selected_genes_{args.holdout_batch}_m{m}.txt",
                             index=False, header=False)
    print(f"[DEBUG] Fold {args.holdout_batch}, m={m}, AUROC={metrics['auroc']:.3f}")

# Save all metrics
df_metrics = pd.DataFrame(records)
df_metrics.to_csv(outdir / f"panel_size_metrics_{args.holdout_batch}.csv", index=False)

# ---------------------------------
# Train final model on the best-performing panel size
# ---------------------------------
best = max(records, key=lambda r: r["auroc"])
best_m = best["m"]
best_genes = pd.read_csv(outdir / f"selected_genes_{args.holdout_batch}_m{best_m}.txt",
                         header=None)[0].tolist()
print(f"[INFO] Refitting final model on top {best_m} genes: {best_genes[:5]} ...")

# Standardize train and test
scaler = StandardScaler().fit(X_tr[best_genes])
Xtr_z = scaler.transform(X_tr[best_genes])
Xte_z = scaler.transform(X_te[best_genes])

base = LogisticRegression(penalty=None, max_iter=5000, n_jobs=-1)
clf = CalibratedClassifierCV(base, method="sigmoid", cv=min(5, max(2, int(min(np.bincount(y_tr))))))
clf.fit(Xtr_z, y_tr)

# ---------------------------------
# Save coefficients and intercepts
# ---------------------------------
# Average calibrated models if multiple calibrators
# ---------------------------------
# Save coefficients and intercepts (version-safe)
# ---------------------------------
coef_list, intercept_list = [], []

if hasattr(clf, "calibrated_classifiers_"):  # when calibration is used
    for estimator in clf.calibrated_classifiers_:
        base = getattr(estimator, "base_estimator", None) or getattr(estimator, "estimator", None)
        if base is not None and hasattr(base, "coef_"):
            coef_list.append(base.coef_.ravel())
            intercept_list.append(float(base.intercept_[0]))
elif hasattr(clf, "base_estimator_"):  # fallback for uncalibrated model
    base = getattr(clf, "base_estimator_", clf)
    coef_list.append(base.coef_.ravel())
    intercept_list.append(float(base.intercept_[0]))
elif hasattr(clf, "coef_"):
    coef_list.append(clf.coef_.ravel())
    intercept_list.append(float(clf.intercept_[0]))

# Average across calibrated folds if applicable
coefs = np.mean(np.vstack(coef_list), axis=0)
intercept = float(np.mean(intercept_list))

coef_df = pd.DataFrame({
    "gene": best_genes,
    "coef": coefs
})
coef_df.loc[len(coef_df)] = ["Intercept", intercept]
coef_df.to_csv(outdir / f"coefficients_{args.holdout_batch}.csv", index=False)
print(f"[INFO] Saved coefficients_{args.holdout_batch}.csv")


# ---------------------------------
# Save per-sample predictions
# ---------------------------------
probs = clf.predict_proba(Xte_z)[:, 1]
preds_df = pd.DataFrame({
    "sample": X_te.index,
    "true_label": y_te,
    "pred_prob_tolerant": probs
})
preds_df.to_csv(outdir / f"preds_{args.holdout_batch}.csv", index=False)
print(f"[INFO] Saved preds_{args.holdout_batch}.csv")

# ---------------------------------
# Save JSON summary (final)
# ---------------------------------
# Save JSON summary
best = max(records, key=lambda r: r["auroc"])
best.update({
    "best_m": best_m,
    "n_genes": len(best_genes),
    "coef_file": f"coefficients_{args.holdout_batch}.csv",
    "preds_file": f"preds_{args.holdout_batch}.csv"
})
with open(outdir / f"best_panel_summary_{args.holdout_batch}.json", "w") as f:
    json.dump(best, f, indent=2)

with open(outdir / f"best_panel_summary_{args.holdout_batch}.json", "w") as f:
    json.dump(best, f, indent=2)
print(f"[INFO] Final summary saved to best_panel_summary_{args.holdout_batch}.json")

print(f"[INFO] Results saved to {outdir}")
