#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Fixed 6-gene classifier (revised):
- Per-batch stratified holdout by 'condition' with safety constraints for tiny classes
- Train on pooled remainder with per-gene z-scoring (fit on train only)
- Fit a single logistic regression (class_weight='balanced'), freeze scaler + model
- Evaluate on the held-out test set (overall and per-batch)
- Report AUROC, AUPR, Brier, log-loss, calibration slope/intercept
- Bootstrap 95% CIs for metrics (overall and per-batch)

Inputs:
  VST_PATH   : nf-core DESEQ2_NORM VST matrix (genes x samples; first column = gene_id)
  META_PATH  : original metadata CSV (columns: sample, condition, batch, ...)
  PANEL_PATH : final_panel_gene_list.txt (six genes; one per line)

Outputs (written to OUTDIR):
  scaler.joblib, model.joblib
  model_coefficients.csv
  test_metrics_overall.json
  test_metrics_per_batch.csv
  test_predictions.csv
"""

import json
from pathlib import Path
import numpy as np
import pandas as pd

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss, log_loss
from sklearn.utils import check_random_state
import joblib


# -------------------------------
# User-editable paths & knobs
# -------------------------------
OUTDIR     = "final_panel_fixed_model"
VST_PATH = "DESEQ2_NORM_all.vst.tsv"
META_PATH = "../../data/differential_abundance_sheets/rnaseq_diffabundance_study1and5_samplesheet_filled.csv"
PANEL_PATH = "final_panel_gene_list.txt"


# Holdout fraction per batch (target; constraints apply per class within batch)
HOLDOUT_FRAC = 0.25  # suggest 0.20–0.30

# Safety constraints for tiny classes (per batch, per class)
MIN_TEST_PER_CLASS  = 1   # aim to put >=1 in test if feasible
MIN_TRAIN_PER_CLASS = 3   # keep >=3 for training if feasible

# Bootstrap settings for CIs on test metrics
N_BOOT   = 1000
CI_LOW   = 2.5
CI_HIGH  = 97.5

# Random seed for reproducibility
RANDOM_SEED = 12345

# Optional: map verbose batch labels to A/B1/B2 for reporting; leave {} to use raw labels
BATCH_MAP = {
    # "P&S 2023": "A",
    # "P&S 2020 2015": "B1",
    # "P&S 2020 2017": "B2",
}


# -------------------------------
# I/O helpers
# -------------------------------
def load_vst(vst_path: str) -> pd.DataFrame:
    """Read genes x samples VST; return samples x genes (index='sample')."""
    df = pd.read_csv(vst_path, sep="\t")
    df = df.rename(columns={df.columns[0]: "gene_id"}).set_index("gene_id")
    sxg = df.T
    sxg.index.name = "sample"
    return sxg


def load_meta(meta_path: str) -> pd.DataFrame:
    """Read raw metadata; normalize headers/values; index by 'sample'."""
    meta = pd.read_csv(meta_path)
    meta.columns = [c.replace("\ufeff", "").strip() for c in meta.columns]
    # normalize key names
    rename = {}
    for c in meta.columns:
        cl = c.lower()
        if cl == "sample":    rename[c] = "sample"
        if cl == "condition": rename[c] = "condition"
        if cl == "batch":     rename[c] = "batch"
    meta = meta.rename(columns=rename)
    req = {"sample", "condition", "batch"}
    miss = req - set(meta.columns)
    if miss:
        raise ValueError(f"Metadata missing required columns: {miss}")

    meta["sample"]    = meta["sample"].astype(str).str.strip()
    meta["condition"] = meta["condition"].astype(str).str.strip().str.lower()
    meta["batch"]     = meta["batch"].astype(str).str.strip()
    if BATCH_MAP:
        meta["batch"] = meta["batch"].map(lambda x: BATCH_MAP.get(x, x))
    return meta.set_index("sample")


# -------------------------------
# Splitting & metrics
# -------------------------------
def stratified_holdout_by_batch(meta: pd.DataFrame,
                                y: np.ndarray,
                                hold_frac=0.25,
                                seed=123,
                                min_test_per_class=1,
                                min_train_per_class=3):
    """
    For each batch, sample a test set stratified by class, enforcing:
      - at least min_test_per_class in test per class (if feasible)
      - at least min_train_per_class in train per class
    If constraints cannot be met for a class in a batch, that class in that batch contributes 0 to test.
    Returns train_mask, test_mask (bool arrays aligned to meta.index order).
    """
    rng = check_random_state(seed)
    n = len(meta)
    test_flags = np.zeros(n, dtype=bool)
    batches = meta["batch"].values

    for b in np.unique(batches):
        idx_b = np.where(batches == b)[0]
        y_b   = y[idx_b]
        for label in np.unique(y_b):
            idx_bc = idx_b[y_b == label]
            n_bc = len(idx_bc)
            if n_bc == 0:
                continue

            # proposed test size from fraction
            k = int(round(n_bc * hold_frac))
            # if not enough to satisfy both minima → all go to train
            if n_bc < (min_test_per_class + min_train_per_class):
                k = 0
            else:
                k = max(k, min_test_per_class)                   # ensure min test
                k = min(k, n_bc - min_train_per_class)           # leave min train

            if k > 0:
                take = rng.choice(idx_bc, size=k, replace=False)
                test_flags[take] = True

    train_flags = ~test_flags
    # guardrails
    assert train_flags.any() and test_flags.any(), "Empty train or test split after holdout."
    return train_flags, test_flags


def calibration_slope_intercept(y_true, y_prob, eps=1e-12):
    """Logistic calibration: regress y on logit(p); return slope, intercept."""
    from sklearn.linear_model import LogisticRegression
    p = np.clip(y_prob, eps, 1 - eps)
    logit = np.log(p / (1 - p))
    lr = LogisticRegression(penalty=None, solver="lbfgs", max_iter=1000)
    lr.fit(logit.reshape(-1, 1), y_true)
    return float(lr.coef_.ravel()[0]), float(lr.intercept_.ravel()[0])


def overall_metrics(y, p):
    return {
        "auroc": float(roc_auc_score(y, p)),
        "aupr":  float(average_precision_score(y, p)),
        "brier": float(brier_score_loss(y, p)),
        "log_loss": float(log_loss(y, np.clip(p, 1e-9, 1 - 1e-9))),
    }


def bootstrap_cis(y, p, n_boot=1000, seed=123, metrics_fn=overall_metrics, stratify_by=None):
    """
    Bootstrap percentile CIs for metrics. If stratify_by is provided (array of group labels),
    resample within each group to preserve composition (useful for per-batch metrics).
    Returns dict: {'auroc_ci': [low, high], 'aupr_ci': [...], ...}
    """
    rng = check_random_state(seed)
    stats = {"auroc": [], "aupr": [], "brier": [], "log_loss": []}

    idx_all = np.arange(len(y))
    groups = None
    if stratify_by is not None:
        groups = {g: np.where(stratify_by == g)[0] for g in np.unique(stratify_by)}

    for _ in range(n_boot):
        if groups is None:
            boot_idx = rng.choice(idx_all, size=len(idx_all), replace=True)
        else:
            # group-wise resample to keep composition
            boot_idx = np.concatenate([
                rng.choice(ix, size=len(ix), replace=True) for ix in groups.values()
            ])
        yb = y[boot_idx]; pb = p[boot_idx]
        try:
            m = metrics_fn(yb, pb)
            for k in stats.keys():
                stats[k].append(m[k])
        except Exception:
            # in rare degenerate resamples (e.g., only one class), skip
            continue

    def ci(arr):
        return [float(np.percentile(arr, CI_LOW)), float(np.percentile(arr, CI_HIGH))]

    return {f"{k}_ci": ci(v) for k, v in stats.items() if len(v) > 0}


# -------------------------------
# Main
# -------------------------------
def main():
    Path(OUTDIR).mkdir(parents=True, exist_ok=True)

    # Load inputs
    X_all = load_vst(VST_PATH)       # samples x genes
    meta  = load_meta(META_PATH)     # indexed by sample
    genes = pd.read_csv(PANEL_PATH, header=None)[0].astype(str).str.strip().tolist()

    # Align & subset to 6 genes
    common = X_all.index.intersection(meta.index)
    if len(common) == 0:
        raise ValueError("No overlapping sample IDs between VST and metadata.")
    X_all = X_all.loc[common]
    meta  = meta.loc[common]

    missing_genes = [g for g in genes if g not in X_all.columns]
    if missing_genes:
        raise ValueError(f"Panel genes missing in expression matrix: {missing_genes}")

    X = X_all[genes].copy()

    # Encode labels
    cond_map = {"tolerant": 1, "sensitive": 0}
    unknown = set(meta["condition"].unique()) - set(cond_map.keys())
    if unknown:
        raise ValueError(f"Unexpected condition values in metadata: {unknown}")
    y = meta["condition"].map(cond_map).astype(int).values
    batch = meta["batch"].values

    # Per-batch stratified holdout with safety constraints
    train_mask, test_mask = stratified_holdout_by_batch(
        meta, y,
        hold_frac=HOLDOUT_FRAC,
        seed=RANDOM_SEED,
        min_test_per_class=MIN_TEST_PER_CLASS,
        min_train_per_class=MIN_TRAIN_PER_CLASS
    )

    X_train, y_train = X.iloc[train_mask], y[train_mask]
    X_test,  y_test  = X.iloc[test_mask],  y[test_mask]
    batch_test       = batch[test_mask]

    # Fit scaler on TRAIN only; transform train & test
    scaler = StandardScaler().fit(X_train.values)
    Xz_train = scaler.transform(X_train.values)
    Xz_test  = scaler.transform(X_test.values)

    # Fit single logistic regression on TRAIN
    clf = LogisticRegression(
        penalty=None, solver="lbfgs", max_iter=5000,
        class_weight="balanced", n_jobs=-1
    ).fit(Xz_train, y_train)

    # Freeze artifacts
    joblib.dump(scaler, f"{OUTDIR}/scaler.joblib")
    joblib.dump(clf,    f"{OUTDIR}/model.joblib")

    coef = pd.Series(clf.coef_.ravel(), index=genes)
    coef_df = pd.DataFrame({
        "feature": list(coef.index) + ["INTERCEPT"],
        "coefficient": list(coef.values) + [float(clf.intercept_.ravel()[0])]
    })
    coef_df.to_csv(f"{OUTDIR}/model_coefficients.csv", index=False)

    # Inference on TEST
    prob_test = clf.predict_proba(Xz_test)[:, 1]
    pred_test = (prob_test >= 0.5).astype(int)

    # Overall metrics + CIs
    metrics_overall = {
        "n_test": int(len(y_test)),
        **overall_metrics(y_test, prob_test)
    }
    slope, intercept = calibration_slope_intercept(y_test, prob_test)
    metrics_overall.update({
        "calibration_slope": slope,
        "calibration_intercept": intercept
    })
    # Bootstrap CIs (overall)
    cis_overall = bootstrap_cis(
        y_test, prob_test,
        n_boot=N_BOOT, seed=RANDOM_SEED + 1
    )
    metrics_overall.update(cis_overall)

    with open(f"{OUTDIR}/test_metrics_overall.json", "w") as f:
        json.dump(metrics_overall, f, indent=2)

    # Per-batch metrics + CIs
    rows = []
    for b in np.unique(batch_test):
        m = (batch_test == b)
        if m.sum() == 0:
            continue
        y_b = y_test[m]
        p_b = prob_test[m]
        base = {
            "batch": str(b),
            "n_test": int(m.sum()),
            **overall_metrics(y_b, p_b)
        }
        s_b, i_b = calibration_slope_intercept(y_b, p_b)
        base.update({
            "calibration_slope": s_b,
            "calibration_intercept": i_b
        })
        # CIs via bootstrap, stratified by class within this batch
        strat = y_b  # preserve class composition within batch
        cis_b = bootstrap_cis(
            y_b, p_b,
            n_boot=N_BOOT, seed=RANDOM_SEED + hash(str(b)) % 10000,
            stratify_by=strat
        )
        base.update(cis_b)
        rows.append(base)

    pd.DataFrame(rows).to_csv(f"{OUTDIR}/test_metrics_per_batch.csv", index=False)

    # Save per-sample predictions
    out_pred = pd.DataFrame({
        "sample": X_test.index,
        "batch":  batch_test,
        "y_true": y_test,
        "p_tolerant": prob_test,
        "y_pred_0.5": pred_test
    })
    out_pred.to_csv(f"{OUTDIR}/test_predictions.csv", index=False)

    # Console summary
    print("[SUMMARY] Overall test metrics:")
    print(json.dumps(metrics_overall, indent=2))
    print(f"[INFO] Wrote: {OUTDIR}/scaler.joblib, {OUTDIR}/model.joblib")
    print(f"[INFO] Wrote: {OUTDIR}/model_coefficients.csv, {OUTDIR}/test_metrics_overall.json")
    print(f"[INFO] Wrote: {OUTDIR}/test_metrics_per_batch.csv, {OUTDIR}/test_predictions.csv")


if __name__ == "__main__":
    main()
