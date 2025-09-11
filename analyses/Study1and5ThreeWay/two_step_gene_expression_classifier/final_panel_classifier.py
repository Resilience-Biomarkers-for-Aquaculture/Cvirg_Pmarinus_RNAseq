#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Repeated seeds: frozen 6-gene (or N-gene) classifier
- For each seed:
  * constrained per-batch stratified holdout by 'condition'
  * train pooled model with per-gene z-scoring (fit on train only)
  * evaluate on held-out test set (overall + per-batch)
  * save coefficients, metrics, predictions

- After all seeds:
  * aggregate metrics distributions (overall + per-batch)
  * summarize coefficient variability and sign stability

Inputs:
  VST_PATH   : nf-core DESEQ2_NORM VST (genes x samples; first col = gene_id)
  META_PATH  : original metadata CSV (must have: sample, condition, batch)
  PANEL_PATH : final_panel_gene_list.txt (one gene per line; any length)

Outputs:
  runs/<seed>/...                         # per-seed artifacts
  summary/metrics_overall_all_seeds.csv   # one row per seed
  summary/metrics_per_batch_all_seeds.csv # one row per seed x batch
  summary/coefficients_all_seeds.csv      # one row per seed; columns=genes+INTERCEPT
  summary/coefficients_stability.csv      # per-gene sign stability and stats
"""

from pathlib import Path
import json
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
VST_PATH = "DESEQ2_NORM_all.vst.tsv"
META_PATH = "../../data/differential_abundance_sheets/rnaseq_diffabundance_study1and5_samplesheet_filled.csv"
PANEL_PATH = "final_panel_gene_list.txt"

OUTROOT    = "repeat_seeds_fixed_model"
SEEDS      = [12345 + i for i in range(30)]   # 30 runs; adjust as needed

# Holdout fraction per batch (target; constraints apply per class within batch)
HOLDOUT_FRAC = 0.25  # 0.20–0.30 is typical

# Safety constraints for tiny classes (per batch, per class)
MIN_TEST_PER_CLASS  = 1
MIN_TRAIN_PER_CLASS = 3

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
    df = pd.read_csv(vst_path, sep="\t")
    df = df.rename(columns={df.columns[0]: "gene_id"}).set_index("gene_id")
    sxg = df.T
    sxg.index.name = "sample"
    return sxg

def load_meta(meta_path: str) -> pd.DataFrame:
    meta = pd.read_csv(meta_path)
    meta.columns = [c.replace("\ufeff", "").strip() for c in meta.columns]
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
# Split & metrics helpers
# -------------------------------
def stratified_holdout_by_batch(meta: pd.DataFrame,
                                y: np.ndarray,
                                hold_frac=0.25,
                                seed=123,
                                min_test_per_class=1,
                                min_train_per_class=3):
    """
    For each batch, sample a test set stratified by class, enforcing:
      - >= min_test_per_class in test per class (if feasible)
      - >= min_train_per_class in train per class
    If constraints cannot be met for a class in a batch, that class in that batch contributes 0 to test.
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
            k = int(round(n_bc * hold_frac))
            if n_bc < (min_test_per_class + min_train_per_class):
                k = 0
            else:
                k = max(k, min_test_per_class)
                k = min(k, n_bc - min_train_per_class)
            if k > 0:
                take = rng.choice(idx_bc, size=k, replace=False)
                test_flags[take] = True
    train_flags = ~test_flags
    assert train_flags.any() and test_flags.any(), "Empty train or test split after holdout."
    return train_flags, test_flags

def overall_metrics(y, p):
    return {
        "auroc": float(roc_auc_score(y, p)),
        "aupr":  float(average_precision_score(y, p)),
        "brier": float(brier_score_loss(y, p)),
        "log_loss": float(log_loss(y, np.clip(p, 1e-9, 1 - 1e-9))),
    }

def calibration_slope_intercept(y_true, y_prob, eps=1e-12):
    from sklearn.linear_model import LogisticRegression
    p = np.clip(y_prob, eps, 1 - eps)
    logit = np.log(p / (1 - p))
    lr = LogisticRegression(penalty=None, solver="lbfgs", max_iter=1000)
    lr.fit(logit.reshape(-1, 1), y_true)
    return float(lr.coef_.ravel()[0]), float(lr.intercept_.ravel()[0])

# -------------------------------
# Main multi-seed experiment
# -------------------------------
def main():
    outroot = Path(OUTROOT)
    (outroot / "runs").mkdir(parents=True, exist_ok=True)
    (outroot / "summary").mkdir(parents=True, exist_ok=True)

    # Load once
    X_all = load_vst(VST_PATH)        # samples x genes
    meta  = load_meta(META_PATH)      # indexed by sample
    genes = pd.read_csv(PANEL_PATH, header=None)[0].astype(str).str.strip().tolist()

    # Align and subset
    common = X_all.index.intersection(meta.index)
    if len(common) == 0:
        raise ValueError("No overlapping sample IDs between VST and metadata.")
    X_all = X_all.loc[common]
    meta  = meta.loc[common]

    missing_genes = [g for g in genes if g not in X_all.columns]
    if missing_genes:
        raise ValueError(f"Panel genes missing in expression matrix: {missing_genes}")

    X = X_all[genes].copy()
    cond_map = {"tolerant": 1, "sensitive": 0}
    unknown = set(meta["condition"].unique()) - set(cond_map.keys())
    if unknown:
        raise ValueError(f"Unexpected condition values: {unknown}")
    y = meta["condition"].map(cond_map).astype(int).values
    batch = meta["batch"].values

    # Accumulators across seeds
    rows_overall = []
    rows_perbatch = []
    rows_coefs = []

    for seed in SEEDS:
        run_dir = outroot / "runs" / f"seed_{seed}"
        run_dir.mkdir(parents=True, exist_ok=True)

        # Split
        train_mask, test_mask = stratified_holdout_by_batch(
            meta, y,
            hold_frac=HOLDOUT_FRAC,
            seed=seed,
            min_test_per_class=MIN_TEST_PER_CLASS,
            min_train_per_class=MIN_TRAIN_PER_CLASS
        )
        X_train, y_train = X.iloc[train_mask], y[train_mask]
        X_test,  y_test  = X.iloc[test_mask],  y[test_mask]
        batch_test       = batch[test_mask]

        # Scale on train
        scaler = StandardScaler().fit(X_train.values)
        Xz_train = scaler.transform(X_train.values)
        Xz_test  = scaler.transform(X_test.values)

        # Train logistic regression (unpenalized; balanced)
        clf = LogisticRegression(
            penalty=None, solver="lbfgs", max_iter=5000,
            class_weight="balanced", n_jobs=-1
        ).fit(Xz_train, y_train)

        # Save artifacts for this run (optional; comment out if not needed)
        joblib.dump(scaler, run_dir / "scaler.joblib")
        joblib.dump(clf,    run_dir / "model.joblib")

        # Save readable coefs for this run
        coef = pd.Series(clf.coef_.ravel(), index=genes)
        coef_df = pd.DataFrame({
            "feature": list(coef.index) + ["INTERCEPT"],
            "coefficient": list(coef.values) + [float(clf.intercept_.ravel()[0])]
        })
        coef_df.to_csv(run_dir / "model_coefficients.csv", index=False)

        # Inference
        prob_test = clf.predict_proba(Xz_test)[:, 1]
        pred_test = (prob_test >= 0.5).astype(int)

        # Overall metrics
        met = {"seed": seed, "n_test": int(len(y_test))}
        met.update(overall_metrics(y_test, prob_test))
        slope, intercept = calibration_slope_intercept(y_test, prob_test)
        met.update({"cal_slope": slope, "cal_intercept": intercept})
        rows_overall.append(met)

        # Per-batch metrics
        for b in np.unique(batch_test):
            m = (batch_test == b)
            if m.sum() == 0: 
                continue
            y_b = y_test[m]; p_b = prob_test[m]
            base = {"seed": seed, "batch": str(b), "n_test": int(m.sum())}
            base.update(overall_metrics(y_b, p_b))
            s_b, i_b = calibration_slope_intercept(y_b, p_b)
            base.update({"cal_slope": s_b, "cal_intercept": i_b})
            rows_perbatch.append(base)

        # Collect coefficients row (wide)
        wide = {g: float(coef[g]) for g in genes}
        wide["INTERCEPT"] = float(clf.intercept_.ravel()[0])
        wide["seed"] = seed
        rows_coefs.append(wide)

        # Save per-sample predictions for this run (optional)
        pd.DataFrame({
            "sample": X_test.index,
            "batch":  batch_test,
            "y_true": y_test,
            "p_tolerant": prob_test,
            "y_pred_0.5": pred_test
        }).to_csv(run_dir / "test_predictions.csv", index=False)

        # Save quick per-run summary JSON
        with open(run_dir / "test_metrics_overall.json", "w") as f:
            json.dump(met, f, indent=2)

    # --- Aggregate across seeds ---
    summary_dir = outroot / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    overall_df   = pd.DataFrame(rows_overall).sort_values("seed")
    perbatch_df  = pd.DataFrame(rows_perbatch).sort_values(["batch","seed"])
    coefs_df     = pd.DataFrame(rows_coefs).sort_values("seed")

    overall_df.to_csv(summary_dir / "metrics_overall_all_seeds.csv", index=False)
    perbatch_df.to_csv(summary_dir / "metrics_per_batch_all_seeds.csv", index=False)
    coefs_df.to_csv(summary_dir / "coefficients_all_seeds.csv", index=False)

    # Coefficient stability summary
    coef_cols = [c for c in coefs_df.columns if c not in ("seed", "INTERCEPT")]
    stab_rows = []
    for g in coef_cols:
        vals = coefs_df[g].values
        sign_pos = np.mean(vals > 0.0)   # proportion of runs with positive sign
        sign_neg = np.mean(vals < 0.0)
        sign_zero = np.mean(vals == 0.0)
        stab_rows.append({
            "gene": g,
            "mean_coef": float(np.mean(vals)),
            "median_coef": float(np.median(vals)),
            "sd_coef": float(np.std(vals, ddof=1)),
            "sign_pos_frac": float(sign_pos),
            "sign_neg_frac": float(sign_neg),
            "sign_zero_frac": float(sign_zero),
            "min_coef": float(np.min(vals)),
            "max_coef": float(np.max(vals)),
        })
    # Intercept variability too
    ival = coefs_df["INTERCEPT"].values
    stab_rows.append({
        "gene": "INTERCEPT",
        "mean_coef": float(np.mean(ival)),
        "median_coef": float(np.median(ival)),
        "sd_coef": float(np.std(ival, ddof=1)),
        "sign_pos_frac": float(np.mean(ival > 0)),
        "sign_neg_frac": float(np.mean(ival < 0)),
        "sign_zero_frac": float(np.mean(ival == 0)),
        "min_coef": float(np.min(ival)),
        "max_coef": float(np.max(ival)),
    })

    stab_df = pd.DataFrame(stab_rows)
    stab_df.to_csv(summary_dir / "coefficients_stability.csv", index=False)

    # Quick console summary
    def q(x):  # 2.5%, 50%, 97.5%
        return np.percentile(x, [2.5, 50, 97.5]).tolist()

    print("\n[SUMMARY across seeds] Overall AUROC / AUPR (2.5%, 50%, 97.5%):")
    print("AUROC:", q(overall_df["auroc"].values))
    print("AUPR :", q(overall_df["aupr"].values))
    print("\n[Coefficient sign stability] (first few rows):")
    print(stab_df.head())

if __name__ == "__main__":
    main()
