#!/usr/bin/env python3
import json
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import (
    roc_curve, roc_auc_score,
    average_precision_score, brier_score_loss
)
from sklearn.calibration import calibration_curve


BASE = Path("results_08_Dec_2025")
FOLDS = sorted(p for p in BASE.glob("loso_*") if p.is_dir())

all_metrics = []
all_coefs = []
all_preds = []

for d in FOLDS:
    holdout = d.name.replace("loso_", "")
    print(f"[INFO] Fold: {holdout}")

    # ---- best_panel_summary_<holdout>.json (required for metrics paths)
    best_json = next(d.glob("best_panel_summary_*.json"), None)
    if not best_json:
        print(f"[WARN] Missing best_panel_summary_*.json in {d}")
        continue
    with open(best_json) as f:
        best = json.load(f)
    # Metrics row
    metrics_row = {
        "batch": best.get("heldout", holdout),
        "n_test": best.get("n_test"),
        "auroc": best.get("auroc"),
        "aupr": best.get("aupr"),
        "brier": best.get("brier"),
        "best_m": best.get("best_m"),
        "genes_kept": best.get("genes_kept"),
    }
    # Optional if present in JSON
    if "log_loss" in best:  # keep column name consistent if you later add it
        metrics_row["log_loss"] = best["log_loss"]
    all_metrics.append(metrics_row)

    # ---- coefficients_<holdout>.csv
    coef_path = d / best.get("coef_file", f"coefficients_{holdout}.csv")
    if coef_path.exists():
        dfc = pd.read_csv(coef_path)
        dfc["batch"] = holdout
        all_coefs.append(dfc)
    else:
        print(f"[WARN] Missing coefficients CSV: {coef_path}")

    # ---- preds_<holdout>.csv
    preds_path = d / best.get("preds_file", f"preds_{holdout}.csv")
    if preds_path.exists():
        dfp = pd.read_csv(preds_path)
        dfp["batch"] = holdout
        all_preds.append(dfp)
    else:
        print(f"[WARN] Missing preds CSV: {preds_path}")

# ---------------------------
# Write summary metrics
# ---------------------------
if all_metrics:
    dfm = pd.DataFrame(all_metrics)
    # Order columns nicely
    cols = ["batch", "n_test", "auroc", "aupr", "brier", "best_m", "genes_kept"]
    if "log_loss" in dfm.columns:
        cols.insert(5, "log_loss")
    dfm = dfm[cols]
    dfm.to_csv(BASE / "loso_summary_metrics.csv", index=False)
    # Weighted means (by n_test)
    w = dfm["n_test"].astype(float).values
    weighted = {
        k: float(np.average(dfm[k].astype(float).values, weights=w))
        for k in dfm.columns if k in {"auroc","aupr","brier","log_loss"}
    }
    print("[INFO] Weighted means:", weighted)
    print(f"[INFO] Wrote {BASE/'loso_summary_metrics.csv'}")
else:
    print("[WARN] No metrics collected.")

# ---------------------------
# Coefficient stability summary
# ---------------------------
if all_coefs:
    dfc = pd.concat(all_coefs, ignore_index=True)
    # Drop intercept for per-gene stats; keep a separate file if you want intercepts
    dfc_no_intercept = dfc[dfc["gene"] != "Intercept"].copy()
    coef_summary = (
        dfc_no_intercept.groupby("gene")["coef"]
        .agg(mean="mean", sd="std", n_folds="count")
        .reset_index()
    )
    # Directional stability across folds
    signs = (
        dfc_no_intercept.assign(sign=lambda x: np.sign(x["coef"]))
        .groupby("gene")["sign"]
        .agg(
            sign_pos_frac=lambda s: (s > 0).mean(),
            sign_neg_frac=lambda s: (s < 0).mean(),
            sign_zero_frac=lambda s: (s == 0).mean(),
        )
        .reset_index()
    )
    coef_summary = coef_summary.merge(signs, on="gene", how="left")
    coef_summary.to_csv(BASE / "loso_coefficients_summary.csv", index=False)

    # Also save intercepts per fold (optional)
    intercepts = dfc[dfc["gene"] == "Intercept"][["batch","coef"]].rename(columns={"coef":"intercept"})
    intercepts.to_csv(BASE / "loso_intercepts_by_fold.csv", index=False)

    print(f"[INFO] Wrote {BASE/'loso_coefficients_summary.csv'} and {BASE/'loso_intercepts_by_fold.csv'}")
else:
    print("[WARN] No coefficients collected.")

# ---------------------------
# Combined predictions
# ---------------------------
if all_preds:
    dfp = pd.concat(all_preds, ignore_index=True)
    # Harmonize expected columns if needed
    rename_map = {
        "true": "true_label",
        "true_y": "true_label",
        "prob": "pred_prob_tolerant",
        "pred": "pred_prob_tolerant",
    }
    for k, v in rename_map.items():
        if k in dfp.columns and v not in dfp.columns:
            dfp = dfp.rename(columns={k: v})

    needed = {"sample", "true_label", "pred_prob_tolerant", "batch"}
    missing = needed - set(dfp.columns)
    if missing:
        print(f"[WARN] Predictions missing columns {missing}; writing what we have.")
    dfp.to_csv(BASE / "loso_all_predictions.csv", index=False)
    print(f"[INFO] Wrote {BASE/'loso_all_predictions.csv'}")

    # pooled AUROC/AUPR quick print (ignores per-fold calibration)
    try:
        y_true = dfp["true_label"].astype(int).values
        y_prob = dfp["pred_prob_tolerant"].astype(float).values
        auroc = roc_auc_score(y_true, y_prob)
        aupr = average_precision_score(y_true, y_prob)
        brier = brier_score_loss(y_true, y_prob)
        print("[INFO] Pooled metrics:\n",
              f"AUROC={auroc:.3f}\n",
              f"AUPR={aupr:.3f}\n",
              f"Brier={brier:.3f}")
    except Exception as e:
        print(f"[WARN] Could not compute pooled AUROC/AUPR: {e}")
        exit(1)
    # -------------------------------
    # ROC curve (overall + per-batch)
    # -------------------------------
    plt.figure(figsize=(6, 6))
    fpr, tpr, _ = roc_curve(y_true, y_prob)
    plt.plot(fpr, tpr, color='black', lw=2, label=f"Overall (AUROC={auroc:.2f})")

    for batch, sub in dfp.groupby("batch"):
        if len(sub["true_label"].unique()) < 2:
            continue  # skip if all one class
        auc = roc_auc_score(sub["true_label"], sub["pred_prob_tolerant"])
        fpr, tpr, _ = roc_curve(sub["true_label"], sub["pred_prob_tolerant"])
        plt.plot(fpr, tpr, lw=1.5, label=f"{batch} (AUROC={auc:.2f})")

    plt.plot([0, 1], [0, 1], 'k--', lw=0.8)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curve by LOSO Fold (Batch)")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(BASE / "roc_by_batch.png", dpi=200)
    plt.close()
    print(f"[INFO] Wrote {BASE/'roc_by_batch.png'}")

    # -------------------------------
    # Calibration plot (overall)
    # -------------------------------
    plt.figure(figsize=(6, 6))
    prob_true, prob_pred = calibration_curve(y_true, y_prob, n_bins=10,
                                             strategy="quantile")

    plt.plot(prob_pred, prob_true, marker="o", label="Observed")
    plt.plot([0, 1], [0, 1], "k--", lw=1, label="Perfect calibration")
    plt.xlabel("Predicted probability")
    plt.ylabel("Observed fraction tolerant")
    plt.title("Calibration Plot (Overall)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(BASE / "calibration_overall.png", dpi=200)
    print(f"[INFO] Wrote {BASE/'calibration_overall.png'}")
    plt.close()

    # -------------------------------
    # Optional: per-batch calibration
    # -------------------------------
    for batch, sub in dfp.groupby("batch"):
        if len(sub["true_label"].unique()) < 2:
            continue
        plt.figure(figsize=(5, 5))
        prob_true, prob_pred = calibration_curve(sub["true_label"],
                                                 sub["pred_prob_tolerant"],
                                                 n_bins=8, strategy="quantile")
        plt.plot(prob_pred, prob_true, marker="o")
        plt.plot([0, 1], [0, 1], "k--", lw=1)
        plt.xlabel("Predicted probability")
        plt.ylabel("Observed fraction tolerant")
        plt.title(f"Calibration: {batch}")
        plt.tight_layout()
        plt.savefig(BASE / f"calibration_{batch.replace(' ', '_')}.png", dpi=200)
        print(f"[INFO] Wrote {BASE/('calibration_'+batch.replace(' ', '_')+'.png')}")
        plt.close()

else:
    print("[WARN] No predictions collected.")

