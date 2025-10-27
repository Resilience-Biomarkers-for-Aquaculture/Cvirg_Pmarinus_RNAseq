#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate ROC and calibration plots from LOSO predictions.

Input: loso_all_predictions.csv with columns:
  sample,true_label,pred_prob_tolerant,batch
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import (
    roc_curve, roc_auc_score,
    average_precision_score, brier_score_loss
)
from sklearn.calibration import calibration_curve

# -------------------------------
# Load data
# -------------------------------
df = pd.read_csv("results/loso_all_predictions.csv")

# Sanity check
print("Shape:", df.shape)
print(df.head())

# -------------------------------
# Compute overall metrics
# -------------------------------
y_true = df["true_label"]
y_prob = df["pred_prob_tolerant"]

auroc = roc_auc_score(y_true, y_prob)
aupr = average_precision_score(y_true, y_prob)
brier = brier_score_loss(y_true, y_prob)

print(f"\nOverall metrics:")
print(f"  AUROC = {auroc:.3f}")
print(f"  AUPR  = {aupr:.3f}")
print(f"  Brier = {brier:.3f}")

# -------------------------------
# ROC curve (overall + per-batch)
# -------------------------------
plt.figure(figsize=(6, 6))
fpr, tpr, _ = roc_curve(y_true, y_prob)
plt.plot(fpr, tpr, color='black', lw=2, label=f"Overall (AUROC={auroc:.2f})")

for batch, sub in df.groupby("batch"):
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
plt.savefig("roc_by_batch.png", dpi=200)
plt.close()

# -------------------------------
# Calibration plot (overall)
# -------------------------------
plt.figure(figsize=(6, 6))
prob_true, prob_pred = calibration_curve(y_true, y_prob, n_bins=10, strategy="quantile")

plt.plot(prob_pred, prob_true, marker="o", label="Observed")
plt.plot([0, 1], [0, 1], "k--", lw=1, label="Perfect calibration")
plt.xlabel("Predicted probability")
plt.ylabel("Observed fraction tolerant")
plt.title("Calibration Plot (Overall)")
plt.legend()
plt.tight_layout()
plt.savefig("calibration_overall.png", dpi=200)
plt.close()

# -------------------------------
# Optional: per-batch calibration
# -------------------------------
for batch, sub in df.groupby("batch"):
    if len(sub["true_label"].unique()) < 2:
        continue
    plt.figure(figsize=(5, 5))
    prob_true, prob_pred = calibration_curve(sub["true_label"], sub["pred_prob_tolerant"],
                                             n_bins=8, strategy="quantile")
    plt.plot(prob_pred, prob_true, marker="o")
    plt.plot([0, 1], [0, 1], "k--", lw=1)
    plt.xlabel("Predicted probability")
    plt.ylabel("Observed fraction tolerant")
    plt.title(f"Calibration: {batch}")
    plt.tight_layout()
    plt.savefig(f"calibration_{batch.replace(' ', '_')}.png", dpi=200)
    plt.close()
