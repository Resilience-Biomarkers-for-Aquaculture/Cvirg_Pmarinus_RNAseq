#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Consensus minimal classifier genes from LOSO outputs.

Usage:
  python loso_consensus_genes.py --metric auroc
  python loso_consensus_genes.py --metric aupr --panel_size 12 --results_dir results

Inputs expected under --results_dir (default: results):
  loso_*/best_panel_summary_*.json  (per-fold metrics incl. 'heldout', 'auroc','aupr','brier','best_m')
  loso_*/selected_genes_*_m{m}.txt  (per-fold selected genes at panel size m)
  loso_summary_metrics.csv          (batch,n_test,auroc,aupr,brier,best_m,genes_kept)
  loso_coefficients_summary.csv     (gene,mean,sd,n_folds,sign_pos_frac,sign_neg_frac,sign_zero_frac)
"""

import argparse
import json
from pathlib import Path
import pandas as pd

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results_dir", default="results")
    ap.add_argument("--metric", choices=["auroc", "aupr"], default="auroc",
                    help="Fold weight for voting")
    ap.add_argument("--panel_size", type=int, default=None,
                    help="If given, force this m; otherwise auto-select")
    ap.add_argument("--output", default="final_consensus_genes.csv")
    ap.add_argument("--min_sign_stable_frac", type=float, default=2/3)
    args = ap.parse_args()

    root = Path(args.results_dir)

    # ---- Step 1: per-fold metrics from JSONs
    metrics_rows = []
    for jf in root.glob("loso_*/best_panel_summary_*.json"):
        with open(jf) as f:
            d = json.load(f)
        metrics_rows.append({
            "fold": d.get("heldout"),
            "auroc": d.get("auroc"),
            "aupr":  d.get("aupr"),
            "brier": d.get("brier"),
            "best_m_json": d.get("best_m"),
            "json_path": str(jf),
        })
    metrics_df = pd.DataFrame(metrics_rows)
    if metrics_df.empty:
        raise FileNotFoundError("No loso_*/best_panel_summary_*.json files found.")
    print("[shape-check] Folds:", metrics_df["fold"].tolist())
    print(metrics_df[["fold","auroc","aupr","brier","best_m_json"]])

    # ---- Step 2: decide panel size m
    m = args.panel_size
    if m is None:
        summary_path = root / "loso_summary_metrics.csv"
        if summary_path.exists():
            sm = pd.read_csv(summary_path)
            needed = {"batch","n_test","auroc","aupr","brier","best_m","genes_kept"}
            if needed.issubset(sm.columns):
                # choose mode of best_m across folds (ties → take smallest)
                m = int(sm["best_m"].mode().sort_values().iloc[0])
                print(f"[panel-size] Using mode(best_m) from loso_summary_metrics.csv → m={m}")
            else:
                print("[panel-size] loso_summary_metrics.csv found but headers differ; "
                      "falling back to per-fold JSON best_m.")
        if m is None:
            m = int(metrics_df["best_m_json"].mode().sort_values().iloc[0])
            print(f"[panel-size] Using mode(best_m) from per-fold JSON → m={m}")
    else:
        print(f"[panel-size] Using user-provided m={m}")

    # ---- Step 3: weighted voting using chosen metric
    # Normalize weights to sum to 1 for robustness
    wcol = args.metric
    metrics_df["weight_raw"] = metrics_df[wcol].fillna(0.0)
    total = metrics_df["weight_raw"].sum()
    if total <= 0:
        metrics_df["weight"] = 1.0 / len(metrics_df)
        print(f"[weights] {wcol} missing/nonpositive → using uniform weights.")
    else:
        metrics_df["weight"] = metrics_df["weight_raw"] / total
        print(f"[weights] Using normalized {wcol} weights (sum=1).")

    # Gather per-fold gene lists at size m
    weighted_counts = {}
    missing = []
    for _, r in metrics_df.iterrows():
        fold = r["fold"]
        w = float(r["weight"])
        gene_file = root / f"loso_{fold}" / f"selected_genes_{fold}_m{m}.txt"
        if not gene_file.exists():
            missing.append(str(gene_file))
            continue
        with open(gene_file) as f:
            genes = [ln.strip() for ln in f if ln.strip()]
        for g in genes:
            weighted_counts[g] = weighted_counts.get(g, 0.0) + w
    if missing:
        print("[warn] Missing per-fold gene files for m=", m)
        for p in missing:
            print("  -", p)

    weighted_df = (pd.DataFrame(list(weighted_counts.items()), columns=["gene","weighted_vote"])
                     .sort_values("weighted_vote", ascending=False)
                     .reset_index(drop=True))
    weighted_df["rank"] = range(1, len(weighted_df)+1)
    weighted_df.insert(1, "panel_size_used", m)
    print("[shape-check] Unique genes tallied:", len(weighted_df))
    print(weighted_df.head(10))

    # ---- Step 4: join sign stability from loso_coefficients_summary.csv
    coef_path = root / "loso_coefficients_summary.csv"
    if coef_path.exists():
        coef = pd.read_csv(coef_path)
        need_coef = {"gene","mean","sd","n_folds","sign_pos_frac","sign_neg_frac","sign_zero_frac"}
        if need_coef.issubset(coef.columns):
            out = weighted_df.merge(coef[list(need_coef)], on="gene", how="left")
            pos = out["sign_pos_frac"].fillna(0.0)
            neg = out["sign_neg_frac"].fillna(0.0)
            zero = out["sign_zero_frac"].fillna(0.0)
            thr = float(args.min_sign_stable_frac)
            out["sign_stable"] = ((pos >= thr) | (neg >= thr)) & (zero <= (1.0 - thr))
            weighted_df = out
        else:
            print("[warn] Unexpected headers in loso_coefficients_summary.csv; skipping join.")
    else:
        print("[info] loso_coefficients_summary.csv not found; skipping sign stability.")

    # ---- Step 5: write outputs
    # Set column order
    weighted_df = weighted_df[["gene", "rank", "weighted_vote", "mean", "sd", "n_folds",
                               "sign_pos_frac", "sign_neg_frac", "sign_zero_frac", "sign_stable"]]
    weighted_df.to_csv(args.output, index=False)
    print(f"[done] Wrote consensus ranking → {args.output}")

    # ---- Step 6: optional quick plot
    try:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(weighted_df["rank"], weighted_df["weighted_vote"])
        plt.title(f"Consensus weighted votes ({args.metric.upper()}), m={m}")
        plt.xlabel("Rank")
        plt.ylabel("Weighted vote (normalized)")
        plt.tight_layout()
        png = Path(args.output).with_suffix(".png")
        plt.savefig(png, dpi=200)
        print(f"[done] Vote curve saved → {png}")
    except Exception as e:
        print("[plot] Skipped:", e)

if __name__ == "__main__":
    main()


