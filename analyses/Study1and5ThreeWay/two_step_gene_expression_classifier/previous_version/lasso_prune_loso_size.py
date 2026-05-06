
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
OBSOLETE:
Note this was originally the second script in the two-step gene
expression classifier pipeline, and has been superseded by
run_loso_pipeline.py (an orchestrator) and lasso_prune_onefold.py.py
(the per-fold classifier).
End-to-end classifier pipeline using:
  - nf-core/differentialabundance DESEQ2_NORM VST matrix (genes x samples)
  - Original metadata CSV (raw; cleaned inline)
  - Tier1+Tier2 gene list restriction (to assist logistic regression convergence)
Implements: stability lasso + redundancy pruning + LOSO/LOGO evaluation.

Output files:
  stability_selection.csv (gene, selection_prob, mean_abs_coef_when_selected)
  panel_size_sweep_metrics.csv (m, auroc_mean, auroc_min, aupr_mean, brier_mean)
  final_panel_gene_list.txt (one gene per line)
  loso_by_heldout.csv (heldout, n_te, auroc, aupr, brier per held-out group)
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Tuple

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedShuffleSplit
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss
from sklearn.calibration import CalibratedClassifierCV


# -------------------------------
# User-editable paths for input data
# -------------------------------
VST_PATH = "DESEQ2_NORM_all.vst.tsv"  # nf-core DESEQ2_NORM VST matrix (genes x samples; first col = gene_id)
META_PATH = "../../../data/differential_abundance_sheets/rnaseq_diffabundance_study1and5_samplesheet_filled.csv"
GENE_LIST_PATH = "panel_candidates_tier12.txt"

# Optional: remap free-form batch labels to abstract groups (A,B1,B2) for LOGO clarity
# Leave empty to use your raw 'batch' labels as-is (e.g., "P&S 2023").
BATCH_MAP = {
    # "P&S 2023": "A",
    # "P&S 2020 2015": "B1",
    # "P&S 2020 2017": "B2",
}

# Optional: mitigate dominance if one batch (e.g., "A") is much larger.
CAP_A_MULTIPLIER = None  # e.g., 2 caps "A" to ≤ 2x the size of non-A in each training fold; None disables.

# Panel sizes to evaluate and selection tolerance
TARGET_SIZES = [6, 8, 10, 12, 16, 24]
DELTA_AUROC = 0.02  # choose smallest m within 0.02 AUROC of best


# -------------------------------
# I/O helpers
# -------------------------------
def load_expr_from_vst(vst_path: str) -> pd.DataFrame:
    """
    Load genes x samples VST (first column = gene_id) and return samples x genes with 'sample' column.
    """
    df = pd.read_csv(vst_path, sep="\t")
    # Normalize header name for gene id (first column)
    df = df.rename(columns={df.columns[0]: "gene_id"})
    # Transpose: make rows = samples, cols = genes
    expr = df.set_index("gene_id").T
    expr.index.name = "sample"
    expr.reset_index(inplace=True)
    return expr  # samples x genes, first column 'sample'


def load_metadata_raw(meta_path: str) -> pd.DataFrame:
    """
    Load original metadata CSV and normalize in-memory:
      - strip BOM/spaces from headers
      - standardize key column names
      - strip whitespace, lower-case condition
      - optional batch remap
    Returns a DataFrame indexed by 'sample'.
    """
    meta = pd.read_csv(meta_path)

    # Clean header names
    meta.columns = [c.replace("\ufeff", "").strip() for c in meta.columns]

    # Standardize column names
    rename_map = {}
    for c in meta.columns:
        cl = c.lower()
        if cl == "sample": rename_map[c] = "sample"
        elif cl == "condition": rename_map[c] = "condition"
        elif cl == "treatment": rename_map[c] = "treatment"
        elif cl == "batch": rename_map[c] = "batch"
        elif cl == "collection_interval_days": rename_map[c] = "Collection_Interval_Days"
    meta = meta.rename(columns=rename_map)

    required = {"sample", "condition", "batch"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"Metadata missing required columns: {missing}")

    # Value normalization
    meta["sample"] = meta["sample"].astype(str).str.strip()
    meta["condition"] = meta["condition"].astype(str).str.strip().str.lower()
    meta["batch"] = meta["batch"].astype(str).str.strip()

    if BATCH_MAP:
        meta["batch"] = meta["batch"].map(lambda x: BATCH_MAP.get(x, x))

    return meta.set_index("sample")


def align_expr_meta(expr_sxg: pd.DataFrame, meta: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Align expression (samples x genes with 'sample' col) and metadata (indexed by sample).
    Returns X_df indexed by sample, and aligned meta.
    """
    if "sample" not in expr_sxg.columns:
        raise ValueError("Expression frame must contain a 'sample' column.")
    expr_sxg["sample"] = expr_sxg["sample"].astype(str).str.strip()
    expr_sxg = expr_sxg.set_index("sample")

    common = meta.index.intersection(expr_sxg.index)
    if len(common) == 0:
        # Diagnostics to help locate mismatch causes
        raise ValueError(
            "No overlapping sample IDs between expression and metadata. "
            "Check SRX/SRR mismatches, trailing spaces, or differing ID conventions."
        )
    expr_sxg = expr_sxg.loc[common]
    meta = meta.loc[common]
    return expr_sxg, meta


# -------------------------------
# Modeling helpers
# -------------------------------
def composite_strata(y, batch):
    """Combine batch and class for balanced resampling in stability selection."""
    return np.array([f"{b}__{c}" for b, c in zip(batch, y)])


def fit_l1_logreg(X, y, C=0.1, l1_ratio=None, seed=0):
    """Logistic regression with L1 (lasso) or elastic net (if l1_ratio provided)."""
    penalty = "l1" if l1_ratio is None else "elasticnet"
    print(f"[DEBUG] fit_li_logreg")
    clf = LogisticRegression(
        penalty=penalty,
        l1_ratio=l1_ratio,
        C=C,
        solver="saga",
        max_iter=5000,
        n_jobs=-1,
        class_weight="balanced",
        random_state=seed,
    )
    clf.fit(X, y)
    return clf


def stability_selection(X_df, y, batch, n_resamples=300, subsample_frac=0.7,
                        C_grid=(0.05, 0.1, 0.2), l1_ratio=None, random_state=123) -> pd.DataFrame:
    """
    Repeated subsampling (stratified by batch×class) + sparse model;
    returns per-gene selection probability and mean |coef| when selected.
    """
    rng = np.random.RandomState(random_state)
    genes = X_df.columns.to_numpy()
    sel_counts = pd.Series(0.0, index=genes)
    coef_sums = pd.Series(0.0, index=genes)

    strata = composite_strata(y, batch)
    sss = StratifiedShuffleSplit(n_splits=n_resamples, train_size=subsample_frac,
                                 random_state=random_state)

    X = X_df.to_numpy()
    for i, (tr, _) in enumerate(sss.split(np.zeros_like(strata), strata), 1):
        scaler = StandardScaler().fit(X[tr])     # fit transforms on training only
        Xt = scaler.transform(X[tr])

        best_coef, best_nz = None, 10**9
        for C in C_grid:
            coef = fit_l1_logreg(Xt, y[tr], C=C, l1_ratio=l1_ratio,
                                 seed=rng.randint(1_000_000_000)).coef_.ravel()
            nz = np.flatnonzero(coef != 0)
            if len(nz) < best_nz:
                best_nz = len(nz)
                best_coef = coef

        nz = np.flatnonzero(best_coef != 0)
        if len(nz):
            sel_counts.iloc[nz] += 1.0
            coef_sums.iloc[nz]  += np.abs(best_coef[nz])

    stab = pd.DataFrame({
        "gene": genes,
        "selection_prob": (sel_counts / n_resamples).values,
        "mean_abs_coef_when_selected": (coef_sums / np.maximum(sel_counts, 1.0)).values
    }).sort_values(["selection_prob","mean_abs_coef_when_selected"], ascending=False)
    return stab


def redundancy_prune(X_df, genes, corr_threshold=0.9):
    """
    Greedy removal of highly correlated genes (|r| > threshold).
    Keeps earlier genes in priority order; drops later near-duplicates.
    """
    kept = []
    corr = np.corrcoef(X_df[genes].to_numpy(), rowvar=False)
    name_to_idx = {g: i for i, g in enumerate(genes)}
    for g in genes:
        gi = name_to_idx[g]
        redundant = any(abs(corr[gi, name_to_idx[h]]) > corr_threshold for h in kept)
        if not redundant:
            kept.append(g)
    return kept


def loso_eval(X_df, y, batch, genes, calibrate=True, seed=123, cap_A_multiplier=None):
    """
    Leave-One-Group-Out where 'batch' defines groups.
    Equal-weight averaging across held-out groups.
    Optional cap of group named "A" in training folds (see CAP_A_MULTIPLIER).
    """
    heldouts = np.unique(batch)
    per_held = []

    for held in heldouts:
        tr = (batch != held)
        te = (batch == held)

        # Optional: cap overly large batch "A" in the training set
        if cap_A_multiplier is not None:
            tr_idx = X_df.index[tr]
            train_batches = pd.Series(batch[tr], index=tr_idx)
            is_A = (train_batches == "A")
            if is_A.any():
                nonA_n = (~is_A).sum()
                cap = min(is_A.sum(), max(1, int(cap_A_multiplier * nonA_n)))
                if is_A.sum() > cap:
                    A_idx = tr_idx[is_A.values]
                    nonA_idx = tr_idx[~is_A.values]
                    A_cond = pd.Series(y[tr], index=tr_idx).loc[A_idx]
                    rng = np.random.RandomState(123)
                    # proportional sample within A by class
                    A_keep = (A_cond.groupby(A_cond)
                               .apply(lambda s: s.sample(max(1, int(round(len(s)*cap/len(A_idx)))), random_state=rng))
                               .index.get_level_values(1))
                    keep_idx = pd.Index(nonA_idx).append(A_keep)
                    tr = X_df.index.isin(keep_idx)

        X_tr = X_df.loc[tr, genes].to_numpy(); y_tr = y[tr]
        X_te = X_df.loc[te, genes].to_numpy(); y_te = y[te]

        scaler = StandardScaler().fit(X_tr)
        X_tr = scaler.transform(X_tr)
        X_te = scaler.transform(X_te)

        base = LogisticRegression(penalty=None, max_iter=5000, n_jobs=-1)
        clf = CalibratedClassifierCV(base, method="sigmoid", cv=5) if calibrate else base
        clf.fit(X_tr, y_tr)
        p = clf.predict_proba(X_te)[:, 1]

        per_held.append({
            "heldout": held,
            "n_te": int(te.sum()),
            "auroc": roc_auc_score(y_te, p),
            "aupr":  average_precision_score(y_te, p),
            "brier": brier_score_loss(y_te, p)
        })

    df = pd.DataFrame(per_held)
    summary = {
        "auroc_mean": df["auroc"].mean(),   # equal-weight across held-out groups
        "auroc_min":  df["auroc"].min(),
        "aupr_mean":  df["aupr"].mean(),
        "brier_mean": df["brier"].mean()
    }
    return df, summary


# -------------------------------
# Main
# -------------------------------
if __name__ == "__main__":
    # 1) Load expression (genes x samples) and transpose to samples x genes
    expr = load_expr_from_vst(VST_PATH)
    print(f"[DEBUG] expr (samples x genes) shape before alignment: {expr.shape}")

    # 2) Load raw metadata and clean inline
    meta = load_metadata_raw(META_PATH)

    # 3) Align by 'sample'
    X_df, meta = align_expr_meta(expr, meta)
    print(f"[DEBUG] After alignment: X_df shape={X_df.shape}, meta shape={meta.shape}")
    print(f"[DEBUG] Meta columns: {list(meta.columns)}")

    # 4) Encode labels and grouping
    cond_map = {"tolerant": 1, "sensitive": 0}
    unknown = set(meta["condition"].unique()) - set(cond_map.keys())
    if unknown:
        raise ValueError(f"Unexpected condition values: {unknown}")
    y = meta["condition"].map(cond_map).astype(int).values
    batch = meta["batch"].astype(str).values

    print(f"[DEBUG] Condition counts: {meta['condition'].value_counts().to_dict()}")
    print(f"[DEBUG] Unique batches: {np.unique(batch)}")

    # 5) Feature restriction (Tier1+Tier2)
    if GENE_LIST_PATH and Path(GENE_LIST_PATH).exists():
        keep_genes = pd.read_csv(GENE_LIST_PATH, header=None, names=["gene"])["gene"].tolist()
        keep_genes = [g for g in keep_genes if g in X_df.columns]
        X_df = X_df[keep_genes]
        print(f"[DEBUG] Restricted to Tier1+Tier2 features present in X: {X_df.shape[1]} genes")
    else:
        print(f"Missing or empty GENE_LIST_PATH {GENE_LIST_PATH}; cannot proceed.")
        exit(1)

    # Sanity: keep only numeric columns (should already be numeric for VST)
    numeric_cols = X_df.columns[X_df.dtypes.apply(lambda t: np.issubdtype(t, np.number))]
    if len(numeric_cols) != X_df.shape[1]:
        X_df = X_df[numeric_cols]
        print(f"[DEBUG] Dropped non-numeric columns; new shape: {X_df.shape}")

    # 6) Stability selection to rank features
    stab = stability_selection(
        X_df=X_df, y=y, batch=batch,
        n_resamples=300, subsample_frac=0.7,
        C_grid=(0.05, 0.1, 0.2), l1_ratio=None, random_state=123
    )
    stab.to_csv("stability_selection.csv", index=False)
    print("[DEBUG] Wrote stability_selection.csv")

    # 7) Panel size sweep with redundancy pruning + LOSO
    ordered_genes = stab["gene"].tolist()
    records = []
    perheld_store = {}

    for m in TARGET_SIZES:
        topm = ordered_genes[:max(m*2, m)]      # oversample then prune
        pruned = redundancy_prune(X_df, topm, corr_threshold=0.9)
        if len(pruned) > m:
            pruned = pruned[:m]

        per_held, summary = loso_eval(
            X_df, y, batch, pruned, calibrate=True, seed=123,
            cap_A_multiplier=CAP_A_MULTIPLIER
        )

        records.append({
            "m": m,
            "genes_kept": len(pruned),
            "auroc_mean": summary["auroc_mean"],
            "auroc_min":  summary["auroc_min"],
            "aupr_mean":  summary["aupr_mean"],
            "brier_mean": summary["brier_mean"]
        })
        perheld_store[m] = (pruned, per_held)

    sweep_df = pd.DataFrame(records).sort_values("m")
    sweep_df.to_csv("panel_size_sweep_metrics.csv", index=False)
    print("[DEBUG] Wrote panel_size_sweep_metrics.csv")
    print(sweep_df)

    # 8) Choose smallest m within ΔAUROC ≤ best; export final panel & LOSO breakdown
    best_auroc = sweep_df["auroc_mean"].max()
    eligible = sweep_df[sweep_df["auroc_mean"] >= best_auroc - DELTA_AUROC].sort_values("m")
    chosen_m = int(eligible.iloc[0]["m"])
    final_genes, per_held = perheld_store[chosen_m]

    pd.Series(final_genes).to_csv("final_panel_gene_list.txt", index=False, header=False)
    per_held.to_csv("loso_by_heldout.csv", index=False)
    print(f"[DEBUG] Chosen m={chosen_m} with mean AUROC≈{eligible.iloc[0]['auroc_mean']:.3f}")
    print("[DEBUG] Wrote final_panel_gene_list.txt and loso_by_heldout.csv")
