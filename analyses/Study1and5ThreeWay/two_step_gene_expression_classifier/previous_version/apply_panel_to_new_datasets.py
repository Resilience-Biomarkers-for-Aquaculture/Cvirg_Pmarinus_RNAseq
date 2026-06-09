#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
apply_panel_to_new_datasets.py
==============================
Leave-one-study-out (LOSO) evaluation of a fixed gene panel classifier
across the three independent datasets that were NOT used during panel
selection (Studies 2, 3, and 4).

Background
----------
The 6-gene panel was selected using the pipeline in this directory
(lasso_prune_loso_size.py), working exclusively from datasets 1 and 5
(P&S lab challenge experiments).  The panel is now FIXED; this script
evaluates how well those six genes generalise to three completely
independent studies with different experimental designs and labs.

This script can also be run with alternative gene panels (e.g. the
consensus panel derived from the parent-directory LOSO pipeline) by
supplying a different --panel argument.

LOSO design
-----------
Because we have three independent test studies (A=Study 2, B=Study 3,
C=Study 4) we use leave-one-study-out cross-validation:

    Fold 1 : Train on (B + C)  →  Test on A  (Study 2 held out)
    Fold 2 : Train on (A + C)  →  Test on B  (Study 3 held out)
    Fold 3 : Train on (A + B)  →  Test on C  (Study 4 held out)

In each fold, a fresh StandardScaler and logistic regression are fit on
the two training studies using the fixed gene panel, then applied to
the held-out study.

Datasets
--------
  Study 2  – Louisiana field survey; tolerant/sensitive labeled by
             field hypnospore load; 40 samples, gill tissue.
  Study 3  – Dermo injection challenge (Chinese Academy of Sciences);
             10 samples, gill tissue.
  Study 4  – VIMS injection challenge; 44 samples, mantle tissue;
             two sub-batches (2015 and 2017).

Normalisation note
------------------
The panel was selected using DESeq2 VST-normalised expression values
(DESEQ2_NORM_all.vst.tsv), which is only available for datasets 1 and
5.  Datasets 2–4 exist only as raw Salmon count matrices.  To keep
training and test normalisation consistent within each LOSO fold, this
script applies log2(count + 1) to all three datasets.  Across the 82
training samples that DO have both, log2(count+1) has Pearson r =
0.93–0.97 with VST for each of the six panel genes, confirming the
approximation is adequate.

Outputs (written to --outdir, default apply_panel_results/)
-----------------------------------------------------------
  <outdir>/
    fold1_heldout_study2/
      predictions.csv         – per-sample probabilities + true labels
      metrics.json            – AUROC, AUPR, Brier for this fold
      model_coefficients.csv  – gene coefficients for this fold's model
      plots/
        probabilities.png     – ranked bar chart of P(tolerant)
        roc.png               – ROC curve
    fold2_heldout_study3/   (same structure)
    fold3_heldout_study4/   (same structure)
    metrics_summary.csv     – one row per fold (AUROC / AUPR / Brier)

Usage
-----
  # Default: use the 6-gene panel (final_panel_gene_list.txt)
  python apply_panel_to_new_datasets.py

  # Alternate panel, explicit output directory
  python apply_panel_to_new_datasets.py \\
      --panel ../consensus_final_panel_gene_list.txt \\
      --outdir apply_consensus_panel_results

Command-line arguments
----------------------
  --panel PANEL_PATH
      Path to a plain-text file listing gene IDs (one per line, no
      header).  Resolved relative to this script's directory if not
      absolute.  Default: ../final_panel_gene_list.txt

  --outdir OUTDIR
      Output directory for all results.  Resolved relative to this
      script's directory if not absolute.  Created automatically if
      it does not exist.  Default: apply_panel_results
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")       # non-interactive backend (no display needed)
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    brier_score_loss,
    roc_curve,
)

# ============================================================
# FIXED PATHS  (not user-editable; derived from script location)
# ============================================================

# Directory containing this script
SCRIPT_DIR = Path(__file__).resolve().parent

# Root of the repository (four levels up from this script)
REPO_ROOT = SCRIPT_DIR.parents[3]

# Raw Salmon count matrices for the three test datasets
# (genes × samples; first col = gene_id, second col = gene_name, dropped automatically)
COUNTS_DS2 = REPO_ROOT / "data/rnaseq_gene_counts/salmon.merged.gene_counts_dataset2.tsv"
COUNTS_DS3 = REPO_ROOT / "data/rnaseq_gene_counts/salmon.merged.gene_counts_dataset3.tsv"
COUNTS_DS4 = REPO_ROOT / "data/rnaseq_gene_counts/salmon.merged.gene_counts_dataset4.tsv"

# SRA run tables for each test dataset
# (Experiment column = SRX sample ID matching count-matrix headers; Trait = condition)
META_DS2 = REPO_ROOT / "data/SraRunTable_study2.csv"
META_DS3 = REPO_ROOT / "data/SraRunTable_study3.csv"
META_DS4 = REPO_ROOT / "data/SraRunTable_study4.csv"

# Logistic regression seed (for reproducibility)
RANDOM_SEED = 12345


# ============================================================
# ARGUMENT PARSING
# ============================================================

def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    argparse.Namespace
        Parsed arguments with attributes:
        - panel : Path  – gene list file (one gene ID per line)
        - outdir : Path – output directory
    """
    parser = argparse.ArgumentParser(
        description=(
            "Leave-one-study-out (LOSO) evaluation of a fixed gene panel "
            "across Studies 2, 3, and 4."
        )
    )
    parser.add_argument(
        "--panel",
        type=str,
        default=str(SCRIPT_DIR.parent / "final_panel_gene_list.txt"),
        help=(
            "Path to a plain-text file listing gene IDs, one per line, "
            "no header.  Relative paths are resolved from this script's "
            "directory.  Default: ../final_panel_gene_list.txt"
        ),
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=str(SCRIPT_DIR / "apply_panel_results"),
        help=(
            "Output directory for all results.  Relative paths are "
            "resolved from this script's directory.  Created automatically "
            "if it does not exist.  Default: apply_panel_results"
        ),
    )
    args = parser.parse_args()

    # Resolve paths relative to the script directory so the script can be
    # called from any working directory
    args.panel  = _resolve(args.panel)
    args.outdir = _resolve(args.outdir)
    return args


def _resolve(path_str: str) -> Path:
    """
    Resolve a path string to an absolute Path.

    If the path is already absolute it is returned unchanged; otherwise it
    is resolved relative to SCRIPT_DIR so the script can be called from
    any working directory.

    Parameters
    ----------
    path_str : str
        Raw path string from command-line argument.

    Returns
    -------
    Path
        Absolute Path object.
    """
    p = Path(path_str)
    if not p.is_absolute():
        p = (SCRIPT_DIR / p).resolve()
    return p

# ============================================================
# HELPER FUNCTIONS
# ============================================================

def load_counts(path: Path) -> pd.DataFrame:
    """
    Load a Salmon merged gene count TSV (genes × samples).

    The file has gene_id as the row index (first column) and an optional
    gene_name column (second column) that is dropped.  Returns a DataFrame
    indexed by gene_id with one column per sample (SRX accession).

    Parameters
    ----------
    path : Path
        Path to the count TSV file.

    Returns
    -------
    pd.DataFrame
        Shape: (n_genes, n_samples), indexed by gene_id.
    """
    df = pd.read_csv(path, sep="\t", index_col=0)
    # Drop the gene_name column if present (always the first data column
    # in nf-core/rnaseq Salmon output)
    if "gene_name" in df.columns:
        df = df.drop(columns=["gene_name"])
    return df


def log2_plus1(counts: pd.DataFrame) -> pd.DataFrame:
    """
    Apply log2(count + 1) transformation element-wise.

    This variance-stabilising pseudo-log transform is used in place of
    DESeq2 VST because VST is only available for the training datasets
    (1 and 5) and not for 2, 3, or 4.  Across the 82 training samples
    that have both, log2(count+1) has Pearson r = 0.93–0.97 with VST
    for each of the six panel genes, confirming the approximation is
    adequate.

    Parameters
    ----------
    counts : pd.DataFrame
        Raw count matrix (genes × samples).

    Returns
    -------
    pd.DataFrame
        log2-transformed matrix with the same shape and index.
    """
    return np.log2(counts + 1)


def load_sra_metadata(path: Path, study_label: str) -> pd.DataFrame:
    """
    Load and clean an SRA run-table CSV for one of the test datasets.

    SRA run tables use 'Experiment' for the SRX accession (which matches
    the column headers in the count matrices) and 'Trait' for the
    tolerant/sensitive label.

    Parameters
    ----------
    path : Path
        Path to the SRA run-table CSV.
    study_label : str
        Human-readable label (e.g. "Study 2") used in log messages.

    Returns
    -------
    pd.DataFrame
        Indexed by SRX accession, containing 'condition' (mapped from
        'Trait') and 'batch' (if a batch column is present).
    """
    meta = pd.read_csv(path)
    meta.columns = [c.replace("\ufeff", "").strip() for c in meta.columns]

    if "Experiment" not in meta.columns:
        raise ValueError(f"{study_label}: SRA run table missing 'Experiment' column")
    if "Trait" not in meta.columns:
        raise ValueError(f"{study_label}: SRA run table missing 'Trait' column")

    meta["sample"]    = meta["Experiment"].astype(str).str.strip()
    # Standardise to lowercase to match the encoding used for training
    meta["condition"] = meta["Trait"].astype(str).str.strip().str.lower()

    # Carry a batch column through if present (used for annotation only)
    batch_col = next(
        (c for c in meta.columns if c.lower() == "batch"), None
    )
    meta["batch"] = meta[batch_col].astype(str).str.strip() if batch_col else "unknown"

    return meta.set_index("sample")[["condition", "batch"]]


def align_and_subset(
    expr_log: pd.DataFrame,
    meta: pd.DataFrame,
    panel_genes: list,
    dataset_label: str,
) -> tuple:
    """
    Align expression matrix with metadata and extract panel genes.

    Transposes the expression matrix from (genes × samples) to
    (samples × genes), intersects sample IDs with the metadata index,
    and subsets to the panel genes.

    Parameters
    ----------
    expr_log : pd.DataFrame
        log2(count+1) expression matrix, shape (n_genes, n_samples).
    meta : pd.DataFrame
        Metadata indexed by sample ID.
    panel_genes : list
        List of gene IDs to extract.
    dataset_label : str
        Used in diagnostic messages.

    Returns
    -------
    X : pd.DataFrame
        Expression matrix (n_samples × n_panel_genes), indexed by sample.
    aligned_meta : pd.DataFrame
        Metadata aligned to the same samples as X.
    """
    expr_t = expr_log.T
    expr_t.index.name = "sample"
    expr_t.index = expr_t.index.astype(str).str.strip()

    common = expr_t.index.intersection(meta.index)
    if len(common) == 0:
        raise ValueError(
            f"{dataset_label}: No overlapping sample IDs between expression matrix "
            "and metadata.  Check that 'Experiment' IDs in the SRA run table match "
            "column headers in the count file."
        )

    n_meta_only = len(meta.index.difference(expr_t.index))
    n_expr_only = len(expr_t.index.difference(meta.index))
    if n_meta_only or n_expr_only:
        print(
            f"[INFO] {dataset_label}: {n_meta_only} samples in metadata but not "
            f"in counts; {n_expr_only} in counts but not in metadata – both dropped."
        )

    missing_genes = [g for g in panel_genes if g not in expr_t.columns]
    if missing_genes:
        raise ValueError(f"{dataset_label}: panel genes missing in expression matrix: {missing_genes}")
    X    = expr_t.loc[common, panel_genes]
    meta = meta.loc[common]
    print(
        f"[INFO] {dataset_label}: {len(common)} samples × "
        f"{len(panel_genes)} genes after alignment."
    )
    return X, meta


# ============================================================
# MODELLING FUNCTIONS
# ============================================================

def fit_classifier(X_train: pd.DataFrame, y_train: np.ndarray):
    """
    Fit a StandardScaler and an unpenalised logistic regression classifier.

    The unpenalised model (penalty=None) mirrors the approach used in
    ../final_panel_classifier.py and train_eval.py. With only 6 features the model
    is well-determined without regularisation, and class_weight='balanced' accounts
    for any class imbalance.

    Parameters
    ----------
    X_train : pd.DataFrame
        Training expression matrix (n_train × n_genes).
    y_train : np.ndarray
        Binary labels (1 = tolerant, 0 = sensitive).

    Returns
    -------
    scaler : StandardScaler
        Fitted scaler (μ and σ estimated from training data only).
    clf : LogisticRegression
        Fitted logistic regression model.
    """
    scaler = StandardScaler().fit(X_train.values)
    Xz = scaler.transform(X_train.values)
    clf = LogisticRegression(
        penalty=None,
        solver="lbfgs",
        max_iter=5000,
        class_weight="balanced",   # compensates for class imbalance
        random_state=RANDOM_SEED,
        n_jobs=-1,
    ).fit(Xz, y_train)
    return scaler, clf


def predict(
    scaler: StandardScaler,
    clf: LogisticRegression,
    X_test: pd.DataFrame,
) -> np.ndarray:
    """
    Apply the trained scaler and classifier to a test expression matrix.

    Parameters
    ----------
    scaler : StandardScaler
        Scaler fitted on training data.
    clf : LogisticRegression
        Classifier fitted on training data.
    X_test : pd.DataFrame
        Test expression matrix (n_test × n_genes); must have the same
        gene columns as the training matrix.

    Returns
    -------
    np.ndarray
        Predicted probability of being 'tolerant' (class 1) for each sample.
    """
    Xz = scaler.transform(X_test.values)
    return clf.predict_proba(Xz)[:, 1]


def compute_metrics(y_true: np.ndarray, y_prob: np.ndarray, label: str) -> dict:
    """
    Compute classification metrics for a set of predictions.

    Parameters
    ----------
    y_true : np.ndarray
        Ground-truth binary labels (1 = tolerant, 0 = sensitive).
    y_prob : np.ndarray
        Predicted probability of class 1.
    label : str
        Dataset label for display.

    Returns
    -------
    dict
        Keys: heldout_study, n_samples, n_tolerant, n_sensitive,
              auroc, aupr, brier.
    """
    n     = len(y_true)
    n_tol = int(y_true.sum())
    n_sen = n - n_tol

    # AUROC / AUPR require at least one sample from each class
    if n_tol == 0 or n_sen == 0:
        print(f"[WARN] {label}: only one class present – AUROC/AUPR undefined.")
        return {
            "heldout_study": label,
            "n_samples":     n,
            "n_tolerant":    n_tol,
            "n_sensitive":   n_sen,
            "auroc":         float("nan"),
            "aupr":          float("nan"),
            "brier":         float(brier_score_loss(y_true, y_prob)),
        }

    return {
        "heldout_study": label,
        "n_samples":     n,
        "n_tolerant":    n_tol,
        "n_sensitive":   n_sen,
        "auroc":         float(roc_auc_score(y_true, y_prob)),
        "aupr":          float(average_precision_score(y_true, y_prob)),
        "brier":         float(brier_score_loss(y_true, y_prob)),
    }


# ============================================================
# PLOTTING FUNCTIONS
# ============================================================

def plot_probabilities(
    samples: pd.Index,
    y_prob: np.ndarray,
    y_true: np.ndarray,
    title: str,
    outpath: Path,
):
    """
    Horizontal bar chart of per-sample predicted P(tolerant), sorted
    highest → lowest.  Bars are coloured by the true Trait label and a
    dashed vertical line marks the 0.5 decision boundary.

    Parameters
    ----------
    samples : pd.Index
        Sample identifiers (used as y-axis tick labels).
    y_prob : np.ndarray
        Predicted probability of 'tolerant' for each sample.
    y_true : np.ndarray
        Ground-truth binary labels (1=tolerant, 0=sensitive).
    title : str
        Plot title.
    outpath : Path
        File path for the saved PNG.
    """
    order = np.argsort(y_prob)   # ascending so highest prob is at the top

    fig, ax = plt.subplots(figsize=(7, max(4, len(samples) * 0.25)))
    colors = ["#d73027" if y_true[i] == 0 else "#4575b4" for i in order]
    ax.barh(range(len(order)), y_prob[order], color=colors, edgecolor="none")
    ax.axvline(0.5, color="black", linestyle="--", linewidth=1,
               label="0.5 threshold")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(samples[order], fontsize=7)
    ax.set_xlabel("P(tolerant)")
    ax.set_title(title)
    ax.set_xlim(0, 1)

    legend_elements = [
        Patch(facecolor="#4575b4", label="Trait: tolerant"),
        Patch(facecolor="#d73027", label="Trait: sensitive"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", fontsize=8)

    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"[INFO] Saved {outpath}")


def plot_roc(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    auroc: float,
    title: str,
    outpath: Path,
):
    """
    ROC curve with AUROC annotation.

    Parameters
    ----------
    y_true : np.ndarray
        Ground-truth binary labels.
    y_prob : np.ndarray
        Predicted probability of class 1.
    auroc : float
        Pre-computed AUROC (avoids recomputing; may be NaN if single class).
    title : str
        Plot title.
    outpath : Path
        File path for the saved PNG.
    """
    if np.isnan(auroc):
        print(f"[INFO] Skipping ROC plot for '{title}' (single class present).")
        return

    fpr, tpr, _ = roc_curve(y_true, y_prob)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.plot(fpr, tpr, color="#2166ac", linewidth=2,
            label=f"AUROC = {auroc:.3f}")
    ax.plot([0, 1], [0, 1], "k--", linewidth=1, label="random")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(title)
    ax.legend(loc="lower right", fontsize=9)
    plt.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"[INFO] Saved {outpath}")


# ============================================================
# MAIN
# ============================================================

def main():
    # -------------------------------------------------------
    # 0. Parse arguments and set up output directory
    # -------------------------------------------------------
    args   = parse_args()
    PANEL_PATH = args.panel
    OUTDIR     = args.outdir

    if not PANEL_PATH.exists():
        raise FileNotFoundError(
            f"Panel gene list not found: {PANEL_PATH}\n"
            "Pass --panel with the path to a plain-text file of gene IDs."
        )

    OUTDIR.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Panel file : {PANEL_PATH}")
    print(f"[INFO] Output dir : {OUTDIR}")

    # -------------------------------------------------------
    # 1. Load the fixed gene panel
    # -------------------------------------------------------
    panel_genes = (
        pd.read_csv(PANEL_PATH, header=None)[0]
        .astype(str).str.strip().tolist()
    )
    print(f"[INFO] Panel genes ({len(panel_genes)}): {panel_genes}")

    # -------------------------------------------------------
    # 2. Load and normalise all three datasets once
    # -------------------------------------------------------
    # Each entry: (count_path, metadata_path, short_key, human_label)
    dataset_specs = [
        (COUNTS_DS2, META_DS2, "study2", "Study 2 (Louisiana field)"),
        (COUNTS_DS3, META_DS3, "study3", "Study 3 (Dermo challenge)"),
        (COUNTS_DS4, META_DS4, "study4", "Study 4 (VIMS challenge)"),
    ]

    datasets = {}  # key → {"X": pd.DataFrame, "meta": pd.DataFrame, "y": np.ndarray}

    for counts_path, meta_path, key, label in dataset_specs:
        print(f"\n[INFO] Loading {label} …")
        counts_raw = load_counts(counts_path)
        counts_log = log2_plus1(counts_raw)
        meta       = load_sra_metadata(meta_path, label)
        X, meta_al = align_and_subset(counts_log, meta, panel_genes, label)

        # Encode condition labels
        cond_map = {"tolerant": 1, "sensitive": 0}
        unknown = set(meta_al["condition"].unique()) - set(cond_map.keys())
        if unknown:
            raise ValueError(
                f"{label}: unexpected condition values {unknown}.  "
                "Expected 'tolerant' or 'sensitive' (case-insensitive)."
            )
        y = meta_al["condition"].map(cond_map).astype(int).values

        print(f"[INFO] {label}: {y.sum()} tolerant, {(y == 0).sum()} sensitive")
        datasets[key] = {"X": X, "meta": meta_al, "y": y, "label": label}

    # -------------------------------------------------------
    # 3. LOSO cross-validation across studies
    #
    #    For each held-out study, the other two serve as training data.
    #    The gene panel is FIXED; only the scaler and logistic regression
    #    weights are re-estimated per fold.
    # -------------------------------------------------------
    all_keys   = list(datasets.keys())      # ["study2", "study3", "study4"]
    all_metrics = []

    for fold_idx, test_key in enumerate(all_keys, start=1):
        train_keys  = [k for k in all_keys if k != test_key]
        test_info   = datasets[test_key]
        test_label  = test_info["label"]

        # Per-fold output directory
        fold_dir  = OUTDIR / f"fold{fold_idx}_heldout_{test_key}"
        plots_dir = fold_dir / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)

        # Training labels and description
        train_labels = " + ".join(datasets[k]["label"] for k in train_keys)
        print(f"\n{'='*60}")
        print(f"[INFO] Fold {fold_idx}/3: train on [{train_labels}]")
        print(f"[INFO]            test  on [{test_label}]")

        # Concatenate the two training datasets (rows = samples, cols = genes)
        X_train = pd.concat([datasets[k]["X"] for k in train_keys])
        y_train = np.concatenate([datasets[k]["y"] for k in train_keys])
        print(
            f"[INFO] Training set: {len(X_train)} samples "
            f"({y_train.sum()} tolerant, {(y_train == 0).sum()} sensitive)"
        )

        # Fit scaler + classifier on the two training datasets
        scaler, clf = fit_classifier(X_train, y_train)

        # Save model coefficients for this fold
        coef_df = pd.DataFrame({
            "feature":     panel_genes + ["INTERCEPT"],
            "coefficient": list(clf.coef_.ravel()) + [float(clf.intercept_[0])],
        })
        coef_df.to_csv(fold_dir / "model_coefficients.csv", index=False)
        print(f"[INFO] Coefficients for fold {fold_idx}:")
        print(coef_df.to_string(index=False))

        # Predict on the held-out study
        X_test    = test_info["X"]
        y_true    = test_info["y"]
        meta_test = test_info["meta"]

        y_prob = predict(scaler, clf, X_test)
        y_pred = (y_prob >= 0.5).astype(int)

        # Save per-sample predictions
        preds_df = pd.DataFrame({
            "sample":      X_test.index,
            "condition":   meta_test["condition"].values,
            "y_true":      y_true,
            "p_tolerant":  y_prob,
            "y_pred_0.5":  y_pred,
            "correct":     (y_pred == y_true).astype(int),
        })
        if "batch" in meta_test.columns:
            preds_df.insert(2, "batch", meta_test["batch"].values)
        preds_df.to_csv(fold_dir / "predictions.csv", index=False)
        print(f"[INFO] Predictions saved to {fold_dir / 'predictions.csv'}")

        # Compute and store metrics
        metrics = compute_metrics(y_true, y_prob, test_label)
        metrics["trained_on"] = train_labels
        all_metrics.append(metrics)
        print(
            f"[INFO] Held-out {test_label}: "
            f"AUROC={metrics['auroc']:.3f}  "
            f"AUPR={metrics['aupr']:.3f}  "
            f"Brier={metrics['brier']:.3f}"
        )

        # Save fold metrics as JSON for easy programmatic access
        with open(fold_dir / "metrics.json", "w") as fh:
            json.dump(metrics, fh, indent=2)

        # Plots
        plot_probabilities(
            samples=X_test.index,
            y_prob=y_prob,
            y_true=y_true,
            title=(
                f"Fold {fold_idx}: held-out {test_label}\n"
                f"Predicted P(tolerant) per sample"
            ),
            outpath=plots_dir / "probabilities.png",
        )
        plot_roc(
            y_true=y_true,
            y_prob=y_prob,
            auroc=metrics["auroc"],
            title=f"Fold {fold_idx}: ROC – {test_label}",
            outpath=plots_dir / "roc.png",
        )

    # -------------------------------------------------------
    # 4. Summary metrics table across all three folds
    # -------------------------------------------------------
    print(f"\n{'='*60}")
    metrics_df = pd.DataFrame(all_metrics)[
        ["heldout_study", "trained_on", "n_samples",
         "n_tolerant", "n_sensitive", "auroc", "aupr", "brier"]
    ]
    metrics_path = OUTDIR / "metrics_summary.csv"
    metrics_df.to_csv(metrics_path, index=False)
    print("[INFO] Summary metrics across all LOSO folds:")
    print(metrics_df.to_string(index=False))
    print(f"\n[INFO] All outputs written to {OUTDIR}")


if __name__ == "__main__":
    main()
