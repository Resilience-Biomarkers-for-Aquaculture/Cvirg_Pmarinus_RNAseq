#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
LOSO Orchestrator for Oyster Condition Classifier

For each LOSO fold (leave one batch out):
  1. Define train vs test sets.
  2. Run R tiering/stability selection *on training data only*
     → produces fold-specific panel_candidates.txt.
  3. Run Python classifier script on training data using that panel file.
  4. Save all artifacts (genes, seeds, configs, metrics, calibration plots).

Outputs:
  results/
    loso_A/
      panel_candidates.txt
      model_coefficients.csv
      test_metrics_overall.json
      test_metrics_per_batch.csv
      calibration_plot.png
      ...
    loso_B1/...
    loso_B2/...
"""

import subprocess, json
from pathlib import Path
from sys import meta_path
import pandas as pd

# -------------------------------
# Paths & configs
# -------------------------------
META_PATH = "../../../data/differential_abundance_sheets/rnaseq_diffabundance_study1and5_samplesheet_filled_with_study5_alldays.csv"
COUNTS_PATH = "../../../data/rnaseq_gene_counts/merged_gene_counts.tsv"
R_TIERING   = "run_tiering.R"           # your R script path
PY_CLASSIFY = "lasso_prune_onefold.py"           # your classifier script path
OUTDIR      = Path("results")

# Which batches to hold out
BATCHES = ["P&S 2023", "P&S 2020 2015", "P&S 2020 2017"]

# Global seed for reproducibility
GLOBAL_SEED = 12345

# -------------------------------
# Helpers
# -------------------------------
def run_r_tiering(train_samples, fold_dir, seed):
    """
    Run the R tiering script on the training samples only.
    Expects run_tiering.R to take args: counts, meta, samples_out, panel_out, seed.
    """
    panel_out = fold_dir / "panel_candidates.txt"
    train_file = fold_dir / "train_samples.txt"
    with open(train_file, "w") as f:
        f.write("\n".join(train_samples))

    cmd = [
        "Rscript", R_TIERING,
        f"--counts={COUNTS_PATH}",
        f"--meta={META_PATH}",
        f"--train_samples={train_file}",
        f"--panel_out={panel_out}",
        f"--seed={seed}"
    ]
    print("[INFO] Running R tiering:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return panel_out

def run_python_classifier(panel_file, fold_dir, seed, test_batch):
    """
    Run the Python classifier script.
    Expects classifier to take args: counts, meta, panel, outdir, seed, holdout_batch.
    """
    cmd = [
        "python", PY_CLASSIFY,
        f"--counts={COUNTS_PATH}",
        f"--meta={META_PATH}",
        f"--panel={panel_file}",
        f"--outdir={fold_dir}",
        f"--seed={seed}",
        f"--holdout_batch={test_batch}"
    ]
    print("[INFO] Running Python classifier:", " ".join(cmd))
    subprocess.run(cmd, check=True)


# -------------------------------
# Main orchestrator
# -------------------------------
def main():
    OUTDIR.mkdir(exist_ok=True)

    # Load metadata to map samples → batch
    meta = pd.read_csv(META_PATH)
    meta.columns = [c.lower() for c in meta.columns]
    meta = meta.rename(columns={"sample": "sample", "batch": "batch"})
    meta = meta.set_index("sample")

    for i, test_batch in enumerate(BATCHES, 1):
        fold_dir = OUTDIR / f"loso_{test_batch}"
        fold_dir.mkdir(parents=True, exist_ok=True)

        # Partition train/test
        train_samples = meta.index[meta["batch"] != test_batch].tolist()
        test_samples  = meta.index[meta["batch"] == test_batch].tolist()

        print(f"\n[INFO] LOSO fold {i}/{len(BATCHES)}: hold out {test_batch} "
              f"({len(test_samples)} samples test, {len(train_samples)} train)")

        # Run tiering on training only
        panel_file = run_r_tiering(train_samples, fold_dir, GLOBAL_SEED + i)
        # --- Diagnostics: summarize DESeq2 success per fold ---
    
        meta_ranked_path = fold_dir / "meta_ranked.csv"
        if meta_ranked_path.exists():
            meta_ranked = pd.read_csv(meta_ranked_path)
            n_total = meta_ranked.shape[0]
            n_A = meta_ranked["padj_A"].notna().sum()
            n_B = meta_ranked["padj_B"].notna().sum()
            n_C = meta_ranked["padj_C"].notna().sum() if "padj_C" in meta_ranked.columns else 0
            print(f"[DIAG] DESeq2 results per fold (genes): A={n_A}, B={n_B}, C={n_C}, total={n_total}")
        else:
            print(f"[WARN] meta_ranked.csv missing for fold {test_batch}")


        # Run classifier with that panel
        run_python_classifier(panel_file, fold_dir, GLOBAL_SEED + i, test_batch)

        # Save config/metadata log
        cfg = {
            "fold": i,
            "holdout_batch": test_batch,
            "train_samples": len(train_samples),
            "test_samples": len(test_samples),
            "seed": GLOBAL_SEED + i
        }
        with open(fold_dir / "config.json", "w") as f:
            json.dump(cfg, f, indent=2)

    print("\n[INFO] All LOSO folds completed. Results in:", OUTDIR)

if __name__ == "__main__":
    main()
