This README refers to an earlier iteration of the pipeline.
The scripts here represent the **current recommended approach** for gene panel
selection and cross-dataset classifier application.  (The parent directory
contains a later structural iteration that introduced LOSO holdouts in the
gene-selection step to prevent data leakage; that design change is not needed
for the cross-dataset application described here.)

# Two-Step Gene Expression Classifier

This directory contains scripts for:

1. **Selecting a minimal gene panel** from training data (Studies 1 and 5).
2. **Applying the selected panel** as a binary classifier to new, independent
   datasets (Studies 2, 3, and 4) that were never used during panel selection
   or model training.

The six-gene panel produced by this pipeline is stored in
`../final_panel_gene_list.txt` (one level up, in the parent directory).

---

## Workflow

### Part A – Panel selection (Studies 1 + 5 only)

1. **Run `metaanalysis_rank_tier.R` first**
   - This R script performs a meta-analysis using I² and meta-p and ranks gene
     candidates into tiers.
   - See this[flow diagram](https://github.com/Resilience-Biomarkers-for-Aquaculture/Cvirg_Pmarinus_RNAseq/blob/main/analyses/Study1and5ThreeWay/two_step_gene_expression_classifier/previous_version/metaanalysis_rank_tier_diagram.svg)
   - **Primary output:** `panel_candidates_tier12.txt`
     - Lists candidate genes from tiers 1 and 2.
     - **Note:** Its use in the next step is optional; `lasso_prune_loso_size.py`
       can run without it but will incorporate it if present.

2. **Run `lasso_prune_loso_size.py`**
   - Performs stability selection to rank features (optionally constrained by
     `panel_candidates_tier12.txt`), followed by a panel-size sweep with
     redundancy pruning + LOSO evaluation.  Selects the smallest panel whose
     AUROC is within ΔAUROC ≤ 0.02 of the best.
   - **Primary output:** `../final_panel_gene_list.txt`
     - One gene ID per line.

### Part B – Visualise the panel (optional)

3. **Run `../plot_final_panel.py`** *(optional)*
   - Visualises VST expression of the panel genes across all training samples,
     stratified by condition and batch.
   - **Outputs** (in `plots/`):
     - `heatmap_panel_genes.png` – genes × samples heatmap
     - `annotation_condition.png` / `annotation_batch.png` – annotation bars
     - `annotation_legends.txt` – code → label mapping
     - `violin_<gene>.png` – per-gene violin plots

### Part C – Stability evaluation of the classifier (optional)

4. **Run `../final_panel_classifier.py`** *(optional)*
   - Trains the classifier over many random seeds on Studies 1 + 5 and reports
     coefficient stability.
   - **Primary output:**
     - `repeat_seeds_fixed_model/summary/coefficients_stability.csv`
     - Per-seed metrics and model artefacts in `repeat_seeds_fixed_model/runs/`

### Part D – Apply classifier to new datasets via leave-one-study-out

5. **Run `apply_panel_to_new_datasets.py`**
   - This is the primary cross-dataset validation script.
   - The 6-gene panel is **fixed** (already selected from Studies 1 + 5).
     This script evaluates how well those six genes generalise to three
     completely independent studies using **leave-one-study-out (LOSO)**
     cross-validation:

     | Fold | Train on | Test on |
     |------|----------|---------|
     | 1 | Study 3 + Study 4 | Study 2 (Louisiana field, 40 samples) |
     | 2 | Study 2 + Study 4 | Study 3 (Dermo challenge, 10 samples) |
     | 3 | Study 2 + Study 3 | Study 4 (VIMS challenge, 44 samples) |

     In each fold a fresh StandardScaler and logistic regression are fit on
     the two training studies using the fixed 6-gene panel.

   - **Normalisation note:** The panel was originally selected using
     DESeq2 VST values, which are only available for datasets 1 and 5.
     For this script all three datasets first undergo simple
     library-size normalisation (divide by sample-total size factors,
     scaled to the median library size within each dataset), then
     `log2(normalised_count + 1)` is applied consistently across all
     three datasets.

   - **Outputs** (in `apply_panel_results/`):
     - `fold{1,2,3}_heldout_study{2,3,4}/predictions.csv`
       – per-sample P(tolerant), hard call, true label
     - `fold{1,2,3}_heldout_study{2,3,4}/metrics.json`
       – AUROC, AUPR, Brier for this fold
     - `fold{1,2,3}_heldout_study{2,3,4}/model_coefficients.csv`
       – gene coefficients for the classifier trained in this fold
     - `fold{1,2,3}_heldout_study{2,3,4}/plots/probabilities.png`
       – ranked bar chart of P(tolerant) per sample
     - `fold{1,2,3}_heldout_study{2,3,4}/plots/roc.png`
       – ROC curve
     - `metrics_summary.csv` – one row per fold (AUROC / AUPR / Brier)

   - **Usage:**
     ```bash
     cd analyses/Study1and5ThreeWay/two_step_gene_expression_classifier/previous_version
     python apply_panel_to_new_datasets.py
     ```
      Optional arguments are available; see `python apply_panel_to_new_datasets.py --help`.
      Use `--panel` to supply an alternative gene list and `--outdir` to change the output directory.

---

## File Overview

| File | Purpose |
|------|---------|
| `metaanalysis_rank_tier.R` | Meta-analysis and tier ranking of gene candidates |
| `panel_candidates_tier12.txt` | Output: tier 1+2 candidates from the R script |
| `lasso_prune_loso_size.py` | LASSO stability selection + LOSO panel size sweep |
| `../final_panel_gene_list.txt` | Output: the final 6-gene panel (parent dir) |
| `combine_panel_candidates.py` | Utility: aggregate panel candidate lists across folds |
| `../plot_final_panel.py` | Optional visualisation of panel gene expression |
| `../final_panel_classifier.py` | Optional repeated-seed evaluation of the classifier |
| `train_eval.py` | *Obsolete* – precursor to `lasso_prune_onefold.py` |
| `apply_panel_to_new_datasets.py` | **Apply 6-gene panel to Studies 2, 3, 4 (LOSO)** |
| `apply_panel_results/` | Output directory created by `apply_panel_to_new_datasets.py` |

---

## Usage Notes

- Ensure required input files exist at the paths specified in each script
  before running.  Outputs from earlier steps may be required by later ones.
- `apply_panel_to_new_datasets.py` is self-contained: it reads raw count files
  and SRA run tables directly from `data/` in the repository root.
- The analysis is flexible; `panel_candidates_tier12.txt` is optional for
  `lasso_prune_loso_size.py`.
- See in-script docstrings and comments for full input/output descriptions.

## Citation

If using these scripts in published work, please cite appropriately or give
credit to the repository authors.
