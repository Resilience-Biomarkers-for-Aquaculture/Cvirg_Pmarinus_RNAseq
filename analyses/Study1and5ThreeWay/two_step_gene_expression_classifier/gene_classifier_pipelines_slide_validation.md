# Validation of `gene_classifier_pipelines (1).pptx` against repository code

This note validates the attached slide deck claims against the two implementations in this repository:

- **Pipeline 1 (previous/global tiering):** `previous_version/metaanalysis_rank_tier.R` + `previous_version/lasso_prune_loso_size.py`
- **Pipeline 2 (nested LOSO/no-leakage):** `1_run_loso_pipeline.py` + `run_tiering.R` + `lasso_prune_onefold.py` + downstream summaries

## What is correct in the slides

- Pipeline 1 does **global** tiering first, then Python stability/pruning/panel sweep.
- Pipeline 1 uses panel sweep sizes including **6, 8, 10, 12, 16, 24**, and picks smallest panel within a ΔAUROC tolerance (`DELTA_AUROC = 0.02`).
- Pipeline 2 runs tiering **per LOSO fold on training-only samples** (`1_run_loso_pipeline.py` writes fold train lists; `run_tiering.R` consumes `--train_samples`).
- Pipeline 2 uses three holdouts corresponding to batches:
  - `P&S 2023`
  - `P&S 2020 2015`
  - `P&S 2020 2017`

## Required corrections to slide content

1. **Pipeline 2 is not VST-first in the implemented code path**
   - Current LOSO pipeline uses raw counts (`merged_gene_counts.tsv`) for both fold-tiering and per-fold classifier scripts (`1_run_loso_pipeline.py`, `run_tiering.R`, `lasso_prune_onefold.py`).
   - VST-based flow applies to earlier/other scripts (e.g., previous version / fixed-panel scripts), not the strict nested LOSO pipeline entrypoint.

2. **“Equal-weight mean AUROC” is inconsistent with implemented summary**
   - `2_summarize_loso_results.py` computes **weighted means by `n_test`**, not equal-weight fold means.
   - Slide text saying equal-weight averaging prevents large-batch dominance should be revised or the code would need to change.

3. **Slide 3 reported sample sizes and AUROC values are inconsistent with saved results**
   - Slide reports approximately: A n≈61 (AUROC 0.59), B1 n≈4 (AUROC 1.00), B2 n≈2 (AUROC 1.00), weighted mean AUROC ~0.84.
   - Repository result artifacts (`results_08_Dec_2025/loso_*/best_panel_summary_*.json`, `loso_summary_metrics.csv`) show:
     - `P&S 2023`: n_test **61**, AUROC **0.7075**
     - `P&S 2020 2015`: n_test **35**, AUROC **1.0**
     - `P&S 2020 2017`: n_test **27**, AUROC **1.0**
   - The weighted mean AUROC from these files is ~**0.855** (not ~0.84).

4. **“Report log-loss, calibration slope/intercept” is not part of the main nested LOSO outputs**
   - `lasso_prune_onefold.py` writes AUROC/AUPR/Brier + predictions + coefficients.
   - `2_summarize_loso_results.py` can include `log_loss` **if present** in JSON, but upstream fold JSON currently does not include it.
   - Calibration slope/intercept are not emitted by the 1→4 LOSO scripts.

5. **“After all folds: retrain one final pooled frozen model” is not implemented in 1→4**
   - The listed main pipeline scripts (`1_...` through `4_...`) do not perform a single pooled retrain/freeze step.
   - Related fixed-panel experimentation exists in `final_panel_classifier.py`, but that is outside the strict 1→4 LOSO chain described in the issue.

## Practical slide-update recommendation

To match repository code, keep the conceptual two-pipeline story, but update Slide 3 metrics/output wording to the actual files in `results_08_Dec_2025`, and replace equal-weight-mean language with **n_test-weighted** summary language unless code is intentionally changed.
