---
layout: post
title: "Lab notebook: 2-stage LOSO classifier pipeline (leakage-aware revision)"
date: 2026-06-09
author: Steve Yost
categories: [lab-notebook, rnaseq, machine-learning]
---

## Context
This note documents how the 2-stage LOSO pipeline was designed and then revised to reduce leakage risk. It is based on:

1. `ChatGPT_Minimal_Biomarker_Longest_Gene_selection_strategy_commented.pdf` (main design thread),
2. `ChatGPT_Minimal_Biomarker_Gene_classifier_explanation_2_commented.pdf` (summary + explicit risk/mitigation framing),
3. current code in this directory, and
4. `previous_version` code representing the pre-pivot implementation.

Primary goal: a **small, reproducible gene panel** and an auditable classifier that generalizes across batches/studies.

## What the initial design proposed
From the main transcript, the intended architecture was:

- **Stage 1 (R):** cross-study ranking/tiering using FE/RE meta-analysis, heterogeneity (`I²`), sign consistency, and tier rules.
- **Stage 2 (Python):** stability selection + redundancy pruning + panel-size sweep + LOSO evaluation.
- select the smallest panel size within `ΔAUROC <= 0.02` of the best-performing size.

This exact structure appears in the implemented scripts:

- `run_tiering.R` (meta + tiering),
- `lasso_prune_onefold.py` (stability selection + pruning + per-fold model fitting),
- `1_run_loso_pipeline.py` (orchestrates fold-wise execution),
- `2_summarize_loso_results.py` (aggregate reporting).

## Why we pivoted
The secondary transcript emphasized key risks (especially leakage and domain imbalance), and that framing drove implementation changes.

### Pre-pivot behavior (`previous_version/lasso_prune_loso_size.py`)
The old script:

- loaded all samples,
- ran stability selection on the full dataset,
- ranked/pruned genes globally,
- then performed LOSO evaluation.

That means held-out groups influenced feature ranking before evaluation. Even with LOSO scoring later, this can inflate apparent generalization.

### Current leakage-aware behavior
The revised pipeline enforces fold isolation:

1. `1_run_loso_pipeline.py` defines train/test samples per held-out batch.
2. It calls `run_tiering.R` with `--train_samples`, so tiering/meta-ranking are fold-specific.
3. It calls `lasso_prune_onefold.py` for that fold only.
4. In `lasso_prune_onefold.py`, stability selection, pruning, and panel-size evaluation use training data; held-out data are only used for fold evaluation outputs.

This converts the process from “global feature discovery + LOSO score” to **nested per-fold feature discovery + fold-held-out testing**, which better matches the intended domain-generalization objective.

## Implemented risk mitigations (from transcript to code)
### 1) Leakage control
- Implemented via fold-specific training sample lists (`train_samples.txt`) and per-fold outputs in `results_08_Dec_2025/loso_*`.
- The Dec-2025 directory name is intentional and reflects when that analysis run was produced, even though this notebook write-up is dated 2026-06-09.
- Tiering and feature ranking are no longer run once on all samples.

### 2) Cross-domain reproducibility emphasis
- `run_tiering.R` computes FE/RE meta quantities and heterogeneity (`I²`), then tier assignments using sign and significance rules.
- Tier 1/2 outputs constrain candidate features passed to modeling (`panel_candidates.txt`).

### 3) Minimal-panel selection discipline
- `lasso_prune_onefold.py` evaluates fixed panel sizes `[6, 8, 10, 12, 16, 24]` and records AUROC/AUPR/Brier.
- Final fold model is refit on the best fold-specific size and exported with coefficients/predictions.

### 4) Probability calibration and reporting
- Calibrated classifiers (`CalibratedClassifierCV`) are used for fold models.
- `2_summarize_loso_results.py` consolidates fold metrics, coefficients, pooled predictions, ROC, and calibration plots.

## Notes on differences vs older docs
`previous_version/README.md` still describes the older single-script workflow (`lasso_prune_loso_size.py`). That documentation is historically useful, but the active pipeline is now the orchestrated per-fold workflow.

## Outcome
The implementation now matches the key methodological intent from both transcripts:

- retain reproducibility-first screening,
- build small interpretable panels,
- and evaluate generalization with stricter fold isolation to reduce leakage risk.

If moved to the website repo, this file is ready to be placed under `_posts` (date-prefixed filename already provided).
