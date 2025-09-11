# Two-Step Gene Expression Classifier

This directory contains scripts for conducting a two-step gene expression classifier analysis, primarily focused on panel gene selection. For further exploration, it includes a script for plotting gene counts for the identified panel genes, and a script for building a concrete classifier using those genes. The recommended workflow and dependencies between scripts are described below.

## Workflow

1. **Run `metaanalysis_rank_tier.R` first**
   - This R script performs a meta-analysis using I^2 and meta-p and ranks gene candidates into tiers. 
   - **Primary output:** `panel_candidates_tier12.txt`
     - This file lists candidate genes from tiers 1 and 2.
     - **Note:** The use of `panel_candidates_tier12.txt` in the next step is optional; `lasso_prune_loso_size.py` can run without it, but will utilize it if present.

2. **Run `lasso_prune_loso_size.py`**
   - This Python script does stability selection to rank features (optionally using `panel_candidates_tier12.txt` output by `metaanalysis_rank_tier.R` to constrain candidates), followed by a panel-size sweep with redundancy pruning + LOSO to select gene panels, selecting for smallest m within ΔAUROC ≤ best, producing a final panel & LOSO breakdown.
   - **Primary output:** `final_panel_gene_list.txt`
     - This file contains the list of genes selected for the final classifier.

3. **Optional: Run `plot_final_panel.py`**
   - Visualizes the gene counts for the identified genes in `final_panel_gene_list.txt`, across all samples, with the samples labeled by study and condition.
   - **Outputs:**
     -    Heatmap (genes × samples, sorted by condition then batch):
     `plots/heatmap_panel_genes.png`

     -   Condition annotation bar (aligned to heatmap width):
     `plots/annotation_condition.png`

     -   Batch annotation bar (aligned to heatmap width):
        `plots/annotation_batch.png`

     -   Mapping of annotation codes to labels:
        `plots/annotation_legends.txt`

     -   Violin plots per gene:
        `plots/violin_<gene>.png`

4. **Optional: Run `final_panel_classifier.py`**
   - This is the (optional) final script in the workflow.
   - It builds and evaluates classifiers using the gene panel defined in `final_panel_gene_list.txt` (produced by `lasso_prune_loso_size.py`). It iterates using many starting seeds to produce output models, which can be examined for stability of their parameters over the iterations.
   - **Primary output:**
     - `summary/coefficients_stability.csv`, and several other files in the `repeat_seeds_fixed_model` subdirectory, includig results from individual runs and other summary data.


## File Overview

- `metaanalysis_rank_tier.R`: R script for meta-analysis and tier ranking of gene candidates.
- `panel_candidates_tier12.txt`: Output from the R script listing tier 1 and 2 gene candidates.
- `lasso_prune_loso_size.py`: Python script for LASSO panel pruning, optionally constrained by `panel_candidates_tier12.txt`.
- `final_panel_gene_list.txt`: Output from `lasso_prune_loso_size.py` listing selected gene panel.
- `final_panel_classifier.py`: Python script to train and evaluate a final classifier using the selected panel.

## Usage Notes

- The analysis is flexible; if you wish to use all genes for panel selection, you may skip `panel_candidates_tier12.txt`.
- Ensure that required input files for each script are available in this directory and relative paths before running them.
- Outputs from earlier steps may be required to run subsequent scripts, as outlined above.

## Citation

If using these scripts in published work, please cite appropriately or give credit to the repository authors.
