# Required packages
library(tidyverse)
library(matrixStats)
library(limma)
library(ggplot2)
library(patchwork)
library(stringr)

setwd("//wsl.localhost/Ubuntu/home/syost/git/Cvirg_Pmarinus_RNAseq/analyses/study4_vs_study4and5")

counts_path <- ("study4and5_genelength_D28_results_all_normalized_counts.tsv")
samplesheet_path <- ("../../data/differential_abundance_sheets/rnaseq_diffabundance_study4and5_D2N_s3_samplesheet.csv")


# ---- Parameters ----
n_features <- 500  # matches nf-core 'exploratory_n_features' default

# ---- Load counts ----
counts <- read_tsv(counts_path, show_col_types = FALSE)
stopifnot("gene_id" %in% colnames(counts))
gene_ids <- counts$gene_id
expr_mat <- counts %>% select(-gene_id) %>% as.matrix()
storage.mode(expr_mat) <- "numeric"
rownames(expr_mat) <- gene_ids

# ---- Load metadata ----
metadata <- read_csv(samplesheet_path, show_col_types = FALSE)
required_cols <- c("sample", "treatment", "batch", "condition")
missing_cols <- setdiff(required_cols, colnames(metadata))
if (length(missing_cols) > 0) {
  stop(sprintf("samplesheet missing columns: %s", paste(missing_cols, collapse=", ")))
}
# Align to counts columns
metadata <- metadata %>%
  filter(sample %in% colnames(expr_mat)) %>%
  arrange(match(sample, colnames(expr_mat)))
stopifnot(all(metadata$sample == colnames(expr_mat)))

# ---- Abbreviate batch labels (Study4_2017 -> S4.7) ----
abbrev_batch <- function(x) {
  out <- str_replace(x, "^Study(\\d+)_[0-9]{3}([0-9])$", "S\\1.\\2")
  out[is.na(out)] <- x[is.na(out)]
  out
}
metadata <- metadata %>%
  mutate(
    batch_label = abbrev_batch(batch),
    treatment = factor(treatment),
    condition = factor(condition),
    batch     = factor(batch)
  )

# ---- Select top-N most variable genes (nf-core exploratory logic) ----
# Note: nf-core typically uses VST/rlog internally; here we use the provided normalized matrix.
vars <- rowVars(expr_mat)
top_idx <- order(vars, decreasing = TRUE)[seq_len(min(n_features, length(vars)))]
expr_top <- expr_mat[top_idx, , drop = FALSE]

# ---- Shape mapping: O = Control, X = Injected (fallback to lexicographic) ----
treat_levels <- levels(metadata$treatment)
shape_map <- if (all(c("Control","Injected") %in% treat_levels)) {
  c("Control" = 1, "Injected" = 4)  # 1 = open circle (O), 4 = X
} else {
  # Fallback: first level -> O, second level -> X
  u <- treat_levels
  setNames(c(1,4, rep(16, max(0,length(u)-2))), u)  # extra levels get filled circles if any
}

# ---- Helpers ----
pca_with_var <- function(mat) {
  p <- prcomp(t(mat), scale. = FALSE)
  imp <- summary(p)$importance
  list(pca = p, var_pc = round(imp[2,1:2]*100, 1))
}
make_plot <- function(pca_list, meta, title) {
  df <- cbind(as.data.frame(pca_list$pca$x), meta)
  ggplot(df, aes(PC1, PC2, shape = treatment, color = condition, label = batch_label)) +
    geom_point(size = 3.5, stroke = 1) +
    geom_text(vjust = -0.9, size = 3) +
    scale_shape_manual(values = shape_map) +
    theme_minimal() +
    labs(
      title = title,
      x = sprintf("PC1 (%.1f%%)", pca_list$var_pc[1]),
      y = sprintf("PC2 (%.1f%%)", pca_list$var_pc[2])
    )
}

# ---- PCA before correction ----
pca_before <- pca_with_var(expr_top)
plot_before <- make_plot(pca_before, metadata, "PCA (Top-500 variance genes) — before batch correction")

# ---- Batch removal (visualization only) ----
expr_top_corrected <- removeBatchEffect(expr_top, batch = metadata$batch)

# ---- PCA after correction ----
pca_after <- pca_with_var(expr_top_corrected)
plot_after <- make_plot(pca_after, metadata, "PCA (Top-500 variance genes) — after batch correction")

# ---- Side-by-side ----
plot_before + plot_after