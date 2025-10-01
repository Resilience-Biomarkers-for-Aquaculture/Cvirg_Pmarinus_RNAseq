# ============================================
# PCA + Dendrogram on Salmon gene-level counts
# ============================================
use_voom <- TRUE   # or FALSE

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(limma)  # for voom
  library(edgeR)  # for voom
  library(tibble)
})

# ---- 1) User Inputs ----
six_genes = readLines("//wsl.localhost/Ubuntu/home/syost/git/Cvirg_Pmarinus_RNAseq/analyses/Study1and5ThreeWay/two_step_gene_expression_classifier/final_panel_gene_list.txt")
six_genes <- trimws(six_genes)  # remove any leading/trailing spaces
counts_fp   <- "//wsl.localhost/Ubuntu/home/syost/git/Cvirg_Pmarinus_RNAseq/data/rnaseq_gene_counts/merged_gene_counts.tsv"   # tab-delim; col1=gene_id, others=samples (raw integer counts)
metadata_fp <- "//wsl.localhost/Ubuntu/home/syost/git/Cvirg_Pmarinus_RNAseq/data/augmented_metadata.csv"        # tab-delim; columns: sample, study, condition
out_dir     <- "//wsl.localhost/Ubuntu/home/syost/git/Cvirg_Pmarinus_RNAseq/analyses/six_gene_biomarker/plots"

if (!dir.exists(out_dir)) dir.create(out_dir)

# ---- 2) Load data ----
counts_df <- read_tsv(counts_fp, col_types = cols(.default = col_double(), gene_id = col_character()))
meta      <- read_csv(metadata_fp, col_types = cols(.default = col_character()))

counts_mat <- counts_df %>% column_to_rownames("gene_id") %>% as.matrix()

# Align samples
common <- intersect(colnames(counts_mat), meta$Experiment)
counts_mat <- counts_mat[, common, drop=FALSE]
meta <- meta %>%
  filter(Experiment %in% common) %>%
  distinct(Experiment, .keep_all=TRUE) %>%
  arrange(match(Experiment, common))
stopifnot(identical(colnames(counts_mat), meta$Experiment))

if (use_voom) {
  library(edgeR)
  library(limma)
  dge <- DGEList(counts = counts_mat)
  v   <- voom(dge, plot = TRUE)
  expr_tx <- v$E
  trans_label <- "voom (log2-CPM with precision weights)"
} else {
  expr_tx <- log2(counts_mat + 1)
  trans_label <- "log2(x+1)"
}

# ---- 4) Batch correction ----
design_mat <- model.matrix(~ Trait, data=meta)
expr_bc <- removeBatchEffect(expr_log,
                             batch=meta$Study,
                             design=design_mat)

# ---- 5) PCA on six genes ----
avail_genes <- intersect(six_genes, rownames(expr_bc))
missing <- setdiff(six_genes, avail_genes)
if (length(missing) > 0) {
  cat("[WARN] Missing genes:", paste(missing, collapse=", "), "\n")
}
stopifnot(length(avail_genes) >= 2)

expr_sub <- t(expr_bc[avail_genes, , drop=FALSE])
pca <- prcomp(expr_sub, scale.=TRUE)

var_expl <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 2)
pca_df <- data.frame(
  PC1 = pca$x[,1], PC2 = pca$x[,2],
  study = meta$Study, condition = meta$Trait, sample=meta$Experiment
)

g <- ggplot(pca_df, aes(PC1, PC2, color=condition, shape=study)) +
  geom_point(size=3) +
  labs(
    title = paste("PCA of Six Candidate Genes,", trans_label, " + batch-corrected"),
    x = paste0("PC1 (", var_expl[1], "%)"),
    y = paste0("PC2 (", var_expl[2], "%)")
  ) + theme_minimal()

ggsave(file.path(out_dir,
                 paste0("PCA_six_genes_", ifelse(use_voom,"voom","log2"), ".png")),
       g, width=7, height=6, dpi=300)

pheatmap(expr_z,
         annotation_col = ann_col,
         clustering_method = "ward.D2",
         show_colnames = FALSE,
         main = paste("Six-Gene Heatmap (", trans_label, ", batch-corrected)", sep=""),
         filename = file.path(out_dir,
                              paste0("Heatmap_six_genes_", ifelse(use_voom,"voom","log2"), ".png")),
         width = 8, height = 6)

# Per-study plots
studies <- unique(meta$Study)

for (st in studies) {
  sub_meta <- meta[meta$Study == st, , drop = FALSE]
  
  # Require >= 3 samples to compute a stable PCA
  if (nrow(sub_meta) < 3) {
    message(sprintf("[SKIP] study=%s has < 3 samples; skipping PCA.", st))
    next
  }
  
  # Keep only columns (samples) for this study in the same order as sub_meta
  samp_ids <- sub_meta[[samp_col]]
  sub_mat  <- expr_bc[avail_genes, samp_ids, drop = FALSE]
  
  # If fewer than 2 genes are available post-filter, skip
  if (nrow(sub_mat) < 2) {
    message(sprintf("[SKIP] study=%s has < 2 available genes; skipping PCA.", st))
    next
  }
  
  # PCA on samples (rows) with genes as features; scale columns defensively
  pca_st <- prcomp(t(sub_mat), scale. = TRUE)
  var_expl <- round(100 * (pca_st$sdev^2 / sum(pca_st$sdev^2))[1:2], 2)
  
  # Build plotting data
  pca_df <- data.frame(
    PC1 = pca_st$x[, 1],
    PC2 = pca_st$x[, 2],
    condition = sub_meta$Trait,
    sample_id = sub_meta[[samp_col]],
    stringsAsFactors = FALSE
  )
  
  # Plot PCA
  g_st <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample_id)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
    labs(
      title = sprintf("Per-study PCA (%s) — study: %s", trans_label, st),
      subtitle = sprintf("Six-gene panel, PC1: %s%%  PC2: %s%%", var_expl[1], var_expl[2]),
      x = sprintf("PC1 (%s%%)", var_expl[1]),
      y = sprintf("PC2 (%s%%)", var_expl[2])
    ) +
    theme_minimal(base_size = 12)
  
  fn_suffix <- ifelse(use_voom, "voom", "log2")
  ggplot2::ggsave(
    filename = file.path(out_dir, sprintf("PCA_six_genes_perStudy_%s_%s.png", st, fn_suffix)),
    plot     = g_st, width = 7, height = 6, dpi = 300
  )
  
  # -------- Optional: dendrogram per study --------
  # Distance between samples using the same six genes
  dist_st <- dist(t(sub_mat))
  hc_st   <- hclust(dist_st, method = "ward.D2")
  
  png(file.path(out_dir, sprintf("Dendrogram_six_genes_perStudy_%s_%s.png", st, fn_suffix)),
      width = 900, height = 700)
  plot(hc_st, labels = samp_ids,
       main = sprintf("Per-study dendrogram — study: %s (%s)", st, trans_label))
  dev.off()
}
