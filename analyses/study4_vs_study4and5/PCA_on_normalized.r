library(tidyverse)
library(limma)

# Read normalized counts table
# This file came from 
counts <- read_tsv("all.normalised_counts.tsv")

# Keep gene IDs separately
gene_ids <- counts$gene_id

# Expression matrix with genes as rows, samples as columns
expr_mat <- counts %>% 
  select(-gene_id) %>% 
  as.data.frame()

# Ensure numeric
expr_mat <- as.matrix(expr_mat)
storage.mode(expr_mat) <- "numeric"

# Load metadata with batch & treatment
metadata <- read_csv("rnaseq_diffabundance_study4and5_D2N_metadata.csv")

# Reorder metadata rows to match expression matrix column order
metadata <- metadata %>% 
  filter(sample_id %in% colnames(expr_mat)) %>% 
  arrange(match(sample_id, colnames(expr_mat)))

# Apply batch correction
batch <- factor(metadata$batch)

expr_corrected <- removeBatchEffect(expr_mat, batch = batch)

# Perform PCA
# Transpose so samples are rows
pca <- prcomp(t(expr_corrected), scale. = FALSE)

# Put into dataframe for plotting
pca_df <- data.frame(pca$x, metadata)

library(ggplot2)
ggplot(pca_df, aes(x = PC1, y = PC2, color = treatment)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA after batch effect removal",
       x = paste0("PC1 (", round(summary(pca)$importance[2,1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca)$importance[2,2] * 100, 1), "%)"))
ggsave("PCA_plot_normalized.pdf", width = 8, height = 6)



