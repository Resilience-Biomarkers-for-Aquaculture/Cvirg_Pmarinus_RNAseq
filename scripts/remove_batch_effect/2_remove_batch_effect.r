# ==============================
# Batch Effect Correction in R
# Using limma's removeBatchEffect
# ==============================

# 📌 1. Load Required Libraries
library(limma)

# 📌 2. Load Expression Data & Covariates Matrix
expr_matrix <- read.delim("log_transformed_expression.tsv", sep = "\t", row.names = 1, check.names = FALSE)
design_matrix <- read.delim("design_matrix.tsv", sep = "\t", row.names = 1, check.names = FALSE)

# 📌 3. Load Metadata to Extract Batch Variable
metadata <- read.csv("augmented_metadata.csv")

# Convert batch variable (BioProject) to a factor
batch <- as.factor(metadata$BioProject)

# 📌 4. Apply Batch Correction with Biological Covariates
corrected_expr_matrix <- removeBatchEffect(
  as.matrix(expr_matrix), 
  batch = batch, 
  design = as.matrix(design_matrix)
)

# 📌 5. Save Batch-Corrected Expression Data
write.csv(corrected_expr_matrix, "batch_corrected_expression.csv", row.names = TRUE)

# ==============================
# 🎯 Next Steps:
# ✅ Upload "batch_corrected_expression.csv" back
# ✅ Perform PCA to check if batch effects are removed
# ==============================
