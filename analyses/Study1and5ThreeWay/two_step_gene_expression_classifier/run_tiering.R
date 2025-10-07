#!/usr/bin/env Rscript

# =============================================================
# run_tiering.R — LOSO-safe tiering from RAW COUNTS + metadata
#   1) Within-fold DESeq2 per training batch: design ~ condition
#   2) Optional combined DE with blocking: design ~ batch + condition
#   3) FE/RE meta-analysis across the two training batches
#   4) Tiering rules (your original logic)
#   5) Outputs: meta_ranked.csv, tiers.csv, panel_candidates_tier12.txt
# =============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(DESeq2)
  library(metafor)
})

# -------------------------------
# CLI
# -------------------------------
option_list <- list(
  make_option(c("--counts"),        type="character", help="RAW counts matrix (genes x samples); first col = gene_id"),
  make_option(c("--meta"),          type="character", help="Metadata CSV (must contain: sample, condition, batch)"),
  make_option(c("--train_samples"), type="character", help="Text file with training sample IDs (one per line)"),
  make_option(c("--panel_out"),     type="character", help="Output path for Tier1+2 panel (txt)"),
  make_option(c("--seed"),          type="integer", default=12345, help="Random seed")
)
opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)
out_dir <- dirname(opt$panel_out)

if (!file.exists(opt$counts)) stop("--counts not found")
if (!file.exists(opt$meta))   stop("--meta not found")
if (!file.exists(opt$train_samples)) stop("--train_samples not found")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# Helpers
# -------------------------------
guess_delim <- function(path) if (str_detect(tolower(path), "\\.tsv(\\.gz)?$")) "\t" else ","

read_counts <- function(path) {
  df <- read_delim(path, delim = guess_delim(path), show_col_types = FALSE)
  message(sprintf("head of %s:", path))
  print(head(df))

  colnames(df)[1] <- "gene_id"
  as.data.frame(df)
}

read_meta <- function(path) {
  m <- read_csv(path, show_col_types = FALSE)
  names(m) <- tolower(names(m))
  req <- c("sample","condition","batch")
  miss <- setdiff(req, names(m))
  if (length(miss)) stop("Metadata missing columns: ", paste(miss, collapse=", "))
  m$sample    <- as.character(m$sample)
  m$condition <- tolower(as.character(m$condition))
  m$batch     <- as.character(m$batch)
  m
}

ensure_two_batches <- function(batches) {
  u <- unique(batches)
  if (length(u) != 2L) {
    stop(sprintf("Training set must contain exactly TWO batches for A–B meta-analysis; found: %s",
                 paste(u, collapse = ", ")))
  }
  u
}

deseq_per_batch <- function(counts_sxg, meta_s, batch_id,
                            min_row_sum = 10L) {
  message("\n[INFO] ==== deseq_per_batch(", batch_id, ") ====")

  # --- subset samples for this batch ---------------------------------------
  keep_idx <- meta_s$batch == batch_id
  meta_b <- meta_s[keep_idx, , drop = FALSE]
  if (nrow(meta_b) < 2) stop("Batch ", batch_id, " has <2 samples in train")

  # --- check both conditions present ---------------------------------------
  cond_tab <- table(meta_b$condition)
  if (!all(c("tolerant","sensitive") %in% names(cond_tab))) {
    stop("Batch ", batch_id, " missing one of {tolerant,sensitive}")
  }
  message("[DEBUG] n_samples: ", nrow(meta_b),
          " | tolerant=", cond_tab["tolerant"],
          " | sensitive=", cond_tab["sensitive"])

  # --- subset count matrix -------------------------------------------------
  smp <- meta_b$sample
  
  
  
  # --- subset count matrix and ensure gene IDs are rownames ------------------
    if ("gene_id" %in% colnames(counts_sxg)) {
    message("[DEBUG] Using 'gene_id' column as rownames for DESeq2 input.")
    X <- counts_sxg
    rownames(X) <- X$gene_id
    X$gene_id <- NULL
    X <- X[, smp, drop = FALSE]
    } else {
    X <- counts_sxg[, smp, drop = FALSE]
    if (is.null(rownames(X))) {
        stop("[ERROR] No gene_id column or rownames in counts_sxg.")
    }
    }
    message("[DEBUG] Example rownames now: ", paste(head(rownames(X)), collapse = ", "))

  
  
  message("[DEBUG] Subset counts: ", nrow(X), " genes x ", ncol(X), " samples")

  # --- ensure gene IDs are rownames ----------------------------------------
  if (!is.null(counts_sxg$gene_id)) {
    message("[DEBUG] Setting gene_id column as rownames")
    rownames(X) <- counts_sxg$gene_id
  } else if (is.null(rownames(X))) {
    stop("[ERROR] counts_sxg has neither rownames nor gene_id column")
  }
  message("[DEBUG] Example rownames before filtering: ",
          paste(head(rownames(X)), collapse=", "))

  # --- filter low-count genes ----------------------------------------------
  keep_genes <- rowSums(X) >= min_row_sum
  n_keep <- sum(keep_genes)
  message("[DEBUG] Keeping ", n_keep, " / ", nrow(X), " genes with sum >= ", min_row_sum)
  X <- X[keep_genes, , drop = FALSE]

  # --- prepare metadata -----------------------------------------------------
  meta_b$condition <- factor(meta_b$condition, levels = c("sensitive","tolerant"))
  meta_b$batch <- make.names(meta_b$batch)  # safe factor levels
  message("[DEBUG] batch levels: ", paste(levels(meta_b$batch), collapse = ", "))
  message("[DEBUG] condition levels: ", paste(levels(meta_b$condition), collapse = ", "))

  # --- build DESeq2 object --------------------------------------------------
  message("[DEBUG] Creating DESeqDataSet...")
  dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(X)),
                                colData = meta_b,
                                design = ~ condition)
  message("[DEBUG] DESeqDataSet has ", nrow(dds), " genes and ", ncol(dds), " samples")

  # --- run DESeq2 -----------------------------------------------------------
  dds <- estimateSizeFactors(dds)
  dds <- DESeq(dds, quiet = TRUE)
  message("[DEBUG] DESeq run complete")

  # --- get coefficient names ------------------------------------------------
  res_names <- resultsNames(dds)
  message("[DEBUG] Available DESeq2 result names: ", paste(res_names, collapse = ", "))

  # --- extract results for tolerant vs sensitive ----------------------------
  res_name <- grep("^condition_.*tolerant.*", res_names, value = TRUE)
  if (length(res_name) == 0) {
    stop("[ERROR] Could not find condition coefficient in resultsNames(dds)")
  }
  message("[DEBUG] Using coefficient: ", res_name[1])

  res <- results(dds, name = res_name[1])
  message("[DEBUG] results() returned ", nrow(res), " rows")

  # --- convert to data.frame, preserve gene IDs -----------------------------
  out <- as.data.frame(res)
  if (is.null(rownames(out))) {
    stop("[ERROR] DESeq2 result lost rownames — cannot assign gene IDs.")
  }
  out <- tibble::rownames_to_column(out, var = "gene_id")

  # --- diagnostics on gene IDs ---------------------------------------------
  head_ids <- head(out$gene_id, 5)
  message("[DEBUG] head(out$gene_id): ", paste(head_ids, collapse = ", "))
  if (all(grepl("^[0-9]+$", head_ids))) {
    warning("[WARN] gene_id values are numeric; likely lost actual IDs earlier.")
  }

  # --- select standardized columns -----------------------------------------
  expected_cols <- c("gene_id", "log2FoldChange", "lfcSE", "pvalue", "padj", "baseMean")
  missing_cols <- setdiff(expected_cols, colnames(out))
  if (length(missing_cols) > 0) {
    warning("[WARN] Missing columns: ", paste(missing_cols, collapse = ", "))
    out <- out[, intersect(expected_cols, colnames(out)), drop = FALSE]
  } else {
    out <- out[, expected_cols, drop = FALSE]
  }

  # --- final checks ---------------------------------------------------------
  message("[DEBUG] Returning ", nrow(out), " genes with expected columns.")
  message("[DEBUG] Example output rows:")
  print(head(out, 3))

  attr(out, "batch_id") <- batch_id
  message("[INFO] ==== Finished deseq_per_batch(", batch_id, ") ====\n")
  return(out)
}

deseq_combined_block <- function(counts_sxg, meta_s, min_row_sum = 10L) {
  message("[INFO] Running combined DESeq2 (block on batch)")

  # Combined DE across both training batches; block on batch
  smp <- meta_s$sample
  X   <- counts_sxg[, smp, drop = FALSE]
  keep_genes <- rowSums(X) >= min_row_sum
  X <- X[keep_genes, , drop = FALSE]

  meta_c <- meta_s
  meta_c$condition <- factor(meta_c$condition, levels = c("sensitive","tolerant"))
  meta_c$batch     <- factor(make.names(meta_c$batch))
  stopifnot(!anyNA(meta_c$condition), !anyNA(meta_c$batch))

  dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(X)),
                                colData = meta_c,
                                design = ~ batch + condition)
  dds <- estimateSizeFactors(dds)
  dds <- DESeq(dds, quiet = TRUE)

  message("[DEBUG] Available DESeq2 result names:")
  print(resultsNames(dds))

  # Safe retrieval of the condition coefficient
  res_name <- grep("^condition_.*tolerant.*", resultsNames(dds), value = TRUE)
  if (length(res_name) == 0) {
    stop("Could not find coefficient for condition in combined model. Found: ",
         paste(resultsNames(dds), collapse = ", "))
  }
  message("[DEBUG] Using coefficient: ", res_name[1])

  res <- results(dds, name = res_name[1])

  out <- as.data.frame(res)
  if (!"gene_id" %in% colnames(out)) {
    out$gene_id <- rownames(out)
  }

  # Verify expected columns exist before selecting
  expected_cols <- c("gene_id","log2FoldChange","lfcSE","pvalue","padj","baseMean")
  missing_cols <- setdiff(expected_cols, colnames(out))
  if (length(missing_cols)) {
    warning("[WARN] Missing columns in combined DESeq2 results: ",
            paste(missing_cols, collapse=", "),
            " — returning full table instead.")
    return(out)
  }

  # Return selected columns, same as per-batch
  out <- out %>% select(all_of(expected_cols))
  return(out)
}


# -------------------------------
# Load data; restrict to TRAIN
# -------------------------------
train_ids <- scan(opt$train_samples, what = character())
meta <- read_meta(opt$meta) %>% filter(sample %in% train_ids)
if (nrow(meta) == 0) stop("No training samples after filtering with --train_samples")

counts <- read_counts(opt$counts)
if (!"gene_id" %in% names(counts))
  stop("Counts file must have 'gene_id' as first column")

if (!exists("debug_mode")) debug_mode <- TRUE
if (debug_mode) {
  n_keep <- 1000  # choose e.g. 500–2000
  all_genes <- counts$gene_id
  keep <- sample(all_genes, n_keep)
  counts <- counts[counts$gene_id %in% keep, ]
  message(sprintf("[DEBUG] Subsampled %d genes (out of %d)", n_keep, length(all_genes)))
}

rownames(counts) <- counts$gene_id
counts$gene_id <- NULL
print("head of counts after setting rownames:")
print(head(counts))
# Keep only training samples present in counts
common_samp <- intersect(colnames(counts), meta$sample)
if (length(common_samp) == 0)
  stop("No overlap between counts columns and training samples")

counts <- counts[, common_samp, drop = FALSE]
meta   <- meta[match(common_samp, meta$sample), , drop = FALSE]

# Force all columns to numeric and drop any nonfinite entries
rn <- rownames(counts)
counts[] <- lapply(counts, function(x) as.numeric(as.character(x)))
rownames(counts) <- rn
print("head of counts after numeric conversion:")
print(head(counts))

# Verify
stopifnot(all(is.finite(as.matrix(counts))))


message(sprintf("[INFO] Training samples: %d; genes: %d", ncol(counts), nrow(counts)))
batches_train <- ensure_two_batches(meta$batch)
batch_A <- batches_train[1]
batch_B <- batches_train[2]
message(sprintf("[INFO] Training batches: A=%s, B=%s", batch_A, batch_B))

# -------------------------------
# Step 1–2: DE per training batch, plus combined "C"
# -------------------------------
# Resume from previous checkpoint if files exist
if (file.exists(file.path(out_dir, "res_A.rds"))) {
  res_A <- readRDS(file.path(out_dir, "res_A.rds"))
  res_B <- readRDS(file.path(out_dir, "res_B.rds"))
  res_C <- readRDS(file.path(out_dir, "res_C.rds"))
  message("[DEBUG] Loaded cached DESeq2 results; skipping DE computation")
} else {
  # Run DESeq2 normally
  res_A <- deseq_per_batch(counts, meta, batch_A)
  res_B <- deseq_per_batch(counts, meta, batch_B)
  res_C <- tryCatch(deseq_combined_block(counts, meta), error=function(e) NULL)
  saveRDS(res_A, file=file.path(out_dir, "res_A.rds"))
  saveRDS(res_B, file=file.path(out_dir, "res_B.rds"))
  saveRDS(res_C, file=file.path(out_dir, "res_C.rds"))
}
message("head of res_A")
print(head(res_A))
message("head of res_B")
print(head(res_B))


# Harmonize universe to intersection of A and B
genes_common <- intersect(res_A$gene_id, res_B$gene_id)
A2 <- res_A %>% filter(gene_id %in% genes_common) %>% rename_with(~paste0(., ".A"), -gene_id)
B2 <- res_B %>% filter(gene_id %in% genes_common) %>% rename_with(~paste0(., ".B"), -gene_id)
message("head of A2")
print(head(A2))
message("head of B2")
print(head(B2))

# -------------------------------
# Step 3: Meta-analysis (FE/RE) across A & B
# -------------------------------
# A2 and B2 already exist and both have a character gene_id
tmp2 <- A2 %>%
  inner_join(B2, by = "gene_id", suffix = c(".A", ".B")) %>%
  mutate(gene_id = as.character(gene_id))
message("head of tmp2")
print(head(tmp2))


meta_df <- tmp2 %>%
  split(.$gene_id) %>%
  map_dfr(function(.x) {
    yi  <- c(.x$log2FoldChange.A, .x$log2FoldChange.B)
    sei <- c(.x$lfcSE.A,          .x$lfcSE.B)

    fe <- tryCatch(rma.uni(yi = yi, sei = sei, method = "FE"),   error = function(e) NULL)
    re <- tryCatch(rma.uni(yi = yi, sei = sei, method = "REML"), error = function(e) NULL)

    tibble(
      gene_id = as.character(.x$gene_id[1]),
      beta_fe = if (!is.null(fe)) fe$b[1] else NA_real_,
      p_fe    = if (!is.null(fe)) fe$pval else NA_real_,
      I2      = if (!is.null(fe)) fe$I2   else NA_real_,
      beta_re = if (!is.null(re)) re$b[1] else NA_real_,
      p_re    = if (!is.null(re)) re$pval else NA_real_,
      sign_A  = sign(.x$log2FoldChange.A),
      sign_B  = sign(.x$log2FoldChange.B),
      same_sign_AB = sign(.x$log2FoldChange.A) == sign(.x$log2FoldChange.B),
      padj_A = .x$padj.A,
      padj_B = .x$padj.B,
      baseMean_A = .x$baseMean.A,
      baseMean_B = .x$baseMean.B,
      log2FC_A   = .x$log2FoldChange.A,
      log2FC_B   = .x$log2FoldChange.B
    )
  })

# confirm that gene_id survived correctly
stopifnot(is.character(meta_df$gene_id))
message("[CHECK] meta_df gene_id example: ", paste(head(meta_df$gene_id, 5), collapse=", "))

message("head of meta_df at 0")
print(head(meta_df))
meta_df <- as.data.frame(meta_df)
meta_df$gene_id <- as.character(meta_df$gene_id)
message("head of meta_df at 1")
print(head(meta_df))


# Join combined DE results "C" if available 
if (!is.null(res_C)) {
  res_C <- as.data.frame(res_C)
  message("[DEBUG] res_C columns: ", paste(colnames(res_C), collapse = ", "))

  if (all(c("gene_id","padj","log2FoldChange") %in% colnames(res_C))) {
    tmp <- res_C[, c("gene_id","padj","log2FoldChange"), drop = FALSE]
    names(tmp)[names(tmp) == "padj"] <- "padj_C"
    names(tmp)[names(tmp) == "log2FoldChange"] <- "log2FC_C"

    meta_df <- merge(meta_df, tmp, by = "gene_id", all.x = TRUE)
    meta_df$sign_C <- sign(meta_df$log2FC_C)
    meta_df$same_sign_ABC <- meta_df$same_sign_AB &
                             !is.na(meta_df$sign_C) &
                             (meta_df$sign_A == meta_df$sign_C)
  } else {
    warning("[WARN] res_C missing expected columns; skipping C join.")
    meta_df$padj_C <- NA_real_
    meta_df$log2FC_C <- NA_real_
    meta_df$same_sign_ABC <- NA
  }
} else {
  warning("[WARN] res_C is NULL — skipping C join.")
  meta_df$padj_C <- NA_real_
  meta_df$log2FC_C <- NA_real_
  meta_df$same_sign_ABC <- NA
}
message("head of meta_df at 2")
print(head(meta_df))




I2_cut <- 30
meta_df <- meta_df %>%
  mutate(
    p_meta = ifelse(!is.na(I2) & I2 <= I2_cut, p_fe, p_re),
    lambda = 0.02,
    bonus_triplicate = ifelse(!is.na(padj_C) & padj_A < 0.05 & padj_B < 0.05 & padj_C < 0.05, 0.3, 0.0),
    bonus_signmatch = ifelse(isTRUE(same_sign_ABC), 0.2, 0.0),
    score = (-log10(p_meta)) - lambda * ifelse(is.na(I2), 0, I2) + bonus_triplicate + bonus_signmatch
  ) %>%
  arrange(desc(score)) %>%
  mutate(meta_rank = row_number())

message("head of meta_df at 3")
print(head(meta_df))

write_csv(meta_df, file.path(out_dir, "meta_ranked.csv"))
message("[INFO] Wrote meta_ranked.csv")

# -------------------------------
# Step 4: Tiering rules (unchanged logic)
# -------------------------------
tier_df <- meta_df %>%
  mutate(
    sig_A = padj_A < 0.05,
    sig_B = padj_B < 0.05,
    sig_C = ifelse(is.na(padj_C), FALSE, padj_C < 0.05),
    sign_ok_AB = same_sign_AB,
    low_hetero = !is.na(I2) & I2 <= 30,
    moderate_hetero = !is.na(I2) & I2 <= 50,
    effect_large_C = ifelse(is.na(log2FC_C), FALSE, abs(log2FC_C) >= 1),
    expr_ok = ifelse(is.na(baseMean_A) | is.na(baseMean_B), TRUE,
                     pmax(baseMean_A, baseMean_B, na.rm = TRUE) > 10),
    Tier = dplyr::case_when(
      (sig_A & sig_B & sign_ok_AB & (sig_C | (p_meta < 0.01 & low_hetero))) ~ 1L,
      (((sig_C & sign_ok_AB) |
        (sig_A & sign_ok_AB) | (sig_B & sign_ok_AB)) & (p_meta < 0.05) & moderate_hetero) ~ 2L,
      (sig_C & !sig_A & !sig_B & effect_large_C & expr_ok) ~ 3L,
      TRUE ~ NA_integer_
    )
  )

write_csv(tier_df, file.path(out_dir, "tiers.csv"))
message("[INFO] Wrote tiers.csv")

tier12_genes <- tier_df %>%
  filter(Tier %in% c(1L,2L)) %>%
  arrange(Tier, meta_rank) %>%
  mutate(gene_id = as.character(gene_id)) %>%
  pull(gene_id)

writeLines(tier12_genes, opt$panel_out)

message(sprintf(
  "[DEBUG Steps 3/5] Wrote %s (n=%d)\n  Example IDs: %s",
  opt$panel_out,
  length(tier12_genes),
  paste(head(tier12_genes, 5), collapse=", ")
))