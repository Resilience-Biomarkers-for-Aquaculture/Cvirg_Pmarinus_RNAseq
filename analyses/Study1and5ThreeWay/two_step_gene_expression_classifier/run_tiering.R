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

deseq_per_batch <- function(counts_sxg, meta_s, batch_id, min_row_sum = 10L) {
  message("\n[INFO] ==== deseq_per_batch(", batch_id, ") ====")

  # --- subset metadata and samples ---
  meta_b <- meta_s[meta_s$batch == batch_id, , drop = FALSE]
  if (nrow(meta_b) < 2) stop("Batch ", batch_id, " has <2 samples")
  if (!all(c("tolerant","sensitive") %in% unique(meta_b$condition)))
    stop("Batch ", batch_id, " missing tolerant or sensitive samples")

  smp <- meta_b$sample

  # --- prepare count matrix with gene IDs as rownames ---
  if ("gene_id" %in% colnames(counts_sxg)) {
    rownames(counts_sxg) <- counts_sxg$gene_id
    counts_sxg$gene_id <- NULL
  }
  X <- counts_sxg[, smp, drop = FALSE]

  # keep low-count filter
  keep <- rowSums(X) >= min_row_sum
  X <- X[keep, , drop = FALSE]

  # convert safely to numeric without losing IDs
  rn <- rownames(X)
  X[] <- lapply(X, function(x) as.numeric(as.character(x)))
  rownames(X) <- rn

  # --- DESeq2 analysis ---
  meta_b$condition <- factor(meta_b$condition, levels = c("sensitive","tolerant"))
  dds <- DESeqDataSetFromMatrix(countData = round(X),
                                colData = meta_b,
                                design = ~ condition)
  dds <- DESeq(dds, quiet = TRUE)

  res <- results(dds, name = "condition_tolerant_vs_sensitive")
  out <- tibble::rownames_to_column(as.data.frame(res), var = "gene_id")
  out <- out %>% select(gene_id, log2FoldChange, lfcSE, pvalue, padj, baseMean)
  attr(out, "batch_id") <- batch_id
  out
}

deseq_combined_block <- function(counts_sxg, meta_s, min_row_sum = 10L) {
  message("\n[INFO] ==== deseq_combined_block() ====")

  # --- prepare combined matrix with gene IDs ---
  if ("gene_id" %in% colnames(counts_sxg)) {
    rownames(counts_sxg) <- counts_sxg$gene_id
    counts_sxg$gene_id <- NULL
  }

  # filter low-count genes
  keep <- rowSums(counts_sxg) >= min_row_sum
  X <- counts_sxg[keep, , drop = FALSE]

  # safe numeric coercion
  rn <- rownames(X)
  X[] <- lapply(X, function(x) as.numeric(as.character(x)))
  rownames(X) <- rn

  meta_s$condition <- factor(meta_s$condition, levels = c("sensitive","tolerant"))
  meta_s$batch <- factor(make.names(meta_s$batch))

  dds <- DESeqDataSetFromMatrix(countData = round(X),
                                colData = meta_s,
                                design = ~ batch + condition)
  dds <- DESeq(dds, quiet = TRUE)

  res <- results(dds, name = "condition_tolerant_vs_sensitive")
  out <- tibble::rownames_to_column(as.data.frame(res), var = "gene_id")
  out <- out %>% select(gene_id, log2FoldChange, lfcSE, pvalue, padj, baseMean)
  out
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

if (FALSE) {
  message("[DEBUG] Running in debug mode: downsampling genes to speed up testing")
  n_keep <- 1000  # choose e.g. 500–2000
  all_genes <- counts$gene_id
  keep <- sample(all_genes, n_keep)
  counts <- counts[counts$gene_id %in% keep, ]
}

rownames(counts) <- counts$gene_id
counts$gene_id <- NULL
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

# Harmonize universe to intersection of A and B
genes_common <- intersect(res_A$gene_id, res_B$gene_id)
A2 <- res_A %>% filter(gene_id %in% genes_common) %>% rename_with(~paste0(., ".A"), -gene_id)
B2 <- res_B %>% filter(gene_id %in% genes_common) %>% rename_with(~paste0(., ".B"), -gene_id)

# -------------------------------
# Step 3: Meta-analysis (FE/RE) across A & B
# -------------------------------
# A2 and B2 already exist and both have a character gene_id
tmp2 <- A2 %>%
  inner_join(B2, by = "gene_id", suffix = c(".A", ".B")) %>%
  mutate(gene_id = as.character(gene_id))


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

meta_df <- as.data.frame(meta_df)
meta_df$gene_id <- as.character(meta_df$gene_id)


# Join combined DE results "C" if available 
if (!is.null(res_C)) {
  res_C <- as.data.frame(res_C)
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
# --- PATCH: multiple-testing adjustment of meta p-values ---
meta_df <- meta_df %>%
  mutate(q_meta = p.adjust(p_meta, method = "BH"))


write_csv(meta_df, file.path(out_dir, "meta_ranked.csv"))
message("[INFO] Wrote meta_ranked.csv")

# -------------------------------
# Step 4: Tiering rules (unchanged logic)
# -------------------------------
min_abs_lfc <- 0.5   # effect-size floor
I2_t1 <- 30
I2_t2 <- 40

tier_df <- meta_df %>%
  mutate(
    sig_A = padj_A < 0.05,
    sig_B = padj_B < 0.05,
    sig_C = ifelse(is.na(padj_C), FALSE, padj_C < 0.05),
    sign_ok_AB = same_sign_AB,
    low_hetero = !is.na(I2) & I2 <= I2_t1,
    moderate_hetero = !is.na(I2) & I2 <= I2_t2,
    effect_ok_AB = abs(log2FC_A) >= min_abs_lfc & abs(log2FC_B) >= min_abs_lfc,
    expr_ok = ifelse(is.na(baseMean_A) | is.na(baseMean_B), TRUE,
                     pmax(baseMean_A, baseMean_B, na.rm = TRUE) > 10),
    Tier = case_when(
      # Tier 1: both studies FDR<0.05, same sign, decent effect, low I², and meta q<0.05
      (sig_A & sig_B & sign_ok_AB & effect_ok_AB &
         q_meta < 0.05 & low_hetero & expr_ok) ~ 1L,

      # Tier 2: one study significant, other nominal, same sign, meta q<0.10, moderate I²
      (((sig_A & sign_ok_AB) | (sig_B & sign_ok_AB) | (sig_C & sign_ok_AB)) &
         q_meta < 0.10 & moderate_hetero & effect_ok_AB & expr_ok) ~ 2L,

      # Tier 3: C-only significant, large effect, adequate expression
      (sig_C & !sig_A & !sig_B &
         abs(log2FC_C) >= min_abs_lfc & expr_ok) ~ 3L,

      TRUE ~ NA_integer_
    )
  )

write_csv(tier_df, file.path(out_dir, "tiers.csv"))
message("[INFO] Wrote tiers.csv")

topN <- 2000  # adjust as desiredtopN 
tier12_genes <- tier_df %>%
  filter(Tier %in% c(1L,2L)) %>%
  arrange(Tier, meta_rank) %>%
  slice_head(n = topN) %>%
  mutate(gene_id = as.character(gene_id)) %>%
  pull(gene_id)

writeLines(tier12_genes, opt$panel_out)

message(sprintf(
  "[DEBUG Steps 3/5] Wrote %s (n=%d)\n  Example IDs: %s",
  opt$panel_out,
  length(tier12_genes),
  paste(head(tier12_genes, 5), collapse=", ")
))