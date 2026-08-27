#!/usr/bin/env Rscript

# =============================================================
# RE_run_tiering.R — LOSO-safe RE-only tiering from RAW COUNTS + metadata
#   1) Within-fold DESeq2 per training batch: design ~ condition
#   2) RE-only meta-analysis across the two training batches (while storing FE/I^2)
#   3) Tiering rules (your original logic)
#   4) Outputs: RE_meta_ranked.csv, RE_tiers.csv
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
panel_out_re <- file.path(out_dir, paste0("RE_", basename(opt$panel_out)))

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

read_batch_cache <- function(path, batch_id, samples) {
  payload <- readRDS(path)
  required_fields <- c("batch_id", "samples", "result")
  if (!is.list(payload) || !all(required_fields %in% names(payload))) {
    return(NULL)
  }

  if (!identical(as.character(payload$batch_id), as.character(batch_id))) {
    return(NULL)
  }

  expected_samples <- sort(as.character(samples))
  cached_samples <- sort(as.character(payload$samples))
  if (!identical(cached_samples, expected_samples)) {
    return(NULL)
  }

  payload$result
}

write_batch_cache <- function(path, result, batch_id, samples) {
  saveRDS(
    list(
      batch_id = as.character(batch_id),
      samples = sort(as.character(samples)),
      result = result
    ),
    file = path
  )
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
# Step 1–2: DE per training batch
# -------------------------------
# Resume from previous checkpoint if files match the current training batches
cache_A_path <- file.path(out_dir, "RE_res_A.rds")
cache_B_path <- file.path(out_dir, "RE_res_B.rds")
samples_A <- meta$sample[meta$batch == batch_A]
samples_B <- meta$sample[meta$batch == batch_B]

res_A <- if (file.exists(cache_A_path)) read_batch_cache(cache_A_path, batch_A, samples_A) else NULL
res_B <- if (file.exists(cache_B_path)) read_batch_cache(cache_B_path, batch_B, samples_B) else NULL

if (!is.null(res_A) && !is.null(res_B)) {
  message("[DEBUG] Loaded cached DESeq2 results matching current training batches; skipping DE computation")
} else {
  if (file.exists(cache_A_path) || file.exists(cache_B_path)) {
    message("[DEBUG] Cached DESeq2 results did not match current training batches; recomputing")
  }

  # Run DESeq2 normally
  res_A <- deseq_per_batch(counts, meta, batch_A)
  res_B <- deseq_per_batch(counts, meta, batch_B)
  write_batch_cache(cache_A_path, res_A, batch_A, samples_A)
  write_batch_cache(cache_B_path, res_B, batch_B, samples_B)
}

# Harmonize universe to intersection of A and B
genes_common <- intersect(res_A$gene_id, res_B$gene_id)
A2 <- res_A %>% filter(gene_id %in% genes_common) %>% rename_with(~paste0(., ".A"), -gene_id)
B2 <- res_B %>% filter(gene_id %in% genes_common) %>% rename_with(~paste0(., ".B"), -gene_id)

# -------------------------------
# Step 3: Meta-analysis (FE/RE) across A & B
# -------------------------------
#   - Compute FE and RE models per gene with metafor::rma.uni
#   - Extract I^2 for heterogeneity and always use RE p-value for ranking
#   - Add reproducibility score: -log10(p_meta)
#   - Track sign-concordance across studies.
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


meta_df <- meta_df %>%
  mutate(
    p_meta = p_re,
    model_used = "RE",
    score = -log10(p_meta)
  ) %>%
  arrange(desc(score)) %>%
  mutate(meta_rank = row_number())
# --- PATCH: multiple-testing adjustment of meta p-values ---
meta_df <- meta_df %>%
  mutate(q_meta = p.adjust(p_meta, method = "BH"))


# Keep full ranking in RE_meta_ranked.csv; interaction-filtering is applied only for tier assignment below.
write_csv(meta_df, file.path(out_dir, "RE_meta_ranked.csv"))
message("[INFO] Wrote RE_meta_ranked.csv")


# -----------------------------------------------
# Optional: detect genes with significant batch×condition interaction
# -----------------------------------------------
# NOTE: In this RE-only version, batch×condition interaction is calculated
# for diagnostic purposes only and is not used as a hard exclusion criterion.
# The random-effects meta-analysis already allows study-specific effect sizes,
# while tiering separately requires same-sign effects and minimum effect size
# in both training studies.
interaction_genes <- character(0)

# Run only if ≥2 batches present and both conditions represented in each
batches_with_both <- tapply(meta$condition, meta$batch,
                            function(x) length(unique(x)) >= 2)
eligible_batches <- names(batches_with_both)[batches_with_both]

if (length(eligible_batches) >= 2) {
  message("[INFO] Performing LRT for batch×condition interaction on training set")

  # Subset counts/metadata to training batches only
  keep_samples <- meta$sample
  count_mat <- as.matrix(counts[, keep_samples, drop = FALSE])
  meta_sub <- meta

  # Design with interaction
  meta_sub$batch <- factor(make.names(meta_sub$batch))
  meta_sub$condition <- factor(meta_sub$condition, levels = c("sensitive","tolerant"))

  dds_full <- DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData   = meta_sub,
    design    = ~ batch + condition + batch:condition
  )

  # Filter low counts
  keep <- rowSums(counts(dds_full)) >= 10
  dds_full <- dds_full[keep, ]

  # Likelihood-ratio test comparing models with vs without the interaction term
  dds_int <- DESeq(dds_full, test = "LRT",
                   reduced = ~ batch + condition, quiet = TRUE)

  res_int <- results(dds_int)
  res_int <- as.data.frame(res_int)
  res_int$gene_id <- rownames(res_int)

  # Genes with nominally significant interaction (condition effect differs by batch)
  interaction_genes <- res_int$gene_id[which(res_int$pvalue < 0.05 & !is.na(res_int$pvalue))]

  message(sprintf("[INFO] Identified %d genes with nominally significant batch×condition interaction.",
                  length(interaction_genes)))
} else {
  message("[INFO] Interaction LRT skipped (insufficient batches with both conditions).")
}


# -------------------------------
# Step 4: Tiering rules (unchanged logic)
# -------------------------------
min_abs_lfc <- 0.5   # effect-size floor

# Do not filter interaction_genes here: in the RE-only workflow,
# interaction is retained as a diagnostic rather than a tiering veto.
tier_df <- meta_df %>%
  mutate(
    sig_A = padj_A < 0.05,
    sig_B = padj_B < 0.05,
    sign_ok_AB = same_sign_AB,
    effect_ok_AB = abs(log2FC_A) >= min_abs_lfc & abs(log2FC_B) >= min_abs_lfc,
    expr_ok = ifelse(is.na(baseMean_A) | is.na(baseMean_B), TRUE,
                     pmax(baseMean_A, baseMean_B, na.rm = TRUE) > 10),
    Tier = case_when(
      # Tier 1: both studies FDR<0.05, same sign, decent effect, and meta q<0.05
      (sig_A & sig_B & sign_ok_AB & effect_ok_AB &
         q_meta < 0.05 & expr_ok) ~ 1L,

      # Tier 2: one study significant, other nominal, same sign, and meta q<0.10
      (((sig_A & sign_ok_AB) | (sig_B & sign_ok_AB)) &
         q_meta < 0.10 & effect_ok_AB & expr_ok) ~ 2L,

      TRUE ~ NA_integer_
    )
  )


write_csv(tier_df, file.path(out_dir, "RE_tiers.csv"))
message("[INFO] Wrote RE_tiers.csv")

topN <- 2000  # adjust as desired
tier12_genes <- tier_df %>%
  filter(Tier %in% c(1L,2L)) %>%
  arrange(Tier, meta_rank) %>%
  slice_head(n = topN) %>%
  mutate(gene_id = as.character(gene_id)) %>%
  pull(gene_id)

writeLines(tier12_genes, panel_out_re)

message(sprintf(
  "[DEBUG Steps 3/5] Wrote %s (n=%d)\n  Example IDs: %s",
  panel_out_re,
  length(tier12_genes),
  paste(head(tier12_genes, 5), collapse=", ")
))