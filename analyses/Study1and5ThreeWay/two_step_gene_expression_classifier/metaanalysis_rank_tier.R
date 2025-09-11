#!/usr/bin/env Rscript

# ===============================
# Purpose:
#   Implements:
#     Step (1) Harmonize universe via intersection of A and B post-filter genes
#     Step (2) Meta-analysis across A & B (FE/RE, I^2), reproducibility score
#     Step (3) Tiering rules with sign-concordance and C support
#   Also preps outputs for Python modeling (Steps 4–5).
# ===============================

suppressPackageStartupMessages({
  library(dplyr)    # tidy data ops
  library(readr)    # robust CSV/TSV reading
  library(metafor)  # Step 2: meta-analysis (FE/RE & I^2)
  library(stringr)  # small helpers (delimiter guess)
})

# -------------------------------
# [User-editable paths]
#   Replace with your actual files from nf-core/differentialabundance exports.
# -------------------------------
A_path <- "Study1_D7_treatment_resistant_sensitive.deseq2.results.tsv"   # Study A DE table
B_path <- "Study5_D7_treatment_resistant_sensitive.deseq2.results.tsv"   # Study B DE table
C_path <- "Issue44_study1and5_D7_treatment_resistant_sensitive_block.deseq2.results.tsv"   # Combined A+B DE table (optional for bonuses/tiering)

# -------------------------------
# [Step 1] Load results & define common universe (A ∩ B)
#   - Upstream filtering harmonization is assumed done in your pipeline.
#   - Here we standardize column names and types; then intersect gene sets.
# -------------------------------

guess_delim <- function(path) {
  # Tiny helper: guess TSV vs CSV from extension (works for .gz too).
  if (str_detect(tolower(path), "\\.tsv(\\.gz)?$")) "\t" else ","
}

read_res <- function(path) {
  # Read (auto-delim), then coerce to a standard schema used below.
  df <- read_delim(path, delim = guess_delim(path), show_col_types = FALSE)
  df %>%
    rename(
      gene_id        = any_of(c("gene_id","Gene","gene","id")),
      log2FoldChange = any_of(c("log2FoldChange","logFC","log2FC")),
      lfcSE          = any_of(c("lfcSE","SE","se")),
      pvalue         = any_of(c("pvalue","pval","PValue")),
      padj           = any_of(c("padj","FDR","qvalue","qval")),
      baseMean       = any_of(c("baseMean","AveExpr","base_mean"))
    ) %>%
    mutate(
      log2FoldChange = as.numeric(log2FoldChange),
      lfcSE          = as.numeric(lfcSE),
      pvalue         = as.numeric(pvalue),
      padj           = as.numeric(padj),
      baseMean       = as.numeric(baseMean)
    )
}

A <- read_res(A_path)
B <- read_res(B_path)
C <- tryCatch(read_res(C_path), error = function(e) NULL)

message(sprintf("[DEBUG Step 1] Loaded rows: A=%d, B=%d, C=%s",
                nrow(A), nrow(B), ifelse(is.null(C), "None", nrow(C))))

# Intersection of gene IDs = our modeling/meta-analysis universe (prevents optimism).
genes_common <- intersect(A$gene_id, B$gene_id)
message(sprintf("[DEBUG Step 1] Common universe size (A ∩ B): %d", length(genes_common)))

# Keep aligned, core columns only (ensures consistent inputs for Step 2).
A2 <- A %>% filter(gene_id %in% genes_common) %>%
  select(gene_id, log2FoldChange, lfcSE, pvalue, padj, baseMean)
B2 <- B %>% filter(gene_id %in% genes_common) %>%
  select(gene_id, log2FoldChange, lfcSE, pvalue, padj, baseMean)

# -------------------------------
# [Step 2] Meta-analysis across A & B
#   - Compute FE and RE models per gene with metafor::rma.uni
#   - Extract I^2 for heterogeneity; choose p_meta = FE if I^2≤30 else RE
#   - Add reproducibility score: -log10(p_meta) - λ*I^2 + small bonuses
#   - Track sign-concordance across studies.
# -------------------------------

meta_df <- A2 %>%
  inner_join(B2, by = "gene_id", suffix = c(".A",".B")) %>%
  rowwise() %>%
  do({
    # Pull per-study effect sizes (log2FC) and standard errors (lfcSE).
    yi  <- c(.$log2FoldChange.A, .$log2FoldChange.B)
    sei <- c(.$lfcSE.A,          .$lfcSE.B)

    # Fit Fixed-Effect (FE) and Random-Effects (REML) models.
    fe <- tryCatch(rma.uni(yi = yi, sei = sei, method = "FE"),    error = function(e) NULL)
    re <- tryCatch(rma.uni(yi = yi, sei = sei, method = "REML"),  error = function(e) NULL)

    # Safely extract estimates/p-values; NA on failures.
    data.frame(
      gene_id = .$gene_id,
      beta_fe = if (!is.null(fe)) fe$b[1] else NA_real_,
      p_fe    = if (!is.null(fe)) fe$pval  else NA_real_,
      I2      = if (!is.null(fe)) fe$I2    else NA_real_,
      beta_re = if (!is.null(re)) re$b[1] else NA_real_,
      p_re    = if (!is.null(re)) re$pval else NA_real_,
      # Sign consistency flags for A vs B (directional concordance requirement).
      sign_A  = sign(.$log2FoldChange.A),
      sign_B  = sign(.$log2FoldChange.B),
      same_sign_AB = sign(.$log2FoldChange.A) == sign(.$log2FoldChange.B),
      # Carry per-study stats for tiering decisions.
      padj_A  = .$padj.A,
      padj_B  = .$padj.B,
      baseMean_A = .$baseMean.A,
      baseMean_B = .$baseMean.B,
      log2FC_A   = .$log2FoldChange.A,
      log2FC_B   = .$log2FoldChange.B
    )
  }) %>% ungroup()

# Merge C run to enable triplicate/consistency bonuses and tier logic.
if (!is.null(C)) {
  meta_df <- meta_df %>%
    left_join(C %>% select(gene_id, padj, log2FoldChange) %>%
                rename(padj_C = padj, log2FC_C = log2FoldChange),
              by = "gene_id") %>%
    mutate(sign_C = sign(log2FC_C),
           same_sign_ABC = same_sign_AB & !is.na(sign_C) & (sign_A == sign_C))
} else {
  meta_df <- meta_df %>% mutate(padj_C = NA_real_, log2FC_C = NA_real_, same_sign_ABC = NA)
}

# Choose meta p-value by heterogeneity; compute reproducibility score.
I2_cut <- 30  # Step 2 threshold; tune if desired
meta_df <- meta_df %>%
  mutate(
    p_meta = ifelse(!is.na(I2) & I2 <= I2_cut, p_fe, p_re),
    lambda = 0.02,  # heterogeneity penalty weight proposed
    bonus_triplicate = ifelse(!is.na(padj_C) & padj_A < 0.05 & padj_B < 0.05 & padj_C < 0.05, 0.3, 0.0),
    bonus_signmatch = ifelse(isTRUE(same_sign_ABC), 0.2, 0.0),
    score = (-log10(p_meta)) - lambda * ifelse(is.na(I2), 0, I2) + bonus_triplicate + bonus_signmatch
  ) %>%
  arrange(desc(score)) %>%
  mutate(meta_rank = row_number())

write_csv(meta_df, "meta_ranked.csv")
message("[DEBUG Step 2] Wrote meta_ranked.csv")

# -------------------------------
# [Step 3] Tiering rules
#   Tier 1: A & B FDR<0.05 + same sign, and [C FDR<0.05 OR (p_meta<0.01 & I^2≤30)]
#   Tier 2: (C FDR<0.05 & A,B nominal p<0.10 & same sign) OR
#           (one-study FDR<0.05, other nominal, p_meta<0.05, I^2≤50)
#   Tier 3: C-only FDR<0.05 with |log2FC_C|≥1 and adequate expression
#   (Plus a light expression floor; adjust to platform.)
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
    Tier = case_when(
      (sig_A & sig_B & sign_ok_AB & (sig_C | (p_meta < 0.01 & low_hetero))) ~ 1L,
      (((sig_C & sign_ok_AB) |
        (sig_A & sign_ok_AB) | (sig_B & sign_ok_AB)) & (p_meta < 0.05) & moderate_hetero) ~ 2L,
      (sig_C & !sig_A & !sig_B & effect_large_C & expr_ok) ~ 3L,
      TRUE ~ NA_integer_
    )
  )

write_csv(tier_df, "tiers.csv")
message("[DEBUG Step 3] Wrote tiers.csv")

# Export Tier1+2 gene list for modeling (Step 5 prefilter).
tier12_genes <- tier_df %>% filter(Tier %in% c(1L,2L)) %>%
  arrange(Tier, meta_rank) %>% pull(gene_id)
write_lines(tier12_genes, "panel_candidates_tier12.txt")
message(sprintf("[DEBUG Steps 3/5] Wrote panel_candidates_tier12.txt (n=%d)", length(tier12_genes)))

# -------------------------------
# Notes:
# - Step 6: We used C only for bonuses/tiering; the universe remains A∩B.
# - Step 7: Thresholds (I^2 cutoffs, λ, effect size, baseMean) are tuneable.
# ===============================
