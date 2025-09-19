library(tidyverse)
library(limma)

counts <- read_csv("counts_long.csv", show_col_types = FALSE)

# Long rRNAs to drop
long_rrna_ids <- c(
  "rRNA_RNA18SN1_166629","rRNA_RNA18SN2_166610",
  "rRNA_RNA28SN1_166631","rRNA_RNA28SN2_166612","rRNA_RNA28SN3_166621",
  "rRNA_RNA28SN4_172372","rRNA_RNA28SN5_178714",
  "rRNA_RNA45SN1_166625","rRNA_RNA45SN2_166607","rRNA_RNA45SN3_166618",
  "rRNA_RNA45SN4_172369","rRNA_RNA45SN5_178711","rRNA_RNR1_192237","rRNA_RNR2_192239"
)

dat0 <- counts %>%
  filter(condition %in% c("Input","IP")) %>%
  filter(!stringr::str_starts(chr_pos, paste0(long_rrna_ids, collapse="|"))) %>%
  mutate(
    condition = factor(condition, levels = c("Input","IP")),
    pair_id   = factor(pair_id)
  )

# Pre-filter: require min coverage per condition
min_depth <- 20L
min_pairs <- 4L

# Helpers
elogit <- function(k, n) log((k + 0.5) / (n - k + 0.5))
w_binom <- function(k, n) {
  p <- (k + 0.5) / (n + 1)
  pmax(n * p * (1 - p), 1e-6)  # avoid zero weights
}

run_super_group <- function(sg) {
  message("Super-group: ", sg)
  df <- filter(dat0, super_group == sg)
  
  # Keep only pairs that have both Input & IP present at least somewhere
  ok_pairs <- df %>%
    count(pair_id, condition) %>%
    tidyr::pivot_wider(names_from = condition, values_from = n, values_fill = 0L) %>%
    filter(Input > 0, IP > 0) %>%
    pull(pair_id)
  df <- filter(df, pair_id %in% ok_pairs)
  if (nrow(df) == 0) { message("  No valid pairs."); return(invisible(NULL)) }
  
  # ---- Aggregate to per-pair & condition ----
  # (filter a bit early to cut size)
  pc <- df %>%
    select(chr_pos, pair_id, condition, K = total_variant, N = depth_raw) %>%
    group_by(chr_pos, pair_id, condition) %>%
    summarize(K = sum(K, na.rm = TRUE), N = sum(N, na.rm = TRUE), .groups = "drop")
  
  # Make one row per (chr_pos, pair_id): K_Input, K_IP, N_Input, N_IP
  wide <- pc %>%
    tidyr::pivot_wider(names_from = condition, values_from = c(K, N), values_fill = 0) %>%
    # ensure columns exist even if some pairs were sparse
    mutate(
      K_Input = pmin(K_Input, N_Input),
      K_IP    = pmin(K_IP,    N_IP),
      ok_depth = (N_Input >= min_depth) & (N_IP >= min_depth)
    ) %>%
    filter(ok_depth)
  
  # Require enough informative pairs per site
  keep_sites <- wide %>% count(chr_pos, name = "n_pairs") %>%
    filter(n_pairs >= min_pairs) %>% pull(chr_pos)
  wide <- filter(wide, chr_pos %in% keep_sites)
  if (nrow(wide) == 0) { message("  No analyzable rows after pair-level QC."); return(invisible(NULL)) }
  
  # ---- Per-pair Δelogit and weights ----
  wide <- wide %>%
    mutate(
      delta = elogit(K_IP, N_IP) - elogit(K_Input, N_Input),
      w     = (w_binom(K_IP, N_IP) + w_binom(K_Input, N_Input)) / 2
    ) %>%
    filter(is.finite(delta), w > 0)
  
  # Build matrices: rows = sites, cols = pairs
  D <- wide %>% select(chr_pos, pair_id, delta) %>%
    tidyr::pivot_wider(names_from = pair_id, values_from = delta) %>%
    tibble::column_to_rownames("chr_pos") %>% as.matrix()
  
  W <- wide %>% select(chr_pos, pair_id, w) %>%
    tidyr::pivot_wider(names_from = pair_id, values_from = w) %>%
    tibble::column_to_rownames("chr_pos") %>% as.matrix()
  
  # Align rows/cols
  common_rows <- intersect(rownames(D), rownames(W))
  D <- D[common_rows, , drop = FALSE]
  W <- W[common_rows, colnames(D), drop = FALSE]
  
  # 1) Drop pair-columns with zero weight everywhere (uninformative pairs)
  col_keep <- colSums(W > 0, na.rm = TRUE) > 0
  D <- D[, col_keep, drop = FALSE]
  W <- W[, col_keep, drop = FALSE]
  
  # 2) Replace any remaining NA deltas with 0 **and** set their weights to 0
  missing <- !is.finite(D)
  if (any(missing)) {
    W[missing] <- 0
    D[missing] <- 0
  }
  
  # 3) Keep rows with at least min_pairs informative pairs
  row_keep <- rowSums(W > 0) >= min_pairs
  D <- D[row_keep, , drop = FALSE]
  W <- W[row_keep, , drop = FALSE]
  if (nrow(D) == 0) { message("  No analyzable rows after QC."); return(invisible(NULL)) }
  
  # 4) Fit with intercept-only design (P × 1) and trend EB
  design <- matrix(1, nrow = ncol(D), ncol = 1)
  
  fit <- limma::lmFit(D, design, weights = W)
  # Ensure eBayes(trend=TRUE) has a finite covariate (no NAs):
  fit$Amean <- rowMeans(D)   # no NA because we just filled them
  fit <- limma::eBayes(fit, trend = TRUE)
  
  tt <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "none") %>%
    tibble::rownames_to_column("chr_pos") %>%
    transmute(
      chr_pos,
      logOR_total = logFC,
      t_total     = t,
      p_total     = P.Value,
      padj_total  = adj.P.Val
    )
  
  write_csv(tt, paste0("stage1_variant_rate_", sg, "_FAST.csv"))
  message("  -> stage1_variant_rate_", sg, "_FAST.csv")
}

# --- run it ---
groups_to_run <- intersect(c("hematopoietic","non-hematopoietic"),
                           unique(dat0$super_group))
invisible(purrr::map(groups_to_run, run_super_group))
message("✅ Stage 1 (FAST, aggregated) complete.")

