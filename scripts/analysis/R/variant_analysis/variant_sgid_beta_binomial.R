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

# Empirical logit + binomial weight helpers
elogit <- function(k, n) log((k + 0.5) / (n - k + 0.5))
w_binom <- function(k, n) {
  p <- (k + 0.5) / (n + 1)
  pmax(n * p * (1 - p), 1e-6)
}

# QC thresholds
min_depth <- 20L
min_pairs <- 2L

run_super_group <- function(sg) {
  message("Super-group: ", sg)
  df <- filter(dat0, super_group == sg)
  
  ## ---------- SGID share among variants ----------
  pc1 <- df %>%
    select(chr_pos, pair_id, condition, SGID, total_variant) %>%
    group_by(chr_pos, pair_id, condition) %>%
    summarize(K = sum(SGID, na.rm = TRUE),
              N = sum(total_variant, na.rm = TRUE),
              .groups = "drop")
  
  # Wide per pair: K/N for Input & IP
  wide1 <- pc1 %>%
    tidyr::pivot_wider(names_from = condition, values_from = c(K, N), values_fill = 0) %>%
    mutate(
      K_Input = pmin(K_Input, N_Input),
      K_IP    = pmin(K_IP,    N_IP),
      ok      = (N_Input > 0 & N_IP > 0)
    ) %>%
    filter(ok)
  
  # Require enough informative pairs per site
  keep_sites1 <- wide1 %>% count(chr_pos, name = "n_pairs") %>%
    filter(n_pairs >= min_pairs) %>% pull(chr_pos)
  wide1 <- filter(wide1, chr_pos %in% keep_sites1)
  
  if (nrow(wide1) > 0) {
    wide1 <- wide1 %>%
      mutate(
        delta = elogit(K_IP, N_IP) - elogit(K_Input, N_Input),
        w     = (w_binom(K_IP, N_IP) + w_binom(K_Input, N_Input)) / 2
      )
    
    # Build D/W
    D1 <- wide1 %>% select(chr_pos, pair_id, delta) %>%
      tidyr::pivot_wider(names_from = pair_id, values_from = delta) %>%
      tibble::column_to_rownames("chr_pos") %>% as.matrix()
    W1 <- wide1 %>% select(chr_pos, pair_id, w) %>%
      tidyr::pivot_wider(names_from = pair_id, values_from = w) %>%
      tibble::column_to_rownames("chr_pos") %>% as.matrix()
    
    # Align & clean
    common_rows <- intersect(rownames(D1), rownames(W1))
    D1 <- D1[common_rows, , drop = FALSE]
    W1 <- W1[common_rows, colnames(D1), drop = FALSE]
    
    # Drop pair-columns with zero weight everywhere
    col_keep <- colSums(W1 > 0, na.rm = TRUE) > 0
    D1 <- D1[, col_keep, drop = FALSE]
    W1 <- W1[, col_keep, drop = FALSE]
    
    # Zero-fill NA deltas and zero their weights
    missing <- !is.finite(D1)
    if (any(missing)) { W1[missing] <- 0; D1[missing] <- 0 }
    
    # Keep rows with enough informative pairs
    row_keep <- rowSums(W1 > 0) >= min_pairs
    D1 <- D1[row_keep, , drop = FALSE]
    W1 <- W1[row_keep, , drop = FALSE]
    
    if (nrow(D1) > 0) {
      design1 <- matrix(1, nrow = ncol(D1), ncol = 1)
      fit1 <- limma::lmFit(D1, design1, weights = W1)
      fit1$Amean <- rowMeans(D1)  # ensure finite covariate for trend EB
      fit1 <- limma::eBayes(fit1, trend = TRUE)
      
      tt1 <- limma::topTable(fit1, coef = 1, number = Inf, sort.by = "none") %>%
        tibble::rownames_to_column("chr_pos") %>%
        transmute(chr_pos,
                  logOR_sgidshare = logFC,
                  t_sgidshare     = t,
                  p_sgidshare     = P.Value,
                  padj_sgidshare  = adj.P.Val)
      readr::write_csv(tt1, paste0("stage2_sgidshare_", sg, "_FAST.csv"))
      message("  -> stage2_sgidshare_", sg, "_FAST.csv")
    } else {
      message("  SGID: no analyzable rows after pair-level QC.")
    }
  } else {
    message("  SGID: no sites pass pair-level coverage.")
  }
}

# --- run it ---
groups_to_run <- intersect(c("hematopoietic","non-hematopoietic"),
                           unique(dat0$super_group))
invisible(purrr::map(groups_to_run, run_super_group))
message("✅ Stage 2 (FAST, SGID-only) complete.")
