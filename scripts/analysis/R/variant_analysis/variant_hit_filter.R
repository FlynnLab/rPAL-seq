library(tidyverse)

# ---------- user knobs ----------
alpha_total <- 0.05
alpha_sgid  <- 0.05

# Effect-size floors (one-sided: enrichment in IP)
min_logOR_total <- log(2.0)
min_logOR_sgid  <- log(1.25)

super_groups <- c("hematopoietic","non-hematopoietic")
# --------------------------------

filter_stage1 <- function(sg) {
  f <- paste0("stage1_variant_rate_", sg, "_FAST.csv")
  if (!file.exists(f)) { message("Skipping (not found): ", f); return(invisible(NULL)) }
  df <- readr::read_csv(f, show_col_types = FALSE)
  
  hits <- df %>%
    filter(!is.na(padj_total)) %>%
    filter(logOR_total > min_logOR_total, padj_total < alpha_total) %>%
    arrange(padj_total)
  
  out <- paste0("hits_total_", sg, ".csv")
  readr::write_csv(hits, out)
  message("Wrote: ", out, " (n = ", nrow(hits), ")")
  invisible(hits)
}

filter_stage2_sgid <- function(sg) {
  f <- paste0("stage2_sgidshare_", sg, "_FAST.csv")
  if (!file.exists(f)) { message("Skipping (not found): ", f); return(invisible(NULL)) }
  df <- readr::read_csv(f, show_col_types = FALSE)
  
  hits <- df %>%
    filter(!is.na(padj_sgidshare)) %>%
    filter(logOR_sgidshare > min_logOR_sgid, padj_sgidshare < alpha_sgid) %>%
    arrange(padj_sgidshare)
  
  out <- paste0("hits_sgid_", sg, ".csv")
  readr::write_csv(hits, out)
  message("Wrote: ", out, " (n = ", nrow(hits), ")")
  invisible(hits)
}

for (sg in super_groups) {
  message("\n=== ", sg, " ===")
  filter_stage1(sg)
  filter_stage2_sgid(sg)
}

message("\n✅ Filtering complete")
