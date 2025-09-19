library(readr)
library(dplyr)
library(eulerr)
library(grid)

# -----------------------------
# HIT files (must have 'family')
# -----------------------------
file_homo_hits  <- "multi_overlap_all_cell_lines.csv"   # Homo (all cell lines)
file_mus_hits   <- "Family_TP_EMT6.csv"                 # Mus (EMT6)
file_canis_hits <- "Family_TP_MDCK_converted.csv"       # Canis (MDCK; converted earlier)

# --------------------------------------------
# REFERENCE files (must have 'family')
# Universe = union(reference families) under loose matching
# Priority for display names: Homo > Mus > Canis
# --------------------------------------------
ref_paths <- c(
  "/path/to/transcriptome/annotation/homo.csv",   # Homo
  "/path/to/transcriptome/annotation/mus.csv",    # Mus
  "/path/to/transcriptome/annotation/canis.csv"   # Canis
)

# --------------------------------------------
# Family harmonization (LOOSE canonicalizer)
# - lowercase, trim, collapse spaces
# - normalize hyphens/commas/underscores
# - robust SNORD/SNORA -> unified canonical forms
#   ("small nucleolar rna, c/d box X", "small nucleolar rna, h/aca box X")
# --------------------------------------------
canon_family <- function(x) {
  s <- tolower(trimws(x))
  s <- gsub("[\u2013\u2014]", "-", s, perl = TRUE)              # normalize en/em-dashes
  s <- gsub(",", " ", s, perl = TRUE)                           # commas -> space
  s <- gsub("_+", " ", s, perl = TRUE)                          # underscores -> space
  s <- gsub("\\s+", " ", s, perl = TRUE)
  
  # SNORD/SNORA to unified canonical keys (very permissive)
  s <- gsub("\\bsnord\\s*[- ]*([[:alnum:]/.-]+)\\b",
            "small nucleolar rna, c/d box \\1", s, perl = TRUE)
  s <- gsub("\\bsnora\\s*[- ]*([[:alnum:]/.-]+)\\b",
            "small nucleolar rna, h/aca box \\1", s, perl = TRUE)
  
  # Normalize already-converted variants
  s <- gsub("\\bsmall\\s+nucleolar\\s+rna\\s*[, ]*c\\s*/?\\s*d\\s*box\\s*[- ]*",
            "small nucleolar rna, c/d box ", s, perl = TRUE)
  s <- gsub("\\bsmall\\s+nucleolar\\s+rna\\s*[, ]*h\\s*/?\\s*aca\\s*box\\s*[- ]*",
            "small nucleolar rna, h/aca box ", s, perl = TRUE)
  
  s <- gsub("\\s+", " ", s, perl = TRUE)
  trimws(s)
}

# Pretty fallback for display when not present in any ref
beautify_display <- function(canon) {
  s <- canon
  s <- sub("^small nucleolar rna, c/d box\\s+", "Small nucleolar RNA, C/D box ", s)
  s <- sub("^small nucleolar rna, h/aca box\\s+", "Small nucleolar RNA, H/ACA box ", s)
  # simple sentence-case fallback for others
  if (!grepl("^Small nucleolar RNA, (C/D|H/ACA) box ", s)) {
    substr(s, 1, 1) <- toupper(substr(s, 1, 1))
  }
  s
}

# Read a reference CSV and return (display, canonical)
read_ref_df <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  if (!("family" %in% names(df))) stop(sprintf("Reference lacks 'family': %s", path))
  fam_display <- trimws(df$family)
  fam_display <- fam_display[!is.na(fam_display)]
  fam_canon   <- canon_family(fam_display)
  unique(data.frame(family_display = fam_display, family_canon = fam_canon, stringsAsFactors = FALSE))
}

# Build display lookup with Homo priority
build_display_lookup <- function(ref_paths) {
  look <- character(0) # named vector: names = canon, value = preferred display
  for (p in ref_paths) {
    r <- read_ref_df(p)
    for (i in seq_len(nrow(r))) {
      k <- r$family_canon[i]; v <- r$family_display[i]
      if (is.na(k) || k == "") next
      if (!nzchar(v)) next
      if (is.na(look[k]) || is.null(look[k])) {
        look[k] <- v  # first seen wins; order of ref_paths sets priority
      }
    }
  }
  look
}

# Read a HIT CSV, keep original + canonical family and (optional) combined_zscore
read_hits_df <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  if (!("family" %in% names(df))) stop(sprintf("Hit file lacks 'family': %s", path))
  df <- df %>%
    mutate(
      family_orig  = family,
      family_canon = canon_family(family_orig)
    )
  if (!("combined_zscore" %in% names(df))) {
    df$combined_zscore <- NA_real_
  } else {
    df$combined_zscore <- suppressWarnings(as.numeric(df$combined_zscore))
  }
  df
}

# -----------------------------
# Read HIT sets (data frames)
# -----------------------------
homo_df  <- read_hits_df(file_homo_hits)
mus_df   <- read_hits_df(file_mus_hits)
canis_df <- read_hits_df(file_canis_hits)

set_homo  <- unique(homo_df$family_canon)
set_mus   <- unique(mus_df$family_canon)
set_canis <- unique(canis_df$family_canon)

# -----------------------------
# Build universe from references (loose match)
# -----------------------------
# Display lookup with priority: Homo > Mus > Canis
display_lookup <- build_display_lookup(ref_paths)

# Universe = union of all refs (canonical)
ref_union <- unique(unlist(lapply(ref_paths, function(p) read_ref_df(p)$family_canon)))
all_families <- ref_union
N <- length(all_families)
universe_note <- "reference-union (loose-matched)"
cat(sprintf("Universe: %s, N = %d\n", universe_note, N))

# Align each HIT set to universe; drop non-testable families
align_to_universe <- function(s, U) {
  kept <- intersect(unique(s), U)
  dropped <- setdiff(unique(s), U)
  list(kept = kept, dropped = dropped)
}
aligned <- list(
  `Homo (all cell lines)` = align_to_universe(set_homo,  all_families),
  `Mus (EMT6)`            = align_to_universe(set_mus,   all_families),
  `Canis (MDCK)`          = align_to_universe(set_canis, all_families)
)
sets <- lapply(aligned, `[[`, "kept")
for (nm in names(aligned)) {
  dct <- length(aligned[[nm]]$dropped)
  if (dct > 0) cat(sprintf("NOTE: %d families dropped from '%s' (not in reference).\n", dct, nm))
}

# -----------------------------
# Exact multivariate hypergeometric
# -----------------------------
exact_intersection_pmf <- function(N, sizes) {
  sizes <- sort(as.integer(sizes))
  pmf <- c(1.0); names(pmf) <- sizes[1]
  for (j in 2:length(sizes)) {
    mj <- sizes[j]
    new_pmf <- numeric(0)
    for (x_name in names(pmf)) {
      x  <- as.integer(x_name); px <- pmf[[x_name]]
      y_min <- max(0, x + mj - N); y_max <- min(x, mj)
      if (y_min > y_max) next
      y_vals <- y_min:y_max
      probs  <- dhyper(y_vals, m = x, n = N - x, k = mj)
      y_keys <- as.character(y_vals)
      cur <- new_pmf[y_keys]; cur[is.na(cur)] <- 0
      new_pmf[y_keys] <- cur + px * probs
    }
    pmf <- new_pmf / sum(new_pmf)
  }
  pmf
}
exact_pval_intersection <- function(N, sizes, k) {
  pmf  <- exact_intersection_pmf(N, sizes)
  supp <- as.integer(names(pmf))
  sum(pmf[supp >= k])
}
expected_overlap <- function(N, sizes) {
  r <- length(sizes)
  prod(sizes) / (N^(r - 1))
}

super_overlap <- Reduce(intersect, sets)
k <- length(super_overlap)
sizes <- vapply(sets, length, integer(1))
pval_exact <- exact_pval_intersection(N, sizes, k)
lambda <- expected_overlap(N, sizes)
enrichment <- ifelse(lambda > 0, k / lambda, Inf)

cat(sprintf("Set sizes: %s\n", paste(names(sizes), sizes, sep = "=", collapse = ", ")))
cat(sprintf("Observed 3-way overlap k = %d; expected λ = %.3f; enrichment = %.3f; exact p = %.3e\n",
            k, lambda, enrichment, pval_exact))

# -----------------------------
# Euler Diagram (3-way)
# -----------------------------
membership_matrix <- sapply(sets, function(x) all_families %in% x)
rownames(membership_matrix) <- all_families
fit <- euler(membership_matrix, shape = "ellipse", control = list(tol = 1e-6))

qty <- fit$original.values
qty_lab <- ifelse(qty == 0, "", as.character(qty))

pdf_file <- "euler_venn_metaspecies.pdf"
pdf(file = pdf_file, width = 5.5, height = 6)

plot(
  fit,
  fills = c("#F9C22E", "#7B66D2", "#B5C8E2"),
  alpha = 0.6,
  edges = FALSE,
  quantities = list(labels = qty_lab, cex = 0.8, col = "grey30"),
  labels = list(font = 2, cex = 0.9),
  main = NULL
)

grid.text(
  "Cross-species TP overlap",
  x = 0.5, y = 0.965, gp = gpar(fontface = "bold", cex = 1.4)
)
grid.text(
  sprintf("Multivariate hypergeometric p = %.2e", pval_exact),
  x = 0.5, y = 0.035, gp = gpar(fontface = "italic", cex = 1.05)
)

dev.off()
cat(sprintf("Wrote %s\n", pdf_file))

# -----------------------------
# Multi-overlap CSV with Homo-style names
#   - 'count' = membership count (1..3)
#   - 'combined_zscore' = mean across all hit rows (by canonical family)
#   - 'family' column uses Homo display when available, else Mus/Canis, else prettified canonical
# -----------------------------
H <- sets[[1]]; M <- sets[[2]]; Cn <- sets[[3]]

# Overlap counts (canonical keys)
overlap_table <- table(unlist(sets))
base_df <- data.frame(
  family_canon = names(overlap_table),
  count        = as.integer(overlap_table),
  stringsAsFactors = FALSE
)

# Combined z-score from all hit rows, grouped by canonical family
all_tp <- dplyr::bind_rows(
  homo_df  %>% select(family_canon, combined_zscore),
  mus_df   %>% select(family_canon, combined_zscore),
  canis_df %>% select(family_canon, combined_zscore)
)
zscores_df <- all_tp %>%
  group_by(family_canon) %>%
  summarise(combined_zscore = mean(combined_zscore, na.rm = TRUE), .groups = "drop") %>%
  mutate(combined_zscore = ifelse(is.nan(combined_zscore), NA_real_, combined_zscore))

# Category labels (by canonical membership)
cat_of <- function(f) {
  c3 <- (f %in% H) + (f %in% M) + (f %in% Cn)
  if (c3 == 3) return("shared_all_three")
  if (c3 == 2) {
    if (f %in% intersect(H, M))  return("Homo∩Mus")
    if (f %in% intersect(H, Cn)) return("Homo∩Canis")
    return("Mus∩Canis")
  }
  if (c3 == 1) {
    if (f %in% H)  return("Homo_only")
    if (f %in% M)  return("Mus_only")
    return("Canis_only")
  }
  return(NA_character_)
}
base_df$category <- vapply(base_df$family_canon, cat_of, character(1))

# Map canonical -> display with Homo priority; fallback to prettified canonical
pick_display <- function(k) {
  if (!is.null(display_lookup[k]) && !is.na(display_lookup[k])) {
    return(unname(display_lookup[k]))
  }
  beautify_display(k)
}

multi_overlap_df <- base_df %>%
  left_join(zscores_df, by = "family_canon") %>%
  mutate(family = vapply(family_canon, pick_display, character(1))) %>%
  select(family, count, category, combined_zscore, family_canon) %>%  # keep canonical as a helper column
  arrange(
    factor(category, levels = c("shared_all_three","Homo∩Mus","Homo∩Canis","Mus∩Canis",
                                "Homo_only","Mus_only","Canis_only")),
    desc(coalesce(combined_zscore, -Inf)),
    family
  )

out_csv <- "multi_overlap_metaspecies.csv"
readr::write_csv(multi_overlap_df, out_csv)
cat(sprintf("Wrote %s (columns: family [Homo-style when available], count, category, combined_zscore, family_canon)\n", out_csv))