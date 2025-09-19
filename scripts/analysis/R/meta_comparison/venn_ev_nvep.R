library(readr)
library(dplyr)
library(eulerr)
library(grid)

# ------------------------------
# Set group name and load TP data (example = DiFi, change to MDCK or DiFi_dir as needed)
# ------------------------------
group_name <- "DiFi"

tp_0  <- read_csv(paste0("Family_TP_", group_name, ".csv"))
tp_a <- read_csv(paste0("Family_TP_", group_name, "a.csv"))
tp_b  <- read_csv(paste0("Family_TP_", group_name, "b.csv"))
tp_c <- read_csv(paste0("Family_TP_", group_name, "c.csv"))

# Extract unique families
set_0  <- unique(tp_0$family)
set_a <- unique(tp_a$family)
set_b  <- unique(tp_b$family)
set_c <- unique(tp_c$family)

sets <- list(
  `DiFi, small RNA`        = set_0,
  `DiFi sEV`               = set_a,
  `DiFi\n exomere`         = set_b,
  `DiFi\n supermere`       = set_c
)

# ------------------------------
# Load family universe
# ------------------------------
family_map <- read_csv("/path/to/transcriptome/annotation.csv")
all_families <- unique(family_map$family)

# ------------------------------
# Compute observed 4-way overlap
# ------------------------------
super_overlap <- Reduce(intersect, sets)
k <- length(super_overlap)  # Observed overlap count

# ------------------------------
# Multivariate hypergeometric (exact) p-value
# ------------------------------
exact_intersection_pmf <- function(N, sizes) {
  sizes <- sort(as.integer(sizes))
  pmf <- c(1.0); names(pmf) <- sizes[1]  # X1 is degenerate at m1
  for (j in 2:length(sizes)) {
    mj <- sizes[j]
    new_pmf <- numeric(0)
    for (x_name in names(pmf)) {
      x  <- as.integer(x_name)
      px <- pmf[[x_name]]
      y_min <- max(0, x + mj - N)
      y_max <- min(x, mj)
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

N <- length(all_families)
sizes <- vapply(sets, length, integer(1))
pval_exact <- exact_pval_intersection(N, sizes, k)
lambda <- expected_overlap(N, sizes)              # (not displayed, but kept if needed)
enrichment <- ifelse(lambda > 0, k / lambda, Inf) # (not displayed, but kept if needed)

# ------------------------------
# Euler Diagram
# ------------------------------
membership_matrix <- sapply(sets, function(x) all_families %in% x)
rownames(membership_matrix) <- all_families
fit <- euler(membership_matrix, shape = "ellipse", control = list(tol = 1e-6))

# Build custom quantity labels: hide zeros
qty <- fit$original.values                 # same order eulerr uses
qty_lab <- ifelse(qty == 0, "", as.character(qty))

pdf_file <- paste0("euler_venn_", group_name, ".pdf")
pdf(file = pdf_file, width = 5, height = 6)

plot(
  fit,
  fills = c("#F9C22E", "#D098EE", "#DC5FBD", "#7B66D2"),
  alpha = 0.6,
  edges = FALSE,
  quantities = list(labels = qty_lab, cex = 0.8, col = "grey30"),
  labels = list(font = 2, cex = 0.9),
  main = NULL
)

# Add vertical spacing for title and annotation
grid.text(
  sprintf("TP: DiFi & DiFi derived EV/NV"),
  x = 0.5, y = 0.96,
  gp = gpar(fontface = "bold", cex = 1.4)
)

grid.text(
  sprintf("Multivariate hypergeometric p = %.2e", pval_exact),
  x = 0.5, y = 0.035,  # lower than before
  gp = gpar(fontface = "italic", cex = 1.05)
)

dev.off()

# ------------------------------
# Multi-overlap summary
# ------------------------------
overlap_table <- table(unlist(sets))
multi_overlap_df <- data.frame(
  family = names(overlap_table),
  count = as.integer(overlap_table)
) %>%
  filter(count >= 1)

# Merge z-scores
all_tp <- bind_rows(tp_0, tp_a, tp_b, tp_c)

# Compute Shared column
ev_union <- union(union(set_a, set_b), set_c)

# ------------------------------
# Multi-overlap summary (cell-unique / shared / EV-unique)
# ------------------------------

# Merge z-scores
all_tp <- bind_rows(tp_0, tp_a, tp_b, tp_c)

# EV union (same as you currently define it)
ev_union <- union(union(set_a, set_b), set_c)

# Families that appear in either the cell or any EV
universe_subset <- union(set_0, ev_union)

multi_overlap_df <- data.frame(family = universe_subset, stringsAsFactors = FALSE) %>%
  mutate(
    category = dplyr::case_when(
      family %in% set_0 & !(family %in% ev_union) ~ "cell-unique",
      family %in% set_0 &  family %in% ev_union   ~ "shared",
      !(family %in% set_0) & family %in% ev_union ~ "EV-unique",
      TRUE ~ NA_character_
    ),
    # Optional but handy: how many of the 4 sets each family is in
    count = as.integer((family %in% set_0) +
                         (family %in% set_a) +
                         (family %in% set_b) +
                         (family %in% set_c))
  ) %>%
  left_join(
    all_tp %>%
      group_by(family) %>%
      summarise(combined_zscore = mean(combined_zscore, na.rm = TRUE), .groups = "drop"),
    by = "family"
  ) %>%
  arrange(factor(category, levels = c("cell-unique", "shared", "EV-unique")),
          desc(combined_zscore))

# Save categorized families
write_csv(multi_overlap_df, paste0("multi_overlap_", group_name, ".csv"))

