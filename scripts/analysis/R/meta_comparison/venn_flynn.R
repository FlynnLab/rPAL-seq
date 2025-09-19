library(readr)
library(dplyr)
library(eulerr)
library(stringr)
library(grid)

# Load the TP files (switch comparison as needed)
bulk <- read_csv("Family_TP_VA.csv")
flynn <- read_csv("/path/to/flynn/unique/family.csv")

# Identify families in each
set_bulk <- unique(bulk$family)
set_flynn <- unique(flynn$family)

# Load full universe
family_map <- read_csv("/path/to/transcriptome/annotation.csv")
all_families <- unique(family_map$family)

# Hypergeometric test
overlap_count <- length(intersect(set_bulk, set_flynn))
m <- length(set_flynn)  # "successes" in population
n <- length(all_families) - m  # "failures" in population
k <- length(set_bulk)  # sample size
pval <- phyper(q = overlap_count - 1, m = m, n = n, k = k, lower.tail = FALSE)

# --------- Euler section ---------
# Euler diagram input
euler_input <- list(
  `rPAL` = set_bulk,
  `ManNAz` = set_flynn
)

# Build Euler diagram
fit <- euler(euler_input, shape = "ellipse", control = list(tol = 1e-6))

# Save to PDF
pdf("euler_AML_overlap_Flynn.pdf", width = 4.2, height = 5.2)

# Plot
plot(
  fit,
  fills = c("#E01A4F", "#7DCFB6"),
  alpha = 0.6,
  edges = FALSE,
  quantities = list(cex = 0.8, col = "grey30"),
  labels = list(font = 2, cex = 0.9),
  main = NULL
)

# Title (bold, top)
grid.text(
  "Shared Hits: rPAL (H9)\n and ManNAz (AML, Flynn2021)",
  x = 0.5, y = 0.94,
  gp = gpar(fontface = "bold", cex = 1.4)
)

# P-value (italic, bottom)
grid.text(
  sprintf("Hypergeometric p = %.2e", pval),
  x = 0.5, y = 0.04,
  gp = gpar(fontface = "italic", cex = 1.05)
)

dev.off()

# ------------------------------
# Overlap summary (rPAL-unique / shared / ManNAz-unique)
# ------------------------------

# Families present in either set
universe_fams <- union(set_bulk, set_flynn)

overlap_df <- data.frame(family = universe_fams, stringsAsFactors = FALSE) %>%
  mutate(
    category = dplyr::case_when(
      family %in% set_bulk  & !(family %in% set_flynn) ~ "rPAL-unique",
      family %in% set_bulk  &  family %in% set_flynn   ~ "shared",
      !(family %in% set_bulk) & family %in% set_flynn  ~ "ManNAz-unique",
      TRUE ~ NA_character_
    )
  ) %>%
  arrange(
    factor(category, levels = c("rPAL-unique", "shared", "ManNAz-unique")),
    family
  )

# Save categorized families
write_csv(overlap_df, "multi_overlap_Flynn_H9.csv")

