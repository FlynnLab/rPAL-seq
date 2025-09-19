library(readr)
library(dplyr)
library(eulerr)
library(stringr)
library(grid)

# Load the TP files (example = HeLa, change paths as needed)
bulk <- read_csv("Family_TP_VH.csv")
lowinput <- read_csv("Family_TP_dVH.csv")

# Identify families in each
set_bulk <- unique(bulk$family)
set_lowinput <- unique(lowinput$family)

# Load full universe
family_map <- read_csv("/path/to/transcriptome/annotation.csv")
all_families <- unique(family_map$family)

# Hypergeometric test
overlap_count <- length(intersect(set_bulk, set_lowinput))
m <- length(set_lowinput)  # "successes" in population
n <- length(all_families) - m  # "failures" in population
k <- length(set_bulk)  # sample size (HeLa set size)
pval <- phyper(q = overlap_count - 1, m = m, n = n, k = k, lower.tail = FALSE)

# --------- Euler section ---------
# Euler diagram input
euler_input <- list(
  `bulk` = set_bulk,
  `low-\ninput` = set_lowinput
)

# Build Euler diagram
fit <- euler(euler_input, shape = "ellipse", control = list(tol = 1e-6))

# Save to PDF
pdf("euler_downsample_overlap_HeLa.pdf", width = 3.2, height = 4.2)

# Plot
plot(
  fit,
  fills = c("#003D5B", "#7DCFB6"),
  alpha = 0.6,
  edges = FALSE,
  quantities = list(cex = 0.8, col = "grey30"),
  labels = list(font = 2, cex = 0.9),
  main = NULL
)

# Title (bold, top)
grid.text(
  "TP Retention Across Scales: HeLa",
  x = 0.5, y = 0.94,
  gp = gpar(fontface = "bold", cex = 1.1)
)

# P-value (italic, bottom)
grid.text(
  sprintf("Hypergeometric p = %.2e", pval),
  x = 0.5, y = 0.04,
  gp = gpar(fontface = "italic", cex = 0.95)
)

dev.off()
