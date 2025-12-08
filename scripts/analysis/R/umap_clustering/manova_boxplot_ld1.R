library(MASS)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(tidyverse)

# -----------------------------
# Choose ONE metric to run (baseexpr_vst_cosine or z_delta_cosine)
# -----------------------------
metric <- "baseexpr_vst_cosine"

# -----------------------------
# Resolve input file from metric
#   1) Try exact: "umap_coordinates_<metric>.csv"
#   2) If missing, search any umap_coordinates_*.csv that contains <metric>
# -----------------------------
candidate <- sprintf("umap_coordinates_%s.csv", metric)
if (file.exists(candidate)) {
  infile <- candidate
} else {
  all_umap <- list.files(pattern = "^umap_coordinates_.*\\.csv$")
  hits <- all_umap[grepl(metric, all_umap, fixed = TRUE)]
  if (length(hits) == 1) {
    infile <- hits
    message("Exact file not found; using match: ", infile)
  } else if (length(hits) > 1) {
    stop("Multiple files match metric '", metric, "': ",
         paste(hits, collapse = ", "),
         "\nPlease refine 'metric' or rename files.")
  } else {
    stop("No file found for metric '", metric,
         "'. Expected '", candidate, "' or a file containing the metric string.")
  }
}

# Title from metric
plot_title <- paste0("LD1 from\n", metric)

# -----------------------------
# Load data for the chosen metric
# -----------------------------
df <- read.csv(infile, stringsAsFactors = FALSE)

# -----------------------------
# Binary type: use lineage as-is (Hematopoietic / Non-hematopoietic)
# -----------------------------
df$binary_type <- factor(df$lineage, levels = c("Hematopoietic", "Non-hematopoietic"))

# -----------------------------
# Classical MANOVA (no resampling)
# -----------------------------
fit <- manova(cbind(UMAP1, UMAP2) ~ binary_type, data = df)
pval <- summary(fit, test = "Pillai")$stats["binary_type", "Pr(>F)"]

cat("Classical MANOVA (", metric, ") p-value: ", signif(pval, 4), "\n", sep = "")

# Save results (single row for this metric)
manova_pvals <- data.frame(metric = metric, p_value = pval)
write.csv(manova_pvals, paste0("manova_pvals_", metric, ".csv"), row.names = FALSE)

# -----------------------------
# Plot function: LD1 boxplot with MANOVA subtitle
# -----------------------------
plot_lda_box <- function(df_input, title_text, manova_pval) {
  df <- df_input
  df$binary_type <- factor(df$binary_type, levels = c("Hematopoietic", "Non-hematopoietic"))
  
  # Fit LDA on lineage factor (canonical variates for MANOVA)
  lda_model <- lda(binary_type ~ UMAP1 + UMAP2, data = df)
  df$LD1 <- predict(lda_model)$x[, 1]
  
  # Base plot (MANOVA p is reported; no t-test or whisker)
  p <- ggplot(df, aes(x = binary_type, y = LD1)) +
    geom_boxplot(aes(fill = binary_type), alpha = 0.7, width = 0.5,
                 outlier.shape = NA, color = "black") +
    geom_jitter(aes(fill = binary_type), shape = 19, width = 0.25, size = 1.5) +
    scale_fill_manual(values = c("Hematopoietic" = "#E01A4F",
                                 "Non-hematopoietic" = "#7DCFB6")) +
    labs(
      title = title_text,
      subtitle = paste0("MANOVA p = ", signif(manova_pval, 3)),
      y = "LD1 (LDA)",
      x = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid   = element_blank(),
      axis.line    = element_line(color = "black"),
      axis.ticks   = element_line(color = "black"),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, face = "italic"),
      axis.title.y = element_text(size = 9),
      axis.text.x = element_text(face = "bold", size = 9, angle = 45, hjust = 1)
    )
}

# -----------------------------
# Generate and save plot for the chosen metric
# -----------------------------
plot_single <- plot_lda_box(df, plot_title, pval)

outfile_pdf <- paste0("LD1_plot_", metric, ".pdf")
ggsave(outfile_pdf, plot = plot_single, width = 2, height = 3, units = "in", dpi = 300)

cat("Input file: ", infile, "\n", sep = "")
cat("Saved plot to: ", outfile_pdf, "\n", sep = "")
cat("Saved MANOVA p-value to: ", paste0("manova_pvals_", metric, ".csv"), "\n", sep = "")
