library(MASS)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(tidyverse)


# ---------------------------------
# Load data
# ---------------------------------
z_df <- read.csv("umap_coordinates_z_LogExp_Weighted.csv", stringsAsFactors = FALSE)
base_df <- read.csv("umap_coordinates_Log_Expression.csv", stringsAsFactors = FALSE)

# Define binary cell type group
z_df$binary_type <- factor(ifelse(z_df$cell_type == "Hematopoietic", "Hematopoietic", "Non-Hematopoietic"))
base_df$binary_type <- factor(ifelse(base_df$cell_type == "Hematopoietic", "Hematopoietic", "Non-Hematopoietic"))

# ---------------------------------
# Classical MANOVA (no resampling)
# ---------------------------------
fit_z <- manova(cbind(UMAP1, UMAP2) ~ binary_type, data = z_df)
pval_z <- summary(fit_z, test = "Pillai")$stats["binary_type", "Pr(>F)"]

fit_base <- manova(cbind(UMAP1, UMAP2) ~ binary_type, data = base_df)
pval_base <- summary(fit_base, test = "Pillai")$stats["binary_type", "Pr(>F)"]

cat("Classical MANOVA (z_LogExp_Weighted) p-value:", signif(pval_z, 4), "\n")
cat("Classical MANOVA (Log_Expression) p-value:", signif(pval_base, 4), "\n")

# Save results
manova_pvals <- data.frame(
  metric = c("z_LogExp_Weighted", "Log_Expression"),
  p_value = c(pval_z, pval_base)
)
write.csv(manova_pvals, "manova_pvals.csv", row.names = FALSE)

plot_lda_box <- function(df_input, title_text, manova_pval) {
  
  df <- df_input
  df$binary_type <- factor(df$binary_type, levels = c("Hematopoietic", "Non-Hematopoietic"))
  
  # Fit LDA
  lda_model <- lda(binary_type ~ UMAP1 + UMAP2, data = df)
  df$LD1 <- predict(lda_model)$x[, 1]
  
  # T-test on LD1
  t_res <- t.test(LD1 ~ binary_type, data = df)
  t_pval <- t_res$p.value
  t_sig <- case_when(
    t_pval > 0.05     ~ "ns",
    t_pval <= 0.0001  ~ "****",
    t_pval <= 0.001   ~ "***",
    t_pval <= 0.01    ~ "**",
    t_pval <= 0.05    ~ "*"
  )
  
  # p-value for whiskers
  comparisons <- data.frame(
    group1 = "Hematopoietic",
    group2 = "Non-Hematopoietic",
    p.adj.signif = t_sig,
    y.position = max(df$LD1, na.rm = TRUE) * 1.2
  )
  
  # LD1 axis label
  coeffs <- lda_model$scaling[, 1]
  ld_formula <- paste0("LD1 = ", signif(coeffs[1], 3), "×UMAP1 + ", signif(coeffs[2], 3), "×UMAP2")
  
  # Plot
  p <- ggplot(df, aes(x = binary_type, y = LD1)) +
    geom_boxplot(aes(fill = binary_type), alpha = 0.7, width = 0.5, outlier.shape = NA, color = "black") +
    geom_jitter(aes(fill = binary_type), 
                 shape = 19, width = 0.25, size = 1.5) +
    scale_fill_manual(values = c("Hematopoietic" = "#E01A4F", "Non-Hematopoietic" = "#7DCFB6")) +
    scale_color_manual(values = c("Hematopoietic" = "#E01A4F", "Non-Hematopoietic" = "#7DCFB6")) +
    stat_pvalue_manual(comparisons, label = "p.adj.signif", y.position = "y.position") +
    labs(
      title = title_text,
      subtitle = paste0("MANOVA p = ", signif(manova_pval, 3), "  |  * = t-test p"),
      y = ld_formula,
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
      axis.text.x = element_text(face = "bold", size = 9),
      )
  
  return(p)
}


# Generate both plots
plot_z <- plot_lda_box(z_df, "LD1 from z_LogExp_Weighted", pval_z)
plot_base <- plot_lda_box(base_df, "LD1 from Log_Expression", pval_base)

# Save plots
ggsave("LD1_plot_z_LogExp_Weighted.pdf", plot = plot_z, width = 2, height = 3, units = "in", dpi = 300)
ggsave("LD1_plot_Log_Expression.pdf", plot = plot_base, width = 2, height = 3, units = "in", dpi = 300)

