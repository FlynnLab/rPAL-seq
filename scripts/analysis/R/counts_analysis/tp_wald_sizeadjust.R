library(ggplot2)
library(DESeq2)

# Load metadata
sample_info <- read.csv("/path/to/metadata.csv", stringsAsFactors = FALSE)

# Create group_name -> plot_name mapping
plot_name_map <- unique(sample_info[, c("group_name", "plot_name")])
plot_name_lookup <- setNames(plot_name_map$plot_name, plot_name_map$group_name)

# Find DESeq results
group_names <- unique(gsub("deseq2_results_(pI|cI)_(.*)\\.csv", "\\2",
                           list.files(pattern = "^deseq2_results_(pI|cI)_.*\\.csv$")))

for (group_name in group_names) {
  message("Processing group: ", group_name)
  
  file_pI <- paste0("deseq2_results_pI_", group_name, ".csv")
  file_cI <- paste0("deseq2_results_cI_", group_name, ".csv")
  file_dds <- paste0("dds_", group_name, ".rds")
  
  if (!file.exists(file_pI) || !file.exists(file_cI) || !file.exists(file_dds)) {
    message("Missing files for group: ", group_name)
    next
  }
  
  res_pI <- read.csv(file_pI)
  res_cI <- read.csv(file_cI)
  dds <- readRDS(file_dds)
  sf <- sizeFactors(dds)
  coldata <- as.data.frame(colData(dds))
  
  # Match samples by replicate
  ip_sf <- sf[coldata$condition == "IP"]
  ctrl_sf <- sf[coldata$condition == "Control"]
  ip_reps <- coldata$replicate[coldata$condition == "IP"]
  ctrl_reps <- coldata$replicate[coldata$condition == "Control"]
  
  correction_factors <- sapply(ctrl_reps, function(rep_id) {
    ip_sample <- names(ip_sf)[ip_reps == rep_id]
    ctrl_sample <- names(ctrl_sf)[ctrl_reps == rep_id]
    if (length(ip_sample) == 1 && length(ctrl_sample) == 1) {
      -log2(sf[ctrl_sample] / sf[ip_sample])
    } else {
      NA
    }
  })
  
  names(correction_factors) <- paste0("rep", ctrl_reps)
  correction_offset <- median(correction_factors, na.rm = TRUE)
  
  # Merge pI and cI results (retain all rows and NAs)
  merged <- merge(res_pI, res_cI, by = "transcript_id", suffixes = c("_pI", "_cI"), all = TRUE)
  
  # Compute delta stats
  merged$delta_log2FC <- merged$log2FoldChange_pI - merged$log2FoldChange_cI + correction_offset
  merged$se_delta <- sqrt(merged$lfcSE_pI^2 + merged$lfcSE_cI^2)
  merged$z_delta <- merged$delta_log2FC / merged$se_delta
  
  # Estimate empirical null SD
  z_null_region <- merged$z_delta[abs(merged$z_delta) < 0.5]
  if (length(z_null_region) < 10) {
    z_null_region <- merged$z_delta[abs(merged$z_delta) < 1]
  }
  if (length(z_null_region) < 10) {
    warning(sprintf("Group '%s': Too few values for empirical null (n = %d), using fallback sd = 0.3",
                    group_name, length(z_null_region)))
    null_sd <- 0.3
  } else {
    null_sd <- sd(z_null_region)
  }
  
  # Compute p-values and adjusted p-values
  merged$pval_delta <- pnorm(-merged$z_delta, mean = 0, sd = null_sd)
  merged$padj_delta <- p.adjust(merged$pval_delta, method = "BH")
  
  # Directionality
  merged$Directionality <- merged$log2FoldChange_pI > 0.5
  
  # Renamed: DESeq strict filter → Significant_IPvsInput
  merged$Significant_IPvsInput <- ifelse(
    is.na(merged$log2FoldChange_pI) | is.na(merged$padj_pI),
    FALSE,
    merged$log2FoldChange_pI > 0.5 & merged$padj_pI < 0.05
  )
  
  # IP vs Control significance
  merged$Significant_IPvsControl <- merged$z_delta > 0 & merged$padj_delta < 0.05
  
  # Remove old Significance column if present
  merged$Significance <- NULL
  
  # Export full result
  write.csv(merged, paste0("wald_all_", group_name, ".csv"), row.names = FALSE)
  
  # Directionality only
  merged_dir <- merged[merged$Directionality == TRUE & !is.na(merged$Directionality), ]
  write.csv(merged_dir, paste0("wald_dir_", group_name, ".csv"), row.names = FALSE)
  
  # TP: Significant_IPvsInput AND Significant_IPvsControl
  merged_tp <- merged[merged$Significant_IPvsInput & merged$Significant_IPvsControl, ]
  write.csv(merged_tp, paste0("wald_TP_", group_name, ".csv"), row.names = FALSE)

# Plot Z-distribution
z_vals <- seq(min(merged$z_delta, na.rm = TRUE), max(merged$z_delta, na.rm = TRUE), length.out = 512)
f_null <- dnorm(z_vals, mean = 0, sd = null_sd)
null_df <- data.frame(z = z_vals, density = f_null)

# Plotting status color
merged$PlotStatus <- "Not Significant"
merged$PlotStatus[merged$Significant_IPvsControl == TRUE] <- "Significant_IPvsControl"
merged$PlotStatus[merged$Significant_IPvsInput == TRUE] <- "Significant_IPvsInput"

# Priority: Significant_IPvsInput > Significant_IPvsControl
merged$PlotStatus[merged$Significant_IPvsInput == TRUE & merged$Significant_IPvsControl == TRUE] <- "Significant_IPvsInput"

plot_combined <- ggplot(merged, aes(x = z_delta, fill = PlotStatus)) +
  geom_histogram(aes(y = after_stat(density)), bins = 60, color = "black", alpha = 0.7) +
  geom_line(data = null_df, aes(x = z, y = density), color = "black", linewidth = 0.9, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
  scale_fill_manual(values = c(
    "Significant_IPvsInput" = "purple",
    "Significant_IPvsControl" = "red",
    "Not Significant" = "gray70"
  )) +
  labs(
    title = plot_name_lookup[[group_name]],
    subtitle = sprintf("Null: N(0, %.2f)", null_sd),
    x = "Z-score (IP vs Control)",
    y = "Density",
    fill = NULL
  ) +
  theme_minimal(base_size = 8.5) +
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 8, hjust = 0.5, margin = margin(b = 4)),
    axis.title.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text = element_text(size = 7),
    legend.position = "bottom",
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.3, "cm"),
    plot.margin = margin(4, 4, 4, 4)
  )

ggsave(paste0("plot_wald_", group_name, ".pdf"),
       plot = plot_combined, width = 3, height = 3, dpi = 300)

}

message("✅ All groups processed with updated filters and outputs.")
