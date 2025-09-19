library(readr)
library(dplyr)

# Define Stouffer's weighted Z-combination function
stouffer_z <- function(z, weights) {
  if (all(is.na(z)) || all(is.na(weights))) return(NA_real_)
  sum(z * weights, na.rm = TRUE) / sqrt(sum(weights^2, na.rm = TRUE))
}

# Load transcript-to-family annotation
annotation <- read_csv(path.expand("/path/to/transcriptome/annotation.csv"), show_col_types = FALSE) %>%
  mutate(transcript_id = as.character(transcript_id))

# Find all TP_<group>.csv files
tp_files <- list.files(pattern = "^wald_TP_.*\\.csv$")

for (file in tp_files) {
  # Extract group name
  group_name <- gsub("^wald_TP_(.*)\\.csv$", "\\1", file)
  message("Processing group: ", group_name)
  
  # Load transcript-level hits
  res <- read_csv(file, show_col_types = FALSE)
  
  output_file <- paste0("Family_TP_", group_name, ".csv")
  
  # Handle empty or malformed files
  if (!"transcript_id" %in% names(res) || nrow(res) == 0) {
    message("Empty or invalid file: ", file, " → writing empty output")
    
    # Create and save empty output with correct columns
    empty_df <- tibble(
      family = character(),
      avg_log2FC_pI = numeric(),
      avg_log10padj = numeric(),
      combined_zscore = numeric()
    )
    
    write_csv(empty_df, output_file)
    next
  }
  
  res <- res %>%
    mutate(transcript_id = as.character(transcript_id))
  
  # Join annotations
  res_annot <- left_join(res, annotation[, c("transcript_id", "biotype", "product", "family")], by = "transcript_id")
  
  # Compute log10(padj) safely
  if (!"log10padj_pI" %in% colnames(res_annot)) {
    res_annot <- res_annot %>%
      mutate(padj_pI = as.numeric(padj_pI),
             log10padj_pI = -log10(padj_pI))
  }
  
  # Summarize hits by family
  family_summary <- res_annot %>%
    filter(!is.na(family)) %>%
    group_by(family) %>%
    summarise(
      avg_log2FC_pI = weighted.mean(log2FoldChange_pI, baseMean_pI, na.rm = TRUE),
      avg_log10padj = weighted.mean(log10padj_pI, baseMean_pI, na.rm = TRUE),
      combined_zscore = stouffer_z(z_delta, baseMean_pI),
      .groups = "drop"
    ) %>%
    arrange(desc(combined_zscore))
  
  # Save summary
  write_csv(family_summary, output_file)
}