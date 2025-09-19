library(readr)
library(dplyr)

# Find all DESeq2 results files (raw input)
deseq_files <- list.files(pattern = "^deseq2_results_.*\\.csv$")

# Exclude already-annotated files
deseq_files <- deseq_files[!grepl("^Annotated_.*\\.csv$", deseq_files)]

# Load full transcript-to-biotype annotation
annotation <- read_csv(path.expand("/path/to/transcriptome/annotation.csv")) %>%
  mutate(transcript_id = as.character(transcript_id))

# Extract all known biotypes from the annotation reference
known_biotypes <- unique(annotation$biotype)

# Process each DESeq2 file
for (file in deseq_files) {
  message("Annotating: ", file)
  
  # Load DESeq2 result
  res <- read_csv(file) %>%
    mutate(transcript_id = as.character(transcript_id))
  
  # Annotate with biotype, product, family
  res_annot <- left_join(res, annotation[, c("transcript_id", "biotype", "product", "family")], by = "transcript_id")
  
  # Identify biotypes missing from this DESeq2 result
  present_biotypes <- unique(res_annot$biotype)
  missing_biotypes <- setdiff(known_biotypes, present_biotypes)
  
  # Add one dummy row per missing biotype
  if (length(missing_biotypes) > 0) {
    dummy_rows <- data.frame(
      transcript_id = NA,
      baseMean = NA,
      log2FoldChange = NA,
      lfcSE = NA,
      stat = NA,
      pvalue = NA,
      padj = NA,
      biotype = missing_biotypes,
      product = NA,
      family = NA
    )
    
    # Add any additional columns from res that might be missing in dummy
    for (col in setdiff(names(res_annot), names(dummy_rows))) {
      dummy_rows[[col]] <- NA
    }
    
    # Match column order
    dummy_rows <- dummy_rows[, names(res_annot)]
    
    # Combine
    res_annot <- bind_rows(res_annot, dummy_rows)
  }
  
  # Save annotated output
  output_file <- paste0("Annotated_", basename(file))
  write_csv(res_annot, output_file)
}
