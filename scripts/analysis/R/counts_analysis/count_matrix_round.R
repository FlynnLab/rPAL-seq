library(dplyr)

# Set the input directory
input_dir <- "path/to/tsv/files"  # Change this to your actual directory

# List all .tsv files in the current directory
files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)

# Initialize a list to store the counts
count_list <- list()

# Process each file
for (file in files) {
  # Extract sample name from filename (remove path and extension)
  sample_name <- tools::file_path_sans_ext(basename(file))
  
  # Read TSV with expected columns: 'transcript' and 'count'
  df_in <- read.delim(file, header = TRUE, stringsAsFactors = FALSE)
  
  if (!all(c("transcript", "count") %in% names(df_in))) {
    stop(paste("File", basename(file), "must have 'transcript' and 'count' columns."))
  }
  
  # Round float counts to nearest integer
  df_in$count <- round(as.numeric(df_in$count))
  
  # Build per-sample frame with 'transcript_id' for compatibility with original output
  df <- data.frame(
    transcript_id = df_in$transcript,
    tmp = df_in$count,
    stringsAsFactors = FALSE
  )
  
  # Rename the count column to the sample name
  colnames(df)[2] <- sample_name
  
  # Add to list
  count_list[[sample_name]] <- df
}

# Merge all counts by 'transcript_id'
count_matrix <- Reduce(function(x, y) full_join(x, y, by = "transcript_id"), count_list)

# Replace NAs with 0 in numeric columns
num_cols <- setdiff(names(count_matrix), "transcript_id")
count_matrix[num_cols][is.na(count_matrix[num_cols])] <- 0

# Write to CSV (same output format as original)
write.csv(count_matrix, file = "count_matrix.csv", row.names = FALSE)
