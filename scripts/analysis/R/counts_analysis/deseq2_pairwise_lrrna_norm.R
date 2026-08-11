library(DESeq2)
library(dplyr)

# Load count matrix
counts <- read.csv("/path/to/count_matrix.csv", row.names = 1, check.names = FALSE)
all_sample_names <- colnames(counts)

# Load sample metadata
sample_info <- read.csv("/path/to/metadata.csv", stringsAsFactors = FALSE)

# Sample and group names are used as lookup keys against the count matrix columns.
# Stray whitespace makes a key silently fail to match, so trim on read.
sample_info$sample_name <- trimws(sample_info$sample_name)
sample_info$group_name <- trimws(sample_info$group_name)

# Read your mapping file
tx2gene <- read.csv("/path/to/transcriptome/annotation.csv")

# Define families of interest
long_rrna_families <- c("45S ribosomal RNA", "28S ribosomal RNA", "18S ribosomal RNA", "l-rRNA", "s-rRNA")

# Filter and extract transcript IDs
long_rrna_ids <- tx2gene %>%
  filter(family %in% long_rrna_families) %>%
  pull(transcript_id)

# Define pairwise comparisons
comparisons <- list(
  pI = c("IP", "Input"),
  cI = c("Control", "Input")
)

# Helper function to generate sample names like "1i", "1p", "1h"
make_sample_names <- function(ids) {
  c(paste0(ids, "i"), paste0(ids, "p"), paste0(ids, "h"))
}

# Loop through each unique group
for (grp in unique(sample_info$group_name)) {
  message("Processing group: ", grp)
  
  group_samples <- subset(sample_info, group_name == grp)
  
  # Expected sample column names: e.g., 1i, 1p, 1h
  expected_names <- unlist(lapply(group_samples$sample_name, make_sample_names))
  
  # Check if all required samples exist. A group that is short of its declared samples
  # would otherwise be analysed at reduced n without any error, so report which columns
  # are missing and skip rather than proceed silently.
  existing_names <- expected_names[expected_names %in% all_sample_names]
  missing_names <- setdiff(expected_names, all_sample_names)
  if (length(missing_names) > 0) {
    message("  Group ", grp, ": ", length(existing_names), " of ", length(expected_names),
            " declared samples found in the count matrix.")
    message("    Missing: ", paste(missing_names, collapse = ", "))
    next
  }
  
  # Parse sample names
  parse_sample <- function(sname) {
    id <- sub("[iph]$", "", sname)
    suffix <- substr(sname, nchar(sname), nchar(sname))
    list(id = id, suffix = suffix)
  }
  
  parsed_samples <- do.call(rbind, lapply(existing_names, function(sn) {
    parsed <- parse_sample(sn)
    data.frame(
      sample = sn,
      sample_id = parsed$id,
      suffix = parsed$suffix,
      stringsAsFactors = FALSE
    )
  }))
  
  # Merge metadata
  coldata <- merge(parsed_samples, group_samples,
                   by.x = "sample_id", by.y = "sample_name",
                   all.x = TRUE, sort = FALSE)
  
  # Assign condition
  coldata$condition <- factor(dplyr::recode(coldata$suffix,
                                            "i" = "Input",
                                            "p" = "IP",
                                            "h" = "Control"),
                              levels = c("Input", "IP", "Control"))
  
  # Finalize coldata
  coldata <- coldata[match(existing_names, coldata$sample), ]
  rownames(coldata) <- coldata$sample
  coldata$replicate <- factor(coldata$replicates)
  
  # Subset counts
  counts_sub <- counts[, existing_names]
  
  # Filter: keep genes with at least 2 non-zero counts in Input **and** IP; ignore Control
  groups <- coldata$condition
  min_count <- 0
  min_reps <- 2
  
  global_keep <- apply(counts_sub, 1, function(x) {
    sum(x[groups == "Input"] > min_count) >= min_reps &&
      sum(x[groups == "IP"] > min_count) >= min_reps &&
      sum(x[groups == "Control"] > min_count) >= 0
  })
  
  counts_sub <- counts_sub[global_keep, ]
  
  
  # Size factors using rRNAs
  rrna_counts <- counts_sub[rownames(counts_sub) %in% long_rrna_ids, , drop = FALSE]
  umi_totals_rrna <- colSums(rrna_counts)
  sf_rrna <- umi_totals_rrna / mean(umi_totals_rrna)
  
  # Remove rRNAs for DE
  counts_sub <- counts_sub[!rownames(counts_sub) %in% long_rrna_ids, ]
  
  # Build and save full DESeq object
  dds_full <- DESeqDataSetFromMatrix(countData = counts_sub, colData = coldata, design = ~ replicate + condition)
  
  # Skip group if any size factor is zero or NA
  if (any(sf_rrna <= 0 | is.na(sf_rrna))) {
    message("  Skipping group ", grp, ": invalid size factors (0 or NA) from rRNA counts.")
    next
  }
  
  sizeFactors(dds_full) <- sf_rrna
  saveRDS(dds_full, paste0("dds_", grp, ".rds"))
  
  # Run comparisons
  for (comp in names(comparisons)) {
    comp_levels <- comparisons[[comp]]
    message("  Running ", comp_levels[1], " vs ", comp_levels[2])
    
    samples_keep <- rownames(coldata)[coldata$condition %in% comp_levels]
    if (length(samples_keep) < 4) {
      message("    Skipping comparison ", comp, ": not enough samples.")
      next
    }
    
    coldata_comp <- coldata[samples_keep, , drop = FALSE]
    coldata_comp$condition <- droplevels(coldata_comp$condition)
    counts_comp <- counts_sub[, samples_keep]
    sf_sub <- sf_rrna[samples_keep]
    
    dds <- DESeqDataSetFromMatrix(countData = counts_comp, colData = coldata_comp, design = ~ replicate + condition)
    sizeFactors(dds) <- sf_sub
    
    dds <- DESeq(dds, fitType = "local")
    res <- results(dds, contrast = c("condition", comp_levels[1], comp_levels[2]))
    res_df <- as.data.frame(res)
    res_df$transcript_id <- rownames(res_df)
    
    out_file <- paste0("deseq2_results_", comp, "_", grp, ".csv")
    write.csv(res_df, out_file, row.names = FALSE)
  }
  
  message("  Completed group: ", grp)
}

message("✅ All groups processed.")
