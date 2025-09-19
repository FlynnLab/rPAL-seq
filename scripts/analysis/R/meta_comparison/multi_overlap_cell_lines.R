library(tidyverse)
library(viridis)
library(stringr)

# --------- PARAMETERS ---------
hematopoietic <- c("AML2", "AML3", "Molm13", "Jurkat", "Jeko1", "Nalm6")
non_hematopoietic <- c("HeLa", "293T", "A549", "Huh7", "DiFi", "LPS853", "LN308")

# Load all matching files
file_list <- list.files(pattern = "^Family_TP_.*\\.csv$")

# --------- HELPER FUNCTIONS ---------

# Determine group based on filename
get_group <- function(filename) {
  if (any(str_detect(filename, hematopoietic))) {
    return("hematopoietic")
  } else if (any(str_detect(filename, non_hematopoietic))) {
    return("non_hematopoietic")
  } else {
    return("unknown")
  }
}

# Main group processing function
process_group <- function(group_name, group_files) {
  message("Processing group: ", group_name)
  
  # Read and combine data from all files in the group
  family_data <- bind_rows(lapply(group_files, function(file) {
    read_csv(file, show_col_types = FALSE) %>%
      select(family, combined_zscore) %>%
      mutate(file = file)
  }))
  
  # Count how many files each family appears in
  family_counts <- family_data %>%
    distinct(family, file) %>%
    count(family, name = "count")
  
  # Compute mean combined z-score per family
  zscores <- family_data %>%
    group_by(family) %>%
    summarise(combined_zscore = mean(combined_zscore, na.rm = TRUE), .groups = "drop")
  
  # Merge count and zscore
  summary_df <- left_join(family_counts, zscores, by = "family")
  
  # Save full overlap table
  write_csv(summary_df, paste0("multi_overlap_", group_name, ".csv"))
  
  # Take top 20 families by count, then z-score
  filtered_df <- summary_df %>%
    arrange(desc(count), desc(combined_zscore)) %>%
    slice_head(n = 10) %>%
    mutate(family_wrapped = str_wrap(family, width = 40)) %>%
    mutate(family_wrapped = factor(family_wrapped, levels = rev(family_wrapped)))
  
  # Barplot
  p <- ggplot(filtered_df, aes(x = family_wrapped, y = count, fill = combined_zscore)) +
    geom_bar(stat = "identity", color = "black", width = 0.8) +
    scale_fill_viridis_c(option = "B", name = "Avg Z-score") +
    labs(
      title = paste("Top 10 Most Conserved GlycoRNAs:\n", group_name),
      x = NULL,
      y = "Overlap Count"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid   = element_blank(),
      axis.line    = element_line(color = "black"),
      axis.ticks   = element_line(color = "black"),
      plot.title   = element_text(face = "bold", size = 10),
      axis.title = element_text(size = 9),
      axis.text = element_text(face = "bold", size = 9),
      legend.title = element_text(face = "bold", size = 9),
      legend.text  = element_text(size = 8),
      plot.margin = margin(5, 5, 5, 5)
    ) +
    coord_flip()
  
  ggsave(
    filename = paste0("multi_overlap_", group_name, "_barplot.pdf"),
    plot = p,
    width = 4,
    height = max(3, 0.2 * nrow(filtered_df)),
    dpi = 300
  )
  
  # Histogram of transcript counts per overlap count
  hist_color <- case_when(
    group_name == "hematopoietic" ~ "#E01A4F",
    group_name == "non_hematopoietic" ~ "#7DCFB6",
    group_name == "all_cell_lines" ~ "#C6BDE8",
    TRUE ~ "grey70"
  )
  
  p_hist <- ggplot(summary_df, aes(x = count)) +
    geom_histogram(binwidth = 1, fill = hist_color, color = "black", boundary = 0.5, closed = "left") +
    labs(
      title = paste("Conserved GlycoRNA Distribution:\n", group_name),
      x = "Overlapping cell line count",
      y = "Transcript family count"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid   = element_blank(),
      axis.line    = element_line(color = "black"),
      axis.ticks   = element_line(color = "black"),
      plot.title   = element_text(face = "bold", size = 10),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      legend.title = element_text(face = "bold", size = 9),
      legend.text  = element_text(size = 8),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  ggsave(
    filename = paste0("multi_overlap_", group_name, "_histogram.pdf"),
    plot = p_hist,
    width = 2.8,
    height = 2.8,
    dpi = 300
  )
}

# --------- MAIN ---------

# Assign group labels
file_groups <- tibble(
  file = file_list,
  group = sapply(file_list, get_group)
)

# Process groups — no need for min_count_required anymore
heme_files <- file_groups %>% filter(group == "hematopoietic") %>% pull(file)
process_group("hematopoietic", heme_files)

nonheme_files <- file_groups %>% filter(group == "non_hematopoietic") %>% pull(file)
process_group("non_hematopoietic", nonheme_files)

all_files <- file_groups %>% filter(group != "unknown") %>% pull(file)
process_group("all_cell_lines", all_files)
