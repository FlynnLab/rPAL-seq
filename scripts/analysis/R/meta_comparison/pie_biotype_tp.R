library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(paletteer)
library(stringr)
library(tidyr)
library(cowplot)

# Load transcript to biotype mapping
tx2gene <- read_csv(path.expand("/path/to/transcriptome/annotation.csv"), show_col_types = FALSE) %>%
  mutate(
    family = as.character(family),
    biotype = as.character(biotype)
  )

# Use the DiFi-style mapping (keep pairs; dedup later with distinct(family, biotype))
tx2gene_uni <- tx2gene %>% distinct(family, biotype)

# Load sample metadata
sample_info <- read.csv("/path/to/metadata.csv", stringsAsFactors = FALSE)

# Define categories
hematopoietic <- c("AML2", "AML3", "Molm13", "Jurkat", "Jeko1", "Nalm6")
non_hematopoietic <- c("HeLa", "293T", "A549", "Huh7", "DiFi", "LPS853", "LN308")

# Define biotype palette
biotype_levels <- c("miRNA", "vault-RNA", "rRNA", "mt-rRNA", "mt-tRNA", "ribozyme",
                    "Y-RNA", "scaRNA", "snoRNA", "snRNA", "others", "tRNA")
biotype_colors <- paletteer::paletteer_d("PrettyCols::Rainbow")[1:length(biotype_levels)]
names(biotype_colors) <- biotype_levels

# Pie chart function (logic fixed; aesthetics unchanged)
make_count_pie <- function(file, sample_name, show_legend = FALSE) {
  df <- read_csv(file, show_col_types = FALSE) %>%
    mutate(family = as.character(family)) %>%
    left_join(tx2gene_uni, by = "family") %>%
    mutate(biotype = ifelse(is.na(biotype), "others", biotype),
           biotype = factor(biotype, levels = biotype_levels))
  
  pie_df <- df %>%
    distinct(family, biotype) %>%                # key fix: one row per (family, biotype)
    count(biotype, name = "count") %>%
    complete(biotype = biotype_levels, fill = list(count = 0)) %>%
    mutate(
      percent = count / sum(count) * 100,
      label = ifelse(percent >= 5, paste0(round(percent, 1), "%"), "")
    )
  
  ggplot(pie_df, aes(x = "", y = count, fill = biotype)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = biotype_colors, breaks = biotype_levels, drop = FALSE) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
    labs(title = sample_name) +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9),
      legend.position = if (show_legend) "right" else "none"
    )
}

# Helper to generate pie grids
generate_pie_grid <- function(sample_info_subset, output_file) {
  sample_info_subset <- sample_info_subset %>%
    distinct(group_name, .keep_all = TRUE) %>%
    mutate(file = paste0("Family_TP_", group_name, ".csv"))
  
  sample_info_subset$show_legend <- FALSE
  if (nrow(sample_info_subset) > 0) {
    sample_info_subset$show_legend[1] <- TRUE
  }
  
  pie_list <- purrr::pmap(
    list(
      file = sample_info_subset$file,
      sample_name = sample_info_subset$plot_name,
      show_legend = sample_info_subset$show_legend
    ),
    make_count_pie
  )
  
  # Dummy pie for legend extraction
  legend_plot <- make_count_pie(
    file = sample_info_subset$file[1],
    sample_name = "",
    show_legend = TRUE
  )
  legend <- cowplot::get_legend(legend_plot + theme(legend.position = "right"))
  
  pie_list_nolegend <- lapply(pie_list, function(p) p + theme(legend.position = "none"))
  pie_grid <- wrap_plots(pie_list_nolegend, ncol = 4)
  
  final_plot <- cowplot::plot_grid(
    pie_grid,
    legend,
    rel_widths = c(1, 0.15),
    nrow = 1
  )
  
  ggsave(output_file, plot = final_plot, width = 8, height = 4.5, dpi = 300)
}

# Generate for hematopoietic group
hematopoietic_samples <- sample_info %>%
  filter(group_name %in% hematopoietic)
generate_pie_grid(hematopoietic_samples, "Pie_grid_TP_hematopoietic.pdf")

# Generate for non-hematopoietic group
non_hematopoietic_samples <- sample_info %>%
  filter(group_name %in% non_hematopoietic)
generate_pie_grid(non_hematopoietic_samples, "Pie_grid_TP_non_hematopoietic.pdf")
