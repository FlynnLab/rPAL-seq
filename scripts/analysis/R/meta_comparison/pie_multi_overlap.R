library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(paletteer)
library(stringr)
library(tidyr)
library(cowplot)

# Load transcript-to-biotype mapping
tx2gene <- read_csv("/path/to/transcriptome/annotation.csv", show_col_types = FALSE) %>%
  mutate(
    family  = as.character(family),
    biotype = as.character(biotype)
  )

# DiFi-style mapping (keep pairs; dedup later with distinct(family, biotype))
tx2gene_uni <- tx2gene %>% distinct(family, biotype)

# Biotype definitions and colors
biotype_levels <- c("miRNA", "vault-RNA", "rRNA", "mt-rRNA", "mt-tRNA", "ribozyme",
                    "Y-RNA", "scaRNA", "snoRNA", "snRNA", "others", "tRNA")
biotype_colors <- paletteer::paletteer_d("PrettyCols::Rainbow")[1:length(biotype_levels)]
names(biotype_colors) <- biotype_levels

# File paths and labels
base_files <- c(
  "multi_overlap_hematopoietic.csv",
  "multi_overlap_non_hematopoietic.csv",
  "multi_overlap_all_cell_lines.csv"
)
labels_full   <- c("Hematopoietic", "Non-hematopoietic", "All cell lines")
labels_top20  <- paste("Top 20 conserved:\n", labels_full)

# ---- Pie chart function (logic fixed; aesthetics unchanged) ----
make_count_pie <- function(df, sample_name, show_legend = FALSE) {
  df <- df %>%
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
      legend.text  = element_text(size = 9),
      legend.position = if (show_legend) "right" else "none"
    )
}

# ---- Read & prepare pies ----
pie_list <- list()

for (i in seq_along(base_files)) {
  # Read full table
  df_full <- read_csv(base_files[i], show_col_types = FALSE)
  
  # Full data pie
  pie_list[[length(pie_list) + 1]] <- make_count_pie(df_full, labels_full[i], show_legend = FALSE)
  
  # Top 20 data
  df_top20 <- df_full %>%
    arrange(desc(count), desc(combined_zscore)) %>%
    slice_head(n = 20)
  
  pie_list[[length(pie_list) + 1]] <- make_count_pie(df_top20, labels_top20[i], show_legend = (i == length(base_files)))
}

# ---- Build legend and pie grid ----
legend <- cowplot::get_legend(pie_list[[length(pie_list)]])

pie_list_nolegend <- lapply(pie_list, function(p) p + theme(legend.position = "none"))

# Arrange pies: top 3 are full, bottom 3 are top20
pie_matrix <- list(
  pie_list_nolegend[[1]],  # Hematopoietic
  pie_list_nolegend[[3]],  # Non-hematopoietic
  pie_list_nolegend[[5]],  # All cell lines
  pie_list_nolegend[[2]],  # Top 20: Hematopoietic
  pie_list_nolegend[[4]],  # Top 20: Non-hematopoietic
  pie_list_nolegend[[6]]   # Top 20: All cell lines
)

pie_grid <- wrap_plots(pie_matrix, ncol = 3, byrow = TRUE)

final_plot <- cowplot::plot_grid(
  pie_grid,
  legend,
  rel_widths = c(1, 0.25),
  nrow = 1
)

# ---- Save ----
ggsave("Pie_grid_overlap.pdf", plot = final_plot, width = 6.5, height = 4.5, dpi = 300)
