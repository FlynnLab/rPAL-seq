library(ggplot2)
library(ggrepel)
library(readr)
library(dplyr)
library(paletteer)

# Set group name
group_name <- "your_sample_group"

# Load sample info and extract plot name based on group_name
sample_info <- read.csv("/path/to/transcriptome/annotation.csv", stringsAsFactors = FALSE)
plot_name <- sample_info %>%
  filter(group_name == !!group_name) %>%
  pull(plot_name)

if (length(plot_name) == 0) {
  stop(paste0("No matching plot_name found for group_name = ", group_name))
}

# Key parameters
y_cap <- 90
label_logfc_cutoff <- 0.5
label_logpadj_cutoff <- 1.3

# Biotype levels and colors
biotype_levels <- c(
  "miRNA", "vault-RNA", "rRNA", "mt-rRNA", "mt-tRNA", "ribozyme", "Y-RNA", "scaRNA", "snoRNA", "snRNA", "others", "tRNA"
)

biotype_colors <- paletteer::paletteer_d("PrettyCols::Rainbow")[1:length(biotype_levels)]
names(biotype_colors) <- biotype_levels

# Construct input filenames
input_file <- paste0("Annotated_deseq2_results_pI_", group_name, ".csv")
tp_file <- paste0("wald_TP_", group_name, ".csv")

# Load data
res_annot <- read_csv(input_file)
tp_data <- read_csv(tp_file)

# Process biotype
res_annot$biotype <- factor(res_annot$biotype, levels = biotype_levels)

# Annotate data
res_annot <- res_annot %>%
  mutate(
    log10padj = ifelse(padj == 0, 300, -log10(padj)),
    TP = transcript_id %in% tp_data$transcript_id,
    alpha = ifelse(TP, 0.7, 0.2)
  )

# Select top TPs to label
top_labels <- res_annot %>%
  filter(TP, log2FoldChange > label_logfc_cutoff, log10padj > label_logpadj_cutoff) %>%
  mutate(score = abs(log2FoldChange) * log10padj) %>%
  arrange(-score) %>%
  slice_head(n = 4)

# Generate volcano plot
p <- ggplot(res_annot, aes(x = log2FoldChange, y = pmin(log10padj, y_cap), color = biotype)) +
  geom_point(aes(alpha = alpha), size = 2.5) +
  scale_color_manual(name = "Biotype", values = biotype_colors, drop = FALSE) +
  
  geom_point(data = res_annot %>% filter(log10padj > y_cap),
             aes(x = log2FoldChange, y = y_cap, alpha = alpha),
             shape = 21, fill = "white", stroke = 0.5, size = 3, inherit.aes = FALSE, show.legend = FALSE) +
  
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", size = 0.6) +
  geom_vline(xintercept = c(-label_logfc_cutoff, label_logfc_cutoff), linetype = "dashed", color = "grey40", size = 0.6) +
  
  geom_text_repel(
    data = top_labels,
    aes(label = product),
    color = "black",
    size = 3,
    box.padding = 0.6,
    point.padding = 0.4,
    force = 50,
    force_pull = 0.02,
    nudge_x = 0,
    nudge_y = 0.1* y_cap,
    segment.color = "grey40",
    segment.size = 0.3,
    min.segment.length = 0.05,
    max.overlaps = 100,
    show.legend = FALSE
  ) +
  
  annotate("text", x = max(res_annot$log2FoldChange, na.rm = TRUE), y = y_cap + 0.5,
           label = paste0("Y-axis capped at ", y_cap), hjust = 1, size = 3.2, color = "grey50") +
  
  labs(
    title = plot_name,
    x = "log2FC (IP/Input)",
    y = expression(-log[10](padj)),
    color = "Biotype"
  ) +
  
  theme_minimal(base_size = 10) +
  theme(
    panel.grid   = element_blank(),
    axis.line    = element_line(color = "black"),
    axis.ticks   = element_line(color = "black"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    legend.title = element_text(face = "bold", size = 9),
    legend.text  = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  
  scale_alpha_continuous(range = c(0.2, 0.7), guide = 'none')

# Save and show
print(p)
ggsave(filename = paste0("Volcano_transcript_", group_name, "_small.pdf"), plot = p, width = 4, height = 3.2, dpi = 300)
