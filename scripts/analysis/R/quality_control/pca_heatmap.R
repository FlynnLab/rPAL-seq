library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)

# Define your group label
group_name <- "your_group_label" # "A" = AML3 in this example

# Define sample metadata
condition <- c(
  "VA1I" = "Input", "VA2I" = "Input", "VA3I" = "Input", "VA4I" = "Input",
  "VA1p" = "IP", "VA2p" = "IP", "VA3p" = "IP", "VA4p" = "IP",
  "dVA1p" = "IP", "dVA2p" = "IP", "dVA3p" = "IP", "dVA4p" = "IP",
  "VA1h" = "Control", "VA2h" = "Control", "VA3h" = "Control", "VA4h" = "Control",
  "dVA1h" = "Control", "dVA2h" = "Control", "dVA3h" = "Control", "dVA4h" = "Control"
)

replicate <- c(
  "VA1I" = "1", "VA2I" = "2", "VA3I" = "3", "VA4I" = "4",
  "VA1p" = "1", "VA2p" = "2", "VA3p" = "3", "VA4p" = "4",
  "dVA1p" = "1", "dVA2p" = "2", "dVA3p" = "3", "dVA4p" = "4",
  "VA1h" = "1", "VA2h" = "2", "VA3h" = "3", "VA4h" = "4",
  "dVA1h" = "1", "dVA2h" = "2", "dVA3h" = "3", "dVA4h" = "4"
)

samples <- names(condition)

# Load count matrix
counts_all <- read.csv("count_matrix.csv", row.names = 1, check.names = FALSE)
counts <- counts_all[, samples]

# Build colData
coldata <- data.frame(
  sample = samples,
  condition = factor(condition[samples], levels = c("Input", "IP", "Control")),
  replicate = factor(replicate[samples]),
  row.names = samples
)

# Define long rRNA IDs
long_rrna_ids <- c(
  "rRNA_RNA18SN1_166629", "rRNA_RNA18SN2_166610",
  "rRNA_RNA28SN1_166631", "rRNA_RNA28SN2_166612", "rRNA_RNA28SN3_166621",
  "rRNA_RNA28SN4_172372", "rRNA_RNA28SN5_178714",
  "rRNA_RNA45SN1_166625", "rRNA_RNA45SN2_166607", "rRNA_RNA45SN3_166618",
  "rRNA_RNA45SN4_172369", "rRNA_RNA45SN5_178711",
  "rRNA_RNR1_192237", "rRNA_RNR2_192239"
)

# Estimate size factors from rRNA
rrna_counts <- counts_all[rownames(counts_all) %in% long_rrna_ids, samples, drop = FALSE]
rrna_sf <- colSums(rrna_counts)
rrna_sf <- rrna_sf / mean(rrna_sf)  # Normalize

# Force size factors from 'p' samples to corresponding 'h' samples
p_samples <- grep("p$", names(rrna_sf), value = TRUE)
h_samples <- grep("h$", names(rrna_sf), value = TRUE)
p_to_h_map <- setNames(p_samples, gsub("p$", "h", p_samples))  # Map p -> h
rrna_sf[h_samples] <- rrna_sf[p_to_h_map[h_samples]]

# Remove rRNAs for DE analysis
counts <- counts[!rownames(counts) %in% long_rrna_ids, ]

# Apply expression filter: remove genes with no expression in groups
min_count <- 0
min_reps <- 2
groups <- condition[colnames(counts)]

global_keep <- apply(counts, 1, function(x) {
  all(sapply(unique(groups), function(g) {
    sum(x[groups == g] > min_count) >= min_reps
  }))
})
counts <- counts[global_keep, ]

# Build DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ replicate + condition)
sizeFactors(dds) <- rrna_sf[colnames(counts)]

# Run DESeq
dds <- DESeq(dds, fitType = "local")
saveRDS(dds, paste0("dds_", group_name, "_5group.rds"))

# Variance stabilizing transformation
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
vsd_mat <- assay(vsd)

# Assign group for plotting (VA / dVA)
sample_names <- rownames(coldata)

group_labels <- dplyr::case_when(
  grepl("I$", sample_names) ~ "Input",
  grepl("^dVA", sample_names) & grepl("p$", sample_names) ~ "Low-input IP",
  grepl("^VA",  sample_names) & grepl("p$", sample_names) ~ "Bulk IP",
  grepl("^dVA", sample_names) & grepl("h$", sample_names) ~ "Low-input Control",
  grepl("^VA",  sample_names) & grepl("h$", sample_names) ~ "Bulk Control",
  TRUE ~ "Unknown"
)

coldata$plot_group <- factor(
  group_labels,
  levels = c("Input", "Bulk IP", "Low-input IP", "Bulk Control", "Low-input Control")
)

# PCA
pca <- prcomp(t(vsd_mat), scale. = TRUE)
pca_data <- as.data.frame(pca$x[, 1:2])
pca_data$sample <- rownames(pca_data)
pca_data$group <- coldata[rownames(pca_data), "plot_group"]

# Convex hulls
get_hull <- function(df) {
  if (nrow(df) >= 3) df[chull(df$PC1, df$PC2), ] else df
}
hull_data <- pca_data %>%
  group_by(group) %>%
  group_modify(~ get_hull(.x)) %>%
  ungroup()

# Colors to match the 5 groups above
group_colors <- c(
  "Input"              = "#8C2981FF",
  "Bulk IP"            = "#E01A4FFF",
  "Low-input IP"       = "#F9C22EFF",
  "Bulk Control"       = "#53B3CBFF",
  "Low-input Control"  = "#7DCFB6FF"
)

# PCA plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 2) +
  geom_polygon(data = hull_data, aes(fill = group), alpha = 0.15, color = NA, show.legend = FALSE) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(title = "PCA, AML3", x = "PC1", y = "PC2", color = "Group") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid   = element_blank(),
    axis.line    = element_line(color = "black"),
    axis.ticks   = element_line(color = "black"),
    plot.title   = element_text(face = "bold", size = 11),
    legend.title = element_text(face = "bold", size = 9),
    legend.text  = element_text(size = 8)
  )

ggsave(paste0("PCA_", group_name, "_5group.pdf"), plot = pca_plot, width = 3.5, height = 2.5)

# Heatmap of Top 500 Variable Genes
top_var_genes <- names(sort(apply(vsd_mat, 1, var), decreasing = TRUE))[1:500]
annotation_col <- data.frame(Group = coldata$plot_group)
rownames(annotation_col) <- rownames(coldata)
annotation_colors <- list(Group = group_colors)

pheatmap(
  vsd_mat[top_var_genes, ],
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  annotation_names_col = FALSE,
  show_rownames = FALSE,
  scale = "row",
  clustering_method = "ward.D2",
  fontsize_col = 8,
  main = "Top 500 Variable Genes, AML3",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  filename = paste0("heatmap_top500_", group_name, "_5group.pdf"),
  width = 6,
  height = 6
)
