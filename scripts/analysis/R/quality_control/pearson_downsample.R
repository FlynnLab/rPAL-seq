library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(DESeq2)

# Define group name
group_name <- "your_group_label" # "VH" = HeLa in this example

# Construct file names
file_pI    <- paste0("deseq2_results_pI_", group_name, ".csv")       # IP vs Input (bulk)
file_cI    <- paste0("deseq2_results_cI_", group_name, ".csv")       # Control vs Input (bulk)
file_pI_d  <- paste0("deseq2_results_pI_d", group_name, ".csv")      # IP vs Input (low input)
file_cI_d  <- paste0("deseq2_results_cI_d", group_name, ".csv")      # Control vs Input (low input)
dds_file   <- paste0("dds_", group_name, ".rds")                     # DESeq2 object for bulk

# Read input data
df_pI    <- read.csv(file_pI, stringsAsFactors = FALSE)
df_cI    <- read.csv(file_cI, stringsAsFactors = FALSE)
df_pI_d  <- read.csv(file_pI_d, stringsAsFactors = FALSE)
df_cI_d  <- read.csv(file_cI_d, stringsAsFactors = FALSE)

# Load DESeq2 object and extract size factors
dds <- readRDS(dds_file)
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
correction_offset <- mean(correction_factors, na.rm = TRUE)

# ===========================
# Use only shared transcript_ids
# ===========================
shared_ids <- Reduce(intersect, list(
  df_pI$transcript_id,
  df_cI$transcript_id,
  df_pI_d$transcript_id,
  df_cI_d$transcript_id
))
message("Shared transcript IDs: ", length(shared_ids))

# Filter all data frames
df_pI    <- df_pI[df_pI$transcript_id %in% shared_ids, ]
df_cI    <- df_cI[df_cI$transcript_id %in% shared_ids, ]
df_pI_d  <- df_pI_d[df_pI_d$transcript_id %in% shared_ids, ]
df_cI_d  <- df_cI_d[df_cI_d$transcript_id %in% shared_ids, ]

# Rename columns using base R
tmp1 <- df_pI[, c("transcript_id", "log2FoldChange", "lfcSE")]
colnames(tmp1) <- c("transcript_id", "log2FC_pI_bulk", "lfcSE_pI")

tmp2 <- df_cI[, c("transcript_id", "log2FoldChange", "lfcSE")]
colnames(tmp2) <- c("transcript_id", "log2FC_cI_bulk", "lfcSE_cI")

tmp3 <- df_pI_d[, c("transcript_id", "log2FoldChange")]
colnames(tmp3) <- c("transcript_id", "log2FC_pI_lowinput")

tmp4 <- df_cI_d[, c("transcript_id", "log2FoldChange")]
colnames(tmp4) <- c("transcript_id", "log2FC_cI_lowinput")

# Join all together
df <- tmp1 %>%
  inner_join(tmp2, by = "transcript_id") %>%
  inner_join(tmp3, by = "transcript_id") %>%
  inner_join(tmp4, by = "transcript_id")

# ===========================
# Plotting
# ===========================

# Plot 1: IP
plot1 <- ggplot(df, aes(x = log2FC_pI_bulk, y = log2FC_pI_lowinput)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = "#7DCFB6") +
  ggtitle("IP vs Input") +
  xlab("log2FC (bulk)") +
  ylab("log2FC (low input)") +
  annotate("text", x = min(df$log2FC_pI_bulk), y = max(df$log2FC_pI_lowinput),
           label = paste("r =", round(cor(df$log2FC_pI_bulk, df$log2FC_pI_lowinput, use = "complete.obs"), 2)),
           hjust = 0, size = 2.8) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid   = element_blank(),
    axis.line    = element_line(color = "black"),
    axis.ticks   = element_line(color = "black"),
    plot.title   = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 9)
  )

# Plot 2: Control
plot2 <- ggplot(df, aes(x = log2FC_cI_bulk, y = log2FC_cI_lowinput)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5, color = "#7DCFB6") +
  ggtitle("Control vs Input") +
  xlab("log2FC (bulk)") +
  ylab("log2FC (low input)") +
  annotate("text", x = min(df$log2FC_cI_bulk), y = max(df$log2FC_cI_lowinput),
           label = paste("r =", round(cor(df$log2FC_cI_bulk, df$log2FC_cI_lowinput, use = "complete.obs"), 2)),
           hjust = 0, size = 2.8) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid   = element_blank(),
    axis.line    = element_line(color = "black"),
    axis.ticks   = element_line(color = "black"),
    plot.title   = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 9)
  )

# Define combined plot title from group_name
plot_titles <- c(
  "VH"  = "HeLa, Bulk vs Low-input",
  "VH"  = "HeLa, Bulk vs Low-input"
)
combined_title <- plot_titles[[group_name]]

# Save output (only plot1 and plot2)
pdf_filename <- paste0("Pearson_bulk_vs_low_", group_name, ".pdf")
ggsave(pdf_filename,
       arrangeGrob(plot1, plot2, ncol = 2, top = textGrob(combined_title, gp = gpar(fontsize = 11, fontface = "bold"))),
       width = 3.5, height = 2.2, dpi = 300)
