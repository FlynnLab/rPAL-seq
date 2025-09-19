library(tidyverse)
library(viridis)
library(stringr)

# ------------------------------
# Parameters
# ------------------------------
group_name <- "metaspecies"                    # label for title/filename
input_file <- "multi_overlap_metaspecies.csv"  # input CSV
out_pdf    <- paste0("barplot_top20_", group_name, ".pdf")
z_cap      <- 50

# ------------------------------
# Load data
# ------------------------------
multi_overlap_df <- read_csv(input_file, show_col_types = FALSE)

# Validate required columns
required_cols <- c("family", "count", "category", "combined_zscore")
if (!all(required_cols %in% names(multi_overlap_df))) {
  stop("Missing required columns in input file. Expecting: family, count, category, combined_zscore")
}

# ------------------------------
# Prepare data (no category filtering)
# ------------------------------
plot_df <- multi_overlap_df %>%
  mutate(
    combined_zscore = suppressWarnings(as.numeric(combined_zscore)),
    z_capped = pmin(z_cap, combined_zscore),          # cap at 50 for coloring
    family_wrapped = str_wrap(family, width = 50)
  )

# ------------------------------
# Build global Top-20 (by count desc, then raw z-score desc)
# ------------------------------
top_df <- plot_df %>%
  arrange(desc(count), desc(combined_zscore)) %>%
  slice_head(n = 10) %>%
  mutate(family_wrapped = factor(family_wrapped, levels = rev(family_wrapped)))

max_ct <- if (nrow(top_df)) max(top_df$count, na.rm = TRUE) else 0

# ------------------------------
# Plot (same aesthetics/size), fill uses capped z-score
# ------------------------------
p <- ggplot(top_df, aes(x = family_wrapped, y = count, fill = z_capped)) +
  geom_bar(stat = "identity", color = "black", width = 0.8) +
  scale_fill_viridis_c(option = "B", name = "Avg Z-score\n (capped 50)") +
  scale_y_continuous(breaks = 0:max_ct) +
  labs(
    title = paste("Top20 TP Cross-species"),
    x = NULL, y = "Overlap Count"
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
  ) +
  coord_flip()

ggsave(
  filename = out_pdf, plot = p, width = 4,
  height = max(3, 0.2 * nrow(top_df)), dpi = 300
)

message(sprintf("Wrote %s (z-score coloring capped at %g)", out_pdf, z_cap))
