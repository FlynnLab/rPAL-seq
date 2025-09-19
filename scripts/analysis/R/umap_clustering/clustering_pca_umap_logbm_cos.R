library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(umap)
library(stringr)
library(cluster)
library(RColorBrewer)
library(proxy)
library(stats)

# -------------------------------
# Load and prepare input
# -------------------------------

ignore_list <- c("sample_names_to_exclude")
files <- list.files(pattern = "^wald_all_.*\\.csv$")
if (length(files) == 0) stop("No wald_all_*.csv files found.")

extract_group <- function(filename) {
  str_match(filename, "^wald_all_(.*)\\.csv$")[, 2]
}

metric <- "Log_Expression"
group_vectors <- list()
transcript_sets <- list()

for (f in files) {
  group <- extract_group(f)
  
  if (group %in% ignore_list) {
    message("⏭️ Skipping ignored group: ", group)
    next
  }
  
  df <- read.csv(f, stringsAsFactors = FALSE)
  
  required_cols <- c("z_delta", "transcript_id", "baseMean_pI", "baseMean_cI")
  if (!all(required_cols %in% colnames(df))) {
    warning(paste("Skipping", f, "- missing required columns"))
    next
  }
  
  df <- df[complete.cases(df[, required_cols]), ]
  if (nrow(df) == 0) {
    warning(paste("Skipping", group, "- no valid data"))
    next
  }
  
  # Compute log-expression (original metric calculation from your first script)
  df$baseMean_avg <- (df$baseMean_pI + df$baseMean_cI) / 2
  df[[metric]] <- log2(1 + df$baseMean_avg)
  
  group_vectors[[group]] <- setNames(df[[metric]], df$transcript_id)
  transcript_sets[[group]] <- df$transcript_id
}

# -------------------------------
# Construct group matrix
# -------------------------------

all_transcripts <- Reduce(union, transcript_sets)
group_matrix <- do.call(rbind, lapply(group_vectors, function(vec) {
  out <- rep(NA, length(all_transcripts))
  names(out) <- all_transcripts
  out[names(vec)] <- vec
  return(out)
}))
rownames(group_matrix) <- names(group_vectors)
colnames(group_matrix) <- all_transcripts

# Filter low-coverage transcripts
min_presence <- 4
keep_cols <- colSums(!is.na(group_matrix)) >= min_presence
group_matrix <- group_matrix[, keep_cols]

# -------------------------------
# Clean and scale matrix  (NA-safe; keep NA)
# -------------------------------
cm <- colMeans(group_matrix, na.rm = TRUE)
cs <- apply(group_matrix, 2, sd, na.rm = TRUE)
cs[!is.finite(cs) | cs == 0] <- 1
group_matrix_scaled <- sweep(sweep(group_matrix, 2, cm, "-"), 2, cs, "/")

# -------------------------------
# Pairwise-complete cosine
# -------------------------------
cosine_pairwise <- function(A) {
  n <- nrow(A)
  D <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in i:n) {
      ok <- is.finite(A[i, ]) & is.finite(A[j, ])
      if (!any(ok)) {
        d <- 1
      } else {
        ai <- A[i, ok]; aj <- A[j, ok]
        num <- sum(ai * aj)
        denom <- sqrt(sum(ai^2)) * sqrt(sum(aj^2))
        d <- if (denom == 0) 1 else 1 - (num / denom)
      }
      D[i, j] <- D[j, i] <- d
    }
  }
  as.dist(D)
}

cosine_dist <- cosine_pairwise(group_matrix_scaled)
mds_result <- cmdscale(cosine_dist, k = 9)

# -------------------------------
# UMAP
# -------------------------------
umap_config <- umap.defaults
umap_config$min_dist <- 0.22
umap_config$n_neighbors <- min(10, nrow(group_matrix_scaled) - 1)

set.seed(42)
umap_result <- umap(mds_result, config = umap_config)

umap_df <- as.data.frame(umap_result$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$group <- rownames(group_matrix_scaled)

# -------------------------------
# K-means clustering (k = 2)
# -------------------------------
set.seed(42)
k_fixed <- 2
kmeans_result <- kmeans(umap_df[, c("UMAP1", "UMAP2")], centers = k_fixed, nstart = 25)
umap_df$cluster <- factor(kmeans_result$cluster)

# -------------------------------
# Classify cell types
# -------------------------------
hematopoietic <- c("AML2", "AML3", "Molm13", "Jurkat", "Jeko1", "Nalm6")
epithelial <- c("HeLa", "293T", "A549", "Huh7", "DiFi")
mesenchymal_neural <- c("LPS853", "LN308")

umap_df$cell_type <- case_when(
  umap_df$group %in% hematopoietic ~ "Hematopoietic",
  umap_df$group %in% epithelial ~ "Epithelial",
  umap_df$group %in% mesenchymal_neural ~ "Mesenchymal/Neural",
  TRUE ~ "Other"
)

# -------------------------------
# Convex hulls for clusters
# -------------------------------
get_hulls <- function(df, cluster_col, x_col, y_col) {
  df %>%
    group_by(!!sym(cluster_col)) %>%
    group_split() %>%
    lapply(function(group) {
      if (nrow(group) >= 3) {
        hull_idx <- chull(group[[x_col]], group[[y_col]])
        group[hull_idx, ]
      } else {
        NULL
      }
    }) %>%
    bind_rows()
}
hull_df <- get_hulls(umap_df, "cluster", "UMAP1", "UMAP2")

# -------------------------------
# Define custom colors
# -------------------------------
cluster_levels <- levels(umap_df$cluster)
cluster_colors <- c("1" = "#F7A1B2", "2" = "#BCE8DD")
names(cluster_colors) <- cluster_levels

cell_type_colors <- c(
  "Hematopoietic" = "#E01A4F",
  "Epithelial" = "#7DCFB6",
  "Mesenchymal/Neural" = "#F9C22E"
)

# -------------------------------
# UMAP Plot
# -------------------------------
umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_polygon(
    data = hull_df,
    aes(x = UMAP1, y = UMAP2, fill = cluster, group = cluster),
    alpha = 0.2,
    color = NA,
    inherit.aes = FALSE,
    show.legend = FALSE
  ) +
  geom_point(aes(color = cell_type), size = 3) +
  geom_text(
    aes(label = group, color = cell_type),
    vjust = 1.5,
    size = 3,
    show.legend = FALSE
  ) +
  scale_fill_manual(values = cluster_colors) +
  scale_color_manual(values = cell_type_colors) +
  coord_cartesian(clip = "off") +
  labs(
    title = paste0("UMAP of Cell Lines (", metric, ") with Clusters"),
    color = "Cell Type",
    fill = "Cluster"
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
    legend.text  = element_text(size = 8)
  )

ggsave(paste0("umap_cluster_", metric, ".pdf"), plot = umap_plot, width = 4, height = 3)

# -------------------------------
# Heatmap
# -------------------------------
pdf(paste0("heatmap_", metric, ".pdf"), width = 6, height = 6)

# Use cosine for BOTH axes
d_rows <- cosine_pairwise(group_matrix_scaled)
d_cols <- cosine_pairwise(t(group_matrix_scaled))

# Drop any rows/cols that still produce NA/Inf distances
Dr <- as.matrix(d_rows); Dc <- as.matrix(d_cols)
bad_r <- which(!is.finite(rowSums(Dr)))
bad_c <- which(!is.finite(rowSums(Dc)))
if (length(bad_r) > 0 || length(bad_c) > 0) {
  keep_r <- setdiff(seq_len(nrow(group_matrix_scaled)), bad_r)
  keep_c <- setdiff(seq_len(ncol(group_matrix_scaled)), bad_c)
  group_matrix_scaled <- group_matrix_scaled[keep_r, keep_c, drop = FALSE]
  d_rows <- as.dist(Dr[keep_r, keep_r, drop = FALSE])
  d_cols <- as.dist(Dc[keep_c, keep_c, drop = FALSE])
}

pheatmap(group_matrix_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = d_rows,
         clustering_distance_cols = d_cols,
         na_col = "#EEEEEE",
         main = paste("Heatmap of", metric),
         show_rownames = TRUE,
         show_colnames = FALSE)
dev.off()

# -------------------------------
# Save outputs
# -------------------------------
write.csv(umap_df, paste0("umap_coordinates_", metric, ".csv"))

message("✅ UMAP and heatmap saved using Log_Expression metric")
