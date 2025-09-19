 library(tidyverse)
 library(readr)
 library(stringr)
 library(ggrepel)

# ======================= USER SETTINGS =======================
hit_files <- c(
  "Annotated_hits_sgid_hematopoietic.csv",
  "Annotated_hits_sgid_non-hematopoietic.csv",
  "Annotated_hits_total_hematopoietic.csv",
  "Annotated_hits_total_non-hematopoietic.csv"
)

# ---- selection mode ----
selection_mode <- "family"           # "family" or "product"
family_query   <- "tRNA-Ser-CGA"     # case-insensitive substring (used when selection_mode=="family")
product_query  <- "tRNA-Ile-AAT-1"                 # case-insensitive substring (used when selection_mode=="product")

# score = logOR * -log10(padj)
min_p       <- 1e-300
score_floor <- 0                     # clamp to >= 0; use -Inf to allow negatives
winsor_pct  <- 0                     # e.g. 0.01 to tame outliers; 0 disables

# colors (metrics)
metric_colors <- c(total = "#E01A4F", sgidshare = "#4B3F72")

# labeling
label_all   <- FALSE                 # now default to top-N only
label_top_n_total   <- 6           # top-6 for Variant (total)
label_top_n_sgid    <- 3           # top-3 for SGID
label_size  <- 3.2

# geometry
dodge_width <- 0.28                  # horizontal offset between metrics at same position

# x-axis (full-length control)
length_override <- NA_integer_       # set to an integer to force length; NA => auto-detect

# output (PDFs; one per super-group)
out_prefix  <- "family_spikes"
# =============================================================

# ---------------- helpers ----------------
metric_from_file <- function(f) {
  b <- basename(f)
  if (str_detect(b, "^Annotated_hits_total_"))      return("total")
  if (str_detect(b, "^Annotated_hits_sgid_"))       return("sgidshare")
  cols <- names(suppressMessages(readr::read_csv(f, n_max = 0, show_col_types = FALSE)))
  if (all(c("logOR_total","padj_total") %in% cols))         return("total")
  if (all(c("logOR_sgidshare","padj_sgidshare") %in% cols)) return("sgidshare")
  stop("Can't infer metric for file: ", f)
}
group_from_file <- function(f) {
  b <- basename(f)
  if (str_detect(b, "_hematopoietic\\.csv$"))     return("hematopoietic")
  if (str_detect(b, "_non-hematopoietic\\.csv$")) return("non-hematopoietic")
  NA_character_
}
winsorize <- function(x, p = 0) {
  if (p <= 0) return(x)
  qs <- stats::quantile(x, c(p, 1-p), na.rm = TRUE)
  pmin(pmax(x, qs[1]), qs[2])
}
nice_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

choose_ref_base <- function(x) {
  u <- toupper(na.omit(x))
  if (!length(u)) return(NA_character_)
  tb <- sort(table(u), decreasing = TRUE)
  names(tb)[1]
}

read_hits_scored <- function(file, min_p = 1e-300, score_floor = 0) {
  met <- metric_from_file(file)
  grp <- group_from_file(file)
  df  <- readr::read_csv(file, show_col_types = FALSE)
  
  # columns for position & family/product & optional length
  pos_col <- intersect(c("pos","pos_in_tx","position","tx_pos"), names(df))
  len_col <- intersect(c("seq_length","tx_length","length","len","tx_len","product_length"), names(df))
  if (length(pos_col) == 0) stop("No position column (pos/pos_in_tx/position/tx_pos) in: ", file)
  if (!"family" %in% names(df)) stop("No 'family' column in: ", file)
  if (!"ref_base" %in% names(df)) stop("No 'ref_base' column in: ", file)
  
  # metric columns
  if (met == "total") {
    lor_col <- "logOR_total";      p_col <- "padj_total"
  } else {
    lor_col <- "logOR_sgidshare";  p_col <- "padj_sgidshare"
  }
  if (!all(c(lor_col, p_col) %in% names(df)))
    stop("Missing metric columns in ", file, " (need ", lor_col, ", ", p_col, ").")
  
  lor  <- suppressWarnings(as.numeric(df[[lor_col]]))
  padj <- pmax(suppressWarnings(as.numeric(df[[p_col]])), min_p)
  score_raw <- lor * (-log10(padj))
  score <- if (is.finite(score_floor)) pmax(score_raw, score_floor) else score_raw
  
  tx_len <- if (length(len_col)) suppressWarnings(as.integer(df[[len_col[1]]])) else NA_integer_
  
  tibble(
    group     = grp,
    metric    = met,
    family    = as.character(df$family),
    product   = dplyr::coalesce(df$product, df$transcript_id),
    ref_base  = toupper(as.character(df$ref_base)),
    pos       = suppressWarnings(as.integer(df[[pos_col[1]]])),
    tx_len    = tx_len,
    score     = score
  ) %>%
    filter(is.finite(score), !is.na(pos), !is.na(family), !is.na(ref_base))
}

detect_length <- function(df, override = NA_integer_) {
  if (is.finite(override)) return(as.integer(override))
  # use the most common non-NA length if available; else max(pos)
  if ("tx_len" %in% names(df) && any(is.finite(df$tx_len))) {
    m <- stats::na.omit(df$tx_len)
    # take the median of the mode(s) to be robust
    tab <- sort(table(m), decreasing = TRUE)
    len <- as.integer(names(tab)[1])
    return(max(len, max(df$pos, na.rm = TRUE)))
  }
  max(df$pos, na.rm = TRUE)
}
# -----------------------------------------

# -------- load & combine all four files --------
all_hits <- purrr::map_dfr(hit_files, read_hits_scored,
                           min_p = min_p, score_floor = score_floor)

# ---- selection filter ----
if (selection_mode == "family") {
  hits_sel <- all_hits %>%
    filter(str_detect(tolower(family), tolower(family_query)))
  target_label <- family_query
} else if (selection_mode == "product") {
  if (!nzchar(product_query)) stop("selection_mode=='product' but product_query is empty.")
  hits_sel <- all_hits %>%
    filter(str_detect(tolower(product), tolower(product_query)))
  target_label <- product_query
} else {
  stop("selection_mode must be 'family' or 'product'.")
}

if (!nrow(hits_sel)) stop("No rows matched the selection (", selection_mode, ").")

# optional winsorization
if (winsor_pct > 0) hits_sel <- hits_sel %>% mutate(score = winsorize(score, winsor_pct))

# detect full length for x-axis (shared across both groups)
full_len <- detect_length(hits_sel, override = length_override)
x_min <- 1L; x_max <- as.integer(full_len)

# average across members at the SAME position; pick majority ref_base
spikes <- hits_sel %>%
  group_by(group, metric, pos) %>%
  summarise(
    score     = mean(score, na.rm = TRUE),
    ref_base  = choose_ref_base(ref_base),
    n_members = dplyr::n(),
    .groups   = "drop"
  ) %>%
  mutate(label = ifelse(is.na(ref_base), paste0("N", pos), paste0(ref_base, pos))) %>%
  arrange(group, metric, pos)

# ------------- plotting (PDF; one per group) -------------
plot_group <- function(gid, df, label_txt) {
  dd <- df %>% filter(group == gid)
  if (!nrow(dd)) { message("No data for group: ", gid); return(invisible(NULL)) }
  
  # ensure both legend levels exist, even if one metric is absent
  metric_levels <- names(metric_colors)
  dd <- dd %>%
    mutate(
      metric = factor(metric, levels = metric_levels),
      # offset but clamp into [x_min, x_max] so axis remains full length
      x_raw = pos + ifelse(metric == "total", -dodge_width/2, dodge_width/2),
      x = pmin(pmax(x_raw, x_min), x_max)
    )
  
  # seed layer so both legend entries appear even if one metric is absent
  legend_seed <- tibble(
    metric = factor(metric_levels, levels = metric_levels),
    x = x_min, xend = x_min, y = 0, yend = 0
  )
  
  # choose which labels to draw: top-6 for total, top-3 for sgidshare (per group)
  if (label_all) {
    lab_df <- dd
  } else {
    lab_total <- dd %>% filter(metric == "total")     %>%
      slice_max(score, n = label_top_n_total, with_ties = FALSE)
    lab_sgid  <- dd %>% filter(metric == "sgidshare") %>%
      slice_max(score, n = label_top_n_sgid,  with_ties = FALSE)
    lab_df <- bind_rows(lab_total, lab_sgid)
  }
  
  
  y_max <- max(dd$score, na.rm = TRUE); if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  lab_nudge <- 0.02 * y_max
  
  p <- ggplot() +
    # main spikes
    geom_segment(
      data = dd,
      aes(x = x, xend = x, y = 0, yend = score, color = metric),
      linewidth = 0.7, lineend = "butt", alpha = 0.95
    ) +
    # seed layer so both legend entries appear
    geom_segment(
      data = legend_seed,
      aes(x = x, xend = xend, y = y, yend = yend, color = metric),
      linewidth = 1.2, inherit.aes = FALSE, show.legend = TRUE
    ) +
    scale_color_manual(
      values = metric_colors,
      breaks = metric_levels,
      limits = metric_levels,
      drop   = FALSE,
      labels = c(total = "Variant ratio hit", sgidshare = "SGID share hit")
    ) +
    scale_x_continuous(limits = c(x_min, x_max), expand = expansion(mult = c(0, 0.01))) +
    coord_cartesian(ylim = c(0, y_max * 1.08)) +
    labs(
      title  = paste0(label_txt, " :\n ", gid),
      x = "Position",
      y = "Score = logOR × -log10(padj)",
      color = "Metric"
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
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.margin = margin(6, 6, 6, 6)
    ) +
    ggrepel::geom_text_repel(
      data = lab_df,
      aes(x = x, y = score, label = label, color = metric),
      size = label_size, show.legend = FALSE,
      min.segment.length = 0, segment.size = 0.25,
      box.padding = 0.25, point.padding = 0.15, max.overlaps = Inf,
      nudge_y = lab_nudge
    )
  
  outfile <- paste0(out_prefix, "_", nice_name(label_txt), "_", gid, ".pdf")
  ggsave(outfile, p, width = 2.5, height = 2.2, units = "in", dpi = 300)
  message("Saved: ", outfile)
}

target_pretty <- target_label
purrr::walk(c("hematopoietic", "non-hematopoietic"),
            ~ plot_group(.x, spikes, target_pretty))

# also write the collapsed numbers (handy for QC)
out_csv <- paste0(out_prefix, "_", nice_name(target_pretty), "_spikes.csv")
readr::write_csv(spikes %>% arrange(group, metric, pos), out_csv)
message("Wrote: ", out_csv)
