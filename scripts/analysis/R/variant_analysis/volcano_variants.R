library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggrastr)

# ---------- user knobs ----------
alpha_total <- 0.05
alpha_sgid  <- 0.05

min_logOR_total <- log(2.0)
min_logOR_sgid  <- log(1.25)

label_top_n <- 3
super_groups <- c("hematopoietic","non-hematopoietic")

# Colors for ref bases (hits). Non-hits use gray via na.value.
base_colors <- c(A = "#7DCFB6", C = "#53B3CB", G = "#F9C22E", T = "#E01A4F")
# --------------------------------

safe_neglog10 <- function(p) -log10(pmax(p, 1e-300))

volcano_with_annotation <- function(df, effect_col, padj_col,
                                    annotated_hits_file, outfile,
                                    title, min_logOR, alpha) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  
  # Core metrics + hit flag
  df <- df %>%
    filter(!is.na(.data[[padj_col]])) %>%
    mutate(
      neglogq = safe_neglog10(.data[[padj_col]]),
      is_hit  = .data[[padj_col]] < alpha & .data[[effect_col]] > min_logOR
    )
  
  # Bring in ref_base + label (product_pos) for HITS only
  ann <- NULL
  if (!is.null(annotated_hits_file) && file.exists(annotated_hits_file)) {
    ann <- readr::read_csv(annotated_hits_file, show_col_types = FALSE) %>%
      mutate(
        ref_base = toupper(ref_base),
        product  = dplyr::coalesce(product, transcript_id),
        label    = paste0(product, "_", pos)
      ) %>%
      select(chr_pos, ref_base, label, any_of(c(padj_col, effect_col)))
  } else {
    message("Annotated hits file not found: ", annotated_hits_file,
            " — proceeding without base coloring or labels.")
  }
  
  df <- df %>%
    left_join(ann %>% select(chr_pos, ref_base, label), by = "chr_pos") %>%
    mutate(
      plot_color = ifelse(is_hit, ref_base, NA_character_)
    )
  
  df <- df %>% mutate(plot_color = ifelse(is_hit, ref_base, NA_character_))
  
  p <- ggplot(df, aes(x = .data[[effect_col]], y = neglogq)) +
    ggrastr::geom_point_rast(
      aes(color = plot_color, alpha = is_hit),
      size = 2.5,
      raster.dpi = 300,
      dev = "ragg_png",   # or dev = "ragg"
      na.rm = TRUE
    ) +
    scale_color_manual(
      values = base_colors,           # A/C/G/T colors
      na.value = "gray80",            # grey for non-hits
      breaks = names(base_colors),    # legend shows only A/C/G/T
      drop = TRUE                     # drop unused bases from legend
    ) +
    scale_alpha_manual(values = c(`TRUE` = 0.9, `FALSE` = 0.35), guide = "none") +
    geom_vline(xintercept = 0, linetype = "solid", color = "grey40", size = 0.6) +
    { if (min_logOR > 0) geom_vline(xintercept = c(-min_logOR, min_logOR),
                                    linetype = "dashed", color = "grey40", size = 0.6) } +
    geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "grey40", size = 0.6) +
    labs(title = title,
         x = "log-odds (IP - Input)",
         y = expression(-log[10](padj)),
         color = "Ref base") +
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
    )
  
  # Label top N hits by product_pos (rank by padj)
  if (!is.null(ann) && nrow(ann) > 0 && label_top_n > 0) {
    top_lab <- ann %>%
      filter(!is.na(.data[[padj_col]])) %>%
      arrange(.data[[padj_col]]) %>%
      slice_head(n = label_top_n) %>%
      select(chr_pos, label)
    
    if (nrow(top_lab) > 0) {
      lab_points <- df %>% semi_join(top_lab, by = "chr_pos")
      p <- p + ggrepel::geom_text_repel(
        data = lab_points,
        aes(label = label),
        size = 3, max.overlaps = 10, box.padding = 0.8,
        point.padding = 0.4, force = 50, force_pull = 0.02,
        #nudge_x = -1, nudge_y = 1,
        segment.color = "grey40", segment.size = 0.3,
        min.segment.length = 0.05, show.legend = FALSE
      )
    }
  }
  
  ggsave(outfile, p, width = 4, height = 3.2, units = "in", dpi = 300)
  message("Saved: ", outfile)
  invisible(NULL)
}

# ---- stage 1: total variants (REF vs non-REF) ----
plot_stage1 <- function(sg) {
  in_df <- paste0("stage1_variant_rate_", sg, "_FAST.csv")
  if (!file.exists(in_df)) { message("Skipping (not found): ", in_df); return(invisible(NULL)) }
  df <- readr::read_csv(in_df, show_col_types = FALSE)
  
  ann_file <- paste0("Annotated_hits_total_", sg, ".csv")
  out_pdf  <- paste0("volcano_total_", sg, "_annot.pdf")
  
  volcano_with_annotation(
    df,
    effect_col = "logOR_total",
    padj_col   = "padj_total",
    annotated_hits_file = ann_file,
    outfile = out_pdf,
    title   = paste0("Total variants: ", sg),
    min_logOR = min_logOR_total,
    alpha     = alpha_total
  )
}

# ---- stage 2: SGID share among variants ----
plot_stage2_sgid <- function(sg) {
  in_df <- paste0("stage2_sgidshare_", sg, "_FAST.csv")
  if (!file.exists(in_df)) { message("Skipping (not found): ", in_df); return(invisible(NULL)) }
  df <- readr::read_csv(in_df, show_col_types = FALSE)
  
  ann_file <- paste0("Annotated_hits_sgid_", sg, ".csv")
  out_pdf  <- paste0("volcano_sgid_", sg, "_annot.pdf")
  
  volcano_with_annotation(
    df,
    effect_col = "logOR_sgidshare",
    padj_col   = "padj_sgidshare",
    annotated_hits_file = ann_file,
    outfile = out_pdf,
    title   = paste0("Skip/Gap/Indel%: ", sg),
    min_logOR = min_logOR_sgid,
    alpha     = alpha_sgid
  )
}

# ---------------- run for each super-group ----------------
for (sg in super_groups) {
  message("\n=== ", sg, " ===")
  plot_stage1(sg)
  plot_stage2_sgid(sg)
}

message("\n✅ Volcano plotting complete (hits colored by ref base; top 3 labeled).")
