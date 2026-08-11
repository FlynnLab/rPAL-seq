library(tidyverse)

# ---------- user knobs ----------
annot_files <- c(
  "Annotated_hits_sgid_hematopoietic.csv",
  "Annotated_hits_sgid_non-hematopoietic.csv",
  "Annotated_hits_total_hematopoietic.csv",
  "Annotated_hits_total_non-hematopoietic.csv"
)
counts_file      <- "/path/to/counts_long.csv"                # MUST contain 'sample'
sample_info_file <- "/path/to/metadata.csv"
out_prefix       <- "conversions_delta"

# Group membership (match your step0)
hematopoietic     <- c("AML2","AML3","Molm13","Jurkat","Jeko1","Nalm6")
non_hematopoietic <- c("HeLa","293T","A549","Huh7","DiFi","LPS853","LN308")

# Analysis knobs
include_sgid_in_denominator <- TRUE  # denom = A+C+G+T (+ SGID if TRUE)
clip_negative_delta         <- TRUE  # show only enrichment if TRUE

# Aesthetics — match your volcano colors for bases; SGID in grey
base_colors <- c(A = "#7DCFB6", C = "#53B3CB", G = "#F9C22E", T = "#E01A4F")
fill_colors <- c(base_colors, SGID = "grey60")

# ---------------------------------

# --- helpers ---
num0 <- function(x) replace_na(suppressWarnings(as.numeric(x)), 0)
parse_condition <- function(x) {
  suf <- stringr::str_sub(x, -1)
  dplyr::case_when(
    suf %in% c("i","I") ~ "Input",
    suf %in% c("p","P") ~ "IP",
    suf %in% c("h","H") ~ "Control",
    TRUE ~ NA_character_
  )
}
base_id <- function(x) stringr::str_replace(x, "[iIhHpP]$", "")

# 1) Annotated hits (chr_pos, ref_base)
ann <- lapply(annot_files, function(f) {
  if (!file.exists(f)) stop("File not found: ", f)
  readr::read_csv(f, show_col_types = FALSE) |>
    transmute(chr_pos = as.character(chr_pos),
              ref_base = toupper(as.character(ref_base)))
}) |> bind_rows() |> distinct(chr_pos, .keep_all = TRUE)

# 2) Sample metadata → super_group
if (!file.exists(sample_info_file)) stop("File not found: ", sample_info_file)
sample_info <- readr::read_csv(sample_info_file, show_col_types = FALSE) |>
  mutate(across(c(sample_name, group_name, replicates), as.character)) |>
  mutate(super_group = case_when(
    group_name %in% hematopoietic     ~ "hematopoietic",
    group_name %in% non_hematopoietic ~ "non-hematopoietic",
    TRUE ~ NA_character_
  )) |>
  select(sample_name, super_group)

# 3) counts_long (per-sample)
counts_needed <- c("chr_pos","sample","REF","SGID","A_star","C_star","G_star","T_star")
counts <- readr::read_csv(counts_file, show_col_types = FALSE)
missing <- setdiff(counts_needed, names(counts))
if (length(missing)) stop("counts_long.csv missing: ", paste(missing, collapse = ", "))

counts <- counts |>
  transmute(chr_pos = as.character(chr_pos),
            sample  = as.character(sample),
            REF     = num0(REF),
            SGID    = num0(SGID),
            A_star  = num0(A_star),
            C_star  = num0(C_star),
            G_star  = num0(G_star),
            T_star  = num0(T_star))

# 4) Attach condition + super_group; keep IP & Input
counts_meta <- counts |>
  mutate(condition = parse_condition(sample),
         sample_id = base_id(sample)) |>
  left_join(sample_info |> rename(sample_id = sample_name), by = "sample_id") |>
  filter(condition %in% c("IP","Input"), !is.na(super_group))
if (!nrow(counts_meta)) stop("No IP/Input samples with known super_group.")

# 5) Keep annotated hits
merged <- inner_join(ann, counts_meta, by = "chr_pos")
if (!nrow(merged)) stop("No overlap between annotated hits and counts_long.csv")

# 6) Off-diagonal-only table (drop REF)
obs <- merged %>%
  mutate(
    A    = if_else(ref_base == "A", 0, A_star),
    T    = if_else(ref_base == "T", 0, T_star),
    C    = if_else(ref_base == "C", 0, C_star),
    G    = if_else(ref_base == "G", 0, G_star),
    SGID = SGID
  ) %>%
  select(super_group, condition, ref_base, A, T, C, G, SGID)

# 7) Aggregate by condition (pool sites & reps per super_group × ref_base)
agg_cond <- obs %>%
  group_by(super_group, ref_base, condition) %>%
  summarise(across(c(A,T,C,G,SGID), ~sum(.x, na.rm = TRUE)), .groups = "drop")

# 8) Fractions within each (super_group × ref_base × condition)
frac_cond <- agg_cond %>%
  mutate(total = if (include_sgid_in_denominator) A+T+C+G+SGID else A+T+C+G) %>%
  mutate(across(c(A,T,C,G,SGID), ~ if_else(total > 0, .x/total, NA_real_))) %>%
  select(super_group, ref_base, condition, A, T, C, G, SGID)

# 9) Δfraction = IP − Input
delta_tbl <- frac_cond %>%
  pivot_longer(c(A,T,C,G,SGID), names_to = "observed", values_to = "fraction") %>%
  pivot_wider(names_from = condition, values_from = fraction) %>%
  mutate(delta = IP - Input,
         delta = if (clip_negative_delta) pmax(delta, 0) else delta) %>%
  select(super_group, ref_base, observed,
         fraction_IP = IP, fraction_Input = Input, delta)

# Save numbers
readr::write_csv(delta_tbl, paste0(out_prefix, "_delta_aggregated.csv"))
message("Wrote: ", paste0(out_prefix, "_delta_aggregated.csv"))

# ---------------- aesthetics & plotting ----------------
plot_df <- delta_tbl %>%
  mutate(
    ref_base = factor(ref_base, levels = c("A","T","C","G")),
    observed = factor(observed,
                      levels = c("A","C","G","T","SGID"),
                      labels = c("A","C","G","T","Skip/Gap/Indel"))
  )

# legend order/levels (include SGID even if absent in data)
fill_levels  <- levels(plot_df$observed)
fill_values  <- c("A"="#7DCFB6", "C"="#53B3CB", "G"="#F9C22E", "T"="#E01A4F",
                  "Skip/Gap/Indel"="grey60")

# knob: share y-axis between panels? (FALSE = each panel auto scales)
share_y_limits <- FALSE

# optional shared ymax if requested
if (isTRUE(share_y_limits)) {
  ymax_global <- plot_df %>%
    group_by(super_group, ref_base) %>%
    summarise(total = sum(delta, na.rm = TRUE), .groups = "drop_last") %>%
    summarise(max_total = max(total, na.rm = TRUE), .groups = "drop") %>%
    summarise(global_max = max(max_total, na.rm = TRUE)) %>% pull(global_max)
  if (!is.finite(ymax_global) || ymax_global <= 0) ymax_global <- 0.1
}

make_panel <- function(sg_label, df, outfile_pdf) {
  dd <- df %>% filter(super_group == sg_label)
  if (!nrow(dd)) { message("No data for ", sg_label); return(invisible(NULL)) }
  
  # per-panel ymax (used when share_y_limits == FALSE)
  if (!isTRUE(share_y_limits)) {
    ymax_local <- dd %>%
      group_by(ref_base) %>%
      summarise(total = sum(delta, na.rm = TRUE), .groups = "drop") %>%
      summarise(max_total = max(total, na.rm = TRUE)) %>% pull(max_total)
    if (!is.finite(ymax_local) || ymax_local <= 0) ymax_local <- 0.1
  }
  
  p <- ggplot(dd, aes(x = ref_base, y = delta, fill = observed)) +
    geom_col(width = 0.85) +
    scale_fill_manual(
      values = fill_values,
      breaks = fill_levels,
      limits = fill_levels,
      drop   = FALSE        # <-- keep SGID in legend even if not present
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = if (isTRUE(share_y_limits))
        c(0, ymax_global * 1.05) else c(0, ymax_local * 1.05),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = paste0("Delta-conversion (IP - Input):\n ", sg_label),
      x = "Reference base",
      y = if (clip_negative_delta) "Delta-fraction enriched in IP" else "Delta-fraction",
      fill = "Observed base"
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
  
  ggsave(outfile_pdf, p, width = 3, height = 2.5, units = "in", dpi = 300)
  message("Saved: ", outfile_pdf)
}

make_panel("hematopoietic",     plot_df, paste0(out_prefix, "_hematopoietic.pdf"))
make_panel("non-hematopoietic", plot_df, paste0(out_prefix, "_non-hematopoietic.pdf"))

# Console quick-look: dominant enriched mismatch per ref base
top1 <- delta_tbl %>%
  group_by(super_group, ref_base) %>%
  slice_max(delta, n = 1, with_ties = FALSE) %>%
  ungroup() %>% arrange(super_group, ref_base)
print(top1)

