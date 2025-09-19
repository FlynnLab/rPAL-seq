library(tidyverse)

# ------------------------- SETTINGS -------------------------
k <- 11
B_boot <- 500               # bootstrap resamples for CI (non-hit resampling)
B_perm <- 500               # permutations for null (label shuffling within pos_key)
alpha_fdr <- 0.05
position_window <- 0        # if you want binning of pos → pos_key, set > 0
set.seed(42)

# Labels/colors (we'll use the "total_*" colors for group-only plots)
group_labels <- c(
  "hematopoietic"      = "Hematopoietic",
  "non-hematopoietic"  = "Non-hematopoietic"
)
group_colors <- c(
  "hematopoietic"      = "#E01A4F",
  "non-hematopoietic"  = "#7DCFB6"
)
ci_alpha_plot <- 0.20

if (k %% 2 == 0) stop("Please use an odd k for symmetric centering; got k = ", k)
pos_axis <- -((k-1)/2):((k-1)/2)

# ------------------------- HELPERS --------------------------
valid_kmer <- function(x, k) {
  x <- as.character(x)
  !is.na(x) & nchar(x) == k & !str_detect(x, "[^ACGTacgt]")
}

gc_curve <- function(seqs, k) {
  if (length(seqs) == 0) return(rep(NA_real_, k))
  ss <- strsplit(toupper(seqs), "", fixed = TRUE)
  ok <- lengths(ss) == k
  if (!any(ok)) return(rep(NA_real_, k))
  mat <- matrix(unlist(ss[ok], use.names = FALSE), ncol = k, byrow = TRUE)
  colMeans(mat == "G" | mat == "C")
}

make_pos_key <- function(pos, w = 0) {
  pos <- as.numeric(pos)
  if (w <= 0) return(pos)
  # e.g., floor(pos / w) * w for binning
  pos
}

# Center-trim any sequence to length k (if >= k). Returns NA if too short.
center_kmer <- function(seqs, k) {
  s <- toupper(as.character(seqs))
  n <- nchar(s)
  out <- rep(NA_character_, length(s))
  ok <- !is.na(s) & n >= k
  if (any(ok)) {
    out[ok] <- mapply(function(str, L) {
      start <- floor((L - k)/2) + 1L
      substr(str, start, start + k - 1L)
    }, s[ok], n[ok], USE.NAMES = FALSE)
  }
  out
}

# --------- OBSERVED ΔGC (pos-matched, hit-weighted background) + BOOTSTRAP CI ---------
bootstrap_ci_pos_matched <- function(hits_df, stage_df, k = 15, B = 500, position_window = 0) {
  # Clean
  hits_df  <- hits_df  %>%
    mutate(sequence_context_k = center_kmer(sequence_context, k)) %>%
    filter(valid_kmer(sequence_context_k, k))
  
  stage_df <- stage_df %>%
    mutate(sequence_context_k = center_kmer(sequence_context, k)) %>%
    filter(valid_kmer(sequence_context_k, k))
  
  stopifnot(nrow(hits_df) > 0, nrow(stage_df) > 0)
  
  # pos_key
  hits_df  <- hits_df  %>%
    mutate(pos = as.numeric(pos), pos_key = make_pos_key(pos, w = position_window))
  stage_df <- stage_df %>%
    mutate(pos = as.numeric(pos), pos_key = make_pos_key(pos, w = position_window))
  
  # Non-hits pool = STAGE \ HITS by locus
  non_hits <- stage_df %>% anti_join(hits_df %>% select(chr_pos), by = "chr_pos")
  
  need <- c("chr_pos","sequence_context","pos_key","pos")
  stopifnot(all(need %in% names(hits_df)), all(need %in% names(non_hits)))
  
  # Hit counts & weights by pos_key
  hit_counts <- hits_df %>% count(pos_key, name = "n_hit_pos")
  pools <- non_hits %>% mutate(.rowid = row_number()) %>%
    group_by(pos_key) %>%
    summarise(idx = list(.rowid), .groups = "drop")
  available <- hit_counts %>% inner_join(pools, by = "pos_key")
  hits_keep <- hits_df %>% inner_join(available %>% select(pos_key), by = "pos_key")
  
  n_drop <- nrow(hits_df) - nrow(hits_keep)
  if (n_drop > 0) message("Dropped ", n_drop, " hit(s) from positions with no non-hit background.")
  
  # Observed line: GC(hits) − sum_s w_hit(s) * GC(non-hit | s)
  gc_hit_obs <- gc_curve(hits_keep$sequence_context_k, k)
  
  bg_per_pos <- non_hits %>%
    group_by(pos_key) %>%
    summarise(gc = list(gc_curve(sequence_context_k, k)), .groups = "drop")
  
  w_hit <- hits_keep %>% count(pos_key, name = "n") %>% mutate(w = n / sum(n))
  bg_join <- w_hit %>% inner_join(bg_per_pos, by = "pos_key")
  if (nrow(bg_join) == 0) stop("No overlapping pos strata for background.")
  w_norm <- bg_join$w / sum(bg_join$w)
  gc_bg_obs <- Reduce(`+`, Map(`*`, bg_join$gc, as.list(w_norm)))
  
  dgc_obs <- gc_hit_obs - gc_bg_obs
  
  # Bootstrap CI: resample non-hits within each pos_key, keep hit weights fixed
  avail_map <- available %>% mutate(n_hit_pos = as.integer(n_hit_pos)) %>% filter(n_hit_pos > 0L)
  seq_vec_bg <- non_hits$sequence_context_k
  
  boot_mat <- matrix(NA_real_, nrow = k, ncol = B)
  for (b in seq_len(B)) {
    chosen <- vector("list", nrow(avail_map))
    for (i in seq_len(nrow(avail_map))) {
      ids <- avail_map$idx[[i]]
      n_h <- avail_map$n_hit_pos[i]
      chosen[[i]] <- if (length(ids) >= n_h) sample(ids, n_h, FALSE) else sample(ids, n_h, TRUE)
    }
    gc_bg_list <- vector("list", nrow(avail_map))
    for (i in seq_len(nrow(avail_map))) {
      bg_idx  <- chosen[[i]]
      bg_seqs <- if (length(bg_idx)) seq_vec_bg[bg_idx] else character(0)
      gc_bg_list[[i]] <- gc_curve(bg_seqs, k)
    }
    gc_bg_b <- Reduce(`+`, Map(`*`, gc_bg_list, as.list(w_norm)))
    boot_mat[, b] <- gc_hit_obs - gc_bg_b
  }
  
  ci_lo <- apply(boot_mat, 1, quantile, probs = 0.025, na.rm = TRUE)
  ci_hi <- apply(boot_mat, 1, quantile, probs = 0.975, na.rm = TRUE)
  
  # exact pos_key ↔ weight table actually used (order-independent)
  w_hit_tbl <- bg_join %>% transmute(pos_key, w = w_norm)
  
  list(
    dgc_obs        = dgc_obs,
    ci_lo          = ci_lo,
    ci_hi          = ci_hi,
    n_hits_used    = nrow(hits_keep),
    n_hits_dropped = n_drop,
    w_hit_tbl      = w_hit_tbl,
    avail_map      = avail_map,
    non_hits       = non_hits
  )
}

# --------- PERMUTATION NULL (labels shuffled within pos_key across STAGE) ---------
perm_test_pos <- function(hits_df, stage_df, k = 15, B = 1000,
                          position_window = 0, w_hit_tbl = NULL) {
  # Clean & keys (create trimmed seqs here too)
  hits_df  <- hits_df %>%
    mutate(sequence_context_k = center_kmer(sequence_context, k)) %>%
    filter(valid_kmer(sequence_context_k, k)) %>%
    mutate(pos = as.numeric(pos), pos_key = make_pos_key(pos, w = position_window))
  stage_df <- stage_df %>%
    mutate(sequence_context_k = center_kmer(sequence_context, k)) %>%
    filter(valid_kmer(sequence_context_k, k)) %>%
    mutate(pos = as.numeric(pos), pos_key = make_pos_key(pos, w = position_window))
  stopifnot(nrow(hits_df) > 0, nrow(stage_df) > 0)
  
  # Hit counts by pos_key (observed composition)
  hit_counts <- hits_df %>% count(pos_key, name = "n_hit_pos")
  
  # If bootstrap provided exact pos weights, restrict to them
  if (!is.null(w_hit_tbl)) {
    hit_counts <- inner_join(hit_counts, w_hit_tbl, by = "pos_key")
  } else {
    hit_counts <- hit_counts %>% mutate(w = n_hit_pos / sum(n_hit_pos))
  }
  
  # Stage pools by pos_key
  stage_df <- stage_df %>% mutate(.rowid = row_number())
  pools <- stage_df %>% group_by(pos_key) %>%
    summarise(idx = list(.rowid), .groups = "drop")
  
  # Keep only pos_key that exist in both (and have enough rows for a proper split)
  pools <- inner_join(pools, hit_counts, by = "pos_key") %>%
    mutate(n_pool = lengths(idx), ok = n_pool > n_hit_pos)  # > so complement isn't empty
  
  if (any(!pools$ok)) {
    dropped <- pools %>% filter(!ok) %>% pull(pos_key)
    warning("Dropping pos_key with insufficient pool size in permutation: ",
            paste(head(dropped, 5), collapse = ", "),
            if (length(dropped) > 5) " ..." else "")
    pools <- pools %>% filter(ok)
  }
  stopifnot(nrow(pools) > 0)
  
  # Renormalize weights to the final pos_key set
  pools <- pools %>% mutate(w = w / sum(w)) %>% arrange(pos_key)
  
  seq_vec_all <- stage_df$sequence_context_k
  
  perm_mat <- matrix(NA_real_, nrow = k, ncol = B)
  for (b in seq_len(B)) {
    pseudo_hit_ids <- vector("list", nrow(pools))
    bg_gc_list     <- vector("list", nrow(pools))
    for (i in seq_len(nrow(pools))) {
      ids <- pools$idx[[i]]
      n_h <- pools$n_hit_pos[i]
      h_idx <- sample(ids, n_h, replace = FALSE)
      pseudo_hit_ids[[i]] <- h_idx
      bg_idx <- setdiff(ids, h_idx)  # non-empty by construction
      bg_gc_list[[i]] <- gc_curve(seq_vec_all[bg_idx], k)
    }
    hit_idx <- unlist(pseudo_hit_ids, use.names = FALSE)
    gc_hits_b <- gc_curve(seq_vec_all[hit_idx], k)
    gc_bg_b   <- Reduce(`+`, Map(`*`, bg_gc_list, as.list(pools$w)))
    perm_mat[, b] <- gc_hits_b - gc_bg_b
  }
  
  null_mean <- rowMeans(perm_mat, na.rm = TRUE)
  list(perm_mat = perm_mat, null_mean = null_mean)
}

# ------------------------- PLOTTING (q-value dots) -------------------------
plot_dgc_ci_qdots_pdf <- function(df_sub, file_pdf, pretty_title,
                                  line_color, alpha_fdr,
                                  show_legend = TRUE) {
  x_breaks <- seq(min(df_sub$position, na.rm = TRUE),
                  max(df_sub$position, na.rm = TRUE), by = 1)
  
  # vertical placement for dots
  y_top  <- max(df_sub$ci_hi, na.rm = TRUE)
  y_bot  <- min(df_sub$ci_lo, na.rm = TRUE)
  y_span <- y_top - y_bot
  y_up   <- y_top + 0.08 * y_span
  
  p <- ggplot(df_sub, aes(x = position)) +
    geom_ribbon(aes(ymin = ci_lo, ymax = ci_hi),
                fill = line_color, alpha = ci_alpha_plot, color = NA) +
    geom_line(aes(y = dgc_obs), color = line_color, linewidth = 1.2) +
    geom_point(aes(y = y_up, size = neglog10q, fill = q < alpha_fdr),
               shape = 21, stroke = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray55") +
    scale_x_continuous(breaks = x_breaks) +
    scale_size_continuous(name = expression(-log[10](padj)),
                          range = c(0.8, 3.2)) +
    scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white"),
                      name = paste0("padj < ", alpha_fdr)) +
    labs(
      title    = paste0("Delta-GC (position-matched), ", pretty_title),
      subtitle = paste0("95% Bootstrap CI; Permutation-derived padj (BH)"),
      x        = "Position",
      y        = "Delta-GC (Hit - Background)"
    ) +
    guides(fill = if (show_legend) "legend" else "none",
           size = if (show_legend) "legend" else "none") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid   = element_blank(),
      axis.line    = element_line(color = "black"),
      axis.ticks   = element_line(color = "black"),
      plot.title = element_text(face = "bold", size = 11, margin = margin(b = 6)),
      plot.title.position = "plot",
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold", size = 9),
      axis.text  = element_text(size = 8),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    ) +
    expand_limits(y = y_up + 0.06 * y_span)
  
  ggsave(file_pdf, p, width = 5, height = 3)
}

# ------------------------- IO: READ & COMBINE ---------------------
hits_files  <- sort(list.files(pattern = "^Annotated_hits_(sgid|total)_.+\\.csv$"))
stage_files <- sort(list.files(pattern = "^Annotated_stage[12]_(variant_rate|sgidshare)_.+\\.csv$"))

map_stage_metric <- function(x) recode(x, sgidshare = "sgid", variant_rate = "total")

hits_info <- tibble(file_hits = hits_files) %>%
  mutate(metric = str_match(file_hits, "^Annotated_hits_(sgid|total)_")[,2],
         group  = str_match(file_hits, "^Annotated_hits_(?:sgid|total)_(.+)\\.csv$")[,2])

stage_info <- tibble(file_stage = stage_files) %>%
  mutate(stage_metric = str_match(file_stage, "^Annotated_stage[12]_(variant_rate|sgidshare)_")[,2],
         metric = map_stage_metric(stage_metric),
         group  = file_stage %>%
           str_remove("^Annotated_stage[12]_(variant_rate|sgidshare)_") %>%
           str_remove("_FAST\\.csv$") %>%
           str_remove("\\.csv$"))

# Groups present in both hits and stage
groups <- sort(intersect(unique(hits_info$group), unique(stage_info$group)))
stopifnot(length(groups) > 0)

# ------------------------- RUN (per GROUP) --------------------------
all_out <- list()

for (grp in groups) {
  pretty <- if (grp %in% names(group_labels)) group_labels[[grp]] else grp
  col    <- if (grp %in% names(group_colors)) group_colors[[grp]] else "#333333"
  
  message("ΔGC (pos-matched across ALL non-hits) for group: ", grp)
  
  # --- Join sgid + total HITS for this group, dedupe by chr_pos
  hits_files_grp <- hits_info %>% filter(group == grp) %>% pull(file_hits)
  stopifnot(length(hits_files_grp) > 0)
  hits_list <- lapply(hits_files_grp, function(f) readr::read_csv(f, show_col_types = FALSE))
  hits <- bind_rows(hits_list) %>%
    mutate(pos = as.numeric(pos)) %>%
    distinct(chr_pos, .keep_all = TRUE)
  
  # --- Join all STAGE files for this group, dedupe by chr_pos
  stage_files_grp <- stage_info %>% filter(group == grp) %>% pull(file_stage)
  stopifnot(length(stage_files_grp) > 0)
  stage_list <- lapply(stage_files_grp, function(f) readr::read_csv(f, show_col_types = FALSE))
  stage <- bind_rows(stage_list) %>%
    mutate(pos = as.numeric(pos)) %>%
    distinct(chr_pos, .keep_all = TRUE)
  
  # raw CSVs should have sequence_context, not sequence_context_k yet
  need <- c("chr_pos","sequence_context","biotype","family","ref_base","pos")
  stopifnot(all(need %in% names(hits)), all(need %in% names(stage)))
  
  # 1) Observed + bootstrap CI
  ci <- bootstrap_ci_pos_matched(
    hits_df = hits, stage_df = stage,
    k = k, B = B_boot, position_window = position_window
  )
  
  # 2) Permutation null uses the same pos set + weights as the bootstrap
  perm <- perm_test_pos(
    hits_df = hits, stage_df = stage,
    k = k, B = B_perm, position_window = position_window,
    w_hit_tbl = ci$w_hit_tbl
  )
  
  # per-position p (for BH) → q
  dgc_obs <- ci$dgc_obs
  dev_obs  <- abs(dgc_obs - perm$null_mean)
  dev_perm <- abs(sweep(perm$perm_mat, 1, perm$null_mean, FUN = "-"))
  
  pvals <- vapply(seq_len(k), function(j) {
    (1 + sum(dev_perm[j, ] >= dev_obs[j], na.rm = TRUE)) / (B_perm + 1)
  }, numeric(1))
  
  q <- p.adjust(pvals, method = "BH")
  neglog10q <- -log10(pmax(q, .Machine$double.eps))
  
  df <- tibble(
    position   = pos_axis,
    dgc_obs    = dgc_obs,
    ci_lo      = ci$ci_lo,
    ci_hi      = ci$ci_hi,
    q          = q,
    neglog10q  = neglog10q,
    group      = grp
  )
  
  out_pdf <- paste0("delta_gc_posmatched_", grp, "_joinedHits_withQ.pdf")
  plot_dgc_ci_qdots_pdf(df, out_pdf, pretty, col, alpha_fdr, show_legend = TRUE)
  
  message("  wrote: ", out_pdf, " (hits used: ", ci$n_hits_used,
          if (ci$n_hits_dropped > 0) paste0(", dropped: ", ci$n_hits_dropped) else "", ")")
  
  all_out[[length(all_out) + 1]] <- df
}

# ------------------------- WRITE TABLE ----------------------
all_res <- bind_rows(all_out)
readr::write_csv(all_res, "delta_gc_posmatched_joinedHits_withQ_results.csv")
message("Saved: delta_gc_posmatched_joinedHits_withQ_results.csv")
