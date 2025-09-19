library(tidyverse)
library(ggseqlogo)

# ------------------------- SETTINGS -------------------------
k <- 11                     # odd k for centered kmers
epsilon <- 1e-6             # pseudocount for log-odds stability
position_window <- 0        # if you later want pos binning, set > 0
set.seed(42)

# Pretty labels (merged by GROUP only)
group_labels <- c(
  "hematopoietic"      = "Hematopoietic",
  "non-hematopoietic"  = "Non-hematopoietic"
)

if (k %% 2 == 0) stop("Please use an odd k for symmetric centering; got k = ", k)
pos_axis <- -((k-1)/2):((k-1)/2)

# ------------------------- HELPERS --------------------------
valid_kmer <- function(x, k) {
  x <- as.character(x)
  !is.na(x) & nchar(x) == k & !str_detect(x, "[^ACGTacgt]")
}

# center-trim any sequence to length k (if >= k). Returns NA if too short.
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

make_pos_key <- function(pos, w = 0) {
  pos <- as.numeric(pos)
  if (w <= 0) return(pos)
  # e.g., floor(pos / w) * w  (left here in case you want binning later)
  pos
}

pwm_from_seqs <- function(seqs, k) {
  # returns 4 x k matrix of base probabilities
  if (length(seqs) == 0) return(matrix(0,4,k, dimnames=list(c("A","C","G","T"), NULL)))
  mat <- matrix(0,4,k, dimnames=list(c("A","C","G","T"), NULL))
  for (i in seq_len(k)) {
    b <- substr(seqs, i, i)
    counts <- table(factor(toupper(b), levels = rownames(mat)))
    mat[, i] <- counts / sum(counts)
  }
  mat
}

logodds_pwm <- function(p_hit, p_bg, eps=1e-6) log2((p_hit + eps)/(p_bg + eps))

gc_curve <- function(seqs, k) {
  if (length(seqs) == 0) return(rep(NA_real_, k))
  ss <- strsplit(toupper(seqs), "", fixed = TRUE)
  ok <- lengths(ss) == k
  if (!any(ok)) return(rep(NA_real_, k))
  m <- matrix(unlist(ss[ok], use.names = FALSE), ncol = k, byrow = TRUE)
  colMeans(m == "G" | m == "C")
}

sanitize_logo_mat <- function(mat) {
  stopifnot(is.matrix(mat), ncol(mat) > 0)
  want_rows <- c("A","C","G","T")
  if (is.null(rownames(mat))) rownames(mat) <- want_rows
  mat <- mat[want_rows, , drop = FALSE]
  if (is.null(colnames(mat))) {
    colnames(mat) <- as.character(seq_len(ncol(mat)))
  } else {
    posnum <- suppressWarnings(as.numeric(colnames(mat)))
    if (any(is.na(posnum)) || length(posnum) != ncol(mat)) {
      colnames(mat) <- as.character(seq_len(ncol(mat)))
    }
  }
  mat[!is.finite(mat)] <- 0
  if (max(abs(mat)) == 0) mat["A", 1] <- 1e-12  # tiny spike so ggseqlogo doesn’t error
  
  ## relabel T as U for plotting
  rownames(mat)[rownames(mat) == "T"] <- "U"
  
  mat
}


plot_logo_pdf <- function(mat, file_pdf, title_main, subtitle = NULL) {
  mat <- sanitize_logo_mat(mat)
  p <- ggseqlogo(mat, method = "custom") +
    scale_x_continuous(breaks = 1:ncol(mat), labels = pos_axis[seq_len(ncol(mat))]) +
    labs(
      title    = title_main,
      x = "Position",
      y = "Log2 enrichment (hits vs pos-matched bg)"
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
  
  ggsave(file_pdf, p, width = 4, height = 3)
}

# ------------------------- IO: READ & MERGE -------------------------
hits_files  <- sort(list.files(pattern = "^Annotated_hits_(sgid|total)_.+\\.csv$"))
stage_files <- sort(list.files(pattern = "^Annotated_stage[12]_(variant_rate|sgidshare)_.+\\.csv$"))

hits_info <- tibble(file_hits = hits_files) %>%
  mutate(metric = str_match(file_hits, "^Annotated_hits_(sgid|total)_")[,2],
         group  = str_match(file_hits, "^Annotated_hits_(?:sgid|total)_(.+)\\.csv$")[,2])

stage_info <- tibble(file_stage = stage_files) %>%
  mutate(stage_metric = str_match(file_stage, "^Annotated_stage[12]_(variant_rate|sgidshare)_")[,2],
         group  = file_stage %>%
           str_remove("^Annotated_stage[12]_(variant_rate|sgidshare)_") %>%
           str_remove("_FAST\\.csv$") %>%
           str_remove("\\.csv$"))

groups <- sort(intersect(unique(hits_info$group), unique(stage_info$group)))
stopifnot(length(groups) > 0)

# --------------------- CORE: POS-MATCHED BACKGROUND ---------------------
build_pos_matched_bg <- function(hits_df, stage_df, k, position_window = 0) {
  # Clean + derive centered kmers and pos_key
  hits_df <- hits_df %>%
    mutate(sequence_context_k = center_kmer(sequence_context, k)) %>%
    filter(valid_kmer(sequence_context_k, k)) %>%
    mutate(pos = as.numeric(pos),
           pos_key = make_pos_key(pos, w = position_window))
  stage_df <- stage_df %>%
    mutate(sequence_context_k = center_kmer(sequence_context, k)) %>%
    filter(valid_kmer(sequence_context_k, k)) %>%
    mutate(pos = as.numeric(pos),
           pos_key = make_pos_key(pos, w = position_window))
  stopifnot(nrow(hits_df) > 0, nrow(stage_df) > 0)
  
  # Non-hits pool = STAGE \ HITS by exact locus
  non_hits <- stage_df %>% anti_join(hits_df %>% select(chr_pos), by = "chr_pos")
  
  # hit counts per pos_key
  hit_counts <- hits_df %>% count(pos_key, name = "n_hit_pos")
  
  # available background pools by pos_key
  pools <- non_hits %>% mutate(.rowid = row_number()) %>%
    group_by(pos_key) %>%
    summarise(idx = list(.rowid), .groups = "drop")
  
  available <- hit_counts %>% inner_join(pools, by = "pos_key")
  
  # Keep hits only where there is background
  hits_keep <- hits_df %>% inner_join(available %>% select(pos_key), by = "pos_key")
  n_drop <- nrow(hits_df) - nrow(hits_keep)
  if (n_drop > 0) message("Dropped ", n_drop, " hit(s) from positions with no non-hit background.")
  
  # Sample background sequences to *match hit counts per pos_key*
  seq_bg <- non_hits$sequence_context_k
  chosen_bg <- vector("list", nrow(available))
  for (i in seq_len(nrow(available))) {
    ids <- available$idx[[i]]
    n_h <- available$n_hit_pos[i]
    # with replacement if not enough in pool
    chosen_bg[[i]] <- if (length(ids) >= n_h) sample(ids, n_h, replace = FALSE)
    else                       sample(ids, n_h, replace = TRUE)
  }
  bg_ids <- unlist(chosen_bg, use.names = FALSE)
  bg_seqs <- if (length(bg_ids)) seq_bg[bg_ids] else character(0)
  
  list(
    hits_keep = hits_keep,
    bg_seqs   = bg_seqs
  )
}

# ------------------------- RUN (per GROUP) -------------------------
for (grp in groups) {
  pretty <- if (grp %in% names(group_labels)) group_labels[[grp]] else grp
  message("Motif log-odds (pos-matched bg), merged hits for group: ", grp)
  
  # Merge: all hits (sgid + total), dedupe by chr_pos
  hits_files_grp <- hits_info %>% filter(group == grp) %>% pull(file_hits)
  stopifnot(length(hits_files_grp) > 0)
  hits <- lapply(hits_files_grp, function(f) readr::read_csv(f, show_col_types = FALSE)) %>%
    bind_rows() %>%
    mutate(pos = as.numeric(pos)) %>%
    distinct(chr_pos, .keep_all = TRUE)
  
  # Merge: all stage files for this group, dedupe by chr_pos
  stage_files_grp <- stage_info %>% filter(group == grp) %>% pull(file_stage)
  stopifnot(length(stage_files_grp) > 0)
  stage <- lapply(stage_files_grp, function(f) readr::read_csv(f, show_col_types = FALSE)) %>%
    bind_rows() %>%
    mutate(pos = as.numeric(pos)) %>%
    distinct(chr_pos, .keep_all = TRUE)
  
  # Required columns check
  need <- c("chr_pos","sequence_context","biotype","family","ref_base","pos")
  stopifnot(all(need %in% names(hits)), all(need %in% names(stage)))
  
  # Build pos-matched background (and keep only hits with available bg)
  pm <- build_pos_matched_bg(hits, stage, k = k, position_window = position_window)
  hits_keep <- pm$hits_keep
  bg_seqs   <- pm$bg_seqs
  
  # PWMs and log-odds (pos-matched)
  p_hit <- pwm_from_seqs(hits_keep$sequence_context_k, k)
  p_bg  <- pwm_from_seqs(bg_seqs, k)
  logodds <- logodds_pwm(p_hit, p_bg, eps = epsilon)
  
  # Sanity: ΔGC derived from the logo equals direct ΔGC (should match near-exactly)
  dgc_logo   <- colSums(p_hit[c("G","C"), , drop = FALSE] - p_bg[c("G","C"), , drop = FALSE])
  dgc_direct <- gc_curve(hits_keep$sequence_context_k, k) - gc_curve(bg_seqs, k)
  max_diff <- max(abs(dgc_logo - dgc_direct), na.rm = TRUE)
  message(sprintf("  ΔGC consistency check (max |logo - direct|): %.3g", max_diff))
  
  # Plot (one PDF per group)
  out_pdf <- paste0("motif_logodds_posmatched_", grp, "_joinedHits.pdf")
  plot_logo_pdf(logodds, out_pdf,
                paste0("Motif log-odds (pos-matched), ", pretty),)
  message("  wrote: ", out_pdf)
}

message("Done. Produced one logo per group (2 PDFs total if both groups present).")
