library(tidyverse)
library(Biostrings)

# --------- Settings ---------
fasta_path <- path.expand("/path/to/transcriptome.fa")
annot_path <- path.expand("/path/to/transcriptome/annotation.csv")
context_left  <- 7
context_right <- 7

# --------- Load data ---------
annotation <- read_csv(annot_path, show_col_types = FALSE)
fasta <- readDNAStringSet(fasta_path)
names(fasta) <- sub(" .*", "", names(fasta))  # Clean FASTA headers

# --------- Helpers ---------

# Simple, because your format is guaranteed:
# transcript_id = everything before last "_"
# pos           = last "_" field, as integer
parse_chr_pos <- function(x) {
  tibble(
    transcript_id = sub("_[^_]+$", "", x),
    pos           = suppressWarnings(as.integer(sub(".*_", "", x)))
  )
}

get_context_and_ref <- function(tid, pos, fasta,
                                win_left = context_left, win_right = context_right) {
  if (is.na(tid) || is.na(pos) || !(tid %in% names(fasta))) {
    return(tibble(sequence_context = NA_character_, ref_base = NA_character_))
  }
  seq <- fasta[[tid]]
  seq_len <- nchar(as.character(seq))  # avoids width() namespace issues
  
  start_pos <- pos - win_left
  end_pos   <- pos + win_right
  
  pad_left  <- if (start_pos < 1) 1 - start_pos else 0
  pad_right <- if (end_pos > seq_len) end_pos - seq_len else 0
  
  start_pos <- max(1, start_pos)
  end_pos   <- min(seq_len, end_pos)
  
  context <- as.character(Biostrings::subseq(seq, start = start_pos, end = end_pos))
  if (pad_left  > 0) context <- paste0(strrep("N", pad_left),  context)
  if (pad_right > 0) context <- paste0(context, strrep("N", pad_right))
  
  ref_base <- if (pos >= 1 && pos <= seq_len) {
    as.character(Biostrings::subseq(seq, start = pos, end = pos))
  } else {
    NA_character_
  }
  
  tibble(sequence_context = context, ref_base = ref_base)
}

# --------- Annotate hits ---------

annotate_hits <- function(file_in, file_out, annotation, fasta) {
  hits <- read_csv(file_in, show_col_types = FALSE)
  if (!"chr_pos" %in% names(hits)) {
    stop("Expected a 'chr_pos' column in: ", file_in)
  }
  
  parsed <- parse_chr_pos(hits$chr_pos)
  
  hits_annotated <- hits %>%
    bind_cols(parsed) %>%
    left_join(annotation, by = "transcript_id") %>%
    mutate(seq_data = purrr::map2(transcript_id, pos,
                                  ~ get_context_and_ref(.x, .y, fasta))) %>%
    tidyr::unnest_wider(seq_data)
  
  write_csv(hits_annotated, file_out)
}

# --------- Annotate all hit files ---------

hit_files <- list.files(pattern = "^(hits_|stage).*\\.csv$", full.names = TRUE)
hit_files <- hit_files[!grepl("(^|/)Annotated_", hit_files)]
out_files <- file.path(dirname(hit_files), paste0("Annotated_", basename(hit_files)))

purrr::walk2(hit_files, out_files, ~ annotate_hits(.x, .y, annotation, fasta))
