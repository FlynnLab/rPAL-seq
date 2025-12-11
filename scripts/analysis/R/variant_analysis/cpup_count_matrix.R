# step0
library(tidyverse)
library(readr)
library(tools)

# ----- user inputs -----
sample_info_path <- "/path/to/metadata.csv"  # must have: sample_name, group_name, replicates

hematopoietic      <- c("AML2","AML3","Molm13","Jurkat","Jeko1","Nalm6")
non_hematopoietic  <- c("HeLa","293T","A549","Huh7","DiFi","LPS853","LN308")

fixed_colnames <- c("depth_raw","A","C","G","T","N","Skip","Gap","Insert","Delete")
file_list <- list.files(pattern = "\\.txt$", full.names = TRUE)
stopifnot(length(file_list) > 0)

# ----- helpers -----
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

read_one <- function(fp) {
  stem <- file_path_sans_ext(basename(fp))
  df <- read_delim(fp, "\t", col_types = cols(.default = col_character()), quote = "")
  depth_split <- strsplit(df[[4]], ",")
  parsed_matrix <- do.call(rbind, lapply(depth_split, function(x) as.numeric(x[1:10])))
  colnames(parsed_matrix) <- fixed_colnames
  parsed <- as_tibble(parsed_matrix)
  
  tibble(
    chr  = df[[1]], pos = df[[2]], ref_base = df[[3]],
    chr_pos = paste(df[[1]], df[[2]], sep = "_"),
    sample  = stem
  ) %>%
    bind_cols(parsed) %>%
    mutate(across(c(depth_raw, A, C, G, T, N, Skip, Gap, Insert, Delete),
                  ~ as.integer(replace_na(., 0L))))
}

# ----- read & join -----
sample_info <- read_csv(sample_info_path, show_col_types = FALSE) %>%
  mutate(across(c(sample_name, group_name, replicates), as.character)) %>%
  mutate(super_group = case_when(
    group_name %in% hematopoietic     ~ "hematopoietic",
    group_name %in% non_hematopoietic ~ "non-hematopoietic",
    TRUE ~ NA_character_
  ))

raw_counts <- purrr::map_dfr(file_list, read_one)

counts <- raw_counts %>%
  mutate(condition = parse_condition(sample),
         sample_id = base_id(sample)) %>%
  left_join(sample_info, by = c("sample_id" = "sample_name")) %>%
  # keep only samples we can classify into a super_group
  filter(!is.na(super_group), !is.na(condition)) %>%
  mutate(
    SGID = Skip + Gap + Insert + Delete,
    REF  = case_when(ref_base == "A" ~ A,
                     ref_base == "C" ~ C,
                     ref_base == "G" ~ G,
                     ref_base == "T" ~ T,
                     TRUE ~ 0L),
    A_star = if_else(ref_base == "A", 0L, A),
    C_star = if_else(ref_base == "C", 0L, C),
    G_star = if_else(ref_base == "G", 0L, G),
    T_star = if_else(ref_base == "T", 0L, T),
    mismatches    = A_star + C_star + G_star + T_star,
    total_variant = mismatches + SGID,
    nonvariant    = pmax(depth_raw - total_variant, 0L),
    pair_id       = paste(group_name, replicates, sep = "_")  # pairing block
  ) %>%
  select(chr_pos, chr, pos, ref_base,
         super_group, group_name, replicates, pair_id,
         sample_id, sample, condition,
         depth_raw, REF, SGID, A_star, C_star, G_star, T_star,
         mismatches, total_variant, nonvariant)

write_csv(counts, "counts_long.csv")

# Optional: depth matrix
counts %>% select(chr_pos, sample, depth_raw) %>%
  pivot_wider(names_from = sample, values_from = depth_raw) %>%
  write_csv("depth_matrix.csv")

message("✅ Wrote counts_long.csv (aggregated to super_groups).")
