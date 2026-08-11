load_template_config <- function(config_file = "config.yaml") {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("The 'yaml' package is required to read ", config_file,
         ". Install it with install.packages('yaml').", call. = FALSE)
  }
  if (!file.exists(config_file)) {
    stop("Config file not found: ", config_file,
         ". Copy config.template.yaml to config.yaml in this run directory and fill in run-specific values.",
         call. = FALSE)
  }
  yaml::read_yaml(config_file)
}

cfg_get <- function(config, path, default = NULL) {
  parts <- strsplit(path, "\\.", fixed = FALSE)[[1]]
  value <- config
  for (part in parts) {
    if (is.null(value) || is.null(value[[part]])) {
      return(default)
    }
    value <- value[[part]]
  }
  value
}

cfg_path <- function(config, path, default = NULL) {
  value <- cfg_get(config, path, default)
  if (is.null(value)) {
    return(NULL)
  }
  path.expand(value)
}

cfg_named_vector <- function(x) {
  if (is.null(x)) {
    return(character())
  }
  out <- unlist(x, use.names = TRUE)
  # Preserve names: as.character() alone drops them, which breaks name-based
  # lookups (suffix_to_condition) and named colour scales (plots.colors.*).
  setNames(as.character(out), names(out))
}

cfg_number <- function(config, path, default) {
  value <- cfg_get(config, path, default)
  as.numeric(value)
}

cfg_integer <- function(config, path, default) {
  value <- cfg_get(config, path, default)
  as.integer(value)
}

cfg_bool <- function(config, path, default = FALSE) {
  value <- cfg_get(config, path, default)
  if (is.logical(value)) {
    return(value)
  }
  tolower(as.character(value)) %in% c("true", "1", "yes", "y")
}

cfg_template <- function(pattern, ...) {
  values <- list(...)
  for (name in names(values)) {
    pattern <- gsub(paste0("\\{", name, "\\}"), values[[name]], pattern)
  }
  pattern
}

cfg_file <- function(config, key, ...) {
  cfg_template(cfg_get(config, paste0("naming.", key)), ...)
}

load_theme_from_config <- function(config) {
  theme_files <- c(cfg_path(config, "theme.file"), unlist(cfg_get(config, "theme.fallback_files", list())))
  theme_files <- path.expand(theme_files[!is.na(theme_files) & nzchar(theme_files)])
  theme_file <- theme_files[file.exists(theme_files)][1]
  if (is.na(theme_file)) {
    stop("None of the configured theme files exist: ", paste(theme_files, collapse = ", "), call. = FALSE)
  }
  source(theme_file)
  invisible(theme_file)
}

read_sample_matrix <- function(config) {
  sample_matrix_file <- cfg_path(config, "metadata.sample_matrix")
  if (!file.exists(sample_matrix_file)) {
    stop("Sample matrix not found: ", sample_matrix_file, call. = FALSE)
  }
  sample_info <- read.csv(sample_matrix_file, stringsAsFactors = FALSE)
  cols <- cfg_get(config, "metadata.columns")
  required <- unlist(cols[c("sample_name", "group_name", "replicate", "plot_name")], use.names = FALSE)
  missing <- setdiff(required, colnames(sample_info))
  if (length(missing) > 0) {
    stop("Sample matrix is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  sample_info
}

sample_col <- function(config, key) {
  cfg_get(config, paste0("metadata.columns.", key))
}

sample_groups <- function(config, sample_info) {
  g <- unique(sample_info[[sample_col(config, "group_name")]])
  # Optional: restrict the main pipeline to groups matching this regex (e.g.
  # "^5-gly_" for a single-construct run). Unset -> all groups (default).
  rx <- cfg_get(config, "analysis.groups_include_regex", NULL)
  if (!is.null(rx) && nzchar(rx)) g <- g[grepl(rx, g)]
  g
}

group_plot_lookup <- function(config, sample_info) {
  group_col <- sample_col(config, "group_name")
  plot_col <- sample_col(config, "plot_name")
  plot_name_map <- unique(sample_info[, c(group_col, plot_col)])
  setNames(plot_name_map[[plot_col]], plot_name_map[[group_col]])
}

plot_name_for_group <- function(config, sample_info, group_name) {
  lookup <- group_plot_lookup(config, sample_info)
  value <- lookup[[group_name]]
  if (is.null(value) || is.na(value) || value == "") {
    return(group_name)
  }
  value
}

sample_suffix_map <- function(config) {
  cfg_named_vector(cfg_get(config, "counts.sample_suffixes"))
}

suffix_condition_map <- function(config) {
  cfg_named_vector(cfg_get(config, "analysis.suffix_to_condition"))
}

make_sample_names <- function(config, ids) {
  suffixes <- unname(sample_suffix_map(config))
  as.vector(vapply(ids, function(id) paste0(id, suffixes), character(length(suffixes))))
}

parse_sample_name <- function(config, sname) {
  suffixes <- unname(sample_suffix_map(config))
  suffix_pattern <- paste0("(", paste(suffixes, collapse = "|"), ")$")
  suffix <- sub(paste0("^.*", suffix_pattern), "\\1", sname)
  sample_id <- sub(suffix_pattern, "", sname)
  list(id = sample_id, suffix = suffix)
}

suffix_to_condition <- function(config, suffix_vec) {
  mapping <- suffix_condition_map(config)
  unname(mapping[as.character(suffix_vec)])
}

make_safe_name <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

biotype_levels_from_config <- function(config) {
  unlist(cfg_get(config, "biotypes.levels"), use.names = FALSE)
}

biotype_colors_from_config <- function(config, levels) {
  palette_name <- cfg_get(config, "biotypes.palette", "PrettyCols::Rainbow")
  colors <- as.character(paletteer::paletteer_d(palette_name, n = length(levels)))
  names(colors) <- levels
  colors
}
