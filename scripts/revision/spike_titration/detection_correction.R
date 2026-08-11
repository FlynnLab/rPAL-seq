# Post-DE detection-correction offsets for the summed 5' glyco spike (paper revision).
# See config titration.detection_correction. The 5' glycan blocks TSO detection of the
# glyco species (glyco detected at e_in in Input, beta in IP; unglyco Fluc unaffected).
# Un-suppressing the glyco component of the SUMMED spike (Input /e_in, IP /beta) and
# re-summing gives a per-dose ADDITIVE log2 offset to the summed IP/Input LFC:
#   delta(dose) = log2[(gp/beta + fp)/(gi/e_in + fi)] - log2[(gp+fp)/(gi+fi)]
# (g=glyco, f=unglyco Fluc; i=Input, p=IP; rRNA-size-factor-normalized, per replicate).
# Reduces to +log2(e_in/beta)=+0.20 at 100% glyco, ->0 at 0%. Deterministic (SE-preserving).
# FIGURES ONLY -- never applied to count_matrix or the DE. c=0.966 (HPLC UV peak area).
#
# Returns data.frame(group, construct, nominal, delta) or NULL if disabled.
detection_offsets <- function(cfg) {
  dc <- cfg_get(cfg, "titration.detection_correction", NULL)
  if (is.null(dc) || !isTRUE(as.logical(dc$enabled))) return(NULL)
  e_in <- as.numeric(dc$e_in); beta <- as.numeric(dc$beta)
  glyco <- dc$spike_glyco; unglyco <- dc$spike_unglyco
  constructs <- as.character(unlist(dc$constructs))
  cm <- read.csv(cfg_path(cfg, "counts.count_matrix"), row.names = 1, check.names = FALSE)
  meta <- read_sample_matrix(cfg); gc <- sample_col(cfg, "group_name"); sn <- sample_col(cfg, "sample_name")
  rrna <- as.character(unlist(cfg_get(cfg, "analysis.size_factor_feature_ids")))
  lv <- cfg_get(cfg, "titration.levels"); levmap <- setNames(as.numeric(unlist(lv)), names(lv))
  rows <- list()
  for (g in unique(meta[[gc]])) {
    constr <- sub("_[^_]+$", "", g); lev <- sub("^.*_", "", g)
    if (!constr %in% constructs) next
    sub <- meta[meta[[gc]] == g, ]; ic <- paste0(sub[[sn]], "i"); pc <- paste0(sub[[sn]], "p")
    ic <- ic[ic %in% colnames(cm)]; pc <- pc[pc %in% colnames(cm)]
    if (!length(ic) || !length(pc)) next
    sf <- colSums(cm[rrna, c(ic, pc), drop = FALSE]); sf <- sf / mean(sf)
    gi <- as.numeric(cm[glyco, ic]) / sf[ic]; gp <- as.numeric(cm[glyco, pc]) / sf[pc]
    fi <- as.numeric(cm[unglyco, ic]) / sf[ic]; fp <- as.numeric(cm[unglyco, pc]) / sf[pc]
    meas <- log2((gp + fp) / (gi + fi)); corr <- log2((gp / beta + fp) / (gi / e_in + fi))
    rows[[g]] <- data.frame(group = g, construct = constr,
                            nominal = as.numeric(levmap[lev]), delta = mean(corr - meas),
                            stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}
