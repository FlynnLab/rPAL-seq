#!/usr/bin/env Rscript
# De-novo: two-stage dose-response for the summed spike (Fluc_spike_total).
#   stage 1  read per-level DESeq2 log2FC +/- lfcSE for the summed spike straight
#            from the canonical deseq2_results_pI_<group>.csv (main pipeline runs
#            with counts.merge_for_de; NO DESeq re-run here).
#   stage 2  inverse-variance-weighted regression of log2FC on nominal glyco%
# Plots the 5 LFC+/-SE points + the weighted fit; slope/p in the subheading.
# Construct(s) plotted follow titration.constructs_include (default all).
#
# Run from this directory:  Rscript denovo_doseresponse_metareg.R

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
source("template_config.R"); cfg <- load_template_config(); load_theme_from_config(cfg)

meta <- read_sample_matrix(cfg)
gc <- sample_col(cfg,"group_name")
levmap <- c("0"=0,"1"=25,"2"=50,"3"=75,"4"=100)
cc <- cfg_get(cfg,"titration.constructs")
labmap <- setNames(vapply(names(cc), function(k) cc[[k]]$label, character(1)), names(cc))
colmap <- setNames(vapply(names(cc), function(k) cc[[k]]$color, character(1)), names(cc))
include_keys <- as.character(unlist(cfg_get(cfg, "titration.constructs_include", names(cc))))
include_keys <- names(cc)[names(cc) %in% include_keys]
de_spike <- names(cfg_get(cfg, "counts.merge_for_de"))[1]   # summed spike feature, e.g. Fluc_spike_total

## stage 1: read summed-spike LFC + lfcSE from canonical deseq2_results_pI (no re-run)
per_level <- function(construct) {
  groups <- unique(meta[[gc]][grepl(paste0("^", construct, "_"), meta[[gc]])])
  do.call(rbind, lapply(groups, function(g) {
    f <- cfg_file(cfg, "deseq", comparison = "pI", group = g)
    if (!file.exists(f)) { warning("missing ", f); return(NULL) }
    x <- read.csv(f); r <- x[x$transcript_id == de_spike, ]
    if (nrow(r) == 0) { warning(de_spike, " not in ", f); return(NULL) }
    data.frame(construct = construct, nominal = as.numeric(levmap[sub("^.*_", "", g)]),
               LFC = r$log2FoldChange[1], SE = r$lfcSE[1], stringsAsFactors = FALSE)
  }))
}
d <- do.call(rbind, lapply(include_keys, per_level))
d$construct_label <- factor(unname(labmap[d$construct]), levels = unname(labmap[include_keys]))

## stage 2: inverse-variance-weighted regression (nominal dose), per construct
fit_lines <- list(); subt_parts <- character()
for (k in include_keys) {
  dk <- d[d$construct == k, ]
  m <- lm(LFC ~ nominal, weights = 1/SE^2, data = dk)
  s <- summary(m)$coefficients["nominal", ]
  nd <- data.frame(nominal = seq(0, 100, 1))
  pr <- as.data.frame(predict(m, nd, interval = "confidence"))
  fit_lines[[k]] <- data.frame(construct = k, construct_label = unname(labmap[k]),
                               nominal = nd$nominal, fit = pr$fit, lwr = pr$lwr, upr = pr$upr)
  subt_parts <- c(subt_parts, sprintf("slope = %+.3f / %%,  p = %.3f", s[1], s[4]))
}
fitdf <- do.call(rbind, fit_lines)
fitdf$construct_label <- factor(fitdf$construct_label, levels = unname(labmap[include_keys]))
subt <- if (length(include_keys) == 1)
  paste0("IVW regression, ", subt_parts[1]) else
  "IVW regression"

# Dot color. Single construct -> the "IP" condition color (#E01A4F), matching
# denovo_spike_coverage.R (this plot is log2 IP/Input). Multiple constructs ->
# per-construct colors (colmap) so the facets stay visually distinct.
if (length(include_keys) == 1) {
  ip_col <- cfg_get(cfg, "plots.colors.conditions.IP", "#E01A4F")
  pcols <- setNames(ip_col, unname(labmap[include_keys]))
} else {
  pcols <- setNames(unname(colmap[include_keys]), unname(labmap[include_keys]))
}

# Per-point SEs are used to weight the fit (inverse-variance) but not drawn:
# the trend's uncertainty is the fit CI ribbon (tighter, pools all levels).
# Set show_se: true in config to overlay light per-level error bars.
# Fig 1l geometry, 2026-08-03. It sits under panel g in the same column, so it is sized to g's
# block width rather than to h/j/k -- and, more importantly, it is now exported at the size it will
# actually OCCUPY on the page.
#
# WHY. Measured from the author's 2026-08-03 export: the 2.0 x 2.0 in panels were being placed at
# ~0.68 scale in Illustrator (j's x-axis rule is 216 px in the source PDF and 147 px in the figure;
# the identical "Nominal glyco-Fluc (%)" string is 195-196 px in the sources and 130 px in the
# figure -- two independent measures agreeing at 0.66-0.68). That silently drops theme_nature's
# 8 pt axis text to ~5.4 pt on the page, and to ~5.2 pt once the 190 mm figure is reduced to
# 180 mm -- Nature's floor is 5 pt. Exporting at the placed size and inserting at 100 % keeps the
# house type size instead of scaling it away.
# Sized from the 2026-08-03 14:26 export, measured: panel g's block above it is **2.33 in** wide,
# and the slot l currently occupies is **1.31 in** tall. So l must fill 2.33 in of width WITHOUT
# growing past 1.31 in of height -- otherwise fixing the aspect ratio makes column 2 (already the
# tallest, hanging 0.40 in below the rest) worse. 1.30 in leaves ~0.95 in of plot area, and the
# rotated y title "log2(IP/Input)" is ~0.88 in at 9 pt, so it still fits.
PANEL_W <- 2.33  # = panel g's block width, so l fills the column instead of 63 % of it
PANEL_H <- 1.30  # = the height l already occupies, so column 2 does not grow

show_se <- cfg_bool(cfg, "titration.metareg_show_se", FALSE)
p <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_ribbon(data = fitdf, aes(nominal, ymin = lwr, ymax = upr), fill = "grey85", alpha = 0.6) +
  geom_line(data = fitdf, aes(nominal, fit), color = "grey20", linewidth = 0.5) +
  {if (show_se) geom_errorbar(data = d, aes(nominal, ymin = LFC - SE, ymax = LFC + SE, color = construct_label),
                              width = 4, linewidth = 0.35, alpha = 0.5, show.legend = FALSE)} +
  geom_point(data = d, aes(nominal, LFC, color = construct_label), size = 2.1, show.legend = FALSE) +
  scale_color_manual(values = pcols) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
  # Fig 1k, 2026-07-31: title/subtitle dropped (panel of a multi-panel figure) and the y
  # label shortened for a 2-inch panel. "DESeq2" and the summed-spike detail go to the legend.
  labs(x = "Nominal glyco-Fluc (%)", y = expression(log[2]*"(IP/Input)")) +
  theme_nature +
  theme(plot.title = element_blank())
if (length(include_keys) > 1) p <- p + facet_wrap(~construct_label)

ggsave("denovo_doseresponse_metareg.pdf", p,
       width = PANEL_W, height = PANEL_H, dpi = 300)
write.csv(d[order(d$construct_label, d$nominal), ], "denovo_doseresponse_metareg.csv", row.names = FALSE)
message("wrote denovo_doseresponse_metareg.{pdf,csv}")

## ---- detection-corrected figure (paper revision; FIGURES ONLY, not counts/DE) ----
## Corrected-only: summed-spike LFC shifted by the per-dose detection offset; the
## correction is disclosed in the subheading (no measured overlay). Mirrors the
## original figure's single-series style.
source("detection_correction.R")
off <- detection_offsets(cfg)
if (!is.null(off)) {
  dc <- merge(d, off[, c("construct", "nominal", "delta")], by = c("construct", "nominal"), all.x = TRUE)
  dc$delta[is.na(dc$delta)] <- 0
  dc$LFC_corr <- dc$LFC + dc$delta
  dc$construct_label <- factor(unname(labmap[dc$construct]), levels = unname(labmap[include_keys]))
  cfit <- list(); csub <- character()
  for (k in include_keys) {
    dk <- dc[dc$construct == k, ]
    m <- lm(LFC_corr ~ nominal, weights = 1 / SE^2, data = dk)
    s <- summary(m)$coefficients["nominal", ]
    nd <- data.frame(nominal = 0:100); pr <- as.data.frame(predict(m, nd, interval = "confidence"))
    cfit[[k]] <- data.frame(construct_label = unname(labmap[k]), nominal = nd$nominal,
                            fit = pr$fit, lwr = pr$lwr, upr = pr$upr)
    csub <- c(csub, sprintf("slope = %+.3f / %%,  p = %.3f", s[1], s[4]))
  }
  cfitdf <- do.call(rbind, cfit); cfitdf$construct_label <- factor(cfitdf$construct_label, levels = unname(labmap[include_keys]))
  c_val <- as.numeric(cfg_get(cfg, "titration.detection_correction.c"))
  subt_c <- if (length(include_keys) == 1)
    sprintf("IVW regression, detection-corrected (c=%.3f)\n%s", c_val, csub[1]) else
    sprintf("IVW regression, detection-corrected (c=%.3f)", c_val)
  pc <- ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_ribbon(data = cfitdf, aes(nominal, ymin = lwr, ymax = upr), fill = "grey85", alpha = 0.6) +
    geom_line(data = cfitdf, aes(nominal, fit), color = "grey20", linewidth = 0.5) +
    {if (show_se) geom_errorbar(data = dc, aes(nominal, ymin = LFC_corr - SE, ymax = LFC_corr + SE, color = construct_label),
                                width = 4, linewidth = 0.35, alpha = 0.5, show.legend = FALSE)} +
    geom_point(data = dc, aes(nominal, LFC_corr, color = construct_label), size = 2.1, show.legend = FALSE) +
    scale_color_manual(values = pcols) +
    scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
    labs(x = "Nominal glyco-Fluc (%)",
         y = expression(log[2]*"(IP/Input)"),
         title = "Spike enrichment dose-response", subtitle = subt_c) +
    theme_nature +
    theme(plot.title = element_blank())
  if (length(include_keys) > 1) pc <- pc + facet_wrap(~construct_label)
  ggsave("denovo_doseresponse_metareg_corrected.pdf", pc,
         width = PANEL_W, height = PANEL_H, dpi = 300)
  write.csv(dc[order(dc$nominal), c("construct", "nominal", "LFC", "SE", "delta", "LFC_corr")],
            "denovo_doseresponse_metareg_corrected.csv", row.names = FALSE)
  message("wrote denovo_doseresponse_metareg_corrected.{pdf,csv}")
}
