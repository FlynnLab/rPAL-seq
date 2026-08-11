#!/usr/bin/env Rscript
# De-novo #3: spike_total IP-enrichment dose-response, using MEASURED input
# glyco% (per replicate) as the dose axis instead of nominal.
# Rationale: 5'-glyco is under-detected in Input (compressed), so the measured
# input glycan fraction is the truer "dose actually present" than the nominal
# mixing ratio. Fit is replicate-blocked (the dominant variance is a per-rep
# offset). Compares nominal vs measured as predictors.
#
# Enrichment metric (paired, per replicate):
#   E = log2[ (spike_total/rRNA)_IP / (spike_total/rRNA)_Input ]
# where spike_total = Fluc_spike + glyco_Fluc_spike.

suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
source("template_config.R"); cfg <- load_template_config(); load_theme_from_config(cfg)

rrna  <- unlist(cfg_get(cfg, "analysis.size_factor_feature_ids"), use.names = FALSE)
glyco <- cfg_get(cfg, "titration.spike.glyco"); fluc <- cfg_get(cfg, "titration.spike.total_ref")
cm <- read.csv(cfg_path(cfg, "counts.count_matrix"), row.names = 1, check.names = FALSE)
meta <- read_sample_matrix(cfg)
sname <- sample_col(cfg,"sample_name"); gcol <- sample_col(cfg,"group_name"); rcol <- sample_col(cfg,"replicate")
levmap <- c("0"=0,"1"=25,"2"=50,"3"=75,"4"=100)
cc <- cfg_get(cfg,"titration.constructs")
labmap <- setNames(vapply(names(cc), function(k) cc[[k]]$label, character(1)), names(cc))
colmap <- setNames(vapply(names(cc), function(k) cc[[k]]$color, character(1)), names(cc))
# Optional filter: only PLOT these constructs (default all). CSV stays full.
include_keys <- as.character(unlist(cfg_get(cfg, "titration.constructs_include", names(cc))))
include_keys <- names(cc)[names(cc) %in% include_keys]

sp  <- function(col) cm[glyco, col] + cm[fluc, col]
rt  <- function(col) sum(cm[rrna, col])

# Detection correction (paper revision; figures only): un-suppress the glyco component
# of the summed spike (Input /e_in, IP /beta) for the blocked constructs. See config.
dc <- cfg_get(cfg, "titration.detection_correction", NULL)
corr_on <- !is.null(dc) && isTRUE(as.logical(dc$enabled))
e_in <- if (corr_on) as.numeric(dc$e_in) else 1
beta <- if (corr_on) as.numeric(dc$beta) else 1
corr_constructs <- if (corr_on) as.character(unlist(dc$constructs)) else character(0)
sp_c <- function(col, frac, this) { f <- if (this && frac == "i") e_in else if (this && frac == "p") beta else 1
                                    cm[glyco, col] / f + cm[fluc, col] }

dat <- do.call(rbind, lapply(seq_len(nrow(meta)), function(i) {
  s <- meta[[sname]][i]; grp <- meta[[gcol]][i]
  construct <- sub("_[^_]+$", "", grp); lvl <- sub("^.*_", "", grp)
  gi <- cm[glyco, paste0(s,"i")]; fi <- cm[fluc, paste0(s,"i")]
  data.frame(
    construct = construct, rep = factor(meta[[rcol]][i]),
    nominal = as.numeric(levmap[lvl]),
    measured = 100 * gi / (gi + fi),
    E = log2((sp(paste0(s,"p"))/rt(paste0(s,"p"))) / (sp(paste0(s,"i"))/rt(paste0(s,"i")))),
    E_corr = { this <- construct %in% corr_constructs
      log2((sp_c(paste0(s,"p"),"p",this)/rt(paste0(s,"p"))) / (sp_c(paste0(s,"i"),"i",this)/rt(paste0(s,"i")))) },
    stringsAsFactors = FALSE
  )
}))
dat$construct_label <- factor(unname(labmap[dat$construct]), levels = unname(labmap))
datp <- dat %>% filter(construct %in% include_keys)
datp$construct_label <- droplevels(datp$construct_label)

## ---- nominal vs measured, replicate-blocked ------------------------------
fit_report <- function(d) {
  f <- function(x) {
    m <- lm(reformulate(c(x, "rep"), "E"), data = d)
    s <- summary(m)$coefficients[x, ]
    c(slope = unname(s[1]), p = unname(s[4]), adjR2 = summary(m)$adj.r.squared)
  }
  list(nominal = f("nominal"), measured = f("measured"))
}
cat(sprintf("%-16s %-9s %10s %8s %7s\n","construct","predictor","slope/%","p","adjR2"))
for (cn in include_keys) {
  r <- fit_report(dat[dat$construct == cn, ])
  for (pr in c("nominal","measured"))
    cat(sprintf("%-16s %-9s %+10.4f %8.4f %7.2f\n", cn, pr, r[[pr]]["slope"], r[[pr]]["p"], r[[pr]]["adjR2"]))
}

## ---- figure: E vs dose, replicate-colored + fitted slope --
# dose axis: "nominal" (evenly spaced; avoids the measured-axis leverage at 100%)
# or "measured" (input glycan %). Configurable via titration.dose_axis.
axis_col <- cfg_get(cfg, "titration.dose_axis", "measured")
axis_lab <- if (axis_col == "nominal") "Nominal glyco-Fluc (%)" else "Measured input glycan (%)"
ttl      <- if (axis_col == "nominal") "Spike enrichment vs nominal dose" else "Spike enrichment vs measured input dose"

# slope + p for the subheading (replicate-adjusted linear model). With a single
# construct we annotate the one fit; with >1 the fit is per-construct (faceted,
# geom_smooth fits per panel) so we use a generic subheading and rely on the
# per-construct table printed above / the CSV for exact slopes.
multi <- length(include_keys) > 1
if (multi) {
  subt <- "Replicate-adjusted linear model (per construct)"
} else {
  m_sub <- lm(reformulate(c(axis_col, "rep"), "E"), data = datp)
  sl <- unname(coef(m_sub)[axis_col]); pv <- summary(m_sub)$coefficients[axis_col, "Pr(>|t|)"]
  subt <- sprintf("Replicate-adjusted linear model\nslope = %+.3f / %%,   p = %.3f", sl, pv)
}

p <- ggplot(datp, aes(x = .data[[axis_col]], y = E)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_smooth(method = "lm", se = TRUE, color = "grey20", fill = "grey80", linewidth = 0.5) +
  geom_point(aes(color = rep), size = 1.5) +
  scale_color_manual(values = c("1"="#8C2981","2"="#E01A4F","3"="#53B3CB"), name = "rep") +
  labs(x = axis_lab, y = expression(log[2]~"spike"[total]~"(IP/Input)"),
       title = ttl, subtitle = subt) +
  theme_nature +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 9),
        plot.subtitle = element_text(size = 7, color = "grey30", lineheight = 0.95))
if (multi) p <- p + facet_wrap(~construct_label)

ggsave("denovo_enrichment_doseresponse.pdf", p,
       width = if (multi) 5.2 else 3.2, height = 3.0, dpi = 300)
write.csv(datp[order(datp$construct_label, datp$nominal, datp$rep), ], "denovo_enrichment_doseresponse.csv", row.names = FALSE)
message("done")

## ---- detection-corrected figure (paper revision; FIGURES ONLY, not counts/DE) ----
## Corrected-only: per-rep enrichment recomputed with the glyco component un-suppressed;
## correction disclosed in the subheading (no measured overlay). Mirrors original style.
if (corr_on && !multi) {
  sc <- summary(lm(reformulate(c(axis_col, "rep"), "E_corr"), data = datp))$coefficients[axis_col, ]
  subt_c <- sprintf("Replicate-adjusted LM, detection-corrected\n(c=%.3f)   slope = %+.3f / %%,  p = %.3f",
                    as.numeric(dc$c), sc[1], sc[4])
  pc <- ggplot(datp, aes(x = .data[[axis_col]], y = E_corr)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_smooth(method = "lm", se = TRUE, color = "grey20", fill = "grey80", linewidth = 0.5) +
    geom_point(aes(color = rep), size = 1.5) +
    scale_color_manual(values = c("1" = "#8C2981", "2" = "#E01A4F", "3" = "#53B3CB"), name = "rep") +
    labs(x = axis_lab, y = expression(log[2]~"spike"[total]~"(IP/Input)"), title = ttl, subtitle = subt_c) +
    theme_nature +
    theme(legend.position = "bottom", plot.title = element_text(size = 9),
          plot.subtitle = element_text(size = 6.4, color = "grey30", lineheight = 0.95))
  ggsave("denovo_enrichment_doseresponse_corrected.pdf", pc, width = 3.3, height = 3.0, dpi = 300)
  write.csv(datp[order(datp$nominal, datp$rep), c("construct", "nominal", "measured", "rep", "E", "E_corr")],
            "denovo_enrichment_doseresponse_corrected.csv", row.names = FALSE)
  message("wrote denovo_enrichment_doseresponse_corrected.{pdf,csv}")
}
