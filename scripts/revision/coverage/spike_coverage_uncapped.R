#!/usr/bin/env Rscript
# Spike 5' coverage, recomputed from UNCAPPED pileups (2026-08-03).
#
# WHY THIS EXISTS. The pipeline's cpup step runs `samtools mpileup -d $mpileup_max_depth` with
# mpileup_max_depth: 10000. On the 51-nt glyco_Fluc_spike the true depth is up to 1,141,186, so
# EVERY position sat at the ceiling and the coverage profile was flat by construction -- the
# published Input value of 1.02 was the cap, not the molecule. Re-ran with `-d 0`
# (SLURM 22142991, spike reads only, otherwise the pipeline's own flags) and the effect is real,
# glycan-specific and far stronger: Input 0.742 -> IP 0.616, paired t = -8.00, P = 6.5e-06,
# 12/12 libraries, Cohen dz = -2.31.
#
# The metric is UNCHANGED from spike_metric.r:30 / meta_depth_spike.r:37 --
#   cov_norm = depth / mean(depth) per library x reference ; five_retention = mean(pos <= 5).
# Only the input changes. Reading the pileups directly rather than through counts_long.csv,
# because the parse_pileup step would have to be re-run to regenerate that.
#
# The UNCONJUGATED partner (Fluc_spike, 89 nt, same tube, same RT, no glycan) is the
# within-library control the capped data could not support. It is the point of the figure.
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
tp <- c("../theme_nature.R")
source(path.expand(tp[file.exists(path.expand(tp))][1]))

PILE   <- "uncapped_pileups"
FIVE_N <- 5                                   # == spike_coverage.five_nt
COND   <- c(Input = "#8C2981", IP = "#E01A4F")   # == plots.colors.conditions
# Short strip labels: the full names clipped at panel width. The legend spells them out.
REFS   <- c(glyco_Fluc_spike = "glyco-Fluc", Fluc_spike = "Fluc, unconjugated")
DOSE <- c(`1`=100,`2`=75,`3`=50,`4`=25,`6`=100,`7`=75,`8`=50,`9`=25,
          `11`=100,`12`=75,`13`=50,`14`=25)
REP  <- c(`1`=1,`2`=1,`3`=1,`4`=1,`6`=2,`7`=2,`8`=2,`9`=2,`11`=3,`12`=3,`13`=3,`14`=3)

files <- list.files(PILE, pattern = "\\.d0\\.txt$", full.names = TRUE)
stopifnot(length(files) == 24)

rows <- do.call(rbind, lapply(files, function(f) {
  lib <- sub("\\.d0\\.txt$", "", basename(f))
  num <- gsub("[^0-9]", "", sub("^spk_", "", lib))
  x <- read.delim(f, stringsAsFactors = FALSE)
  x$depth <- as.numeric(sub(",.*$", "", x[[4]]))
  x <- x[x$chr %in% names(REFS), c("chr", "pos", "depth")]
  x$lib <- lib
  x$condition <- ifelse(grepl("i$", lib), "Input", "IP")
  x$dose <- DOSE[[num]]; x$rep <- REP[[num]]
  x
}))
rows <- rows %>% group_by(lib, chr) %>% mutate(cov_norm = depth / mean(depth)) %>% ungroup()
rows$condition <- factor(rows$condition, levels = c("Input", "IP"))
rows$ref <- factor(unname(REFS[rows$chr]), levels = unname(REFS))

## ---- panel 1: per-position ribbon, both references ---------------------------------------
rib <- rows %>% group_by(ref, condition, pos) %>%
  summarise(m = mean(cov_norm), se = sd(cov_norm)/sqrt(n()), .groups = "drop")
p_rib <- ggplot(rib, aes(pos, m, colour = condition, fill = condition)) +
  geom_ribbon(aes(ymin = m - se, ymax = m + se), alpha = 0.25, colour = NA) +
  geom_line(linewidth = 0.5) +
  facet_wrap(~ ref, scales = "free_x") +
  scale_colour_manual(values = COND, name = NULL) +
  scale_fill_manual(values = COND, guide = "none") +
  scale_x_continuous(breaks = c(1, 20, 40, 60, 80)) +
  labs(x = "Position (nt)", y = "Normalised coverage") +
  theme_nature +
  theme(legend.position = "bottom", legend.key.width = unit(11, "pt"),
        legend.margin = margin(t = -4, b = 0), legend.box.spacing = unit(1, "pt"),
        strip.text = element_text(size = 8))

## ---- panel 2: 5'-retention box, both references ------------------------------------------
met <- rows %>% filter(pos <= FIVE_N) %>%
  group_by(ref, lib, condition, dose, rep) %>%
  summarise(five_retention = mean(cov_norm), .groups = "drop")
# paired test per reference: the libraries are matched Input/IP pairs, so pair them.
stats <- met %>% tidyr::pivot_wider(id_cols = c(ref, dose, rep),
                                    names_from = condition, values_from = five_retention) %>%
  group_by(ref) %>%
  summarise(p = t.test(IP, Input, paired = TRUE)$p.value,
            d = mean(IP - Input), n = dplyr::n(), .groups = "drop")
print(as.data.frame(stats))

# House style (2026-08-04): significance ASTERISKS on the panel, exact P in the Results text.
# Thresholds match the manuscript's own key -- "**P < 0.01, ****P < 0.0001" -- so the tiers are
# the standard * <0.05, ** <0.01, *** <0.001, **** <0.0001, and "n.s." is the wording already
# used in the legends ("n.s., not significant").
stars <- function(p) ifelse(p < 1e-4, "****", ifelse(p < 1e-3, "***",
                     ifelse(p < 1e-2, "**", ifelse(p < 5e-2, "*", "n.s."))))
# Bracket sits above each facet's OWN maximum, not a shared height: the two references differ by
# ~0.2 in retention, so one common height would float far above the glyco-Fluc boxes.
top <- met %>% group_by(ref) %>% summarise(ymax = max(five_retention), .groups = "drop")
lab <- stats %>% left_join(top, by = "ref") %>%
  mutate(txt = stars(p),
         yb  = ymax + 0.035,                  # bracket bar
         yt  = yb + 0.012,                    # asterisks / n.s.
         tick = 0.014)                        # downward end ticks
p_box <- ggplot(met, aes(condition, five_retention, colour = condition)) +
  # staplewidth draws the T-caps on the whiskers, matching the other box panels in the paper.
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.35, staplewidth = 0.35) +
  geom_jitter(width = 0.12, size = 0.9) +
  geom_segment(data = lab, aes(x = 1, xend = 2, y = yb, yend = yb),
               inherit.aes = FALSE, linewidth = 0.3, colour = "grey30") +
  geom_segment(data = lab, aes(x = 1, xend = 1, y = yb, yend = yb - tick),
               inherit.aes = FALSE, linewidth = 0.3, colour = "grey30") +
  geom_segment(data = lab, aes(x = 2, xend = 2, y = yb, yend = yb - tick),
               inherit.aes = FALSE, linewidth = 0.3, colour = "grey30") +
  geom_text(data = lab, aes(x = 1.5, y = yt, label = txt), inherit.aes = FALSE,
            size = 8/.pt, colour = "grey20", vjust = 0) +
  facet_wrap(~ ref) +
  scale_colour_manual(values = COND, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL, y = sprintf("5' coverage retention (%d nt)", FIVE_N)) +
  theme_nature + theme(strip.text = element_text(size = 8))

for (nm in c("spike_coverage_ribbon_uncapped", "spike_coverage_metric_box_uncapped")) {
  p <- if (grepl("ribbon", nm)) p_rib else p_box
  w <- if (grepl("ribbon", nm)) 3.4 else 3.0
  ggsave(paste0(nm, ".pdf"), p, width = w, height = 2.0)
  ggsave(paste0(nm, ".png"), p, width = w, height = 2.0, dpi = 300)
}
write.csv(met, "spike_coverage_metric_uncapped.csv", row.names = FALSE)
cat("wrote spike_coverage_{ribbon,metric_box}_uncapped.{pdf,png} + metric CSV\n")
print(met %>% group_by(ref, condition) %>%
      summarise(mean = round(mean(five_retention), 4), sd = round(sd(five_retention), 4),
                .groups = "drop") %>% as.data.frame())
