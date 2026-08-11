#!/usr/bin/env Rscript
# Supplementary: the WORST CASE for EM, and its boundary.
# "Decoy" = a modification-induced mismatch at a discriminating base that converts it to the
# PARALOG's base. When the two references differ at only ONE base this makes the read a
# perfect, full-length match to the wrong reference: the evidence is destroyed, not merely
# weakened, and no estimator can recover it. EM does worst because its abundance prior also
# drags the non-anchor tied reads the same way, so a linearly corrupted likelihood becomes a
# super-linearly corrupted estimate.
# The boundary: with k divergent bases, one lesion removes only 1/k of the evidence. Sweeping
# k = 1, 2, 3 shows the effect collapsing (EM slope -1.29 -> -0.57 -> -0.13), so the failure is
# confined to near-identical (1-2 base) pairs at high occupancy.
suppressMessages({library(ggplot2); library(dplyr)})
tp <- c("../theme_nature.R",  # repo copy
        "~/git/rPAL-seq-R/theme_nature/theme_nature.R",
        "~/Seqdata/git/rPAL_seq_R/theme_nature/theme_nature.R")
source(path.expand(tp[file.exists(path.expand(tp))][1]))
suppressMessages(library(paletteer))
luc <- as.character(paletteer_d("PrettyCols::Lucent"))
rel <- as.character(paletteer_d("PrettyCols::Relax"))

d <- read.delim("results_mod_divergence_prod.tsv", check.names = FALSE)
d$q <- as.numeric(d$q)
d$recovered_ref1_frac <- as.numeric(d$recovered_ref1_frac)
d$ref1_frac_true <- as.numeric(d$ref1_frac_true)

lv  <- c("EM+UMI", "default", "random", "fractional_top", "fractional_all", "unique")
lab <- c("EM (UMI)", "default", "random",
         "fractional (top-AS)", "fractional (all)", "unique (discards)")
pal <- c(luc[1], luc[4], rel[5], luc[5], luc[3], rel[1])
names(pal) <- lv; names(lab) <- lv
d <- d[d$method %in% lv, ]
d$method <- factor(d$method, levels = lv)
d$div <- factor(sprintf("%d divergent base%s", d$block_size,
                        ifelse(d$block_size == 1, "", "s")),
                levels = sprintf("%d divergent base%s", 1:3, c("", "s", "s")))
tru <- d %>% group_by(div) %>% summarise(p = mean(ref1_frac_true), .groups = "drop")

p1 <- ggplot(d, aes(q, recovered_ref1_frac, colour = method)) +
  geom_hline(data = tru, aes(yintercept = p), linetype = 2, linewidth = 0.3,
             colour = "grey45", inherit.aes = FALSE) +
  geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.25, colour = "grey80") +
  geom_line(linewidth = 0.45) + geom_point(size = 1.0) +
  facet_wrap(~ div, nrow = 1) +
  scale_colour_manual(values = pal, labels = lab, name = NULL) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75),
                     labels = c("0", ".25", ".5", ".75"), expand = expansion(mult = 0.09)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Worst case: mismatch converts the base to the paralog's",
       subtitle = "production EM (delta_as=2); full mimicry needs 1-base divergence and collapses by 3 (p 0.75, 35 nt)",
       x = "site occupancy q", y = "recovered fraction") +
  theme_nature +
  theme(legend.position = "bottom", legend.key.width = unit(10, "pt"),
        legend.text = element_text(size = 6), legend.box.spacing = unit(2, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        panel.spacing = unit(3.5, "pt"), strip.text = element_text(size = 7.5),
        plot.subtitle = element_text(size = 7, colour = "grey35")) +
  guides(colour = guide_legend(nrow = 2))

ggsave("mod_decoy_divergence.pdf", p1, width = 5.8, height = 3.2)
ggsave("mod_decoy_divergence.png", p1, width = 5.8, height = 3.2, dpi = 190)

sl <- d %>% group_by(block_size, method) %>%
  summarise(slope = round(coef(lm(recovered_ref1_frac ~ q))[2], 4), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = block_size, values_from = slope,
                     names_prefix = "div")
cat("\n== decoy slope d(recovered)/dq by divergence ==\n")
print(as.data.frame(sl), row.names = FALSE)
cat("\nwrote mod_decoy_divergence.pdf / .png\n")
