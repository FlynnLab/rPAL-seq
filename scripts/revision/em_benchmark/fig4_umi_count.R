#!/usr/bin/env Rscript
# Fig 4: the UMI contribution. Recovered count / TRUE molecule count under PCR duplication
# (geometric, mean 4x). Only UMI-aware EM returns the molecule count (~1x); every non-UMI
# method -- including eLife-style EM without UMIs -- counts PCR duplicates as molecules (~4x).
# Uses the synthetic grid (all anchor sizes / ratios / read lengths) as replicate cells.
suppressMessages({library(ggplot2); library(dplyr)})
tp <- c("../theme_nature.R",  # repo copy
        "~/git/rPAL-seq-R/theme_nature/theme_nature.R",
        "~/Seqdata/git/rPAL_seq_R/theme_nature/theme_nature.R")
source(path.expand(tp[file.exists(path.expand(tp))][1]))
suppressMessages(library(paletteer))
luc <- as.character(paletteer_d("PrettyCols::Lucent"))

d <- read.delim("results_synth_grid_prod.tsv", check.names = FALSE)
# drop the no-pre-filter EM arm: it is deduplicated like EM+UMI, so bucketing it with the
# read-level methods would pull their mean away from the true PCR factor.
d <- d[d$method != "EM_unfiltered", ]
# EM+UMI vs EM_noUMI carry the real distinction; the read-level baselines all equal the raw
# read count (no dedup), so collapse them into one "read-level methods" series.
d$grp <- ifelse(d$method == "EM+UMI", "EM + UMI",
         ifelse(d$method == "EM_noUMI", "EM, no UMI", "read-level (no dedup)"))
d <- d %>% filter(!is.na(count_over_molecules)) %>%
  group_by(grp, read_len) %>%
  summarise(m = mean(count_over_molecules), sd = sd(count_over_molecules), .groups = "drop")
lv  <- c("EM + UMI", "EM, no UMI", "read-level (no dedup)")
pal <- c(luc[1], luc[2], luc[4]); names(pal) <- lv
d$grp <- factor(d$grp, levels = lv)

p <- ggplot(d, aes(factor(read_len), m, fill = grp)) +
  geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3, colour = "grey45") +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_errorbar(aes(ymin = pmax(0, m - sd), ymax = m + sd),
                position = position_dodge(width = 0.75), width = 0.18, linewidth = 0.25) +
  annotate("text", x = 0.55, y = 1, label = "truth", vjust = -0.5, hjust = 0,
           size = 2.3, colour = "grey45") +
  scale_fill_manual(values = pal, name = NULL) +
  labs(title = "UMI dedup recovers molecule counts", subtitle = "PCR mean 4x; mean +/- SD over grid cells",
       x = "read length (nt)", y = "count / true molecules") +
  theme_nature +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.key.height = unit(7, "pt"), legend.key.width = unit(7, "pt"),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
        legend.text = element_text(size = 6.5),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"))

ggsave("umi_count_recovery.pdf", p, width = 3.4, height = 3.1)
ggsave("umi_count_recovery.png", p, width = 3.4, height = 3.1, dpi = 200)
cat("wrote umi_count_recovery\n")
print(as.data.frame(d), row.names = FALSE)
