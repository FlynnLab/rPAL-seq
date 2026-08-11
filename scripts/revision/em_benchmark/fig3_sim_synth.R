#!/usr/bin/env Rscript
# Fig 3: simulated SYNTHETIC symmetric paralogs (200 nt, identical except one mid-molecule
# divergent "anchor" block). Recovered vs TRUE fraction, faceted by anchor size, read 35 nt.
# EM tracks truth; default and the eLife-faithful "fractional (top-AS)" coincide and compress
# toward 50:50; the naive "fractional (all)" collapses fully.
suppressMessages(library(ggplot2))
tp <- c("../theme_nature.R",  # repo copy
        "~/git/rPAL-seq-R/theme_nature/theme_nature.R",
        "~/Seqdata/git/rPAL_seq_R/theme_nature/theme_nature.R")
source(path.expand(tp[file.exists(path.expand(tp))][1]))
suppressMessages(library(paletteer))
luc <- as.character(paletteer_d("PrettyCols::Lucent"))

em <- read.delim("results_synth_grid_prod.tsv", check.names = FALSE)
em <- em[em$method == "EM+UMI" & em$read_len == 35,
         c("block_size", "ref1_frac_true", "recovered_ref1_frac")]
em$method <- "EM"
as <- read.delim("synth_assignment_rescored.tsv", check.names = FALSE)
as$block_size    <- as.integer(sub("^bs([0-9]+)_.*$", "\\1", as$tag))
as$ref1_frac_true <- as.numeric(sub("^.*_f([0-9.]+)_L.*$", "\\1", as$tag))
as$read_len      <- as.integer(sub("^.*_L([0-9]+)$", "\\1", as$tag))
as <- as[as$read_len == 35, c("block_size", "ref1_frac_true", "recovered_ref1_frac", "method")]
d <- rbind(em, as)

lv  <- c("EM", "default", "random", "fractional_top", "fractional_all")
lab <- c("EM", "default", "random", "fractional (top-AS)", "fractional (all)")
pal <- c(luc[1], luc[4], luc[2], luc[5], luc[3]); names(pal) <- lv; names(lab) <- lv
d <- d[d$method %in% lv, ]; d$method <- factor(d$method, levels = lv)
d$anchor <- factor(paste0(d$block_size, ifelse(d$block_size == 1, " nt (SNP)", " nt")),
                   levels = paste0(sort(unique(d$block_size)),
                                   ifelse(sort(unique(d$block_size)) == 1, " nt (SNP)", " nt")))

p <- ggplot(d, aes(ref1_frac_true, recovered_ref1_frac, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3, colour = "grey55") +
  geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.3, colour = "grey75") +
  geom_line(linewidth = 0.5) + geom_point(size = 1.1) +
  facet_wrap(~ anchor, nrow = 1) +
  scale_colour_manual(values = pal, labels = lab, name = NULL) +
  scale_x_continuous(breaks = c(0.5, 0.75, 0.9)) +
  coord_fixed(ratio = 1, xlim = c(0.48, 0.92), ylim = c(0.45, 1.0)) +
  labs(title = "Synthetic paralogs (anchor sweep)", subtitle = "read 35 nt; production EM (delta_as=2); dashed = truth",
       x = "true fraction", y = "recovered fraction") +
  theme_nature +
  theme(legend.position = "bottom", legend.key.width = unit(10, "pt"),
        legend.text = element_text(size = 6.5),
        legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
        panel.spacing = unit(4, "pt"), strip.text = element_text(size = 8),
        plot.subtitle = element_text(size = 8, colour = "grey35"))

ggsave("sim_synth_anchor.pdf", p, width = 5.8, height = 3.0)
ggsave("sim_synth_anchor.png", p, width = 5.8, height = 3.0, dpi = 200)
cat("wrote sim_synth_anchor\n")
