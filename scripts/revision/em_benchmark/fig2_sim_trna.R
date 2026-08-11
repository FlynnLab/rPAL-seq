#!/usr/bin/env Rscript
# Fig 2: simulated REAL tRNA pair (tRNA-Tyr-GTA1 vs GTA2, single SNP), true ratio 0.75.
# Recovered fraction vs read length. Baselines follow Smith et al. (eLife 2024): "fractional
# (top-AS)" splits only the read's co-best alignments (= their Fractional, filtered) and ~=
# default; "fractional (all)" is the naive split-every-alignment strawman. EM uses abundance.
suppressMessages(library(ggplot2))
tp <- c("../theme_nature.R",  # repo copy
        "~/git/rPAL-seq-R/theme_nature/theme_nature.R",
        "~/Seqdata/git/rPAL_seq_R/theme_nature/theme_nature.R")
source(path.expand(tp[file.exists(path.expand(tp))][1]))
suppressMessages(library(paletteer))
luc <- as.character(paletteer_d("PrettyCols::Lucent"))   # crimson orange yellow cyan teal

# production EM (delta_as=2, max_iters 40, tol 1e-3), rescored on the same BAMs
em <- read.delim("prodEM_em_paralog_grid.tsv", check.names = FALSE)
em <- data.frame(read_len = as.integer(sub(".*_L", "", em$tag)),
                 recovered_ref1_frac = as.numeric(em$recovered_ref1_frac),
                 method = "EM")
as <- read.delim("paralog_assignment_rescored.tsv", check.names = FALSE)
as$read_len <- as.integer(sub(".*_L", "", as$tag))
d <- rbind(em, as[, c("read_len", "recovered_ref1_frac", "method")])

lv  <- c("EM", "default", "random", "fractional_top", "fractional_all")
lab <- c("EM", "default", "random", "fractional (top-AS)", "fractional (all)")
pal <- c(luc[1], luc[4], luc[2], luc[5], luc[3]); names(pal) <- lv; names(lab) <- lv
d <- d[d$method %in% lv, ]; d$method <- factor(d$method, levels = lv)
truth <- 0.75

p <- ggplot(d, aes(read_len, recovered_ref1_frac, colour = method)) +
  geom_hline(yintercept = truth, linetype = 2, linewidth = 0.3, colour = "grey45") +
  geom_line(linewidth = 0.5) + geom_point(size = 1.3) +
  annotate("text", x = min(d$read_len), y = truth, label = "truth", vjust = -0.5,
           hjust = 0, size = 2.3, colour = "grey45") +
  scale_colour_manual(values = pal, labels = lab, name = NULL) +
  scale_x_continuous(breaks = sort(unique(d$read_len))) +
  labs(title = "Simulated tRNA pair (1 SNP)", subtitle = "true ratio 0.75; production EM (delta_as=2)",
       x = "read length (nt)", y = "recovered fraction") +
  theme_nature +
  theme(legend.position = c(0.99, 0.02), legend.justification = c(1, 0),
        legend.key.height = unit(8, "pt"), legend.background = element_blank(),
        legend.text = element_text(size = 6.5),
        plot.subtitle = element_text(size = 8, colour = "grey35"))

ggsave("sim_trna_1snp.pdf", p, width = 3.3, height = 3.0)
ggsave("sim_trna_1snp.png", p, width = 3.3, height = 3.0, dpi = 200)
cat("wrote sim_trna_1snp\n")
