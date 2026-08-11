#!/usr/bin/env Rscript
# CANDIDATE main-figure panel (2026-08-02, author's RE2): the 1-nt (SNP) facet of the
# synthetic anchor sweep, rebuilt at Fig 1 size.
#
# WHY THIS RATHER THAN umi_count_recovery. The author's point is the right one: a reader
# without a computational background does not doubt that UMI dedup recovers molecule counts
# -- that result is expected. What they doubt is the EM. And Fig 1j, despite being the
# real-data arm, does NOT show EM beating hard assignment: there EM and hard assignment
# track each other (3.43/10.15/25.89/99.95 vs 2.63/8.07/22.44/99.84) and both separate from
# fractional (all). So no main figure currently answers the reviewer's actual question --
# how much better is the EM than what people normally do. This panel answers it.
#
# WHY THE 1-nt FACET AND NOT THE 3-FACET SWEEP. The sweep barely moves: hard assignment at
# true 0.9 reads 0.579 / 0.587 / 0.605 at 1 / 4 / 12 nt and EM reads 0.888 / 0.891 / 0.887.
# All three facets tell one story, so the sweep is a completeness check (keep it in ED) and
# a single facet carries the main figure. 1 nt is the right one: it is the hardest case, it
# has the largest EM-vs-hard gap, and it is what real tRNA isodecoders look like.
#
# SERIES COLLAPSED TO MATCH Fig 1j. default / random / fractional (top-AS) coincide (max
# spread 0.0040 over the 12 cells at 35 nt), so they are drawn as one "hard assignment" line
# -- the same three-entry key as Fig 1j, so the pair reads as one unit and the key is
# learned once.
suppressMessages(library(ggplot2))
tp <- c("../theme_nature.R")
source(path.expand(tp[file.exists(path.expand(tp))][1]))
suppressMessages(library(paletteer))
luc <- as.character(paletteer_d("PrettyCols::Lucent"))

BS  <- 1       # divergent block size (nt)
RL  <- 35      # read length

em <- read.delim("results_synth_grid_prod.tsv", check.names = FALSE)
em <- em[em$method == "EM+UMI" & em$read_len == RL & em$block_size == BS,
         c("ref1_frac_true", "recovered_ref1_frac")]
em$method <- "EM"

as <- read.delim("synth_assignment_rescored.tsv", check.names = FALSE)
as$block_size     <- as.integer(sub("^bs([0-9]+)_.*$", "\\1", as$tag))
as$ref1_frac_true <- as.numeric(sub("^.*_f([0-9.]+)_L.*$", "\\1", as$tag))
as$read_len       <- as.integer(sub("^.*_L([0-9]+)$", "\\1", as$tag))
as <- as[as$read_len == RL & as$block_size == BS, ]

# Collapse the three coinciding hard-assignment methods into one line, and assert that they
# really do coincide before doing it -- if a future re-run separates them this must fail.
hardset <- c("default", "random", "fractional_top")
h <- as[as$method %in% hardset, c("ref1_frac_true", "recovered_ref1_frac", "method")]
sp <- tapply(h$recovered_ref1_frac, h$ref1_frac_true, function(v) max(v) - min(v))
stopifnot(all(sp < 0.01))
hard <- aggregate(recovered_ref1_frac ~ ref1_frac_true, h, mean)
hard$method <- "hard"

fa <- as[as$method == "fractional_all", c("ref1_frac_true", "recovered_ref1_frac")]
fa$method <- "fracall"

d <- rbind(em, hard, fa)
lv  <- c("EM", "hard", "fracall")
lab <- c("EM", "hard assignment", "fractional (all)")
pal <- c(luc[1], luc[5], luc[3])          # identical to Fig 1j
names(pal) <- lv; names(lab) <- lv
d$method <- factor(d$method, levels = lv)

LIM <- c(0.47, 0.93)                       # equal spans so coord_fixed gives a true 45 deg
stopifnot(all(d$recovered_ref1_frac > LIM[1]), all(d$recovered_ref1_frac < LIM[2]))

base <- ggplot(d, aes(ref1_frac_true, recovered_ref1_frac, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3, colour = "grey55") +
  geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.3, colour = "grey75") +
  geom_line(linewidth = 0.5) + geom_point(size = 1.1) +
  scale_colour_manual(values = pal, labels = lab, name = NULL, breaks = lv) +
  # Both axes forced to 2 dp: x carries 0.75, so leaving y to default gave "0.9" against
  # "0.90" on the same panel.
  scale_x_continuous(breaks = c(0.5, 0.75, 0.9), labels = function(v) sprintf("%.2f", v)) +
  scale_y_continuous(breaks = c(0.5, 0.7, 0.9), labels = function(v) sprintf("%.2f", v)) +
  coord_fixed(ratio = 1, xlim = LIM, ylim = LIM, clip = "off") +
  # Both annotations sit above the diagonal, which is empty by construction: no method can
  # recover more than the true fraction.
  annotate("text", x = 0.475, y = 0.925, hjust = 0, vjust = 1, size = 8/.pt,
           colour = "grey45", label = "1-nt divergence") +
  # Hugging the dashed line, or it reads as a floating word rather than its label.
  annotate("text", x = 0.665, y = 0.700, hjust = 0, vjust = 0, size = 8/.pt,
           colour = "grey45", label = "truth") +
  labs(x = "True fraction", y = "Recovered fraction") +
  theme_nature

key <- theme(legend.position = "bottom",
             legend.key.height = unit(7, "pt"), legend.key.width = unit(11, "pt"),
             legend.key.spacing.y = unit(0, "pt"),
             legend.margin = margin(t = -4, b = 0), legend.box.spacing = unit(1, "pt"))

# Two versions. The nokey one is the intended main-figure form -- placed beside Fig 1j it
# shares that panel's key, which is identical, and saves the 14 pt the key costs.
p_key <- base + guides(colour = guide_legend(ncol = 1)) + key
p_no  <- base + theme(legend.position = "none")

ggsave("synth_snp_1nt.pdf",       p_key, width = 2.0, height = 2.2)
ggsave("synth_snp_1nt.png",       p_key, width = 2.0, height = 2.2, dpi = 300)
ggsave("synth_snp_1nt_nokey.pdf", p_no,  width = 2.0, height = 1.85)
ggsave("synth_snp_1nt_nokey.png", p_no,  width = 2.0, height = 1.85, dpi = 300)
cat("wrote synth_snp_1nt (+ _nokey)\n")
print(reshape(d, idvar = "ref1_frac_true", timevar = "method", direction = "wide"),
      row.names = FALSE)

## ---- placed-size export (2026-08-03) --------------------------------------------------------
# Measured from the author's 2026-08-03 figure: the 2.0 in panels were being inserted at ~0.68
# scale (j's x-axis rule is 216 px in this PDF and 147 px in the figure; the identical
# "Nominal glyco-Fluc (%)" string is 196 px here and 130 px there -- two independent measures
# that agree). That drops theme_nature's 8 pt axis text to ~5.4 pt on the page, and ~5.2 pt once
# the 190 mm figure is reduced to 180 mm, against Nature's 5 pt floor. Exporting at the size the
# panel will OCCUPY, and inserting at 100 %, keeps the house type size instead of scaling it away.
# 1.45 in each so the j,k pair spans the 2.97 in of columns 3+4. Change these two numbers if the
# pair is re-flowed; nothing else depends on them.
FINAL_W <- 1.45
FINAL_H <- 1.45   # 1.35 clipped the rotated y title, which is ~1.06 in at 9 pt
# Annotations dropped rather than moved -- "1-nt divergence" alone is ~1.0 in at 8 pt against a
# ~0.95 in plot area here, and both facts are already in the drafted legend ("simulated pairs of
# equal-length 200-nt references differing at a single position"; "Dashed line, truth; dotted line,
# an even split").
p_final <- ggplot(d, aes(ref1_frac_true, recovered_ref1_frac, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3, colour = "grey55") +
  geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.3, colour = "grey75") +
  geom_line(linewidth = 0.5) + geom_point(size = 1.1) +
  scale_colour_manual(values = pal, labels = lab, name = NULL, breaks = lv) +
  scale_x_continuous(breaks = c(0.5, 0.7, 0.9), labels = function(v) sprintf("%.2f", v)) +
  scale_y_continuous(breaks = c(0.5, 0.7, 0.9), labels = function(v) sprintf("%.2f", v)) +
  coord_cartesian(xlim = LIM, ylim = LIM) +
  labs(x = "True fraction", y = "Recovered fraction") +
  theme_nature + theme(legend.position = "none")
ggsave("synth_snp_1nt_final.pdf", p_final, width = FINAL_W, height = FINAL_H)
ggsave("synth_snp_1nt_final.png", p_final, width = FINAL_W, height = FINAL_H, dpi = 300)
cat(sprintf("wrote synth_snp_1nt_final: %.2f x %.2f in -- insert at 100%%\n", FINAL_W, FINAL_H))
