#!/usr/bin/env Rscript
# Fig 1j: REAL 5-gly probe spike re-assignment. Recovered glyco fraction of the spike vs the
# nominal glyco dose, over the FULL five-point series (0/25/50/75/100 %), all three series at
# n = 3. EM (rPAL-seq pipeline) vs hard assignment (bowtie2 primary == Smith et al. eLife 2024
# "Fractional", i.e. split only co-best -- identical here) vs "fractional (all)", naive split
# over every reported alignment (abundance-blind, flat ~50 %).
#
# 2026-08-01/02, TWO CORRECTIONS to the 2026-07-31 version:
#   (1) the nominal-ZERO arm was missing. It had been OOM-KILLED in the original benchmark run
#       (job 22009179: ReqMem 6G, MaxRSS 6224292K -- it died on spk_5i, which is the LARGEST
#       library in the series because at nominal 0 every spike read is plain Fluc). The
#       production pipeline was unaffected, which is why Results already quote all five EM
#       values ("0.002, 3.4, 10.2, 25.9 and 99.9 %"). Nominal 0 is the single most informative
#       point: the one place where fractional (all) reports ~50 % against a true value of zero.
#   (2) the EM series was a mean of n = 3 libraries while the baselines were replicate 1 ALONE
#       (run_real_spike_check.sh ships SAMPLES="spk_1i .. spk_5i"). Mixing an n = 3 mean with
#       single-library values on the same axes confounds method with library.
#
# Both fixed by re-running the re-assignment over all 15 5-gly Input libraries as SLURM array
# 22120338 (--mem=64G, 10x the original). Nothing is hardcoded here any more: every plotted
# value is computed from a source table.
suppressMessages(library(ggplot2))
tp <- c("../theme_nature.R")
source(path.expand(tp[file.exists(path.expand(tp))][1]))
suppressMessages(library(paletteer))
luc <- as.character(paletteer_d("PrettyCols::Lucent"))

GLY <- "."   # repo: denovo_input_glycan_titration_summary.csv ships beside this script   # glycoFluc tree ONLY -- counts_analysis/ includes int-gly
nominal <- c(0, 25, 50, 75, 100)

## ---- 1. EM: the production pipeline, n = 3 per occupancy -------------------------------
# Read rather than hardcode, so this panel can never drift from the titration numbers the
# Results and the other Fig 1 spike panels quote.
em <- read.csv(file.path(GLY, "denovo_input_glycan_titration_summary.csv"),
               stringsAsFactors = FALSE)
em <- em[em$construct == "5-gly", ]
em <- em[order(em$nominal_pct), ]
stopifnot(identical(as.numeric(em$nominal_pct), as.numeric(nominal)), all(em$n == 3))

## ---- 2+3. Baselines: re-assignment of the same -k10 BAMs, n = 3 ------------------------
# Long format from analyze_real_spike_assignment.py: label / key / value, labels
# spk_<N>i_<dose>pct_rep<R>. `hard assignment` is computed from the RAW primary-alignment
# counts rather than from default_glyco_frac, which the tool rounds to 4 dp -- at nominal 0
# that rounds to 0.0000 and would lose the point entirely.
ra <- read.delim("real_spike_assignment_n3.tsv", header = FALSE,
                 col.names = c("label", "key", "value"), stringsAsFactors = FALSE)
w <- reshape(ra, idvar = "label", timevar = "key", direction = "wide")
names(w) <- sub("^value\\.", "", names(w))
w$dose <- as.numeric(sub("^spk_[0-9]+i_([0-9]+)pct_rep[0-9]+$", "\\1", w$label))
w$rep  <- as.integer(sub("^.*_rep([0-9]+)$", "\\1", w$label))
for (k in c("glyco_prim", "fluc_prim", "oneN_all_glyco_frac", "frac_top_glyco_frac",
            "spike_reads")) w[[k]] <- as.numeric(w[[k]])
stopifnot(nrow(w) == 15, all(table(w$dose) == 3), all(sort(unique(w$rep)) == 1:3),
          all(w$glyco_prim + w$fluc_prim == w$spike_reads))

w$hard    <- 100 * w$glyco_prim / (w$glyco_prim + w$fluc_prim)
w$fracall <- 100 * w$oneN_all_glyco_frac
# The two hard-assignment definitions must agree, or "plotted as one line" is a false claim.
stopifnot(max(abs(w$hard / 100 - w$frac_top_glyco_frac)) < 0.001)
# And replicate 1 must still reproduce the four published values.
r1 <- w[w$rep == 1, ]; r1 <- r1[order(r1$dose), ]
stopifnot(all.equal(round(r1$hard[r1$dose > 0], 2), c(2.42, 8.33, 22.37, 99.86)))

agg <- function(v) {
  a <- aggregate(v ~ w$dose, FUN = mean); a[order(a[[1]]), 2]
}
sdev <- function(v) {
  a <- aggregate(v ~ w$dose, FUN = sd);   a[order(a[[1]]), 2]
}

lv  <- c("EM", "hard", "fracall")
# Short labels so the key fits a 2-inch panel at house type. "hard assignment" is the
# manuscript's own term for default / fractional (top-AS); the full method names are spelled
# out in the figure legend.
lab <- c("EM", "hard assignment", "fractional (all)")
pal <- c(luc[1], luc[5], luc[3]); names(pal) <- lv; names(lab) <- lv

d <- rbind(
  data.frame(nominal, method = "EM",      val = em$mean_input_glycan_pct,
             sd = em$sd_input_glycan_pct),
  data.frame(nominal, method = "hard",    val = agg(w$hard),    sd = sdev(w$hard)),
  data.frame(nominal, method = "fracall", val = agg(w$fracall), sd = sdev(w$fracall)))
d$method <- factor(d$method, levels = lv)
stopifnot(!anyNA(d$val))

## ---- plot --------------------------------------------------------------------------------
# Fig 1j geometry, 2026-07-31: matched to Fig 1h (2.0 in wide). House type at theme_nature
# sizes throughout -- the space comes from dropping the title/subtitle, not from shrinking text.
PANEL_W <- 2.0
PANEL_H <- 2.2        # a little taller than h/k to pay for the key

p <- ggplot(d, aes(nominal, val, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3, colour = "grey55") +
  geom_line(linewidth = 0.5) + geom_point(size = 1.1) +
  scale_colour_manual(values = pal, labels = lab, name = NULL, breaks = lv) +
  # Axis anchored at 0: the series is a five-point titration and the origin is a measured
  # level, not padding. The additive 2 units keep the 0 % points off the panel edge.
  scale_x_continuous(breaks = nominal, limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.05), add = c(2, 0))) +
  coord_cartesian(ylim = c(0, 100), clip = "off") +
  # The dashed identity line is the reference the whole panel is read against; unlabelled it
  # is just a stray diagonal.
  # x = 85, not 97: at 97 the dashed line runs straight through the word.
  annotate("text", x = 85, y = 99, hjust = 1, vjust = 1, size = 8/.pt, colour = "grey45",
           label = "nominal") +
  labs(x = "Nominal glyco-Fluc (%)", y = "Recovered glyco (%)") +
  theme_nature +
  # Single column: a 2-column key put the right-hand entry against the page edge.
  guides(colour = guide_legend(ncol = 1)) +
  theme(legend.position = "bottom",
        legend.key.height = unit(7, "pt"), legend.key.width = unit(11, "pt"),
        legend.key.spacing.y = unit(0, "pt"),
        legend.margin = margin(t = -4, b = 0), legend.box.spacing = unit(1, "pt"))

ggsave("real_5gly_reassignment.pdf", p, width = PANEL_W, height = PANEL_H)
ggsave("real_5gly_reassignment.png", p, width = PANEL_W, height = PANEL_H, dpi = 300)

## ---- key-sharing variants (2026-08-02) -----------------------------------------------------
# Fig 1j and the candidate Fig 1l carry an IDENTICAL three-entry key. Printing it twice costs
# 2 x 25 pt of column height and reads as two unrelated panels. These two outputs let the pair
# share one horizontal key placed once, underneath both:
#   real_5gly_reassignment_nokey.pdf  (144 x 133 pt)  <- j, key suppressed
#   synth_snp_1nt_nokey.pdf           (144 x 133 pt)  <- l, from fig8_synth_snp.R
#   assignment_key_h.pdf              (one row)       <- the key, drawn once
p_nokey <- p + theme(legend.position = "none")
ggsave("real_5gly_reassignment_nokey.pdf", p_nokey, width = PANEL_W, height = 1.85)
ggsave("real_5gly_reassignment_nokey.png", p_nokey, width = PANEL_W, height = 1.85, dpi = 300)

# The key on its own, horizontal. Built from the SAME `pal`/`lab`/`lv` objects the panels use,
# so it cannot drift from them -- but with a neutral legend theme: the panel version carries
# legend.margin = margin(t = -4) to tuck the key under the plot, and standing alone that
# negative margin collapses the box to ~4 pt.
p_h  <- ggplot(d, aes(nominal, val, colour = method)) + geom_line() +
        scale_colour_manual(values = pal, labels = lab, name = NULL, breaks = lv) +
        theme_nature +
        guides(colour = guide_legend(nrow = 1)) +
        theme(legend.position = "bottom", legend.key.width = unit(11, "pt"))
gt   <- ggplotGrob(p_h)
i    <- grep("guide-box", vapply(gt$grobs, function(z) z$name, character(1)))
stopifnot(length(i) >= 1)
lg   <- gt$grobs[[i[1]]]
# The guide box is a gtable padded with 0.5null rows/cols. sum() on the unit vector, or
# grobHeight(), resolves those nulls to something meaningless (3 pt) -- convert ELEMENTWISE
# first, so the nulls go to 0 and only the real content row is counted.
kw   <- sum(grid::convertWidth(lg$widths,   "in", valueOnly = TRUE))
kh   <- sum(grid::convertHeight(lg$heights, "in", valueOnly = TRUE))
stopifnot(kw > 1.5, kh > 0.08)
pdf("assignment_key_h.pdf", width = kw, height = kh); grid::grid.draw(lg); dev.off()
png("assignment_key_h.png", width = kw, height = kh, units = "in", res = 300)
grid::grid.draw(lg); dev.off()
cat(sprintf("wrote horizontal key: %.2f x %.2f in (%.0f x %.0f pt)\n",
            kw, kh, kw * 72, kh * 72))

cat("wrote real_5gly_reassignment (all series n = 3, five occupancies)\n")
o <- reshape(d[, c("nominal", "method", "val")], idvar = "nominal", timevar = "method",
             direction = "wide")
s <- reshape(d[, c("nominal", "method", "sd")], idvar = "nominal", timevar = "method",
             direction = "wide")
print(merge(o, s, by = "nominal"), row.names = FALSE, digits = 4)

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
FINAL_H <- 1.45   # square, so k's coord_fixed partner matches
# THE BINDING CONSTRAINT AT THIS SIZE IS THE AXIS TITLES, not the data. "Nominal glyco-Fluc (%)"
# is ~1.30 in wide at 9 pt (measured: 195 px at 150 dpi in the 2.0 in export), and it is centred on
# a plot area only ~1.0 in wide, so it overhangs both sides and runs off a 1.45 in canvas. A panel
# this small therefore needs either shorter titles or ~1.8 in of width. Shorter titles is the cheap
# option and costs nothing: the drafted legend defines both axes in its first sentence.
# Also dropped: the in-panel "nominal" annotation (the legend says "Dashed line, nominal"), and two
# of the five tick labels.
p_final <- ggplot(d, aes(nominal, val, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3, colour = "grey55") +
  geom_line(linewidth = 0.5) + geom_point(size = 1.1) +
  scale_colour_manual(values = pal, labels = lab, name = NULL, breaks = lv) +
  scale_x_continuous(breaks = c(0, 50, 100), limits = c(0, 100),
                     expand = expansion(mult = c(0.05, 0.05), add = c(0, 0))) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Nominal (%)", y = "Recovered (%)") +
  theme_nature + theme(legend.position = "none")
ggsave("real_5gly_reassignment_final.pdf", p_final, width = FINAL_W, height = FINAL_H)
ggsave("real_5gly_reassignment_final.png", p_final, width = FINAL_W, height = FINAL_H, dpi = 300)
cat(sprintf("wrote real_5gly_reassignment_final: %.2f x %.2f in -- insert at 100%%\n", FINAL_W, FINAL_H))
