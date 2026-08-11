#!/usr/bin/env Rscript
# Fig 5: modification-induced mismatch arm. RT misincorporation at the discriminating
# position converts cleanly-resolvable reads into either genuine ties (mod-mode "tie":
# substituted base matches NEITHER reference) or actively misleading reads (mod-mode
# "decoy": substituted base matches the PARALOG). x = site occupancy q; y = recovered
# major-paralog fraction, scored against the EMITTED truth (dashed) so that fragment
# dropout is not misattributed to the estimator. All seven methods shown.
suppressMessages({library(ggplot2); library(dplyr)})
tp <- c("../theme_nature.R")
source(path.expand(tp[file.exists(path.expand(tp))][1]))
suppressMessages(library(paletteer))
luc <- as.character(paletteer_d("PrettyCols::Lucent"))   # crimson orange yellow cyan teal
rel <- as.character(paletteer_d("PrettyCols::Relax"))

d <- read.delim("results_mod_grid_prod.tsv", check.names = FALSE)
d$q <- as.numeric(d$q)
d$recovered_ref1_frac <- as.numeric(d$recovered_ref1_frac)
d$ref1_frac_true <- as.numeric(d$ref1_frac_true)

lv  <- c("EM+UMI", "default", "random", "fractional_top", "fractional_all", "unique")
lab <- c("EM (UMI)", "default", "random",
         "fractional (top-AS)", "fractional (all)", "unique (discards)")
pal <- c(luc[1], luc[4], rel[5], luc[5], luc[3], rel[1])
names(pal) <- lv; names(lab) <- lv
d <- d[d$method %in% lv, ]
d <- d[d$mod_mode == "tie", ]   # main figure: symmetric arm only
d$method <- factor(d$method, levels = lv)
d$mode_lab <- factor(ifelse(d$mod_mode == "tie", "tie: anchor destroyed",
                            "decoy: mimics paralog"),
                     levels = c("tie: anchor destroyed", "decoy: mimics paralog"))
# nominal ratio is not a column (only the realized draw): snap to the nearest requested value
d$pn <- vapply(as.numeric(d$ref1_frac_drawn),
               function(x) c(0.75, 0.90)[which.min(abs(c(0.75, 0.90) - x))], numeric(1))
d$panel <- factor(sprintf("p %.2f, %d nt", d$pn, d$read_len))

# truth line per facet (emitted truth shifts slightly with q through dropout)
tru <- d %>% group_by(mode_lab, panel) %>%
  summarise(p = mean(ref1_frac_true), .groups = "drop")

p1 <- ggplot(d, aes(q, recovered_ref1_frac, colour = method)) +
  geom_hline(data = tru, aes(yintercept = p), linetype = 2,
             linewidth = 0.3, colour = "grey45", inherit.aes = FALSE) +
  geom_hline(yintercept = 0.5, linetype = 3, linewidth = 0.25, colour = "grey80") +
  geom_line(linewidth = 0.45) + geom_point(size = 1.0) +
  facet_wrap(~ panel, nrow = 2) +
  scale_colour_manual(values = pal, labels = lab, name = NULL) +
  scale_x_continuous(breaks = sort(unique(d$q))) +
  labs(title = "Modification-induced mismatch destroys the discriminator",
       subtitle = "mismatch matches neither reference; production EM (delta_as=2); dashed = emitted truth",
       x = "site occupancy q", y = "recovered fraction") +
  theme_nature +
  theme(legend.position = "bottom", legend.key.width = unit(10, "pt"),
        legend.text = element_text(size = 6), legend.box.spacing = unit(2, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        panel.spacing = unit(3.5, "pt"), strip.text = element_text(size = 7),
        plot.subtitle = element_text(size = 7.5, colour = "grey35")) +
  guides(colour = guide_legend(nrow = 2))

ggsave("mod_mismatch_tie.pdf", p1, width = 5.0, height = 4.2)
ggsave("mod_mismatch_tie.png", p1, width = 5.0, height = 4.2, dpi = 190)

# --- slope table: the pre-registered endpoint (d recovered / d q) per method and mode ---
sl <- d %>% group_by(mod_mode, method, read_len, pn) %>%
  filter(n() > 2) %>%
  summarise(slope = coef(lm(recovered_ref1_frac ~ q))[2], .groups = "drop") %>%
  group_by(mod_mode, method) %>%
  summarise(mean_slope = round(mean(slope), 4), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = mod_mode, values_from = mean_slope)
cat("\n== mean d(recovered)/dq by method ==\n"); print(as.data.frame(sl), row.names = FALSE)

# --- |assign_err| averaged over the grid, per mode ---
err <- d %>% mutate(ae = abs(as.numeric(assign_err))) %>%
  group_by(mod_mode, method) %>% summarise(MAE = round(mean(ae), 4), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = mod_mode, values_from = MAE)
cat("\n== mean |recovered - emitted truth| ==\n"); print(as.data.frame(err), row.names = FALSE)
cat("\nwrote mod_mismatch_tie.pdf / .png\n")
