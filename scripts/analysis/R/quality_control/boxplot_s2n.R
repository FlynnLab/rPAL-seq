library(ggplot2)
library(dplyr)
library(ggpubr)

# ---- Inputs ----
date_str   <- "date_of_your_tapestation_run"
input_file <- sprintf("tapequant_%s.csv", date_str)

# ---- Load ----
data <- read.csv(input_file, stringsAsFactors = FALSE)

# ---- Keep relevant groups & map display labels ----
type_map <- c(
  "I"   = "Input",
  "Vp"  = "Bulk IP",
  "Vh"  = "Bulk Control",
  "dVp" = "Low-input IP",
  "dVh" = "Low-input Control"
)

order_levels <- c("Input", "Bulk IP", "Bulk Control", "Low-input IP", "Low-input Control")

data <- data %>%
  filter(!is.na(Molarity_pM), !is.na(Type_plot)) %>%
  filter(Type_plot %in% names(type_map)) %>%
  mutate(Type_Display = dplyr::recode(Type_plot, !!!type_map)) %>%
  mutate(Type_Display = factor(Type_Display, levels = order_levels))

# ---- Comparisons (VC vs HI at both dose levels) ----
comparisons <- list(
  c("Bulk IP", "Bulk Control"),
  c("Low-input IP", "Low-input Control")
)

# ---- Aesthetics to resemble your example ----
fill_pal <- c(
  "Input"             = "#8C2981FF",
  "Bulk IP"           = "#E01A4FFF",
  "Low-input IP"      = "#F9C22EFF",
  "Bulk Control"      = "#53B3CBFF",
  "Low-input Control" = "#7DCFB6FF"
)

# ---- Set axis caps ----
y_min <- 1.5
y_max <- 3.8

# ---- Identify capped points ----
data <- data %>%
  mutate(
    logM = log10(Molarity_pM + 1e-6),
    is_low  = logM < y_min,
    is_high = logM > y_max
  )

# ---- Plot with capped axis and custom outlier markers ----
p <- ggplot(
  data,
  aes(x = Type_Display, y = logM)
) +
  geom_boxplot(
    aes(fill = Type_Display),
    width = 0.5,
    alpha = 0.7,
    outlier.shape = NA,
    color = "black"
  ) +
  # Jitter only in-range points
  geom_jitter(
    data = subset(data, !is_low & !is_high),
    shape = 19,
    size = 1.5,
    width = 0.25,
    color = "black",
  ) +
  # ---- CHANGED: Hollow markers for capped low outliers (size like the volcano plot) ----
geom_point(
  data = subset(data, is_low),
  aes(x = Type_Display, y = y_min),
  shape = 21,
  size = 3,        # match example
  stroke = 0.5,    # lighter outline like your volcano plot
  fill = "white",
  color = "black"
) +
  # ---- CHANGED: Hollow markers for capped high outliers (size like the volcano plot) ----
geom_point(
  data = subset(data, is_high),
  aes(x = Type_Display, y = y_max),
  shape = 21,
  size = 3,        # match example
  stroke = 0.5,    # lighter outline like your volcano plot
  fill = "white",
  color = "black"
) +
  scale_fill_manual(values = fill_pal, drop = FALSE) +
  scale_color_manual(values = fill_pal, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
  stat_compare_means(
    comparisons = comparisons,
    method = "t.test",
    p.adjust.method = "BH",
    label = "p.signif",
    size = 4,
    step.increase = 0,              # not used since we hand-place labels
    bracket.nudge.y = -0.005,       # tiny nudge down to stay inside the panel
    label.y = c(y_max - 0.18,       # first comparison (Bulk IP vs Bulk Control)
                y_max - 0.10)       # second comparison (Low-input IP vs Control)
  ) +
  labs(
    title = "Signal-to-noise,\n under equivalent amplification cycles",
    subtitle = "* = t-test p",
    x = NULL,
    y = expression(log[10](Molarity/pM + 10^-6))
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(face = "italic", size = 9),
    axis.title.y       = element_text(size = 9),
    axis.text.x        = element_text(angle = 45, hjust = 1, face = "bold", size = 9),
    panel.grid         = element_blank(),   # remove all grid lines
    axis.line          = element_line(color = "black"), # add axis lines
    axis.ticks         = element_line(color = "black"), # add axis ticks
    # Optional extra headroom above, in case fonts render larger:
    plot.margin        = margin(5, 5, 15, 5)
  )

ggsave(sprintf("Tapequant_%s_S2N.pdf", date_str),
       plot = p, width = 3, height = 4.5, units = "in", dpi = 300)
