library(ggplot2)

theme_nature <- {
  theme_minimal(base_size = 10) +
    theme(
      panel.grid   = element_blank(),
      axis.line    = element_line(color = "black"),
      axis.ticks   = element_line(color = "black"),
      plot.title   = element_text(face = "bold", size = 10),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8),
      legend.title = element_text(face = "bold", size = 9),
      legend.text  = element_text(size = 8),
      plot.margin = margin(5, 5, 5, 5)
    )
}

# theme_minimal(base_size = 10) default text sizes:
# - plot.title       : 12 pt  (rel(1.2) * base_size)
# - axis.title       : 10 pt  (rel(1.0) * base_size)
# - axis.text        :  8 pt  (rel(0.8) * base_size)
# - legend.title     : 10 pt  (rel(1.0) * base_size)
# - legend.text      :  8 pt  (rel(0.8) * base_size)
