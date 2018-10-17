cowplot_title <- function(p, title, rel_heights = c(0.1, 1)) {
  title <- cowplot::ggdraw() +
    cowplot::draw_label(title, fontface = "bold")
  plot_grid(title, p, ncol =1 , rel_heights = rel_heights)
}
