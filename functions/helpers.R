# utilities
u.s <- function() {
  R.utils::sourceDirectory("functions/")
}

# plot
cowplot_title <- function(p, title, rel_heights = c(0.1, 1)) {
  title <- cowplot::ggdraw() +
    cowplot::draw_label(title, fontface = "bold")
  plot_grid(title, p, ncol =1 , rel_heights = rel_heights)
}
rotate_xaxis <- function(angle) {
  theme(axis.text.x = element_text(angle = angle, hjust = 1))
}
no_label_xaxis <- function() {
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
}
