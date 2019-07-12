suppFig_continuousModel <- function(path = "supp_materials/suppFigures/suppFig_continuousModel/network.pdf") {
  g_toy <- igraph::graph_from_literal("Study 1, PC1"-"Study 2, PC2"-"Study 3, PC1"-"Study 1, PC1",
                                      "Study 1, PC2"-"Study 2, PC1"-"Study 3, PC2"-"Study 1, PC2")
  mat_coord <- list(c(-0.5, 0), c(0.5, 0))  %>%
    purrr::map_dfr(function(center) {
                      n_total <- 3
                      interval <- 2/n_total
                      radius <- 0.4
                      (1:n_total) %>%
                        purrr::map_dfr(function(i_vertex) {
                          tibble::tibble(
                            x = center[1] + radius*sinpi((i_vertex - 1)*interval),
                            y = center[2] - sign(center[1]) * radius*cospi((i_vertex - 1)*interval)
                          )
                        })
                    }) %>%
    as.data.frame()
  
  pdf(path, width = 10, height = 5)
  igraph::plot.igraph(g_toy, 
                      layout = mat_coord,
                      xlim = c(-1, 1),
                      ylim = c(-0.6, 0.6),
                      vertex.size = 30,
                      vertex.frame.color = NA,
                      vertex.color = "grey",
                      vertex.label.family = "sans",
                      vertex.label.cex = 1,
                      vertex.label.color = "black",
                      # edge.curved = 0.2,
                      # edge.label = NA,
                      edge.width = 1,
                      rescale = FALSE,
                      margin = 0)
  dev.off()
  
}