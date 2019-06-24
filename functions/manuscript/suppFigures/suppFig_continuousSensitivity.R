# sensitivity plots
l_fitContinuous <- list()
for(i_cutoff in c(0.5, 0.65, 0.8)) {
  load(paste0("results/6-unsupervised_continuous/genera/cor_cutoff_",
              i_cutoff, "/fit_continuous.RData"))
  l_fitContinuous <- c(l_fitContinuous, list(fit_continuous))
}
path <- "supp_materials/suppFigures/suppFig_continuousSensitivity/"
l_fitContinuous %>% purrr::map_dbl(~ .x$membership[[2]] %>% length())
central_nodes <- l_fitContinuous[[2]]$clustered.network$communities[1:2]
colors <- smar::gg_color_hue(c("CS-PRISM_stool_CD",
                               "CS-PRISM_stool_UC",
                               "Pouchitis_biopsy_UC",
                               "PROTECT_biopsy_UC",
                               "PROTECT_stool_UC",
                               "RISK_biopsy_CD",
                               "RISK_stool_CD"))
set.seed(1)
for(i in 1:length(l_fitContinuous)) {
  i.graph <- l_fitContinuous[[i]]$clustered.network$pc.graph
  i.groups.all <- l_fitContinuous[[i]]$clustered.network$communities[
    1:max(l_fitContinuous[[i]]$clustered.network$communities$membership)]
  i.groups.all <- i.groups.all[i.groups.all %>% purrr::map_lgl(~ length(.x) > 2)]
  i.groups.sub <- i.groups.all[1:2]
  i.graph.sub <- igraph::delete.vertices(i.graph, 
                                         igraph::V(i.graph)[setdiff(names(igraph::V(i.graph)),
                                                                    unlist(i.groups.all))])
      
  v.size <- ifelse(igraph::V(i.graph.sub)$name %in% unlist(i.groups.sub),
                   25 / max(igraph::degree(i.graph.sub)) * igraph::degree(i.graph.sub),
                   5)
  v.color <- igraph::V(i.graph.sub)$name %>% 
    stringr::str_replace("\\, PC.*", "") %>% 
    `[`(colors, .)
    
  v.label <- igraph::V(i.graph.sub)$name %>% 
    stringr::str_replace(stringr::fixed("_"), "\n") %>% 
    stringr::str_replace(stringr::fixed("stool_"), "Stool ") %>% 
    stringr::str_replace(stringr::fixed("biopsy_"), "Biopsy ") %>% 
    stringr::str_replace(stringr::fixed(","), "\n") %>% 
    ifelse(igraph::V(i.graph.sub)$name %in% unlist(i.groups.sub),
           .,
           NA)
  v.label.color <- ifelse(igraph::V(i.graph.sub)$name %in% central_nodes[[2]],
                          "red",
                          "black")
  font.size <- ifelse(i == 1, 0.8, 1.5)
  
  if(i == 1) {
    mat_coord <- i.groups.all %>%
      purrr::map2_dfr(list(c(-0.5, 0.25, 0.2), c(0.5, 0.25, 0.25),
                           c(-0.75, -0.25, 0.1), c(-0.25, -0.25, 0.1),
                           c(0.167, -0.2, 0.05), c(0.5, -0.2, 0.05), c(0.833, -0.2, 0.05),
                           c(0.167, -0.4, 0.05), c(0.5, -0.4, 0.05), c(0.833, -0.4, 0.05)),
                      function(vertices, params) {
                        n_total <- length(vertices)
                        interval <- 2/n_total
                        radius <- params[3]
                        # x_center <- sign(x_center) * (1- radius)
                        (1:n_total) %>%
                          purrr::map_dfr(function(i_vertex) {
                            tibble::tibble(
                              vertex = vertices[i_vertex],
                              x = params[1] + radius*sinpi((i_vertex - 1)*interval),
                              y = params[2] + radius*cospi((i_vertex - 1)*interval)
                            )
                          })
                      }) %>%
      tibble::as_data_frame() %>%
      tibble::column_to_rownames("vertex") %>%
      as.matrix()
  }
  if(i == 2) {
    mat_coord <- i.groups.all %>%
      purrr::map2_dfr(list(c(-0.5, 0.2, 0.3), c(0.5, 0.2, 0.3),
                           c(-0.5, -0.4, 0.1), c(0.5, -0.4, 0.1)),
                      function(vertices, params) {
                        n_total <- length(vertices)
                        interval <- 2/n_total
                        radius <- params[3]
                        # x_center <- sign(x_center) * (1- radius)
                        (1:n_total) %>%
                          purrr::map_dfr(function(i_vertex) {
                            tibble::tibble(
                              vertex = vertices[i_vertex],
                              x = params[1] + radius*sinpi((i_vertex - 1)*interval),
                              y = params[2] + radius*cospi((i_vertex - 1)*interval)
                            )
                          })
                      }) %>%
      tibble::as_data_frame() %>%
      tibble::column_to_rownames("vertex") %>%
      as.matrix()
  }
  if(i == 3) {
    mat_coord <- i.groups.all %>%
      purrr::map2_dfr(list(c(-0.5, 0, 0.4), c(0.5, 0, 0.3)),
                      function(vertices, params) {
                        n_total <- length(vertices)
                        interval <- 2/n_total
                        radius <- params[3]
                        # x_center <- sign(x_center) * (1- radius)
                        (1:n_total) %>%
                          purrr::map_dfr(function(i_vertex) {
                            tibble::tibble(
                              vertex = vertices[i_vertex],
                              x = params[1] + radius*sinpi((i_vertex - 1)*interval),
                              y = params[2] + radius*cospi((i_vertex - 1)*interval)
                            )
                          })
                      }) %>%
      tibble::as_data_frame() %>%
      tibble::column_to_rownames("vertex") %>%
      as.matrix()
  }
  
  
  pdf(paste0(path, "cor_", i, "_network.pdf"), width = 15, height = 9)
  igraph::plot.igraph(i.graph.sub, 
                      layout = mat_coord[igraph::V(i.graph.sub)$name, ],
                      xlim = c(-1, 1),
                      ylim = c(-0.6, 0.6),
                      # mark.groups = i.groups.sub,
                      # mark.shape = 1,
                      # mark.col =
                      vertex.size = v.size,
                      vertex.frame.color = NA,
                      vertex.color = v.color,
                      vertex.label = v.label,
                      vertex.label.family = "sans",
                      vertex.label.cex = font.size,
                      vertex.label.color = v.label.color,
                      # edge.curved = 0.2,
                      # edge.label = NA,
                      edge.width = 3,
                      rescale = FALSE,
                      margin = 0)
  dev.off()
}

# legend for labels
graph.legend <- igraph::graph_from_literal("CS-PRISM Stool CD",
                                           "CS-PRISM Stool UC",
                                           "Pouchitis Biopsy UC",
                                           "PROTECT Biopsy UC",
                                           "PROTECT Stool UC",
                                           "RISK Biopsy CD",
                                           "RISK Stool CD")
tb_legend <- tibble::tibble(x = 0,
                            y = seq(1, -1, length.out = 7),
                            label = c("CS-PRISM Stool CD",
                                      "CS-PRISM Stool UC",
                                      "Pouchitis Biopsy UC",
                                      "PROTECT Biopsy UC",
                                      "PROTECT Stool UC",
                                      "RISK Biopsy CD",
                                      "RISK Stool CD"))
colors2 <- colors
names(colors2) <- tb_legend$label
p_legend <- ggplot(tb_legend, aes(x = x, y = y, label = label)) +
  geom_point(aes(color = label), size = 10) +
  geom_text(aes(x = 0.003), hjust = 0) +
  scale_color_manual(values = colors2) +
  scale_x_continuous(limits = c(-0.05, 0.05)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank())
ggsave(paste0(path, "legends.pdf"), p_legend, width = 7, height = 4)