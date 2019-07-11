tb_all_CHSW <- c("bray", "jaccard", "jsd") %>% 
  purrr::map_dfr(
    function(i_dist) {
    c("IBD", "CD", "UC") %>% 
      purrr::map_dfr(
        function(i_disease) {
          readr::read_tsv(paste0("results/5-unsupervised_discrete/",
                                 i_dist, "/",
                                 "results_CHSW_", i_disease, ".tsv")) %>% 
            dplyr::mutate(disease = i_disease)
        }) %>% 
        dplyr::mutate(dist = i_dist)
  }) %>% 
  dplyr::left_join(tb_1, by = c("study" = "study_full_name"))
tb_all_predstr <- c("bray", "jaccard", "jsd") %>% 
  purrr::map_dfr(
    function(i_dist) {
      c("IBD", "CD", "UC") %>% 
        purrr::map_dfr(
          function(i_disease) {
            readr::read_tsv(paste0("results/5-unsupervised_discrete/",
                                   i_dist, "/",
                                   "results_MMUPHin_", i_disease, ".tsv")) %>% 
              dplyr::mutate(disease = i_disease)
          }) %>% 
        dplyr::mutate(dist = i_dist)
    }) %>% 
  dplyr::left_join(tb_1, by = c("study" = "study_full_name")) %>% 
  dplyr::mutate(evaluation = evaluation %>% factor(levels = c("internal", "external")))

suppFig_discreteSummary <- 
  function(tb_all_CHSW, tb_all_predstr,
           path = "supp_materials/suppFigures/suppFig_discreteSummary.pdf") {
    p_list <- c("bray", "jaccard", "jsd") %>% 
      purrr::map2(c("Bray-Curtis", "Jaccard", "root Jensen-Shannon"), function(i_dist, i_dist_full) {
        p_dlist <- c("IBD", "CD", "UC") %>% 
          purrr::map(function(i_disease) {
            p_predStr <- tb_all_predstr %>% 
              dplyr::filter(disease == i_disease, dist == i_dist) %>% 
              ggplot(aes(x = k, y = mean, color = evaluation)) +
              geom_point(size = 3, position = position_dodge(width = 0.75)) +
              geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                            width = 0.5,
                            position = position_dodge(width = 0.75)) +
              facet_wrap(.~Study, scales = "free_y", nrow = 1) +
              scale_color_manual(values = c("internal" = "black", "external" = "red"),
                                 breaks = "external",
                                 name = NULL,
                                 labels = "External") +
              scale_x_continuous(breaks = 2:8) +
              ylab("Prediction strength")
            if(i_dist == "bray" & i_disease == "IBD")
              p_predStr <- p_predStr + 
                theme(legend.position = c(1, 1),
                      legend.justification = c(1, 1),
                      legend.background = element_blank()) +
                xlab("# Clusters")
            else
              p_predStr <- p_predStr + 
                theme(legend.position = "none") +
                xlab("")
            
            p_CH <- tb_all_CHSW %>% 
              dplyr::filter(disease == i_disease, dist == i_dist, measure == "ch") %>% 
              ggplot(aes(x = k, y = value)) +
              geom_point(size = 3) +
              geom_line() +
              facet_wrap(.~Study, scales = "free_y", nrow = 1) +
              scale_x_continuous(breaks = 2:8) +
              xlab("") + ylab("Calinski-Harabasz index")
            
            p_ASW <- tb_all_CHSW %>% 
              dplyr::filter(disease == i_disease, dist == i_dist, measure == "asw") %>% 
              ggplot(aes(x = k, y = value)) +
              geom_point(size = 3) +
              geom_line() +
              facet_wrap(.~Study, scales = "free_y", nrow = 1) +
              scale_x_continuous(breaks = 2:8) +
              xlab("") + ylab("Average silhouette width")
            
            cowplot::plot_grid(p_predStr, p_CH, p_ASW, ncol = 1, rel_heights = c(1, 1, 1)) %>% 
              smar::cowplot_title(title = i_disease, rel_heights = c(0.1, 1))
          })
        cowplot::plot_grid(plotlist = p_dlist,
                           nrow = 1, rel_widths = c(1, 1, 1))
      })
    cowplot::plot_grid(plotlist = p_list, labels = c("a Bray-Curtis", "b Jaccard", "c root Jensen-Shannon"),
                       ncol = 1, rel_widths = c(1, 1, 1), align = "h",
                       label_size = 30, label_fontface = "plain",
                       hjust = 0
                       ) %>% 
      ggsave(path, ., width = 45, height = 40, limitsize = FALSE)
    return(p_list)
  }
plist <- suppFig_discreteSummary(tb_all_CHSW = tb_all_CHSW, tb_all_predstr = tb_all_predstr)
