tb_predstr <- readr::read_tsv("results/5-unsupervised_discrete/bray/results_MMUPHin_IBD.tsv") %>% 
  dplyr::left_join(tb_1, by = c("study" = "study_full_name")) %>% 
  dplyr::mutate(evaluation = evaluation %>% factor(levels = c("internal", "external")))

fig_discrete <- 
  function(tb_predstr) {
    p <- tb_predstr %>% 
      ggplot(aes(x = k, y = mean, color = evaluation)) +
      geom_point(size = 3, position = position_dodge(width = 0.75)) +
      geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                    width = 0.5,
                    position = position_dodge(width = 0.75)) +
      facet_grid(.~Study, scales = "free_y") +
      scale_color_manual(values = c("internal" = "black", "external" = "red"),
                         breaks = "external",
                         name = NULL,
                         labels = "External") +
      scale_x_continuous(breaks = 2:8) +
      xlab("# Clusters") + ylab("Prediction strength") +
      theme(legend.position = c(1, 1),
            legend.justification = c(1, 1),
            legend.background = element_blank())
    
    return(p)
  }
p_discrete <- fig_discrete(tb_predstr)

tb_loading <- readr::read_tsv("results/6-unsupervised_continuous/genera/cor_cutoff_0.65/avg_loading.tsv")
load("data/physeq/genera_adj.RData")
tb_loading <- smar::tax_table2(physeq_genera_adj) %>% 
  as.data.frame(check.names = FALSE, stringAsFactors = FALSE) %>% 
  tibble::rownames_to_column("feature") %>% 
  dplyr::right_join(tb_loading, by = "feature")

fig_loading <- 
  function(tb_loading) {
    p <- tb_loading %>% 
      plyr::mutate(feature_plot = betterGeneraNames(Rank4, Rank5, Rank6),) %>% 
      dplyr::arrange(desc(abs(Cluster_1))) %>% 
      dplyr::slice(1:20) %>% 
      dplyr::arrange(Cluster_1) %>% 
      dplyr::mutate(feature_plot = factor(feature_plot, levels = unique(feature_plot))) %>% 
      ggplot(aes(x = feature_plot, y = Cluster_1)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      ylab("Loading (dysbiosis)") +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_text(face = "italic"))
  }
p_loading <- fig_loading(tb_loading)

tb_scores <- readr::read_tsv("results/6-unsupervised_continuous/genera/scores.tsv")
load("data/ordinate/bray/ord_genera_adj.RData")
metadata <- smar::sample_data2(physeq_genera_adj) %>% 
  dplyr::left_join(tb_scores, by = c("sample_accession_16S" = "sample")) %>% 
  cbind(ord_genera_adj$vectors) %>% 
  dplyr::mutate(Disease = dplyr::case_when(
    disease == "CD" ~ "CD",
    disease == "UC" ~ "UC",
    control == "nonIBD" ~ "Control (non-IBD)",
    control == "HC" ~ "Control (healthy)"
  ) %>% 
    factor(levels = c("CD", "UC", "Control (non-IBD)", "Control (healthy)")))
fig_score <- function(metadata) {
  p <- metadata %>% 
    ggplot(aes(x = Axis.1, y = Axis.2, color = `Score (dysbiosis)`)) +
    geom_point() +
    facet_grid(.~Disease) +
    scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                          breaks = c(-1, -0.5, 0, 0.5), labels = c("-1", "-0.5", "0", "0.5"), 
                          name = "Score\n(dysbiosis)") +
    coord_fixed() +
    xlab(paste0("Axis 1 (", round(ord_genera_adj$values$Relative_eig[1] * 100, digits = 2), "%)")) +
    ylab(paste0("Axis 2 (", round(ord_genera_adj$values$Relative_eig[2] * 100, digits = 2), "%)")) +
    theme(axis.ticks = element_blank(), axis.text = element_blank(),
          legend.position = c(0.99, 0),
          legend.justification = c(1, 0),
          legend.direction = "horizontal",
          legend.background = element_blank())
}
p_scores <- fig_score(metadata)

path <- "figures/figure5/figure5.pdf"
cowplot::plot_grid(p_discrete,
                   cowplot::plot_grid(NULL, p_loading, nrow = 1, rel_widths = c(1, 1),
                                      labels = c("b", "c"),
                                      label_fontface = "plain", label_size = 20),
                   p_scores, 
                   ncol = 1, rel_heights = c(3.2, 4.2, 4), labels = c("a", "", "d"),
                   label_fontface = "plain", label_size = 20) %>% 
  # {cowplot::ggdraw() +
  #     cowplot::draw_image("figures/figure5/assets/network.pdf", x = 0.01, y = 4/(3.2 + 4.2 + 4),
  #                         width = 0.58, height = 4.2/(3.2 + 4.2 + 4)) +
  #     cowplot::draw_plot(.)} %>% 
  ggsave(path, ., width = 15, height = 12.5)

