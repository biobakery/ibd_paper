load("data/physeq/genera_adj.RData")
tb_loading <- readr::read_tsv("results/6-unsupervised_continuous/genera/cor_cutoff_0.65/avg_loading.tsv")
tb_scores <- readr::read_tsv("results/6-unsupervised_continuous/genera/scores.tsv")
path <- "supp_materials/suppFigures/suppFig_continuous2.pdf"

tb_loading <- smar::tax_table2(physeq_genera_adj) %>% 
  as.data.frame(check.names = FALSE, stringAsFactors = FALSE) %>% 
  tibble::rownames_to_column("feature") %>% 
  dplyr::right_join(tb_loading, by = "feature")

colors_mapping <- smar::gg_color_hue(c("Firmicutes", "Actinobacteria",
                                       "Bacteroidetes", "Proteobacteria", "Fusobacteria"))
p_loading <- tb_loading %>% 
  plyr::mutate(feature_plot = betterGeneraNames(Rank4, Rank5, Rank6),) %>% 
  dplyr::arrange(desc(abs(Cluster_2))) %>% 
  dplyr::slice(1:20) %>% 
  dplyr::arrange(Cluster_2) %>% 
  dplyr::mutate(feature_plot = factor(feature_plot, levels = unique(feature_plot))) %>% 
  dplyr::mutate(Phylum = Rank2 %>% 
                  stringr::str_replace(stringr::fixed("p__"), "") %>% 
                  factor(levels = names(colors_mapping))) %>% 
  ggplot(aes(x = feature_plot, y = Cluster_2)) +
  geom_bar(stat = "identity", aes(color = Phylum), fill = "white") +
  scale_color_manual(values = colors_mapping) +
  coord_flip() +
  ylab("Loading (phyla tradeoff)") +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(face = "italic"),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank())

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
p_scores <- metadata %>% 
  dplyr::mutate(Disease = dplyr::case_when(
    disease == "CD" ~ "CD",
    disease == "UC" ~ "UC",
    control == "nonIBD" ~ "Control (non-IBD)",
    TRUE ~ "Control (healthy)"
  ) %>% factor(levels = c("CD", "UC","Control (non-IBD)", "Control (healthy)"))) %>% 
  ggplot(aes(x = Axis.1, y = Axis.2, color = `Score (phyla tradeoff)`)) +
  geom_point() +
  facet_grid(.~Disease) +
  scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, 
                        breaks = c(-0.5, 0, 0.5), name = "Score\n(phyla tradeoff)") +
  coord_fixed() +
  xlab(paste0("Axis 1 (", round(ord_genera_adj$values$Relative_eig[1] * 100, digits = 2), "%)")) +
  ylab(paste0("Axis 2 (", round(ord_genera_adj$values$Relative_eig[2] * 100, digits = 2), "%)")) +
  theme(legend.position = c(1, 0),
        legend.justification = c(1, 0),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

p <- cowplot::plot_grid(p_loading, p_scores, nrow = 1, rel_widths = c(1, 2.2),
                        align = "hv", axis = "tb",
                        labels = c("a", "b"),
                        label_fontface = "plain",
                        label_size = 30)
ggsave(path, p, 
       width = 22, height = 5)
