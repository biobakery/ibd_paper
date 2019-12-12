load("data/ordinate/bray/ord_genera_adj.RData")
load("data/physeq/genera_adj.RData")
df_plot <- smar::sample_data2(physeq_genera_adj) %>% 
  cbind(ord_genera_adj$vectors) %>% 
  dplyr::left_join(tb_1, by = c("dataset_name" = "study_full_name"))
set.seed(1)
df_plot <- df_plot[sample.int(nrow(df_plot), replace = FALSE), ]
x.axis <- paste0("Axis 1 (", round(ord_genera_adj$values$Relative_eig[1] * 100,
                                   digits = 2), "%)")
y.axis <- paste0("Axis 2 (", round(ord_genera_adj$values$Relative_eig[2] * 100,
                                   digits = 2), "%)")

# Studies
colors_study <- c(jcolors::jcolors("default")[c(1:4)], jcolors::jcolors("pal5"))
names(colors_study) <- c(c("CS-PRISM","Pouchitis","PROTECT","RISK"),
                         setdiff(tb_1$Study, c("CS-PRISM","Pouchitis","PROTECT","RISK", NA_character_)))
df_plot_study <- df_plot %>% 
  dplyr::group_by(Study) %>% 
  dplyr::mutate(n_subject = dplyr::n_distinct(subject_accession)) %>% 
  dplyr::ungroup()
p_studies1 <- df_plot_study %>% 
  dplyr::filter(n_subject >= 200) %>% 
  ggplot(aes(x = Axis.1, y = Axis.2, color = Study)) +
  geom_point() +
  scale_color_manual(values = colors_study,
                     name = NULL) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        legend.direction = "horizontal",
        title = element_text(size = 16)) +
  guides(color=guide_legend(nrow=2,byrow=TRUE)) +
  coord_fixed() +
  ggtitle("Study (>= 200 subjects)")

p_studies2 <- df_plot_study %>% 
  dplyr::filter(n_subject < 200) %>% 
  ggplot(aes(x = Axis.1, y = Axis.2, color = Study)) +
  geom_point() +
  scale_color_manual(values = colors_study,
                     name = NULL) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        title = element_text(size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 11),
        legend.background = element_blank(),
        legend.direction = "horizontal") +
  guides(color=guide_legend(nrow=2,byrow=TRUE)) +
  coord_fixed() +
  ggtitle("Study (< 200 subjects)")
ggsave("figures/figure1/assets/ordination_study.pdf",
       cowplot::plot_grid(p_studies1, p_studies2, ncol = 1),
       width = 4,
       height = 8)

colors_disease <- c("CD" = "red", "UC" = "green4", "Control" = "black")
p_disease <- df_plot %>% 
  dplyr::mutate(Disease = ifelse(disease == "control",
                                 "Control", 
                                 disease) %>% 
                  factor(levels = c("CD", "UC", "Control"))) %>% 
  ggplot(aes(x = Axis.1, y = Axis.2, color = Disease)) +
  geom_point() +
  scale_color_manual(values = colors_disease, name = NULL) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        title = element_text(size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.background = element_blank(),
        legend.direction = "horizontal") +
  xlab(x.axis) + ylab(y.axis) +
  coord_fixed() +
  ggtitle("Disease")

# bodysite
colors_bodysite <- smar::gg_color_hue(c("CD", "Biopsy", "UC", "Stool"))[c(2, 4)]
p_bodysite <- df_plot %>% 
  dplyr::mutate(`Sample type` = sample_type %>% 
                  dplyr::recode("biopsy" = "Biopsy",
                                "stool" = "Stool") %>% 
                  factor(levels = c("Biopsy", "Stool"))) %>% 
  ggplot(aes(x = Axis.1, y = Axis.2, color = `Sample type`)) +
  geom_point() +
  scale_color_manual(values = colors_bodysite, name = NULL) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 16),
        title = element_text(size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 16),
        legend.background = element_blank(),
        legend.direction = "horizontal") +
  xlab("") + ylab("") + 
  coord_fixed() +
  ggtitle("Sample type")

# abundance
Firmicutes <- 
  physeq_genera_adj %>% 
  to_relativeAbundance() %>% 
  smar::phyloseq_to_tb() %>% 
  dplyr::filter(Rank2 == "p__Firmicutes") %>% 
  dplyr::group_by(sample_accession_16S) %>% 
  dplyr::summarise(Firmicutes = sum(abundance))
Bacteroidetes <- 
  physeq_genera_adj %>% 
  to_relativeAbundance() %>% 
  smar::phyloseq_to_tb() %>% 
  dplyr::filter(Rank2 == "p__Bacteroidetes") %>% 
  dplyr::group_by(sample_accession_16S) %>% 
  dplyr::summarise(Bacteroidetes = sum(abundance))
p_Bacteroidetes <- df_plot %>% 
  dplyr::left_join(Bacteroidetes, by = "sample_accession_16S") %>% 
  ggplot(aes(x = Axis.1, y = Axis.2, color = Bacteroidetes)) +
  geom_point() +
  scale_color_gradientn(colours = c("black","green3"), 
                        values = scales::rescale(c(0, 1)),
                        limits = c(0, 1),
                        name = "Relative abundance",
                        breaks = c(0, 0.5, 1)) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        title = element_text(size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.background = element_blank(),
        legend.direction = "horizontal") +
  coord_fixed() +
  ggtitle("Bacteroidetes")

p_Firmicutes <- df_plot %>% 
  dplyr::left_join(Firmicutes, by = "sample_accession_16S") %>% 
  ggplot(aes(x = Axis.1, y = Axis.2, color = Firmicutes)) +
  geom_point() +
  scale_color_gradientn(colours = c("black","green3"), 
                        values = scales::rescale(c(0, 1)),
                        limits = c(0, 1),
                        name = "Relative abundance",
                        breaks = c(0, 0.5, 1)) +
  theme_classic() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        title = element_text(size = 16),
        legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.background = element_blank(),
        legend.direction = "horizontal") +
  coord_fixed() +
  ggtitle("Firmicutes")

cowplot::plot_grid(p_disease, p_bodysite, 
                   p_studies1,
                   p_studies2,
                   p_Firmicutes,
                   p_Bacteroidetes,
                   ncol = 2) %>% 
  ggsave(file = "figures/figure1/assets/ordinations.pdf",
         width = 9.2,
         height = 15)



# # horizontal version? -----------------------------------------------------
# 
# p_studies1 <- df_plot_study %>% 
#   dplyr::filter(n_subject >= 200) %>% 
#   ggplot(aes(x = Axis.1, y = Axis.2, color = Study)) +
#   geom_point() +
#   scale_color_manual(values = colors_study,
#                      name = NULL) +
#   theme_classic() +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         legend.position = "bottom",
#         legend.text = element_text(size = 10),
#         legend.background = element_blank(),
#         legend.direction = "horizontal",
#         title = element_text(size = 16)) +
#   guides(color=guide_legend(nrow=2,byrow=TRUE)) +
#   coord_fixed() +
#   ggtitle("Study (>= 200 subjects)")
# 
# p_studies2 <- df_plot_study %>% 
#   dplyr::filter(n_subject < 200) %>% 
#   ggplot(aes(x = Axis.1, y = Axis.2, color = Study)) +
#   geom_point() +
#   scale_color_manual(values = colors_study,
#                      name = NULL) +
#   theme_classic() +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         title = element_text(size = 16),
#         legend.position = "bottom",
#         legend.text = element_text(size = 10),
#         legend.background = element_blank(),
#         legend.direction = "horizontal") +
#   guides(color=guide_legend(nrow=2,byrow=TRUE)) +
#   coord_fixed() +
#   ggtitle("Study (< 200 subjects)")
# ggsave("figures/figure1/assets/ordination_study.pdf",
#        cowplot::plot_grid(p_studies1, p_studies2, ncol = 1),
#        width = 4,
#        height = 8)
# # p_studies_with_legend <- df_plot_study %>% 
# #   dplyr::mutate(Study = factor(Study, levels = names(colors_study))) %>% 
# #   ggplot(aes(x = Axis.1, y = Axis.2, color = Study)) +
# #   geom_point() +
# #   theme_classic() +
# #   theme(axis.text = element_blank(),
# #         axis.ticks = element_blank()) +
# #   scale_color_manual(values = colors_study) +
# #   xlab(x.axis) + ylab(y.axis)
# # ggsave("figures/figure1/assets/ordination_study_legend.pdf",
# #        p_studies_with_legend,
# #        width = 6,
# #        height = 4)
# # disease
# colors_disease <- gg_color_hue(c("CD", "Biopsy", "UC", "Stool"))[c(1, 3)] %>% 
#   c("Control" = "black")
# p_disease <- df_plot %>% 
#   dplyr::mutate(Disease = ifelse(disease == "control",
#                                  "Control", 
#                                  disease) %>% 
#                   factor(levels = c("CD", "UC", "Control"))) %>% 
#   ggplot(aes(x = Axis.1, y = Axis.2, color = Disease)) +
#   geom_point() +
#   scale_color_manual(values = colors_disease, name = NULL) +
#   theme_classic() +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_text(size = 14),
#         title = element_text(size = 16),
#         legend.position = "bottom",
#         legend.text = element_text(size = 12),
#         legend.background = element_blank(),
#         legend.direction = "horizontal") +
#   xlab(x.axis) + ylab(y.axis) +
#   coord_fixed() +
#   ggtitle("Disease")
# # p_disease_with_legend <- df_plot %>% 
# #   dplyr::mutate(Disease = ifelse(disease == "control",
# #                                  "Control", 
# #                                  disease) %>% 
# #                   factor(levels = c("CD", "UC", "Control"))) %>% 
# #   ggplot(aes(x = Axis.1, y = Axis.2, color = Disease)) +
# #   geom_point() +
# #   theme_classic() +
# #   theme(axis.text = element_blank(),
# #         axis.ticks = element_blank()) +
# #   scale_color_manual(values = colors_disease) +
# #   xlab(x.axis) + ylab(y.axis)
# # ggsave("figures/figure1/assets/ordination_disease.pdf", 
# #        p_disease, width = 4, height = 4)
# # ggsave("figures/figure1/assets/ordination_disease_legend.pdf", 
# #        p_disease_with_legend, width = 4, height = 4)
# 
# # bodysite
# colors_bodysite <- gg_color_hue(c("CD", "Biopsy", "UC", "Stool"))[c(2, 4)]
# p_bodysite <- df_plot %>% 
#   dplyr::mutate(`Sample type` = sample_type %>% 
#                   dplyr::recode("biopsy" = "Biopsy",
#                                 "stool" = "Stool") %>% 
#                   factor(levels = c("Biopsy", "Stool"))) %>% 
#   ggplot(aes(x = Axis.1, y = Axis.2, color = `Sample type`)) +
#   geom_point() +
#   scale_color_manual(values = colors_bodysite, name = NULL) +
#   theme_classic() +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_text(size = 14),
#         title = element_text(size = 16),
#         legend.position = "bottom",
#         legend.text = element_text(size = 12),
#         legend.background = element_blank(),
#         legend.direction = "horizontal") +
#   xlab("") + ylab("") + 
#   coord_fixed() +
#   ggtitle("Sample type")
# 
# # p_bodysite <- df_plot %>% 
# #   dplyr::mutate(`Sample type` = sample_type %>% 
# #                   dplyr::recode("biopsy" = "Biopsy",
# #                                 "stool" = "Stool") %>% 
# #                   factor(levels = c("Biopsy", "Stool"))) %>% 
# #   ggplot(aes(x = Axis.1, y = Axis.2, color = `Sample type`)) +
# #   geom_point() +
# #   theme_classic() +
# #   theme(axis.text = element_blank(),
# #         axis.ticks = element_blank(),
# #         axis.title = element_blank()) +
# #   scale_color_manual(values = colors_bodysite, name = "Sample type", guide = FALSE)
# # ggsave("figures/figure1/assets/ordination_sampleType.pdf", 
# #        p_bodysite, width = 4, height = 4)
# # ggsave("figures/figure1/assets/ordination_sampleType_legend.pdf", 
# #        p_bodysite_with_legend, width = 4, height = 4)
# 
# # abundance
# Firmicutes <- 
#   physeq_genera_adj %>% 
#   to_relativeAbundance() %>% 
#   phyloseq_to_tb() %>% 
#   dplyr::filter(Rank2 == "p__Firmicutes") %>% 
#   dplyr::group_by(sample_accession_16S) %>% 
#   dplyr::summarise(Firmicutes = sum(abundance))
# Bacteroidetes <- 
#   physeq_genera_adj %>% 
#   to_relativeAbundance() %>% 
#   phyloseq_to_tb() %>% 
#   dplyr::filter(Rank2 == "p__Bacteroidetes") %>% 
#   dplyr::group_by(sample_accession_16S) %>% 
#   dplyr::summarise(Bacteroidetes = sum(abundance))
# p_Bacteroidetes <- df_plot %>% 
#   dplyr::left_join(Bacteroidetes, by = "sample_accession_16S") %>% 
#   ggplot(aes(x = Axis.1, y = Axis.2, color = Bacteroidetes)) +
#   geom_point() +
#   scale_color_gradientn(colours = c("black","green3"), 
#                         values = scales::rescale(c(0, 1)),
#                         limits = c(0, 1),
#                         name = "Relative abundance",
#                         breaks = c(0, 0.5, 1)) +
#   theme_classic() +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         title = element_text(size = 16),
#         legend.position = "bottom",
#         legend.text = element_text(size = 12),
#         legend.background = element_blank(),
#         legend.direction = "horizontal") +
#   coord_fixed() +
#   ggtitle("Bacteroidetes")
# 
# p_Firmicutes <- df_plot %>% 
#   dplyr::left_join(Firmicutes, by = "sample_accession_16S") %>% 
#   ggplot(aes(x = Axis.1, y = Axis.2, color = Firmicutes)) +
#   geom_point() +
#   scale_color_gradientn(colours = c("black","green3"), 
#                         values = scales::rescale(c(0, 1)),
#                         limits = c(0, 1),
#                         name = "Relative abundance",
#                         breaks = c(0, 0.5, 1)) +
#   theme_classic() +
#   theme(axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         title = element_text(size = 16),
#         legend.position = "bottom",
#         legend.text = element_text(size = 12),
#         legend.background = element_blank(),
#         legend.direction = "horizontal") +
#   coord_fixed() +
#   ggtitle("Firmicutes")
# # p_Bacteroidetes_legend <- df_plot %>% 
# #   dplyr::left_join(Bacteroidetes, by = "sample_accession_16S") %>% 
# #   ggplot(aes(x = Axis.1, y = Axis.2, color = Bacteroidetes)) +
# #   geom_point() +
# #   theme_classic() +
# #   theme(axis.text = element_blank(),
# #         axis.ticks = element_blank(),
# #         axis.title = element_blank(), legend.direction = "horizontal") +
# #   jcolors::scale_colour_jcolors_contin("pal2")
# # ggsave("figures/figure1/assets/ordination_Bacteroidetes.pdf", 
# #        p_Bacteroidetes, width = 4, height = 4)
# # ggsave("figures/figure1/assets/ordination_Firmicutes.pdf", 
# #        p_Firmicutes, width = 4, height = 4)
# # ggsave("figures/figure1/assets/ordination_Bacteroidetes_legend.pdf", 
# #        p_Bacteroidetes_legend, width = 4, height = 4)
# 
# cowplot::plot_grid(p_disease, p_bodysite, 
#                    p_studies1,
#                    p_studies2,
#                    p_Firmicutes,
#                    p_Bacteroidetes,
#                    ncol = 2) %>% 
#   ggsave(file = "figures/figure1/assets/ordinations.pdf",
#          width = 8.5,
#          height = 15)
# 
