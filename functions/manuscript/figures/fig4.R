source("assets/ma_tests.R")
l_results <- list()
tests_CD <- c("B3_vs_B1", "B2_vs_B1")
tests_UC <- c("E3_vs_E1", "E2_vs_E1")
for(test in c(tests_CD, tests_UC, "CD_vs_control", "steroids")) {
  load(paste0("results/4.2-maaslin2/genera/", test, "/result.RData"))
  l_results[[test]] <- result
}
load("data/physeq/genera_adj.RData")
tax_table <- smar::tax_table2(physeq_genera_adj)
fig_MAPhenotypeSeverity <- function(l_results,
                                    tests_CD,
                                    tests_UC,
                                    tax_table,
                                    path = "figures/figure4/assets/") {
  
  colors_subtype <- c("B3_vs_B1" = "red",
                      "B2_vs_B1" = scales::muted("brown"),
                      "E3_vs_E1" = "green4",
                      "E2_vs_E1" = scales::muted("aquamarine"))
  
  list(`CD Behavior` = tests_CD,
       `UC Extent` = tests_UC) %>% 
    purrr::imap(function(tests, name_tests) {
      result_all <- l_results[tests] %>% 
        purrr::imap_dfr(function(result, test) {
          result$tb_summarise %>% 
            tidyr::gather(key = study, value = weight, dplyr::matches("weight_")) %>% 
            dplyr::group_by_at(vars(-weight, -study)) %>% 
            dplyr::arrange(desc(weight)) %>% 
            dplyr::summarise(max_weight = weight[1],
                             study_max_weight = study[1]) %>% 
            dplyr::ungroup() %>% 
            dplyr::mutate(study_max_weight = study_max_weight %>% 
                            stringr::str_replace_all("weight\\_", ""),
                          test = test) %>% 
            dplyr::filter(k > 1) %>% 
            dplyr::mutate(qval.fdr_meta = p.adjust(pval, method = "fdr"))
        }) %>% 
        dplyr::group_by(feature) %>% 
        dplyr::mutate(consistent = sign(coef[1]) == sign(coef[2]),
                      min_q = min(qval.fdr_meta, na.rm = TRUE)) %>% 
        dplyr::ungroup() %>% 
        dplyr::filter(consistent, min_q < 0.1)
      
      load(paste0("results/4.3-montreal/genera/",
                  name_tests, "/result.RData"))
      
      tb_contrast <- result$tb_summarise %>%
        dplyr::select(feature, coef, k, pval, qval.fdr) %>%
        dplyr::mutate(test = "comparison") %>% 
        dplyr::filter(k > 1, feature %in% result_all$feature) %>% 
        dplyr::mutate(qval.fdr_meta = p.adjust(pval, method = "fdr"))

      tb_plot <- result_all %>%
        dplyr::bind_rows(tb_contrast) %>% 
        dplyr::mutate(test = factor(test, levels = c(tests, "comparison"))) %>% 
        dplyr::group_by(feature) %>% 
        dplyr::arrange(test) %>% 
        dplyr::ungroup()
      
      tb_plot <- tb_plot %>% 
        dplyr::mutate(q_levels = 
                        dplyr::case_when(
                          qval.fdr_meta < 0.01 ~ "***",
                          qval.fdr_meta < 0.05 ~ "**",
                          qval.fdr_meta < 0.1 ~ "*",
                          TRUE ~ ""
                        )) %>% 
        dplyr::left_join(tax_table %>% 
                           as.data.frame() %>% 
                           tibble::rownames_to_column("feature"),
                         by = "feature") %>% # create better feature names
        dplyr::filter(Rank5 != "f__") %>%
        dplyr::mutate(feature_plot = betterGeneraNames(Rank4, Rank5, Rank6))
      
      # Order features
      feature_plot_orderedByBeta <- tb_plot %>%
        dplyr::filter(test == tests[1]) %>% 
        dplyr::group_by(sign(coef), Rank5) %>%
        dplyr::mutate(max_coef = sign(coef)*max(abs(coef), na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(sign(coef), max_coef, coef) %>%
        dplyr::pull(feature_plot)
      
      tb_comparison <- tb_plot %>%
        dplyr::group_by(feature) %>% 
        dplyr::arrange(test) %>% 
        dplyr::filter(qval.fdr_meta[3] < 0.1) %>%
        dplyr::summarise(x = factor(feature_plot[3],
                                    levels = feature_plot_orderedByBeta) %>%
                           as.numeric,
                         y = coef[1] + sign(coef[3])*0.008,
                         q_levels = q_levels[3]) 
      
      
      p <- tb_plot %>% 
        dplyr::mutate(feature_plot = factor(feature_plot, 
                                            levels = feature_plot_orderedByBeta)) %>%
        dplyr::filter(test != "comparison") %>% 
        ggplot(aes(x = feature_plot, y = coef)) +
        geom_bar(aes(fill = test, group = test), stat = "identity",
                 position = position_dodge(width = 0.75),
                 width = 0.75) +
        geom_text(aes(y = coef,
                      label = q_levels, 
                      vjust = ifelse(coef > 0, 1, 0.5),
                      group = test,
                      angle = 90), 
                  position = position_dodge(width = 0.75),
                  inherit.aes = TRUE,
                  hjust = "center",
                  size = 7) +
        geom_segment(data = tb_comparison,
                     aes(x = x - 0.375/2, xend = x + 0.375/2,
                         y = y, yend = y)) +
        geom_text(data = tb_comparison,
                  aes(x = x,
                      y = ifelse(y > 0,
                                 y + 0.008,
                                 y - 0.002),
                      label = q_levels,
                      vjust = 0.5),
                  hjust = "center",
                  size = 8,
                  angle = 90) +
      scale_fill_manual(values = colors_subtype,
                        guide = guide_legend(reverse = TRUE),
                        labels = tests %>% 
                          stringr::str_replace(stringr::fixed("_vs_"),
                                               " vs. "),
                        name = "Comparison") +
        coord_flip() +
        ylab("Effect size") +
        theme(title = element_text(size = 12),
              axis.title.y = element_blank(),
              axis.text.y = element_text(face = "italic"),
              legend.position = c(0, 1),
              legend.justification = c(0, 1),
              legend.background = element_blank())
      
      if(name_tests == "CD Behavior") 
        p <- p + ggtitle("CD Behavior\nB1 non-stricturing, non-penetraiting\nB2 stricturing\nB3 penetrating")
      if(name_tests == "UC Extent") 
        p <- p + ggtitle("UC Extent\nE1 Ulcerative proctitis\nE2 Distal UC\nE3 Pancolitis") +
        annotate(geom = "text",
                 x = 1,
                 y = 0.02,
                 label = paste0("* q < 0.1\n",
                                "** q < 0.05\n",
                                "*** q < 0.01\n"),
                 hjust = 0,
                 vjust = 0,
                 size = 6)
      
      ggsave(paste0(path,
                    name_tests,
                    ".pdf"),
             p,
             width = 8, height = 6) 
      return(p)
    })
}
l_p_montreal <- fig_MAPhenotypeSeverity(l_results = l_results,
                                      tests_CD = tests_CD,
                                      tests_UC = tests_UC,
                                      tax_table = tax_table)

# interaction -------------------------------------------------------------
metadata <- physeq_genera_adj %>% smar::sample_data2()
metadata_test <- readxl::read_xlsx("assets/table1_metadata.xlsx") %>% 
  dplyr::right_join(metadata, by = c("study_full_name" = "dataset_name")) %>% 
  dplyr::mutate(sample_type2 = 
                  sample_type %>% 
                  dplyr::recode("biopsy" = "Biopsy",
                                "stool" = "Stool"),
                disease2 = 
                  disease %>% 
                  dplyr::recode("control" = "Control")) %>% 
  dplyr::mutate(study_site = paste(Study,
                                   sample_type2, 
                                   sep = ", "),
                dataset_site = paste(study_full_name,
                                     sample_type, 
                                     sep = "_"),
                study_site_disease = paste(Study,
                                           sample_type2, 
                                           disease2,
                                           sep = ", "),
                dataset_site_disease = paste(study_full_name,
                                             sample_type,
                                             disease,
                                             sep = "_"))
df_moderator_site <- metadata_test %>% 
  dplyr::select(study_site, dataset_site, sample_type, 
                Study, study_full_name) %>% 
  dplyr::group_by(study_site) %>% 
  dplyr::summarise_all(function(x) unique(x)) %>% 
  dplyr::mutate(sample_type = factor(sample_type)) %>% 
  as.data.frame() %>% tibble::column_to_rownames("dataset_site") 
df_moderator_siteDisease <- metadata_test %>% 
  dplyr::select(study_site_disease, dataset_site_disease, sample_type, disease,
                Study, study_full_name) %>% 
  dplyr::filter(disease %in% c("CD", "UC")) %>% 
  dplyr::group_by(study_site_disease) %>% 
  dplyr::summarise_all(function(x) unique(x)) %>% 
  dplyr::mutate(sample_type = factor(sample_type),
                disease = factor(disease)) %>% 
  as.data.frame() %>% tibble::column_to_rownames("dataset_site_disease")

test <- "CD_vs_control"
test_param <- l_tests[[test]]
result <- l_results[[test]]
tb_result_moderator <- readr::read_tsv("results/4.4-interaction/genera/CD_vs_control_moderator_results.tsv") %>% 
  dplyr::filter(moderator_level == "difference") %>% 
  dplyr::group_by(moderator) %>% 
  dplyr::mutate(q = pval %>% p.adjust("fdr")) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(stringr::str_detect(feature, "Dehalobacterium")) %>% 
  dplyr::mutate(moderator_level = moderator_level %>% 
                  dplyr::recode("difference" = "Difference")) %>% 
  dplyr::mutate(class = "MA effect",
                batch = paste0(class, " (difference)"),
                study = "MA effect",
                study_full_name = "MA effect") %>% 
  dplyr::select(feature, coef, stderr, class, batch, moderator_level, 
                study, study_full_name, q)
tb_result_ind <- result$result$ind.results %>% 
  purrr::reduce(rbind) %>% 
  dplyr::mutate(test = test) %>% 
  dplyr::filter(stringr::str_detect(feature, "Dehalobacterium")) %>% 
  dplyr::left_join(df_moderator_site %>% 
                     tibble::rownames_to_column("dataset_site"),
                   by = c("batch" = "dataset_site")) %>% 
  dplyr::mutate(batch = study_site,
                class = "individual",
                moderator_level = sample_type %>% 
                  dplyr::recode("biopsy" = "Biopsy",
                                "stool" = "Stool"),
                study = Study) %>% 
  dplyr::select(feature, coef, stderr, class, batch, moderator_level, 
                study, study_full_name)
tb_result_MA <- result$result$meta.results %>% 
  dplyr::filter(stringr::str_detect(feature, "Dehalobacterium")) %>% 
  dplyr::mutate(batch = "MA effect (combined)",
                class = "MA effect",
                moderator_level = "Combined",
                study = "MA effect",
                study_full_name = "MA effect") %>% 
  dplyr::select(feature, coef, stderr, class, batch, moderator_level, 
                study, study_full_name)
tb_plot <- rbind(tb_result_ind, 
                 tb_result_moderator %>% dplyr::select(-q), 
                 tb_result_MA) %>% 
  dplyr::mutate(study_full_name = factor(study_full_name, 
                                         levels = c(studies, "MA effect")),
                moderator_level = factor(moderator_level, 
                                         levels = c("Biopsy", "Stool", 
                                                    "Combined",
                                                    "Difference"))) %>% 
  dplyr::arrange(moderator_level, study_full_name) %>% 
  dplyr::mutate(batch_plot = factor(batch, levels = rev(batch)),
                moderator_color = moderator_level %>% 
                  dplyr::recode_factor("Biopsy" = "Biopsy",
                                       "Stool" = "Stool",
                                       "Combined" = "MA effect",
                                       "Difference" = "MA effect"))
colors_bodysite <- smar::gg_color_hue(c("CD", "Biopsy", "UC", "Stool"))[c(2, 4)]
p_CD <- tb_plot %>% 
  ggplot(aes(x = batch_plot, y = coef)) +
  geom_point(size = 5, aes(shape = class, color = moderator_color)) +
  geom_errorbar(aes(ymin = coef - 1.96*stderr, 
                    ymax = coef + 1.96*stderr, color = moderator_color), 
                width = 0.5) +
  # scale_size_manual(values = c("individual" = 1, "MA effect" = 2.5),
  #                   guide = FALSE) +
  scale_shape_manual(values = c("individual" = 16,
                                "MA effect" = 15),
                     name = "Effect",
                     labels = c("Individual cohort",
                                "MA effect")) +
  scale_color_manual(values = c(colors_bodysite, "MA effect" = "black"),
                     name = "Sample type") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme(axis.title.y = element_blank(),
        title = element_text(size = 12)) +
  ylab("Effect size") +
  ggtitle(expression(paste("CD vs. control effects on ", 
                           italic("Dehalobacterium")))) +
  geom_text(data = tb_result_moderator,
            aes(x = 1,
                y = coef - 1.96*stderr,
                label = paste0("q = ", 
                               formatC(q, format = "e", digits = 2),
                               " ")),
            size = 5,
            hjust = 1)

test <- "steroids"
test_param <- l_tests[[test]]
result <- l_results[[test]]

tb_result_moderator <- readr::read_tsv("results/4.4-interaction/genera/steroids_moderator_results.tsv") %>% 
  dplyr::filter(moderator_level == "difference") %>% 
  dplyr::group_by(moderator) %>% 
  dplyr::mutate(q = pval %>% p.adjust("fdr")) %>% 
  dplyr::ungroup() %>% 
  dplyr::filter(stringr::str_detect(feature, 
                                    stringr::fixed("f__Enterobacteriaceae|g__|"))) %>% 
  dplyr::mutate(moderator_level = moderator_level %>% 
                  dplyr::recode("difference" = "Difference")) %>% 
  dplyr::mutate(class = "MA effect",
                batch = paste0(class, " (difference)"),
                study = "MA effect",
                study_full_name = "MA effect") %>% 
  dplyr::filter(moderator == "disease") %>%
  dplyr::select(feature, coef, stderr, class, batch, moderator_level, 
                study, study_full_name, q)
tb_result_ind <- result$result$ind.results %>% 
  purrr::reduce(rbind) %>% 
  dplyr::mutate(test = test) %>% 
  dplyr::filter(stringr::str_detect(feature, 
                                    stringr::fixed("f__Enterobacteriaceae|g__|"))) %>% 
  dplyr::left_join(df_moderator_siteDisease %>% 
                     tibble::rownames_to_column("dataset_siteDisease"),
                   by = c("batch" = "dataset_siteDisease")) %>% 
  dplyr::mutate(batch = study_site_disease,
                class = "individual",
                moderator_level = disease,
                study = Study) %>% 
  dplyr::select(feature, coef, stderr, class, batch, moderator_level, 
                study, study_full_name)
tb_result_MA <- result$result$meta.results %>% 
  dplyr::filter(stringr::str_detect(feature, 
                                    stringr::fixed("f__Enterobacteriaceae|g__|"))) %>% 
  dplyr::mutate(batch = "MA effect (combined)",
                class = "MA effect",
                moderator_level = "Combined",
                study = "MA effect",
                study_full_name = "MA effect") %>% 
  dplyr::select(feature, coef, stderr, class, batch, moderator_level, 
                study, study_full_name)
tb_plot <- rbind(tb_result_ind, 
                 tb_result_moderator %>% dplyr::select(-q), 
                 tb_result_MA) %>% 
  dplyr::mutate(study_full_name = factor(study_full_name, 
                                         levels = c(studies, "MA effect")),
                moderator_level = factor(moderator_level, 
                                         levels = c("CD", "UC", 
                                                    "Combined",
                                                    "Difference"))) %>% 
  dplyr::arrange(moderator_level, study_full_name) %>% 
  dplyr::mutate(batch_plot = factor(batch, levels = rev(batch)),
                moderator_color = moderator_level %>% 
                  dplyr::recode_factor("CD" = "CD",
                                       "UC" = "UC",
                                       "Combined" = "MA effect",
                                       "Difference" = "MA effect"))
colors_disease <- c("CD" = "red", "UC" = "green4")
p_steroids <- tb_plot %>% 
  ggplot(aes(x = batch_plot, y = coef)) +
  geom_point(size = 5, aes(shape = class, color = moderator_color)) +
  geom_errorbar(aes(ymin = coef - 1.96*stderr, 
                    ymax = coef + 1.96*stderr, color = moderator_color), 
                width = 0.5) +
  theme(axis.title.y = element_blank()) +
  scale_shape_manual(values = c("individual" = 16,
                                "MA effect" = 15),
                     name = "Effect",
                     labels = c("Individual cohort",
                                "MA effect")) +
  scale_color_manual(values = c(colors_disease, "MA effect" = "black"),
                     name = "Disease") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme(axis.title.y = element_blank(),
        title = element_text(size = 12)) +
  ylab("Effect size") +
  ggtitle(expression(paste("Steroids effects on ", 
                           italic("Enterobacteriaceae(f) unclassified")))) +
  geom_text(data = tb_result_moderator,
            aes(x = 1,
                y = coef + 1.96*stderr,
                label = paste0(" q = ", 
                               formatC(q, format = "e", digits = 2))),
            size = 5,
            hjust = 0)
ggsave(paste0("figures/figure4/assets/MA_",
              "steroids",
              ".pdf"),
       p,
       width = 7.5, height = 4.5)

cowplot::plot_grid(plotlist = l_p_montreal, align = "v", axis = "tb", nrow = 1) %>% 
  cowplot::plot_grid(cowplot::plot_grid(p_CD, p_steroids, align = "h", axis = "tb", nrow = 1),
                     nrow = 2,
                     labels = c("a", "b"),
                     align = "v",
                     axis = "l",
                     label_size = 30,
                     label_fontface = "plain",
                     rel_heights = c(5.5, 4.5)) %>% 
  ggsave("figures/figure4/figure4.pdf",
         .,
         width = 15, height = 10)

