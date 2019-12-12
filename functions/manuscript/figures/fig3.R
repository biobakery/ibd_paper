load("results/4.1-permanova/permanova_config.RData")
l_adonis_fit <- list()
for(i in 1:nrow(tb_permanova)) {
  load(paste0("results/4.1-permanova/adonis_fit_", i, ".RData"))
  l_adonis_fit[[i]] <- list(adonis_fit) # can't just assign adonis_fit because setting an
  # element of a list to NULL removes that element
}

fig_PERMANOVA <- function(l_adonis_fit, tb_permanova, 
                          path = "figures/figure3/assets/") {
  # Extract and format R2 and q values
  tb_permanova <- tb_permanova %>% 
    dplyr::mutate(adonis_fit = l_adonis_fit) %>% 
    dplyr::mutate(R2 = purrr::map2_dbl(variable, adonis_fit,
                                       function(i_variable, i_adonis_fit) {
                                         if(is.null(i_adonis_fit[[1]])) return(NA_real_)
                                         i_adonis_fit[[1]]$aov.tab[i_variable, "R2"]
                                       }),
                  p = purrr::map2_dbl(variable, adonis_fit,
                                      function(i_variable, i_adonis_fit) {
                                        if(is.null(i_adonis_fit[[1]])) return(NA_real_)
                                        i_adonis_fit[[1]]$aov.tab[i_variable, "Pr(>F)"]
                                      })) %>% 
    dplyr::group_by(study) %>% 
    dplyr::mutate(q = p.adjust(p, "fdr")) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(qlevel = dplyr::case_when(
      q < 0.001 ~ "***",
      q < 0.01 ~ "**",
      q < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label_plot = paste0(round(R2*100, digits = 2), "%", "\n",
                        qlevel) %>% 
      dplyr::recode("NA%\n" = ""))
  
  # Format studies
  source("assets/tb_1.R")
  tb_permanova <- tb_permanova %>% 
    dplyr::left_join(tb_1, by = c("study" = "study_full_name")) %>% 
    dplyr::mutate(Study = ifelse(study == "all", 
                                 "All",
                                 Study))
  
  # Format variables
  vars_PERMANOVA_map <- c("dataset_name" = "Study", 
                          "batch" = "Batch",
                          "sample_type" = "Biopsy vs. stool", 
                          "body_site" = "Biopsy location",
                          "IBD" = "IBD vs. control",
                          "disease" = "CD vs. UC",
                          "L.cat" = "CD Location", 
                          "B.cat" = "CD Behavior", 
                          "E.cat" = "UC Extent",
                          "antibiotics" = "Antibiotics", 
                          "immunosuppressants" = "Immunosuppressants", 
                          "steroids" = "Steroids", 
                          "mesalamine_5ASA" = "5-ASA", 
                          "age.cat" = "Age (< 18)", 
                          "age_at_diagnosis.cat" = "Age at diagnosis", 
                          "gender" = "Gender", 
                          "race" = "Race",
                          "subject_accession" = "Subject")
  tb_permanova <- tb_permanova %>% 
    dplyr::mutate(variable = variable %>% 
                    dplyr::recode_factor(!!!vars_PERMANOVA_map)) %>% 
    dplyr::arrange(Study, variable, (physeq != "original")*1) %>% 
    dplyr::mutate(variable = ifelse(physeq != "original",
                                    as.character(variable),
                                    paste0(variable, " (unadjusted)")) %>% 
                    forcats::as_factor() %>% 
                    forcats::fct_rev())
  
  p <- tb_permanova %>%
    ggplot(aes(x = Study, y = variable, fill = R2)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "black", 
                        trans = "log10", na.value = "white") +
    geom_text(aes(label = label_plot,
                  color = R2 < quantile(R2, 0.25, na.rm = TRUE))) +
    scale_color_manual(values = c("TRUE" = "black", "FALSE" = "white"), guide = FALSE) +
    scale_x_discrete(position = "top") +
    theme(line = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          axis.text.y = element_text(size = 14, colour="black"),
          axis.text.x.top = element_text(size = 14, vjust = -0.05, colour="black")) +
    smar::rotate_xaxis(30, vjust = 0, hjust = 0) +
    theme(plot.margin=unit(c(0, 0, 0, 0), "points"))
  
  ggsave(paste0(path, "permanova.pdf"),
         p, width = 9, height = 8)
  return(p)
}

p_permanova <- fig_PERMANOVA(l_adonis_fit = l_adonis_fit, tb_permanova = tb_permanova)

source("assets/ma_tests.R")
l_results <- list()
for(test in c(tests_disease, tests_treatment)) {
  load(paste0("results/4.2-maaslin2/genera/", test, "/result.RData"))
  l_results[[test]] <- result
}
load("data/physeq/genera_adj.RData")
tax_table <- smar::tax_table2(physeq_genera_adj)
fig_MADiseaseTreatment <- function(l_results,
                                   tests_disease,
                                   tests_treatment,
                                   tax_table,
                                   path = "figures/figure3/assets/") {
  list(disease = tests_disease,
       treatment = tests_treatment) %>% 
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
        })
      tb_plot <- result_all %>% 
        dplyr::filter(k >= 2) %>% 
        dplyr::group_by(feature) %>% # calculate minimal q value across tests
        dplyr::mutate(min_qval.fdr = min(qval.fdr_meta, na.rm = TRUE)) %>% 
        dplyr::ungroup() %>% 
        dplyr::filter(min_qval.fdr < 0.05) %>% 
        dplyr::mutate(p_levels = 
                        dplyr::case_when(
                          pval < 0.001 ~ "...",
                          pval < 0.01 ~ "..",
                          pval < 0.05 ~ ".",
                          TRUE ~ ""
                        ),
                      q_levels = 
                        dplyr::case_when(
                          qval.fdr_meta < 0.001 ~ "***",
                          qval.fdr_meta < 0.01 ~ "**",
                          qval.fdr_meta < 0.05 ~ "*",
                          TRUE ~ ""
                        )) %>% 
        dplyr::left_join(tax_table %>% 
                           as.data.frame() %>% 
                           tibble::rownames_to_column("feature"),
                         by = "feature") %>% # create better feature names
        dplyr::filter(Rank5 != "f__") %>%
        dplyr::mutate(feature_plot = betterGeneraNames(Rank4, Rank5, Rank6),
                      test = factor(test, levels = rev(tests)))
      feature_plot_orderedByBeta <- tb_plot %>%
        dplyr::filter(test %in% tests[1]) %>% 
        dplyr::group_by(sign(coef), Rank5) %>%
        dplyr::mutate(max_coef = sign(coef)*max(abs(coef), na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(sign(coef), max_coef, coef) %>%
        dplyr::pull(feature_plot)
      p <- tb_plot %>%  
        dplyr::mutate(feature_plot = factor(feature_plot, levels = rev(feature_plot_orderedByBeta))) %>%
        dplyr::mutate(coef = ifelse(coef > 0.1, 0.1, ifelse(coef < -0.1, -0.1, coef))) %>%
        dplyr::mutate(test = test %>% dplyr::recode_factor("IBD_vs_control" = "IBD vs. control",
                                                           "CD_vs_control" = "CD vs. control",
                                                           "UC_vs_control" = "UC vs. control",
                                                           "CD_vs_UC" = "CD vs. UC",
                                                           "antibiotics" = "Antibiotics",
                                                           "immunosuppressants" = "Immunosuppressants",
                                                           "steroids" = "Steroids",
                                                           "mesalamine_5ASA" = "5-ASA"
                                                           )) %>% 
        ggplot(aes(x = test, y = feature_plot, fill = coef)) +
        geom_tile() +
        geom_text(aes(label = q_levels), size = 5) +
        scale_fill_gradient2(low = "purple", high = "red", midpoint = 0, breaks = c(-0.1, -0.05,
                                                                                    0, 0.05, 0.1),
                             limits = c(-0.1, 0.1),
                             name = "Effect size") +
        smar::rotate_xaxis(30, vjust = 0, hjust = 0) +
        scale_x_discrete(position = "top") +
        scale_y_discrete(position = "right", 
                         labels = paste0(rev(feature_plot_orderedByBeta)," <-")) +
        theme(line = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              panel.border = element_blank(),
              axis.text.y.right = element_text(size = 14, face = "italic", colour="black"),
              axis.text.x.top = element_text(size = 14, vjust = -0.05, colour="black"))
      
      ggsave(paste0(path, "MaaslinPlot_", name_tests, ".pdf"),
             p,
             width = 5, height = 8)
      
      if(name_tests == "treatment") {
        p <-  p + 
          coord_cartesian(xlim = c(0, 4), # This focuses the x-axis on the range of interest
                          clip = 'off') +
          annotate(geom = "text",
                   label = "* q < 0.05\n** q < 0.01\n*** q < 0.001",
                   x = 9.1,
                   y = 5,
                   hjust = 0,
                   vjust = 0,
                   size = 6,
                   fontface = "plain")  +
          theme(plot.margin=unit(c(0, 10, 0, 0), "points"))
      } else {
        p <- p + 
          # cowplot::draw_plot_label(label = "",
          #                          x = 5,
          #                          y = 0.5,
          #                          hjust = 0,
          #                          vjust = 0,
          #                          size = 16,
          #                          fontface = "plain") +
          theme(legend.position = "none",
                plot.margin=unit(c(0, 0, 0, 0), "points"))
        # this is needed so alignment can be done correctly with cowplot::plot_grid
      }
      return(p)
    })
}
l_p_MA <- fig_MADiseaseTreatment(l_results = l_results,
                                 tests_disease = tests_disease,
                                 tests_treatment = tests_treatment,
                                 tax_table = tax_table)
cowplot::plot_grid(plotlist = l_p_MA,
                   nrow = 1,
                   rel_widths = c(2.22, 3),
                   align = "hv",
                   axis = "tb") %>% 
  cowplot::plot_grid(p_permanova, ., 
                     nrow = 1, align = "v", axis = "tb",
                     labels = c("a", "b"),
                     rel_widths = c(2, 3),
                     label_size = 30,
                     label_fontface = "plain") %>% 
  ggsave("figures/figure3/figure3.pdf", .,
         width = 24,
         height = 10)
