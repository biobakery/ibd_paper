is_outlier <- function(x, threshold_upper = -Inf, threshold_lower = -Inf) {
  return((x < quantile(x, 0.25) - 1.5 * IQR(x) & x > threshold_lower) | 
           (x > quantile(x, 0.75) + 1.5 * IQR(x) & x > threshold_upper))
}

is_outlier_label <- function(x, adj = 0) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) - adj | 
         x > quantile(x, 0.75) + 1.5 * IQR(x) + adj)
}

is_outlier_upper_label <- function(x, adj = 0) {
  return(x > quantile(x, 0.75) + 1.5 * IQR(x) + adj)
}

suppFig_meidumHighPrevalent <- 
  function(physeq,
           path = "supp_materials/suppFigures/suppFig_mediumHighPrev.pdf") {
    tb_long <- physeq %>% 
      to_relativeAbundance() %>% 
      smar::phyloseq_to_tb()
    tb_summarise <- tb_long %>% 
      dplyr::group_by(feature) %>% 
      dplyr::mutate(prevalence = mean(abundance > 0)) %>% 
      dplyr::ungroup() %>% 
      dplyr::group_by(feature, dataset_name, prevalence) %>% 
      dplyr::summarise(prevalence_study = mean(abundance > 0),
                       mean_abundance = mean(abundance)) %>% 
      dplyr::group_by(feature) %>% 
      dplyr::mutate(n_absent = sum(prevalence_study == 0),
                    max_mean_abundance = max(mean_abundance)) %>% 
      dplyr::ungroup() %>% 
      dplyr::left_join(tb_1, 
                       by = c("dataset_name" = "study_full_name")) %>% 
      dplyr::mutate(feature = feature %>% 
                      stringr::str_replace_all("\\|NA.*", ""))
    
    tb_readDepth_summary <- physeq %>% 
      smar::sample_data2() %>% 
      dplyr::group_by(dataset_name) %>% 
      dplyr::summarise(median_read_depth = median(log2(read_depth)))
    
    tb_mediumPrev <- tb_summarise %>% 
      dplyr::filter(prevalence <= 0.7, n_absent == 0) %>% 
      dplyr::group_by(feature) %>% 
      dplyr::mutate(spread = IQR(prevalence_study),
                    label = ifelse(is_outlier_label(prevalence_study, adj = 0.05),
                                   Study, 
                                   NA_character_),
                    prevalence_study_not_outlier = ifelse(is_outlier(prevalence_study),
                                                          NA_real_,
                                                          prevalence_study)) %>% 
      dplyr::ungroup() %>% 
      dplyr::arrange(spread) %>% 
      dplyr::mutate(feature = factor(feature, levels = unique(feature))) %>% 
      dplyr::left_join(tb_readDepth_summary, by = "dataset_name") %>% 
      dplyr::group_by(feature) %>% 
      dplyr::mutate(rho.spearman = cor(prevalence_study, median_read_depth, method = "spearman") %>% 
                      ifelse((1:dplyr::n()) == 1, ., NA_real_),
                    p.spearman = 
                      cor.test(prevalence_study, median_read_depth, 
                               method = "spearman", 
                               alternative = "greater",
                               exact = FALSE)$p.value %>% 
                      ifelse((1:dplyr::n()) == 1, ., NA_real_)) %>% 
      dplyr::ungroup() %>% 
      dplyr::mutate(p.adj = p.adjust(p.spearman, method = "fdr")) %>% 
      dplyr::mutate(annotation = 
                      dplyr::case_when(is.na(p.spearman) ~ NA_character_,
                                       p.adj < 0.001 ~ "***",
                                       p.adj < 0.01 ~ "**",
                                       p.adj < 0.05 ~ "*"))
    
    p_mediumPrev <- tb_mediumPrev %>% 
      ggplot(aes(x = feature, y = prevalence_study)) +
      geom_boxplot(outlier.alpha = 0.5,
                   outlier.size = 0.5,
                   color = "grey50") +
      geom_point(aes(y = prevalence_study_not_outlier), 
                 position = position_jitter(width = 0.2),
                 alpha = 0.5,
                 size = 0.5) +
      ggrepel::geom_text_repel(aes(label = label),
                               angle = 270,
                               min.segment.length = 0, 
                               size = 3.5) +
      xlab("Feature") + ylab("Study prevalence") +
      coord_flip() +
      theme(axis.title.y = element_blank(),
            axis.title.x = element_text(size = 12),
            axis.text = element_text(size = 10)) 
    
    p_cor <- tb_mediumPrev %>%
      ggplot(aes(x = feature, y = rho.spearman)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = annotation), size = 6,
                hjust = 0.5,
                angle = 270) +
      coord_flip() +
      annotate(geom = "text", 
               label = "* FDR q < 0.05\n** FDR q < 0.05", 
               x = 1,
               y = 1.25,
               hjust = 1, vjust = 0,
               size = 3.5) +
      scale_y_continuous(breaks = c(-0.5, 0, 0.5, 1)) +
      ylab("Correlation with study median log2 read depth") +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_text(size = 12),
            axis.text = element_text(size = 10)) 
    
    # check Bifido age correlation
    tb_age_summary <- physeq %>% 
      smar::sample_data2() %>% 
      dplyr::group_by(dataset_name) %>% 
      dplyr::summarise(median_age = median(age, na.rm = TRUE))
    tb_summarise %>% 
      dplyr::filter(feature %in% 
                      "k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Bifidobacteriales|f__Bifidobacteriaceae|g__Bifidobacterium|s__bifidum" ) %>% 
      dplyr::left_join(tb_age_summary, by = "dataset_name") %>% 
      dplyr::summarise(rho.spearman = cor(prevalence_study, median_age, method = "spearman",
                                          use = "pairwise.complete.obs"),
                       p.spearman = cor.test(prevalence_study, median_age, 
                                             method = "spearman", 
                                             alternative = "less", 
                                             exact = FALSE)$p.value) %>% 
      print()
    
    p_highPrev <- tb_summarise %>% 
      dplyr::filter(prevalence > 0.7) %>% 
      dplyr::group_by(feature) %>% 
      dplyr::mutate(spread = IQR(mean_abundance),
                    label = ifelse(is_outlier_label(mean_abundance, adj = 0.005),
                                   Study, 
                                   NA_character_),
                    mean_abundance_not_outlier = ifelse(is_outlier(mean_abundance),
                                                        NA_real_,
                                                        mean_abundance)) %>% 
      dplyr::ungroup() %>% 
      dplyr::arrange(spread) %>% 
      dplyr::mutate(feature = factor(feature, levels = unique(feature))) %>% 
      ggplot(aes(x = feature, y = mean_abundance)) +
      geom_boxplot(color = "grey50",
                   outlier.size = 0.5, 
                   outlier.alpha = 0.5) +
      geom_point(aes(x = feature, y = mean_abundance_not_outlier),
                 position = position_jitter(width = 0.2),
                 size = 0.5,
                 alpha = 0.5) +
      ggrepel::geom_text_repel(aes(label = label),
                               angle = 270,
                               size = 3.5,
                               min.segment.length = 0) +
      xlab("Feature") + ylab("Study mean relative abundance") +
      coord_flip() +
      theme(axis.title.y = element_blank(),
            axis.title.x = element_text(size = 12),
            axis.text = element_text(size = 10)) 
    
    p_list <- list(p_mediumPrev = p_mediumPrev,
                   p_cor = p_cor,
                   p_highPrev = p_highPrev)
    if(!is.null(path))
      cowplot::plot_grid(p_mediumPrev, p_cor, rel_widths = c(3, 1), nrow = 1) %>% 
      cowplot::plot_grid(p_highPrev, rel_heights = c(5, 2.5), ncol = 1,
                         labels = c("a", "b"),
                         label_size = 16) %>% 
      ggsave(filename = path, 
             .,
             height = 15,
             width = 16)
    return(p_list)
  }

suppFig_featureMissingFromOne <- 
  function(physeq,
           path = "supp_materials/suppFigures/suppFig_featureMissingFromOne.pdf") {
    tb_long <- physeq %>% 
      to_relativeAbundance() %>% 
      smar::phyloseq_to_tb()
    tb_summarise <- tb_long %>% 
      dplyr::group_by(feature) %>% 
      dplyr::mutate(prevalence = mean(abundance > 0)) %>% 
      dplyr::ungroup() %>% 
      dplyr::group_by(feature, dataset_name, prevalence) %>% 
      dplyr::summarise(prevalence_study = mean(abundance > 0),
                       mean_abundance = mean(abundance)) %>% 
      dplyr::group_by(feature) %>% 
      dplyr::mutate(n_absent = sum(prevalence_study == 0),
                    max_mean_abundance = max(mean_abundance)) %>% 
      dplyr::ungroup() %>% 
      dplyr::left_join(tb_1, 
                       by = c("dataset_name" = "study_full_name")) %>% 
      dplyr::mutate(feature = feature %>% 
                      stringr::str_replace_all("\\|NA.*", ""))
    
    p_featureMissingFromOne <- tb_summarise %>% 
      dplyr::filter(prevalence <= 0.7, n_absent >= 1, prevalence_study > 0) %>% 
      dplyr::group_by(feature, n_absent) %>% 
      dplyr::mutate(n_prevalent = dplyr::n(),
                    include = max(prevalence_study) > 0.05,
                    spread = IQR(mean_abundance),
                    label = ifelse(is_outlier_upper_label(mean_abundance, adj = 0.001),
                                   Study,
                                   NA_character_),
                    abundance_study_not_outlier = ifelse(is_outlier(mean_abundance),
                                                          NA_real_,
                                                         mean_abundance),
                    study_present = ifelse(n_absent[1] == 9,
                                           Study[1],
                                           NA_character_)) %>% 
      dplyr::ungroup() %>% 
      dplyr::filter(include) %>% 
      dplyr::arrange(-n_prevalent) %>% 
      dplyr::mutate(n_prevalent_group = paste0("# Studies prevalent ", n_prevalent) %>% 
                      forcats::as_factor()) %>% 
      dplyr::arrange(n_prevalent, study_present, spread, mean_abundance) %>% 
      dplyr::mutate(feature = factor(feature, levels = unique(feature))) %>% 
      ggplot(aes(x = feature, y = mean_abundance)) +
      geom_boxplot(outlier.size = 0.5,
                   outlier.alpha = 0.5,
                   color = "grey50") +
      geom_point(aes(x = feature, y = abundance_study_not_outlier),
                 position = position_jitter(width = 0.2),
                 size = 0.5,
                 alpha = 0.5) +
      ggrepel::geom_text_repel(aes(label = label),
                               angle = 270,
                               size = 2,
                               min.segment.length = 0) +
      geom_text(aes(label = study_present,
                    y = mean_abundance + 0.001, 
                    hjust = 0), size = 2) +
      facet_wrap(.~n_prevalent_group, scales = "free_y", nrow = 3) +
      theme(axis.text.y = element_text(size = 5)) +
      xlab("Feature") + ylab("Study mean relative abundance") +
      coord_flip() +
      theme(axis.title.y = element_blank(),
            axis.title.x = element_text(size = 10),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            strip.text = element_text(size = 8))
    if(!is.null(path))
      ggsave(filename = path, 
             p_featureMissingFromOne,
             height = 15,
             width = 30)
    return(p_featureMissingFromOne)
  }
