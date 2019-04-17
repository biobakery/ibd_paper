suppFig_readDepth <- function(phylo_readCount, # for plotting read counts
                              phylo_prefilter, # these two are for generating 
                              phylo_postfilter,  # nFeature vs. lib size plots
                              path = "supp_materials/suppFigures/suppFig_readDepth.pdf", 
                              width = 10, height = 10, 
                              rel_heights = c(4, 3), 
                              nrow = 2, 
                              labels = c("a", "b")) {
  library(ggplot2)
  tb_compareReadCounts <- sample_data2(phylo_readCount) %>% 
    tidyr::gather(key = "Type of read count",
                  value = "N",
                  original_read_count,
                  read_depth,
                  reads_mapping_to_OTU_with_taxonomy,
                  reads_mapping_to_unclassified_OTU,
                  factor_key = TRUE) %>% 
    dplyr::mutate(`Type of read count` = `Type of read count` %>% 
                    dplyr::recode_factor("original_read_count" = "Original read count",
                                         "read_depth" = "Read depth",
                                         "reads_mapping_to_OTU_with_taxonomy" =
                                           "Reads mapped to OTUs with taxonomy",
                                         "reads_mapping_to_unclassified_OTU" =
                                           "Reads mapped to unclassified OTUs")) %>% 
    dplyr::left_join(tb_1,
                     by = c("dataset_name" = "study_full_name"))
  p_compareReadCounts <- tb_compareReadCounts %>% 
    dplyr::filter(`Type of read count` != "Read depth") %>% 
    ggplot(aes(x = Study, y = log2(N + 1), color = `Type of read count`)) + 
    geom_boxplot() + rotate_xaxis(30) +
    ylab("log2(read count + 1)") + 
    theme(legend.position = c(0, 1), 
          legend.justification = c(0, 1), 
          legend.background = element_blank(), 
          legend.title = element_blank(),
          axis.title.x = element_blank())
  
  
  tb_nFeatureVsReadDepth <- list(pre_filter = phylo_prefilter,
                                 post_filter = phylo_postfilter) %>% 
    purrr::imap_dfr(function(phylo, filtering) {
      phylo %>% phyloseq_to_tb %>% 
        dplyr::group_by(dataset_name, feature) %>% 
        dplyr::summarise(present = any(abundance > 0),
                         median_read_depth = median(log2(read_depth))) %>% 
        dplyr::group_by(dataset_name) %>% 
        dplyr::summarise(n_present = sum(present),
                         median_read_depth = unique(median_read_depth)) %>% 
        dplyr::ungroup() %>% 
        dplyr::mutate(filtering = filtering)
    }) %>% 
    dplyr::mutate(filtering = filtering %>% 
                    dplyr::recode_factor("pre_filter" = "Pre filtering",
                                         "post_filter" = "Post filtering")) %>% 
    dplyr::group_by(filtering) %>% 
    dplyr::mutate(rho_spearman = cor(median_read_depth, n_present, method = "spearman"),
                  p_spearman = cor.test(median_read_depth, n_present, 
                                        method = "spearman", alternative = "greater",
                                        exact = FALSE)$p.value,
                  label = ifelse((1:n()) == 1,
                                 paste0("Spearman correlation = ", 
                                        round(rho_spearman, digits = 3),
                                        "\nOne-sided p = ", 
                                        round(p_spearman, digits = 3)),
                                 NA_character_),
                  x_pos = min(median_read_depth),
                  y_pos = max(n_present)) %>% 
    dplyr::ungroup()
  p_nFeatureVsReadDepth <- tb_nFeatureVsReadDepth %>% 
    dplyr::left_join(tb_1,
                     by = c("dataset_name" = "study_full_name")) %>% 
    ggplot(aes(x = median_read_depth,
               y = n_present)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = Study)) +
    facet_wrap(~filtering,
               scales = "free_y") +
    geom_text(aes(x = x_pos, y = y_pos, label = label), hjust = 0, vjust = 0.8) +
    xlab("Median log2 read depth") + ylab("# Present species")
  
  p_list <- list(p_compareReadCounts = p_compareReadCounts,
                 p_nFeatureVsReadDepth = p_nFeatureVsReadDepth)
  if(!is.null(path)) {
    ggsave(path, 
           cowplot::plot_grid(plotlist = p_list, 
                              nrow = nrow, rel_heights = rel_heights, 
                              labels = labels,
                              label_size = 16),
           width = width,
           height = height)
  }
  return(p_list)
}
