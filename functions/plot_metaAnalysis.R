plot_metaAnalysis <- function(l_results, tests, name, directory) {
  if(!all(tests %in% names(l_results)))
    stop("All of tests must be in l_results!")
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
                      test = test,
                      q.fdr_meta = ifelse(k > 2,
                                          pval,
                                          NA_real_) %>% 
                        p.adjust(method = "fdr"))
    })
  tb_toplot <- result_all %>% 
    dplyr::mutate(test = factor(test, levels = tests)) %>% 
    dplyr::group_by(feature) %>% 
    dplyr::mutate(min_fdr = min(q.fdr_meta, na.rm = TRUE)) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(min_fdr < 0.05) %>% 
    dplyr::mutate(zval = beta/se) %>% 
    dplyr::mutate(p_levels = 
                    dplyr::case_when(
                      pval < 0.001 ~ "...",
                      pval < 0.01 ~ "..",
                      pval < 0.05 ~ ".",
                      TRUE ~ ""
                    ),
                  q_levels = 
                    dplyr::case_when(
                      q.fdr_meta < 0.001 ~ "***",
                      q.fdr_meta < 0.01 ~ "**",
                      q.fdr_meta < 0.05 ~ "*",
                      TRUE ~ ""
                    ))
  features_orderedByBeta <- tb_toplot %>% 
    dplyr::select(test, feature, beta) %>% 
    dplyr::mutate(beta = ifelse(is.na(beta), 0, beta)) %>% 
    tidyr::spread(key = test, value = beta) %>% 
    as.data.frame() %>% 
    tibble::column_to_rownames("feature") %>% 
    as.matrix() %>% 
    dist() %>% 
    hclust(method = "average") %>% 
    extract2("order") %>% 
    `[`(unique(tb_toplot$feature), .)
  p <- tb_toplot %>% 
    dplyr::mutate(feature = factor(feature, levels = features_orderedByBeta)) %>%
    ggplot(aes(x = test, y = feature, fill = beta)) +
    geom_tile() +
    geom_text(aes(label = paste0(p_levels, "\n", q_levels))) +
    scale_fill_gradient2(low = "blue", high = "red", midpoint = 0) +
    rotate_xaxis(15)
  # features_orderedByZval <- tb_toplot %>% 
  #   dplyr::select(test, feature, zval) %>% 
  #   dplyr::mutate(zval = ifelse(is.na(zval), 0, zval)) %>% 
  #   tidyr::spread(key = test, value = zval) %>% 
  #   as.data.frame() %>% 
  #   tibble::column_to_rownames("feature") %>% 
  #   as.matrix() %>% 
  #   dist() %>% 
  #   hclust(method = "average") %>% 
  #   extract2("order") %>% 
  #   `[`(unique(tb_toplot$feature), .)
  # p2 <- tb_toplot %>% 
  #   dplyr::mutate(feature = factor(feature, levels = features_orderedByZval)) %>%
  #   ggplot(aes(x = test, y = feature, fill = zval)) +
  #   geom_tile() +
  #   geom_text(aes(label = paste0(p_levels, "\n", q_levels))) +
  #   scale_fill_gradient2(low = "blue", high = "red", midpoint = 0)
  ggsave(file = paste0(directory, name, ".pdf"),
         p, 
         width = 8 + length(tests),
         height = length(unique(tb_toplot$feature)) / 3,
         limitsize = FALSE)
  tb_toplot_forest <- l_results[tests] %>% 
    purrr::imap_dfr(function(result, test) {
      result$result$ind.results %>% 
        purrr::reduce(rbind) %>% 
        dplyr::mutate(test = test)
    })
  pdf(paste0(directory, name, "_forests.pdf"),
      width = 6*length(tests), 
      height = tb_toplot_forest %>% 
        dplyr::group_by(feature, test) %>% 
        dplyr::summarise(n_batch = length(unique(batch))) %>% 
        dplyr::pull(n_batch) %>% 
        max %>% divide_by(2))
  for(i_feature in unique(tb_toplot$feature)) {
    l_p <- tests %>% 
      purrr::map(function(i_test) {
        tb_toplot_forest %>% 
          dplyr::filter(feature == i_feature) %>%
          dplyr::filter(test == i_test) %>% 
          dplyr::filter(!is.na(coef), !is.na(stderr)) %>% 
          dplyr::mutate(batch = factor(batch, levels = c("MA Effect", rev(batch)))) %>% 
          dplyr::add_row(batch = "MA Effect", 
                         coef = tb_toplot %>% 
                           dplyr::filter(feature == i_feature, test == i_test) %>% 
                           dplyr::pull(beta),
                         stderr = tb_toplot %>% 
                           dplyr::filter(feature == i_feature, test == i_test) %>% 
                           dplyr::pull(se)) %>% 
          ggplot(aes(x = batch, y = coef, color = batch == "MA Effect")) +
          geom_point() +
          geom_errorbar(aes(ymin = coef - 1.96*stderr, ymax = coef + 1.96*stderr), width = 0.5) +
          scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = FALSE) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          coord_flip() +
          ggtitle(i_test)
      })
    p <- cowplot::plot_grid(plotlist = l_p, nrow = 1) %>% 
      cowplot_title(i_feature)
    print(p)
  }
  dev.off()
  return(NULL)
}
