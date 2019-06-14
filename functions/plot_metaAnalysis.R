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
  tb_plot <- result_all %>% 
    dplyr::mutate(test = factor(test, levels = tests)) %>% 
    dplyr::group_by(feature) %>% 
    dplyr::mutate(min_fdr = min(q.fdr_meta, na.rm = TRUE)) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(min_fdr < 0.05) %>% 
    dplyr::mutate(zval = coef/stderr) %>% 
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
  features_orderedByBeta <- tb_plot %>% 
    dplyr::select(test, feature, coef) %>% 
    dplyr::mutate(coef = ifelse(is.na(coef), 0, coef)) %>% 
    tidyr::spread(key = test, value = coef) %>% 
    as.data.frame() %>% 
    tibble::column_to_rownames("feature") %>% 
    as.matrix() %>% 
    dist() %>% 
    hclust(method = "average") %>% 
    extract2("order") %>% 
    `[`(unique(tb_plot$feature), .)
  p <- tb_plot %>% 
    dplyr::mutate(feature = factor(feature, levels = features_orderedByBeta)) %>%
    ggplot(aes(x = test, y = feature, fill = coef)) +
    geom_tile() +
    geom_text(aes(label = paste0(p_levels, "\n", q_levels))) +
    scale_fill_gradient2(low = "blue", high = "red", midpoint = 0) +
    smar::rotate_xaxis(15)
  # features_orderedByZval <- tb_plot %>% 
  #   dplyr::select(test, feature, zval) %>% 
  #   dplyr::mutate(zval = ifelse(is.na(zval), 0, zval)) %>% 
  #   tidyr::spread(key = test, value = zval) %>% 
  #   as.data.frame() %>% 
  #   tibble::column_to_rownames("feature") %>% 
  #   as.matrix() %>% 
  #   dist() %>% 
  #   hclust(method = "average") %>% 
  #   extract2("order") %>% 
  #   `[`(unique(tb_plot$feature), .)
  # p2 <- tb_plot %>% 
  #   dplyr::mutate(feature = factor(feature, levels = features_orderedByZval)) %>%
  #   ggplot(aes(x = test, y = feature, fill = zval)) +
  #   geom_tile() +
  #   geom_text(aes(label = paste0(p_levels, "\n", q_levels))) +
  #   scale_fill_gradient2(low = "blue", high = "red", midpoint = 0)
  ggsave(file = paste0(directory, name, ".pdf"),
         p, 
         width = 12 + length(tests)*1.5,
         height = length(unique(tb_plot$feature)) / 3,
         limitsize = FALSE)
  tb_plot_forest <- l_results[tests] %>% 
    purrr::imap_dfr(function(result, test) {
      result$result$ind.results %>% 
        purrr::reduce(rbind) %>% 
        dplyr::mutate(test = test)
    })
  pdf(paste0(directory, name, "_forests.pdf"),
      width = 9*length(tests), 
      height = tb_plot_forest %>% 
        dplyr::group_by(feature, test) %>% 
        dplyr::summarise(n_batch = length(unique(batch))) %>% 
        dplyr::pull(n_batch) %>% 
        max %>% magrittr::divide_by(2))
  for(i_feature in unique(tb_plot$feature)) {
    l_p <- tests %>% 
      purrr::map(function(i_test) {
        tb_plot_forest %>% 
          dplyr::filter(feature == i_feature) %>%
          dplyr::filter(test == i_test) %>% 
          dplyr::filter(!is.na(coef), !is.na(stderr)) %>% 
          dplyr::mutate(batch = factor(batch, levels = c("MA Effect", rev(batch)))) %>% 
          dplyr::add_row(batch = "MA Effect", 
                         coef = tb_plot %>% 
                           dplyr::filter(feature == i_feature, test == i_test) %>% 
                           dplyr::pull(coef),
                         stderr = tb_plot %>% 
                           dplyr::filter(feature == i_feature, test == i_test) %>% 
                           dplyr::pull(stderr)) %>% 
          ggplot(aes(x = batch, y = coef, color = batch == "MA Effect")) +
          geom_point() +
          geom_errorbar(aes(ymin = coef - 1.96*stderr, ymax = coef + 1.96*stderr), width = 0.5) +
          scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = FALSE) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          coord_flip() +
          ggtitle(i_test)
      })
    p <- cowplot::plot_grid(plotlist = l_p, nrow = 1) %>% 
      smar::cowplot_title(i_feature)
    print(p)
  }
  dev.off()
  return(NULL)
}

plot_metaAnalysis_moderator <- function(fit.lm.meta, 
                                        fit.rma.mod, 
                                        data.moderator,
                                        test, 
                                        q.cutoff = 0.05,
                                        directory) {
  
  result <- fit.lm.meta$meta.results %>% 
    dplyr::select(feature, exposure, coef, stderr, pval, k) %>% 
    dplyr::mutate(moderator = "",
                  moderator_level = "overall",
                  MA_effect = "overall",
                  test = test)
  
  levels.MA_effect <- lapply(colnames(data.moderator), 
                             function(variable.moderator) {
                               paste0(variable.moderator, ":",
                                      c(levels(data.moderator[, variable.moderator]),
                                        "difference"))
                             }) %>% 
    unlist()
  tb_plot <- result %>% 
    rbind(fit.rma.mod %>% 
            dplyr::mutate(MA_effect = paste0(moderator, ":", moderator_level))) %>% 
    dplyr::mutate(MA_effect = factor(MA_effect, levels = c("overall", levels.MA_effect))) %>% 
    dplyr::filter(k >= 2) %>% 
    dplyr::group_by(MA_effect) %>% 
    dplyr::mutate(q.fdr_meta = p.adjust(pval, method = "fdr")) %>% 
    dplyr::group_by(feature) %>% 
    dplyr::mutate(with_difference = any(moderator_level %in% "difference")) %>% 
    dplyr::ungroup() %>% dplyr::filter(with_difference) %>% 
    dplyr::group_by(feature) %>% 
    dplyr::mutate(min_fdr = min(q.fdr_meta[moderator_level %in% "difference"], na.rm = TRUE)) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(min_fdr < q.cutoff) %>% 
    dplyr::mutate(zval = coef/stderr) %>% 
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
  if(length(unique(tb_plot$feature)) > 4) {
    features_orderedByBeta <- tb_plot %>% 
      dplyr::select(MA_effect, feature, coef) %>% 
      dplyr::mutate(coef = ifelse(is.na(coef), 0, coef)) %>% 
      tidyr::spread(key = MA_effect, value = coef) %>% 
      as.data.frame() %>% 
      tibble::column_to_rownames("feature") %>% 
      as.matrix() %>% 
      dist() %>% 
      hclust(method = "average") %>% 
      extract2("order") %>% 
      `[`(unique(tb_plot$feature), .)
  } else features_orderedByBeta <- unique(tb_plot$feature)
  p <- tb_plot %>% 
    dplyr::mutate(feature = factor(feature, levels = features_orderedByBeta)) %>%
    ggplot(aes(x = MA_effect, y = feature, fill = coef)) +
    geom_tile() +
    geom_text(aes(label = paste0(p_levels, "\n", q_levels))) +
    scale_fill_gradient2(low = "blue", high = "red", midpoint = 0) +
    smar::rotate_xaxis(15)
  # features_orderedByZval <- tb_plot %>% 
  #   dplyr::select(test, feature, zval) %>% 
  #   dplyr::mutate(zval = ifelse(is.na(zval), 0, zval)) %>% 
  #   tidyr::spread(key = test, value = zval) %>% 
  #   as.data.frame() %>% 
  #   tibble::column_to_rownames("feature") %>% 
  #   as.matrix() %>% 
  #   dist() %>% 
  #   hclust(method = "average") %>% 
  #   extract2("order") %>% 
  #   `[`(unique(tb_plot$feature), .)
  # p2 <- tb_plot %>% 
  #   dplyr::mutate(feature = factor(feature, levels = features_orderedByZval)) %>%
  #   ggplot(aes(x = test, y = feature, fill = zval)) +
  #   geom_tile() +
  #   geom_text(aes(label = paste0(p_levels, "\n", q_levels))) +
  #   scale_fill_gradient2(low = "blue", high = "red", midpoint = 0)
  ggsave(file = paste0(directory, test, ".pdf"),
         p, 
         width = 12 + length(levels.MA_effect)*1.5,
         height = max(c(length(unique(tb_plot$feature)) / 3), 6),
         limitsize = FALSE)
  tb_plot_forest <- fit.lm.meta$ind.results %>% 
    purrr::reduce(rbind) %>% 
    dplyr::mutate(test = test) %>% 
    dplyr::filter(feature %in% tb_plot$feature)
  suppressWarnings(tb_plot_forest <- tb_plot %>% 
                     dplyr::transmute(feature = feature,
                                      value = exposure,
                                      coef = coef,
                                      stderr = stderr,
                                      pval = pval,
                                      batch = MA_effect, 
                                      moderator = moderator,
                                      moderator_level = moderator_level,
                                      test = test) %>% 
                     dplyr::bind_rows(tb_plot_forest) %>% 
                     dplyr::left_join(data.moderator %>% tibble::rownames_to_column("batch"), 
                                      by = c("batch")))
    
  pdf(paste0(directory, test, "_forests.pdf"),
      width = 12*ncol(data.moderator), 
      height = tb_plot_forest %>% 
        dplyr::group_by(feature) %>% 
        dplyr::summarise(n_batch = length(unique(batch))) %>% 
        dplyr::pull(n_batch) %>% 
        max %>% divide_by(2) %>% 
        c(6) %>% 
        max())
  for(i_feature in unique(tb_plot$feature)) {
    l_p <- colnames(data.moderator) %>% 
      purrr::map(function(variable.moderator) {
        levels.moderator <- levels(data.moderator[, variable.moderator])
        tb_plot_forest %>% 
          dplyr::filter(feature == i_feature) %>% 
          dplyr::mutate(group_variable = dplyr::case_when(
            !is.na(!!sym(variable.moderator)) ~ as.character(!!sym(variable.moderator)),
            moderator %in% c(variable.moderator, "") ~ "MA Effect",
            TRUE ~ NA_character_
          )) %>% 
          dplyr::filter(!is.na(group_variable), !(moderator_level %in% "difference")) %>% 
          dplyr::mutate(group_variable = factor(group_variable, 
                                                levels = c(levels.moderator, "MA Effect"))) %>% 
          dplyr::arrange(group_variable) %>% 
          dplyr::mutate(batch = factor(batch, 
                                       levels = c("overall",
                                                  rev(paste0(variable.moderator, ":",
                                                             levels.moderator)),
                                                  rev(batch[group_variable != "MA Effect"])))) %>% 
          ggplot(aes(x = batch, y = coef, color = group_variable == "MA Effect")) +
          geom_point() +
          geom_errorbar(aes(ymin = coef - 1.96*stderr, ymax = coef + 1.96*stderr), width = 0.5) +
          scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black"), guide = FALSE) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          coord_flip() +
          ggtitle(tb_plot %>% 
                    dplyr::filter(feature == i_feature, moderator == variable.moderator, 
                                  moderator_level == "difference") %>% 
                    dplyr::select(pval, q.fdr_meta) %>% 
                    unlist() %>% paste(collapse = "\n"))
      })
    p <- cowplot::plot_grid(plotlist = l_p, nrow = 1) %>% 
      smar::cowplot_title(i_feature, rel_heights = c(0.5, 10))
    print(p)
  }
  dev.off()
  return(NULL)
}
