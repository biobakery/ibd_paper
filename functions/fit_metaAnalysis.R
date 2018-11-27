fit_metaAnalysis <- function(mat_otu,
                             data,
                             test_variable,
                             contrasts,
                             batch_variable,
                             covariates = NULL,
                             covariates.random = NULL,
                             moderator_variables = NULL,
                             n_features = 50,
                             normalization = "TSS",
                             transform = "AST",
                             analysis_method = "LM",
                             rma_method = "REML",
                             directory,
                             forest.plots = TRUE,
                             verbose = TRUE) {
  if(!all(data$sample_accession_16S %in% colnames(mat_otu)))
    stop("Sample names of data must be present in mat_otu!")
  if(!all(c(test_variable, batch_variable, covariates, covariates.random,
            moderator_variables) %in% colnames(data)))
    stop("Some variable names aren't present in data!")
  if(!all(names(contrasts) %in% 
          dplyr::pull(data, 
                      !!sym(test_variable))))
    stop("All levels of contrasts must be in test variable values!")
  dir.create(directory, recursive = TRUE, showWarnings = TRUE)
  
  # prepare data
  data <- data %>% 
    dplyr::mutate(test_variable = !!sym(test_variable) %>% 
                    dplyr::recode(!!!contrasts,
                                  .missing = NA_character_))
  data <- data %>% 
    dplyr::group_by(!!sym(batch_variable)) %>% 
    dplyr::mutate(availale = dplyr::n_distinct(test_variable, 
                                               na.rm = TRUE) > 1) %>% 
    dplyr::ungroup() %>% 
    dplyr::filter(availale)
  data <- as.data.frame(data)
  rownames(data) <- data$sample_accession_16S
  mat_otu <- mat_otu[, data$sample_accession_16S]
  
  # fit MA models
  result_no_covariate <- MMUPHin::lm.meta(
    feature.count = mat_otu,
    batch = batch_variable,
    exposure = "test_variable",
    data = data,
    normalization = normalization,
    transform = transform,
    analysis_method = analysis_method,
    rma_method = rma_method,
    forest.plots = forest.plots,
    directory = paste0(directory, "no_covariate/"),
    verbose = verbose
  )
  result_no_random <- MMUPHin::lm.meta(
    feature.count = mat_otu,
    batch = batch_variable,
    exposure = "test_variable",
    covariates = covariates,
    data = data,
    normalization = normalization,
    transform = transform,
    analysis_method = analysis_method,
    rma_method = rma_method,
    forest.plots = forest.plots,
    directory = paste0(directory, "no_random_effect/"),
    verbose = verbose
  )
  result <- MMUPHin::lm.meta(
    feature.count = mat_otu,
    batch = batch_variable,
    exposure = "test_variable",
    covariates = covariates,
    covariates.random = covariates.random,
    data = data,
    normalization = normalization,
    transform = transform,
    analysis_method = analysis_method,
    rma_method = rma_method,
    forest.plots = forest.plots,
    directory = paste0(directory, "with_all/"),
    verbose = verbose
  )
  if(!is.null(moderator_variables)) {
    data.moderator <- data %>% 
      dplyr::select(dplyr::one_of(c(batch_variable, moderator_variables))) %>% 
      dplyr::group_by(!!sym(batch_variable)) %>% 
      dplyr::summarise_all(function(x) unique(x)) %>% 
      dplyr::mutate_at(moderator_variables, factor) %>% 
      as.data.frame() %>% 
      tibble::column_to_rownames(batch_variable)
    df_moderator <- MMUPHin:::rma.mod.wrapper(result$ind.results,
                                              data.moderator = data.moderator,
                                              method = rma_method)
  }
  
  # output & visualization
  # sensitivity analysis of meta-regression coefficients
  p1 <- rbind(
    data.frame(beta = result$meta.results$beta,
               beta_alt = result_no_covariate$meta.results$beta,
               alt_mod = "No covariates"),
    data.frame(beta = result$meta.results$beta,
               beta_alt = result_no_random$meta.results$beta,
               alt_mod = "No random covariates")) %>% 
    dplyr::group_by(alt_mod) %>% 
    dplyr::mutate(
      has_NA = is.na(beta) | is.na(beta_alt),
      beta = dplyr::case_when(
        !is.na(beta) ~ beta,
        is.na(beta) & !is.na(beta_alt) ~ min(beta, na.rm = TRUE) - 
          (max(beta, na.rm = TRUE) - min(beta, na.rm = TRUE)) / 10,
        TRUE ~ NA_real_
      ),
      beta_alt = dplyr::case_when(
        !is.na(beta_alt) ~ beta_alt,
        is.na(beta_alt) & !is.na(beta_alt) ~ min(beta_alt, na.rm = TRUE) - 
          (max(beta_alt, na.rm = TRUE) - min(beta_alt, na.rm = TRUE)) / 10,
        TRUE ~ NA_real_
      )) %>% 
    dplyr::ungroup() %>% 
    ggplot(aes(x = beta, y = beta_alt, color = has_NA)) + geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    facet_grid(.~alt_mod)
  batch_longitudinal <- 
    data[, c("dataset_name", batch_variable)] %>% 
    unique() %>% 
    dplyr::filter(dataset_name %in% studies_longitudinal) %>% 
    dplyr::pull(!!sym(batch_variable)) 
  if(is.null(covariates.random)) {
    ggsave(filename = paste0(directory, "sensitivity.pdf"),
           p1, width = 8, height = 5)
  }
  else {
    p2 <- rbind(
      result$ind.results[batch_longitudinal] %>% 
        purrr::reduce(rbind) %>% 
        dplyr::mutate(model = "with_random"),
      result_no_random$ind.results[batch_longitudinal] %>% 
        purrr::reduce(rbind) %>% 
        dplyr::mutate(model = "no_random")) %>% 
      dplyr::select(batch, feature, model, coef) %>% 
      tidyr::spread(key = model, value = coef) %>% 
      ggplot(aes(x = with_random, y = no_random)) +
      geom_abline(intercept = 0, slope = 1) +
      facet_grid(.~batch) +
      geom_point()
    cowplot::plot_grid(p1, p2, ncol = 1) %>% 
      ggsave(filename = paste0(directory, "sensitivity.pdf"),
             ., width = 8, height = 8)
  }
  
  # significant results
  tb_towrite <- result$meta.results
  if(!is.null(moderator_variables))
    tb_towrite <- tb_towrite %>% 
    dplyr::left_join(df_moderator,
                     by = c("feature", "exposure"),
                     suffix = c("", paste0("_",
                                           paste0(moderator_variables,
                                                  collapse = "+")))) 
  tb_towrite <- tb_towrite %>% 
    dplyr::arrange(exposure, pval) %>% 
    readr::write_tsv(paste0(directory, "results.tsv"))
  tb_toplot <- tb_towrite %>%
    dplyr::filter(k > 1) %>% 
    dplyr::mutate(pval.bonf = p.adjust(pval, method = "bonf"),
                  qval.fdr = p.adjust(pval, method = "fdr")) %>% 
    dplyr::top_n(n_features, -pval) %>% 
    dplyr::arrange(beta) %>% 
    dplyr::mutate(feature = feature %>% 
                    stringr::str_replace_all("\\|NA$", ""),
                  feature = feature %>% 
                    factor(levels = rev(feature)))
  p_effect <- tb_toplot %>% 
    ggplot(aes(x = feature, y = beta)) +
    geom_point() + geom_errorbar(aes(ymin = beta - se, ymax = beta + se)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip()
  p_pval <- tb_toplot %>% 
    ggplot(aes(x = feature, y = -log10(pval))) +
    geom_point() + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    coord_flip() + no_label_yaxis()
  p_qval <- tb_toplot %>% 
    ggplot(aes(x = feature, y = -log10(qval.fdr))) +
    geom_point() + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    coord_flip() + no_label_yaxis()
  cowplot::plot_grid(p_effect, p_pval, p_qval, nrow = 1, rel_widths = c(5, 1, 1)) %>% 
    ggsave(filename = paste0(directory, "effect.pdf"),
           ., width = 21, height = n_features/5)
  # heterogeneity
  if(is.null(moderator_variables)) colors_moderator <- gg_color_hue(n = 1)
  else colors_moderator <- gg_color_hue(n = 2)
  p_forest <- Reduce("rbind", result$ind.results) %>% 
    dplyr::mutate(feature = feature %>% 
                    stringr::str_replace_all("\\|NA$", "")) %>% 
    dplyr::filter(feature %in% tb_toplot$feature) %>% 
    dplyr::select(feature, 
                  exposure = value,
                  coef = coef,
                  batch) %>% 
    dplyr::left_join(tb_toplot %>% 
                       dplyr::select(feature,
                                     exposure,
                                     beta),
                     by = c("feature", "exposure")) %>% 
    dplyr::mutate(feature = factor(feature,
                                   levels = levels(tb_toplot$feature))) %>% 
    ggplot() +
    geom_boxplot(aes(x = feature, y = coef)) +
    geom_point(aes(x = feature, y = beta), color = "red") +
    coord_flip()
  p_I2 <- tb_toplot %>% 
    dplyr::select(feature, dplyr::matches("I2")) %>% 
    tidyr::gather(key = adjustment, value = I2, -feature) %>% 
    dplyr::mutate(adjustment = adjustment %>% 
                    stringr::str_replace_all(stringr::fixed("I2"), "") %>% 
                    stringr::str_replace_all(("^\\_"), "") %>% 
                    dplyr::recode("moderator" = "missing values"),
                  adjustment = ifelse(adjustment == "",
                                      "no adjustment",
                                      adjustment),
                  adjustment = adjustment %>% 
                    factor(levels = unique(adjustment))) %>% 
    ggplot(aes(x = feature, y = I2)) +
    geom_point(aes(color = adjustment)) +
    geom_line(aes(group = feature)) +
    scale_color_manual(values = c("black", colors_moderator), guide = FALSE) +
    coord_flip() + no_label_yaxis()
  p_R2 <- tb_toplot %>% 
    dplyr::select(feature, dplyr::matches("R2")) %>% 
    tidyr::gather(key = adjustment, value = R2, -feature) %>% 
    dplyr::mutate(adjustment = adjustment %>% 
                    stringr::str_replace_all(stringr::fixed("R2"), "") %>% 
                    stringr::str_replace_all(("^\\_"), ""),
                  adjustment = ifelse(adjustment == "",
                                      "missing values",
                                      adjustment),
                  adjustment = adjustment %>% 
                    factor(levels = unique(adjustment))) %>% 
    ggplot(aes(x = feature, y = R2)) +
    geom_point(aes(color = adjustment)) +
    geom_line(aes(group = feature)) +
    scale_color_manual(values = colors_moderator, guide = FALSE) +
    coord_flip() + no_label_yaxis()
  p_p.I2 <- tb_toplot %>% 
    dplyr::select(feature, dplyr::matches("p.tau2")) %>% 
    tidyr::gather(key = adjustment, value = p.I2, -feature) %>% 
    dplyr::mutate(adjustment = adjustment %>% 
                    stringr::str_replace_all(stringr::fixed("p.tau2"), "") %>% 
                    stringr::str_replace_all(("^\\_"), "") %>% 
                    dplyr::recode("moderator" = "missing values"),
                  adjustment = ifelse(adjustment == "",
                                      "no adjustment",
                                      adjustment),
                  adjustment = adjustment %>% 
                    factor(levels = unique(adjustment))) %>% 
    ggplot(aes(x = feature, y = -log10(p.I2))) +
    geom_point(aes(color = adjustment)) +
    geom_line(aes(group = feature)) +
    geom_hline(yintercept = -log10(0.05)) +
    scale_color_manual(values = c("black", colors_moderator)) +
    coord_flip() + no_label_yaxis() 
  p_comp <- tb_toplot %>% 
    dplyr::select(feature, dplyr::matches("weight_")) %>% 
    tidyr::gather(key = batch, value = weight, -feature) %>% 
    dplyr::mutate(batch = batch %>% 
                    stringr::str_replace_all("^weight_", "")) %>% 
    ggplot(aes(x = feature, y = weight, fill = batch)) +
    geom_bar(stat = "identity") +
    coord_flip() + no_label_yaxis()
  cowplot::plot_grid(p_forest, p_I2, p_R2, p_p.I2, p_comp,
                     nrow = 1, rel_widths = c(5, 1, 1, 1.5, 2.5)) %>% 
    ggsave(filename = paste0(directory, "heterogeneity.pdf"),
           ., width = 27.5, height = n_features/5)
  
  # save result
  result <- list(result = result,
                 result_no_covariate = result_no_covariate,
                 result_no_random = result_no_random,
                 tb_summarise = tb_towrite)
  
  return(result)
}

# This is adopted from MMUPHin::rma.mod.wrapper
# To examine effect modification in Maaslin meta-analysis
fit_rma.mod <- function(l.Maaslin.fit, data.moderator,
                        method = "REML"){
  lvl.batch <- names(l.Maaslin.fit)
  if(!all(lvl.batch %in% rownames(data.moderator)))
    stop("data.moderator must have all the batches fitted in Maaslin!")
  data.moderator <- data.moderator[lvl.batch, , drop = FALSE]
  if(!all(sapply(data.moderator, class) %in% "factor"))
    stop("data.moderator must have all columns be factor")
  if(!all(sapply(data.moderator, nlevels) == 2))
    stop("data.moderator must have all columns be two-level factors")
  if(is.null(colnames(data.moderator)))
    stop("data.moderator must have colnames!")
  
  exposure <- unique(l.Maaslin.fit[[1]]$metadata)
  values.exposure <- unique(l.Maaslin.fit[[1]]$value)
  features <- unique(l.Maaslin.fit[[1]]$feature)
  
  l.results <- list()
  for(value.exposure in values.exposure) {
    # sanity check
    if(any(features != l.Maaslin.fit[[2]][l.Maaslin.fit[[2]]$value == value.exposure, "feature"]))
      stop("Feature names don't match between l.Maaslin.fit components!")
    betas <- sapply(l.Maaslin.fit, function(i.Maaslin.fit) {
      i.Maaslin.fit[i.Maaslin.fit$value == value.exposure, "coef"]
    })
    sds <- sapply(l.Maaslin.fit, function(i.Maaslin.fit) {
      i.Maaslin.fit[i.Maaslin.fit$value == value.exposure, "stderr"]
    })
    rownames(betas) <- rownames(sds) <- features
    ind.feature <- !is.na(betas) & !is.na(sds) & (sds != 0)
    count.feature <- apply(ind.feature, 1, sum)
    
    for(variable.moderator in colnames(data.moderator)) {
      levels.moderator <- levels(data.moderator[, variable.moderator])
      i.result <- expand.grid(feature = features, 
                              exposure = value.exposure, 
                              moderator = variable.moderator,
                              moderator_level = c(levels.moderator, "difference"),
                              beta = NA_real_, 
                              se = NA_real_,
                              pval = NA_real_,
                              k = NA_integer_,
                              stringsAsFactors = FALSE)
      i.data.moderator <- data.frame(dummy_var = (data.moderator[, variable.moderator] ==
                                                    levels.moderator[2]) * 1)
      i.data.moderator2 <- data.frame(dummy_var = 1 - i.data.moderator$dummy_var)
      rownames(i.data.moderator) <- rownames(i.data.moderator2) <- rownames(data.moderator)
      
      for(feature in features) {
        if(count.feature[feature] <= 1) next
        suppressWarnings(tmp.rma.fit <-
                           try(metafor::rma.uni(yi = betas[feature, ind.feature[feature, ]],
                                                sei = sds[feature, ind.feature[feature, ]],
                                                mod = ~.,
                                                data = i.data.moderator[ind.feature[feature, ], ,
                                                                        drop = FALSE],
                                                method = method,
                                                control = list(threshold = 1e-10,
                                                               maxiter = 1000)),
                               silent = TRUE)) # FIXME
        if("try-error" %in% class(tmp.rma.fit))
          next
        if(is.null(tmp.rma.fit$R2))
          next
        # Sanity check
        if(nrow(tmp.rma.fit$beta) != 2)
          stop("Something went wrong!")
        
        i.result[i.result$feature == feature &
                   i.result$moderator_level %in% c(levels.moderator[1],"difference"), 
                 c("beta", "se", "pval", "k")] <- cbind(tmp.rma.fit$beta[, 1],
                                                        tmp.rma.fit$se,
                                                        tmp.rma.fit$pval,
                                                        tmp.rma.fit$k)
        
        tmp.rma.fit <-
          metafor::rma.uni(yi = betas[feature, ind.feature[feature, ]],
                           sei = sds[feature, ind.feature[feature, ]],
                           mod = ~.,
                           data = i.data.moderator2[ind.feature[feature, ], ,
                                                    drop = FALSE],
                           method = method,
                           control = list(threshold = 1e-10,
                                          maxiter = 1000))
        i.result[i.result$feature == feature &
                   i.result$moderator_level %in% c(levels.moderator[2]), 
                 c("beta", "se", "pval", "k")] <- c(tmp.rma.fit$beta[1, 1],
                                               tmp.rma.fit$se[1],
                                               tmp.rma.fit$pval[1],
                                               tmp.rma.fit$k)
        
      }
      
      l.results <- c(l.results, list(i.result))
    }
  }
  results <- Reduce("rbind", l.results)
  return(results)
}
