load("results/simulations/discrete/tb_sim.RData")
load("results/simulations/discrete/evaluate_before.RData")
load("results/simulations/discrete/evaluate_after.RData")

suppFig_simulation_discrete <- 
  function(tb_sim, results_before, results_correct,
           path = "supp_materials/suppFigures/suppFig_simulation_discrete.pdf") {
    # Format parameter table
    # effect sizes
    tb_sim <- tb_sim %>% 
      dplyr::mutate(effect_cluster = effectSize %>% 
                      purrr::map_dbl("cluster"),
                    effect_batch = effectSize %>% 
                      purrr::map_dbl("batch") %>% 
                      forcats::as_factor() %>% forcats::fct_inseq(),
                    nCluster = nCluster %>% 
                      forcats::as_factor() %>% forcats::fct_inseq())
    # nSample and nBatch
    tb_sim <- tb_sim %>% 
      dplyr::arrange(nSample_perBatch, nBatch) %>% 
      dplyr::mutate(nSample_nBatch = paste0("nSample per batch = ", nSample_perBatch, ", ",
                                            "nBatch = ", nBatch) %>% 
                      forcats::as_factor())
    # nMicrobe and percMicrobe spiked
    tb_sim <- tb_sim %>% 
      dplyr::arrange(nMicrobe, spikeMicrobes) %>% 
      dplyr::mutate(nMicrobe_percSpike = paste0("nFeature = ", nMicrobe, ", ",
                                                spikeMicrobes*100, "% spiked") %>% 
                      forcats::as_factor())
    
    # summarize results
    tb_results <- list(results, results_correct) %>% 
      purrr::map2_dfr(
        c("Original", "MMUPHin corrected"),
        ~ .x %>% 
          purrr::imap_dfr(
            ~ tibble::tibble(results = list(.x$pred.str),
                             i = .y)
          ) %>% 
          dplyr::mutate(Data = .y)
      ) 
    # format variables
    tb_results <- tb_results %>% 
      dplyr::mutate(Data = factor(Data, levels = c("Original", "MMUPHin corrected")))
    # sanity check
    # all simulation runs were successful
    # mean(tb_results$results %>% purrr::map_lgl(~is.null(.x)))
    # format results
    # calculate successes
    tb_results <- tb_results %>% 
      dplyr::left_join(tb_sim, by = "i") %>% 
      dplyr::mutate(max_k = results %>% 
                      purrr::map_int(~ .x %>% 
                                       dplyr::arrange(-mean_pred) %>% 
                                       magrittr::extract2("k") %>% 
                                       magrittr::extract(1)),
                    success = (max_k == nCluster)) 
    # summarise success rates
    tb_summary <- tb_results %>%
      dplyr::group_by(i_setup, Data) %>% 
      dplyr::summarise(mean_success = mean(success)) %>% 
      dplyr::ungroup()
    
    # colors for the mapped methods
    colors <- smar::gg_color_hue(n=4)
    colors_simulation_discrete <- c("Original" = "black",
                                    "MMUPHin corrected" = colors[1])
    
    
    tb_plot <- tb_summary %>% 
      dplyr::left_join(dplyr::filter(tb_sim, !duplicated(i_setup)),
                       by = "i_setup") %>% 
      dplyr::filter(!(effect_batch == 0 & Data == "MMUPHin corrected"))
    plist <- (1:nlevels(tb_plot$nCluster)) %>% 
      purrr::map(function(i) {
        i.nCluster <- levels(tb_plot$nCluster)[i]
        p <- tb_plot %>% 
          dplyr::filter(nCluster == i.nCluster) %>% 
          ggplot(aes(x = effect_batch, y = mean_success, color = Data)) +
          geom_point(position = position_dodge(width = 0.5)) +
          facet_grid(nMicrobe_percSpike ~ nSample_nBatch) +
          scale_color_manual(values = colors_simulation_discrete) +
          theme(legend.position = c(1, 1),
                legend.justification = c(1, 1),
                legend.background = element_blank(),
                strip.text = element_text(size = 8),
                legend.text = element_text(size = 8)) +
          xlab("Simulatd strength of batch effect") +
          ylab("Success rate") +
          ggtitle(paste0("#Clusters = ", i.nCluster))
        if(i != 2)
          p <- p + theme(legend.position = "none")
        return(p)
      })
      cowplot::plot_grid(plotlist = plist,
                         nrow = 2,
                         labels = c("a", "b", "c", "d"),
                         label_fontface = "plain",
                         label_size = 20) %>% 
      ggsave(path,
             .,
             width = 20, height = 6, limitsize = FALSE)
    return(plist)
  }

ps <- suppFig_simulation_discrete(tb_sim = tb_sim, 
                                  results_before = results_before,
                                  results_correct = results_correct)

