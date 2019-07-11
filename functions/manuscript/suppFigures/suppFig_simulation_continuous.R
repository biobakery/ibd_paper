load("results/simulations/continuous/tb_sim.RData")
load("results/simulations/continuous/evaluate_before.RData")
load("results/simulations/continuous/evaluate_after.RData")

suppFig_simulation_continuous <- 
  function(tb_sim, results_before, results_correct,
           path = "supp_materials/suppFigures/suppFig_simulation_continuous.pdf") {
    # Format parameter table
    # effect sizes
    tb_sim <- tb_sim %>% 
      dplyr::mutate(effect_score = effectSize %>% 
                      purrr::map_dbl("score"),
                    effect_batch = effectSize %>%
                      purrr::map_dbl("batch") %>% 
                      factor(levels = 0:10) %>% 
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
    tb_summary <- list(results, results_correct) %>% 
      purrr::map2_dfr(
        c("Original", "MMUPHin corrected"),
        ~ .x %>% 
          purrr::imap_dfr(
            function(result, i_result) {
              if(.y == "Original")
                tibble::tibble(Cor = result$cor.pca,
                               i = i_result,
                               analysis = "naive analyzed")
              else tibble::tibble(Cor = c(result$cor.pca, 
                                          result$fit.continuous$correlation),
                                  i = i_result,
                                  analysis = c("naive analyzed", 
                                               "MMUPHin analyzed"))
            }
          ) %>% 
          dplyr::mutate(Data = .y)
      ) %>% 
      dplyr::mutate(`Data and model` = paste0(Data, ", ", analysis) %>% 
                      factor(levels = c("Original, naive analyzed",
                                        "MMUPHin corrected, naive analyzed",
                                        "MMUPHin corrected, MMUPHin analyzed")))
    
    # colors for the mapped methods
    tb_plot <- tb_summary %>% 
      dplyr::left_join(tb_sim, by = "i") %>% 
      dplyr::filter(!(effect_batch == 0 & `Data and model` != "Original, naive analyzed"))
    colors <- smar::gg_color_hue(n=4)
    colors_simulation_continuous <- c("Original, naive analyzed" = "black",
                                      "MMUPHin corrected, naive analyzed" = "brown",
                                      "MMUPHin corrected, MMUPHin analyzed" = colors[1])
    p <- tb_plot %>% 
      ggplot(aes(x = effect_batch, y = abs(Cor), color = `Data and model`)) +
      geom_boxplot() +
      facet_grid(nMicrobe_percSpike ~ nSample_nBatch) +
      scale_color_manual(values = colors_simulation_continuous,
                         name = "Data and model") +
      theme(legend.position = c(0, 0),
            legend.justification = c(0, 0),
            legend.background = element_blank(),
            legend.title = element_text(size = 8),
            strip.text = element_text(size = 8),
            legend.text = element_text(size = 8)) +
      xlab("Batch effect") +
      ylab("Correlation")
    ggsave(path,
           p,
           width = 20, height = 4, limitsize = FALSE)
    return(p)
  }

p <- suppFig_simulation_continuous(tb_sim = tb_sim, 
                                    results_before = results_before,
                                    results_correct = results_correct)

