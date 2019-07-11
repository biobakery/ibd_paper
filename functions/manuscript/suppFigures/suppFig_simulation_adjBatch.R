load("results/simulations/adjust.batch/tb_sim.RData")
load("results/simulations/adjust.batch/evaluate_before.RData")
load("results/simulations/adjust.batch/evaluate_after.RData")
load("results/simulations/adjust.batch/evaluate_ComBat.RData")

suppFig_simulation_adjBatch <- 
  function(tb_sim, results_before, results_correct, results_ComBat,
           path = "supp_materials/suppFigures/suppFig_simulation_adjBatch.pdf") {
    # Format parameter table
    # effect sizes
    tb_sim <- tb_sim %>% 
      dplyr::mutate(effect_poscontrol = effectSize %>% 
                      purrr::map_dbl("positive_control1"),
                    effect_batch = effectSize %>% 
                      purrr::map_dbl("batch") %>% 
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
    
    # summarize R2s
    tb_R2 <- list(results, results_ComBat, results_correct) %>% 
      purrr::map2_dfr(
        c("Original", "ComBat corrected", "MMUPHin corrected"),
        ~ .x %>% 
          purrr::imap_dfr(
            ~ tibble::tibble(Variable = names(.x$R2),
                             R2 = .x$R2,
                             i = .y)
          ) %>% 
          dplyr::mutate(Data = .y)
      ) 
    # format variables
    tb_R2 <- tb_R2 %>% 
      dplyr::mutate(
        Variable = Variable %>% 
          dplyr::recode_factor(
            "batch" = "Batch",
            "positive_control1" = "Positive control (binary)",
            "positive_control2" = "Positive control (continuous)",
            "negative_control" = "Negative control"
          ),
        Data = factor(Data, levels = c("Original",
                                       "ComBat corrected",
                                       "MMUPHin corrected"))
      )
    
    # colors for the mapped methods
    colors <- smar::gg_color_hue(n=4)
    colors_simulation_adjust.batch <- c("Original" = "black",
                                        "ComBat corrected" = colors[2],
                                        "MMUPHin corrected" = colors[1])
    
    tb_plot <- tb_R2 %>% dplyr::left_join(tb_sim, by = "i")
    plist <- (1:nlevels(tb_plot$Variable)) %>% 
      purrr::map(function(i) {
        i.Variable <- levels(tb_plot$Variable)[i]
        p <- tb_plot %>% 
          dplyr::filter(Variable == i.Variable) %>% 
          ggplot(aes(x = effect_batch, y = R2, color = Data)) +
          geom_boxplot() +
          facet_grid(nMicrobe_percSpike ~ nSample_nBatch,
                     scales = "free_y") +
          scale_color_manual(values = colors_simulation_adjust.batch) +
          theme(legend.position = c(0, 1),
                legend.justification = c(0, 1),
                legend.background = element_blank(),
                strip.text = element_text(size = 8),
                legend.text = element_text(size = 8)) +
          xlab("Simulatd strength of batch effect") +
          ggtitle(i.Variable)
        if(i >= 2)
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
             width = 50, height = 30, limitsize = FALSE)
    return(plist)
  }

ps <- suppFig_simulation_adjBatch(tb_sim = tb_sim, 
                                  results_before = results_before,
                                  results_correct = results_correct,
                                  results_ComBat = results_ComBat)

