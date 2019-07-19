load("results/simulations/lm.meta/tb_sim.RData")
load("results/simulations/lm.meta/evaluate_before.RData")
load("results/simulations/lm.meta/evaluate_after.RData")
load("results/simulations/lm.meta/evaluate_qnorm.RData")

load("results/simulations/lm.meta/time_MMUPHin.RData")
load("results/simulations/lm.meta/time_BDMMA.RData")
times_qnorm <- readr::read_tsv("results/simulations/lm.meta/qnorm/time_qnorm.txt", col_names = FALSE) %>% 
  dplyr::transmute(`Computation time (min)` = X1 / 60,
                   Class = "Quantile normalization") %>% 
  dplyr::slice(-1)
tb_time <- rbind(
  times_qnorm,
  times_BDMMA %>% 
    purrr::map_dbl(~ as.double(.x, units = "mins")) %>% 
    {tibble::tibble(`Computation time (min)` = .,
                    Class = "BDMMA")},
  times_MMUPHin %>% 
    purrr::map_dbl(~ as.double(.x, units = "mins")) %>% 
    {tibble::tibble(`Computation time (min)` = .,
                    Class = "MMUPHin")}
) %>% 
  dplyr::mutate(Class = Class %>% forcats::as_factor())

suppFig_simulation_lm <- function(tb_sim, 
                                  evaluate_before, evaluate_after, 
                                  evaluate_qnorm, tb_time,
                                  path = "supp_materials/suppFigures/suppFig_simulation_lm.pdf") {
  # Format results
  # effect sizes
  tb_sim <- tb_sim %>% 
    dplyr::mutate(effect_exposure = imbalance %>% 
                    forcats::as_factor() %>% forcats::fct_inseq(),
                  effect_batch = effectSize %>%
                    purrr::map_dbl("batch") %>% 
                    forcats::as_factor() %>% forcats::fct_inseq())
  # nSample and nBatch
  tb_sim <- tb_sim %>% 
    dplyr::arrange(nSample, nBatch) %>% 
    dplyr::mutate(nSample_nBatch = paste0("nSample per batch = ", nSample_perBatch, ", ",
                                          "nBatch = ", nBatch) %>% 
                    forcats::as_factor())
  # nMicrobe and percMicrobe spiked
  tb_sim <- tb_sim %>% 
    dplyr::arrange(nMicrobe, spikeMicrobes) %>% 
    dplyr::mutate(nMicrobe_percSpike = paste0("nFeature = ", nMicrobe, ", ",
                                              spikeMicrobes*100, "% spiked") %>% 
                    forcats::as_factor())
  
  # aggregate results into df format
  tb_results <- list(results, results_correct) %>% 
    purrr::map2_dfr(
      c("Original", "MMUPHin corrected"),
      ~ .x %>% 
        purrr::imap_dfr(
          function(result, i_result) {
            if(.y == "Original")
              return(tibble::tibble(model = "naive",
                                    results = result["naive"],
                                    i = i_result))
            if(.y == "MMUPHin corrected")
              return(tibble::tibble(model = "MMUPHin",
                                    results = result["MMUPHin"],
                                    i = i_result))
          }) %>% 
        dplyr::mutate(Data = .y)
    ) 
  tb_results_qnorm <- results_qnorm %>% 
    purrr::imap_dfr(
      ~ tibble::tibble(results = list(.x),
                       i = .y)
    ) %>% 
    dplyr::mutate(Data = "Quantile corrected",
                  model = "Wilcoxon")
  # format variables
  tb_results <- rbind(tb_results, tb_results_qnorm) %>% 
    dplyr::mutate(Data = factor(Data, levels = c("Original", 
                                                 "Quantile corrected",
                                                 "MMUPHin corrected")),
                  model = factor(model, levels = c("naive", "Wilcoxon", "MMUPHin"))) %>% 
    dplyr::arrange(Data, model) %>% 
    dplyr::mutate(Class = paste0(Data, ", ", model, " analyzed") %>% 
                    dplyr::recode("Original, naive analyzed" = "Naive model",
                                  "Quantile corrected, Wilcoxon analyzed" = "Quantile normalization",
                                  "BDMMA" = "BDMMA",
                                  "MMUPHin corrected, MMUPHin analyzed" = "MMUPHin") %>% 
                    forcats::as_factor())
  
  # calculate FPR
  tb_results <- tb_results %>% 
    dplyr::left_join(tb_sim, by = "i") %>% 
    dplyr::mutate(FPR = results %>% 
                    purrr::map_dbl(~sum(.x$pval < 0.05, na.rm = TRUE)) %>% 
                    magrittr::divide_by(nMicrobe))
  # summarise mean and sd FPR across replicates
  tb_summary <- tb_results %>%
    dplyr::group_by(i_setup, Class) %>% 
    dplyr::summarise(mean_FPR = mean(FPR),
                     sd_FPR = sd(FPR)) %>% 
    dplyr::ungroup()
  
  # Set up colors
  colors <- smar::gg_color_hue(n=4)
  colors_simulation_lm.meta <- c("Naive model" = "black",
                                 "Quantile normalization" = colors[3],
                                 "BDMMA" = colors[4],
                                 "MMUPHin" = colors[1])
  
  tb_plot <- tb_summary %>% 
    dplyr::left_join(dplyr::filter(tb_sim, !duplicated(i_setup)), 
                     by = "i_setup")
  plist <- (1:nlevels(tb_sim$effect_exposure)) %>% 
    purrr::map(function(i) {
      i.effect <- levels(tb_sim$effect_exposure)[i]
      p <- tb_plot %>% 
        dplyr::filter(effect_exposure == i.effect) %>% 
        ggplot(aes(x = effect_batch, y = mean_FPR, color = Class)) +
        geom_point(size = 3, position = position_dodge(width = 1)) +
        geom_errorbar(aes(ymin = mean_FPR - sd_FPR, ymax = mean_FPR + sd_FPR),
                      position = position_dodge(width = 1)) +
        geom_hline(yintercept = 0.05, linetype = "dashed") +
        facet_grid(nMicrobe_percSpike ~ nSample_nBatch) +
        scale_color_manual(values = colors_simulation_lm.meta) +
        theme(legend.position = c(1, 1),
              legend.justification = c(1, 1),
              legend.background = element_blank(),
              legend.title = element_blank(),
              strip.text = element_text(size = 8),
              legend.text = element_text(size = 8)) +
        xlab("Batch effect") +
        ylab("False positive rate") +
        ggtitle(paste0("Exposure imbalance = ", 
                       i.effect %>% 
                         as.numeric %>% 
                         magrittr::multiply_by(200),
                       "%"))
      if(i != 3)
        p <- p + theme(legend.position = "none")
      return(p)
    })
  
  p_time <- tb_time %>% 
    ggplot(aes(x = Class, y = `Computation time (min)`)) +
    geom_boxplot() +
    ylab("Computation time (min)") +
    theme(axis.title.x = element_blank()) +
    ggtitle("Computation performance comparison")
  
  cowplot::plot_grid(plotlist = c(plist,
                                  cowplot::plot_grid(p_time, NULL, nrow = 1,
                                                     rel_widths = c(2, 1)) %>% 
                                    cowplot::plot_grid(NULL, ncol = 1,
                                                       rel_heights = c(1, 3)) %>% 
                                    list()),
                     nrow = 2) %>% 
    ggsave(path,
           .,
           width = 30, height = 30)
}


