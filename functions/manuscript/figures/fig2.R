rm(list = ls())
smar::sourceDir("functions/")
library(magrittr)
library(ggplot2)

# batch adjustment --------------------------------------------------------
dir_data <- "results/simulations/adjust.batch/"
load(paste0(dir_data, "tb_sim.RData"))
load(paste0(dir_data, "evaluate_before.RData"))
load(paste0(dir_data, "evaluate_after.RData"))
load(paste0(dir_data, "evaluate_ComBat.RData"))

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

# aggregate results into df format
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

colors <- smar::gg_color_hue(n=4)
colors_simulation_adjust.batch <- c("Original" = "black",
                                    "ComBat corrected" = colors[3],
                                    "MMUPHin corrected" = colors[1])

tb_R2 <- list(results, results_ComBat, results_correct) %>% 
  purrr::map2_dfr(names(colors_simulation_adjust.batch), 
                  ~ .x %>% 
                    purrr::imap_dfr(~ data.frame(R2_positive_control1 = .x$R2["positive_control1"],
                                                 R2_positive_control2 = .x$R2["positive_control2"],
                                                 R2_negative_control = .x$R2["negative_control"],
                                                 R2_batch = .x$R2["batch"],
                                                 i = .y)) %>% 
                    data.frame(Data = .y)) %>% 
  dplyr::left_join(tb_sim, by = "i") %>% 
  dplyr::filter(nMicrobe == 200,
                spikeMicrobes == 0.05,
                nBatch == 4,
                nSample_perBatch == 500) %>% 
  tidyr::gather(key = "Variable",
                value = "R2",
                R2_batch, R2_positive_control1, R2_positive_control2, R2_negative_control) %>% 
  dplyr::mutate(Variable = Variable %>% 
                  dplyr::recode_factor("R2_batch" = "Batch",
                                       "R2_positive_control1" = "Pos. control (binary)",
                                       "R2_positive_control2" = "Pos. control (continuous)",
                                       "R2_negative_control" = "Neg. control"),
                Data = factor(Data, levels = names(colors_simulation_adjust.batch)))

p_batch <- tb_R2 %>% 
  ggplot(aes(x = as.factor(effect_batch),
             y = R2,
             color = Data)) +
  geom_boxplot() +
  scale_color_manual(values = colors_simulation_adjust.batch) +
  facet_grid(.~Variable) +
  xlab("Simulated strength of batch effect") +
  theme(
    # axis.title = element_text(size = 15),
    #     axis.text = element_text(size = 12),
        # strip.text.x = element_text(size = 15),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        # legend.title = element_text(size = 12),
        # legend.text = element_text(size = 12),
        # legend.key.size = unit(15, "points")
        )

# lm.meta
dir_data <- "results/simulations/lm.meta/"

load(paste0(dir_data, "tb_sim.RData"))
load(paste0(dir_data, "evaluate_before.RData"))
load(paste0(dir_data, "evaluate_after.RData"))
load(paste0(dir_data, "evaluate_qnorm.RData"))
load(paste0(dir_data, "evaluate_BDMMA.RData"))

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
tb_sim_subset <- tb_sim %>%
  dplyr::filter(nSample_perBatch == 100,
                nMicrobe == 200,
                spikeMicrobes == 0.05)
tb_results_BDMMA <- results_BDMMA %>% 
  purrr::map2_dfr(tb_sim_subset$i,
                  ~ tibble::tibble(
                    results = list(
                      .x$parameter_summary %>% 
                        tibble::rownames_to_column("name") %>% 
                        dplyr::filter(
                          stringr::str_detect(name,
                                              stringr::fixed("L_"))
                        ) %>%
                        dplyr::mutate(pval = ifelse(mean > 0.5,
                                                    0, 
                                                    1))
                    ),
                    i = .y)) %>% 
  dplyr::mutate(Data = "Original",
                model = "BDMMA")

# format variables
tb_results <- rbind(tb_results, tb_results_qnorm, tb_results_BDMMA) %>% 
  dplyr::mutate(Class = paste0(Data, ", ", model) %>% 
                  dplyr::recode_factor("Original, naive" = "Naive model",
                                       "Quantile corrected, Wilcoxon" = "Quantile normalization",
                                       "Original, BDMMA" = "BDMMA",
                                       "MMUPHin corrected, MMUPHin" = "MMUPHin"))


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

p_FPR <- tb_plot %>% 
  dplyr::filter(effect_exposure == 0.4, nMicrobe == 200, 
                nSample_perBatch == 500, spikeMicrobes == 0.05) %>% 
  ggplot(aes(x = effect_batch,
             y = mean_FPR,
             color = Class)) +
  geom_point(position = position_dodge(width = 1), size = 5) +
  # geom_line(aes(linetype = as.factor(effect_batch))) +
  geom_errorbar(aes(ymin = mean_FPR - sd_FPR,
                    ymax = mean_FPR + sd_FPR),
                width = 0.5,
                position = position_dodge(width = 1)) +
  scale_color_manual(values = colors_simulation_lm.meta, name = "Data and model") +
  geom_hline(yintercept = 0.05, color = "black", linetype = "dashed") +
  xlab("Simulated strength of batch effect") +
  ylab("FPR") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 15),
        legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.background = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.size = unit(15, "points"))


# discrete clusters -------------------------------------------------------
dir_data <- "../mmuphin/simulation/results/discrete/"
load(paste0(dir_data, "tb_sim.RData"))
tb_sim <- tb_sim %>% 
  dplyr::mutate(effect_cluster = effectSize %>% 
                  purrr::map_dbl("cluster"),
                effect_batch = effectSize %>% 
                  purrr::map_dbl("batch"),
                i = 1:n()) %>% 
  dplyr::select(-`1:n()`)
load(paste0(dir_data, "evaluate_before.RData"))
load(paste0(dir_data, "evaluate_after.RData"))

colors_simulation_unsupervised <- c("Original" = "black",
                                    "MMUPHin corrected" = colors[1])

tb_predStr <- list(results, results_correct) %>% 
  purrr::map2_dfr(names(colors_simulation_unsupervised), 
                  ~ .x %>% 
                    purrr::imap_dfr(function(result, i) {
                      nCluster_true <- tb_sim[i, ]$nCluster
                      tb_result <- result$pred.str
                      tb_result %>% 
                        dplyr::arrange(-(k == nCluster_true),
                                       -mean_pred) %>% 
                        dplyr::slice(1:2) %>% 
                        dplyr::arrange(-mean_pred) %>% 
                        dplyr::mutate(pred_str_rank = c("highest", "second"),
                                      i = i)
                    }) %>% 
                    data.frame(Data = .y)) %>% 
  dplyr::mutate(Data = factor(Data, levels = names(colors_simulation_unsupervised)))

p_discrete <- tb_predStr %>% 
  dplyr::filter(pred_str_rank == "highest") %>% 
  dplyr::left_join(tb_sim, by = "i") %>% 
  dplyr::filter(nSample == 1000, nBatch == 2) %>%
  dplyr::filter(!(effect_batch == 0 & Data == "MMUPHin corrected")) %>% 
  dplyr::filter(effect_cluster == 10) %>% 
  dplyr::mutate(success = (k == nCluster)) %>% 
  dplyr::group_by(effect_batch,
                  nCluster, nBatch,
                  nMicrobe, Data) %>% 
  dplyr::summarise(mean_success = mean(success)) %>% 
  ggplot(aes(x = as.factor(effect_batch), y = mean_success*100, 
             color = Data)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  facet_grid(.~paste0("True #clusters = ", nCluster)) +
  scale_color_manual(values = colors_simulation_unsupervised) +
  xlab("Simulated strength of batch effect") +
  ylab("%Success\n
       identifying clusters") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 15),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.size = unit(15, "points"))


# continuous --------------------------------------------------------------
dir_data <- "../mmuphin/simulation/results/continuous/"
load(paste0(dir_data, "tb_sim.RData"))
tb_sim <- tb_sim %>% 
  dplyr::mutate(effect_score = effectSize %>% 
                  purrr::map_dbl("score"),
                effect_batch = effectSize %>% 
                  purrr::map_dbl("batch"),
                effect_score_plot = paste0("effect score=", effect_score) %>% 
                  factor(levels = paste0("effect score=", c(1, 5, 10))),
                i = 1:n()) %>% 
  dplyr::select(-`1:n()`)
load(paste0(dir_data, "evaluate_before.RData"))
load(paste0(dir_data, "evaluate_after.RData"))

colors_simulation_score <- c("violet", colors_simulation_unsupervised[2])
names(colors_simulation_score) <- c("Original, MMUPHin analyzed", "MMUPHin corrected, MMUPHin analyzed")
tb_continuous <- list(results, results_correct) %>% 
  purrr::map2_dfr(c("Original", "MMUPHin corrected"), 
                  ~ .x %>% 
                    purrr::imap_dfr(~ data.frame(cor_cutoff = .x$fit.continuous$cor.cutoff,
                                                 correlation = .x$fit.continuous$correlation,
                                                 i = .y)) %>% 
                    data.frame(Data = .y)) %>% 
  dplyr::mutate(Data = factor(Data %>% 
                                dplyr::recode("Original" = "Original, MMUPHin analyzed",
                                              "MMUPHin corrected" = "MMUPHin corrected, MMUPHin analyzed"), 
                              levels = names(colors_simulation_score)))
p_continuous <-  tb_continuous %>% 
  dplyr::left_join(tb_sim, by = "i") %>% 
  dplyr::filter(nSample == 1000, effect_score == 10, nBatch == 4) %>%
  dplyr::filter(!(effect_batch == 0 & Data == "MMUPHin corrected, MMUPHin analyzed")) %>% 
  ggplot(aes(x = as.factor(effect_batch), y = abs(correlation), 
             color = Data)) +
  geom_boxplot() +
  # facet_grid(paste0("nBatch=", nBatch)~effect_score_plot) +
  scale_color_manual(values = colors_simulation_score) +
  xlab("Simulated strength of batch effect") +
  ylab("Correlation between true and\n identified scores") +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text.x = element_text(size = 15),
        legend.position = c(1, 1),
        legend.justification = c(1, 1),
        legend.background = element_blank(),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.key.size = unit(15, "points"))

ps_summary <- list(adjust.batch = p_batch, 
                   lm.metaFP = p_FPR,
                   discrete = p_discrete,
                   continuous = p_continuous)
save(ps_summary,
     file = "figures/figure2/assets/summary.RData")



rm(list = ls())
load("figures/figure2/assets/toys.RData")
load("figures/figure2/assets/summary.RData")

cowplot::plot_grid(ps_summary[[1]], ps_toy[[1]],
                   ps_summary[[2]], ps_toy[[2]],
                   ps_summary[[3]], ps_toy[[3]],
                   ps_summary[[4]], ps_toy[[4]],
                   ncol = 2,
                   rel_widths = c(4, 1,
                                  4, 1, 
                                  4, 1,
                                  4, 1),
                   labels = c("a", "b", "c", "d", "e", "f", "g", "h"),
                   label_size = 25) %>% 
  ggsave("figures/figure2/figure2.pdf",
         ., height = 24, width = 3*5)
