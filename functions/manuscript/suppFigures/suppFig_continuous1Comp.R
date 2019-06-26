load("data/physeq/genera_adj.RData")
tb_loading <- readr::read_tsv("results/6-unsupervised_continuous/genera/cor_cutoff_0.65/avg_loading.tsv")
tb_scores <- readr::read_tsv("results/6-unsupervised_continuous/genera/scores.tsv")
path <- "supp_materials/suppFigures/suppFig_continuous1Comp.pdf"


metadata <- smar::sample_data2(physeq_genera_adj) %>% 
  dplyr::left_join(tb_scores, by = c("sample_accession_16S" = "sample")) %>% 
  dplyr::mutate(Disease = dplyr::case_when(
    disease == "CD" ~ "CD",
    disease == "UC" ~ "UC",
    control == "nonIBD" ~ "Control (non-IBD)",
    control == "HC" ~ "Control (healthy)"
  ) %>% 
    factor(levels = c("CD", "UC", "Control (non-IBD)", "Control (healthy)")))

tb_tests <- (1:3) %>% 
  purrr::map_dfr(function(i.lvl) {
    ((i.lvl + 1):4) %>% 
      purrr::map_dfr(function(j.lvl) {
        fit.wilcox <- wilcox.test(metadata$`Score (dysbiosis)`[metadata$Disease == levels(metadata$Disease)[i.lvl]],
                                  metadata$`Score (dysbiosis)`[metadata$Disease == levels(metadata$Disease)[j.lvl]],
                                  alternative = "less") 
        tibble::tibble(i.lvl = levels(metadata$Disease)[i.lvl],
                       j.lvl = levels(metadata$Disease)[j.lvl],
                       p = fit.wilcox$p.value)
      }
      )
  }) %>% 
  dplyr::mutate(q = p.adjust(p, "fdr"))

# q value annotations
tb_q <- tb_tests %>% 
  # dplyr::filter(q < 0.05) %>% 
  dplyr::mutate(x = c(1.5, 2, 2.5, 2.5, 3, 3.5),
                y = c(0.72, 1.12, 1.32, 0.72, 0.92, 0.72),
                q_print = format(q, digits = 2, scientific = TRUE))
tb_segments <- list(tibble::tibble(xstart = 2.05, # UC vs. nonIBD
                                   xend = 2.95,
                                   y = 0.7),
                    tibble::tibble(xstart = 2.05, # UC vs. HC
                                   xend = 3.95,
                                   y = 0.9),
                    tibble::tibble(xstart = 1.05, # CD vs. UC
                                   xend = 1.95,
                                   y = 0.7),
                    tibble::tibble(xstart = 1.05, # CD vs. nonIBD
                                   xend = 2.95,
                                   y = 1.1),
                    tibble::tibble(xstart = 1.05, # CD vs. HC
                                   xend = 3.95,
                                   y = 1.3),
                    tibble::tibble(xstart = 3.05, # nonIBD vs. HC
                                   xend = 3.95,
                                   y = 0.7)
                    ) %>% 
  purrr::reduce(rbind)
tb_segments_vert <- list(tibble::tibble(x = c(2.05, 2.95), # UC vs. nonIBD
                                        ystart = rep(0.7, 2),
                                        yend = rep(0.65, 2)),
                         tibble::tibble(x = c(2.05, 3.95), # UC vs. HC
                                        ystart = rep(0.9, 2),
                                        yend = rep(0.85, 2)),
                         tibble::tibble(x = c(1.05, 1.95), # CD vs. UC
                                        ystart = rep(0.7, 2),
                                        yend = rep(0.65, 2)),
                         tibble::tibble(x = c(1.05, 2.95), # CD vs. nonIBD
                                        ystart = rep(1.1, 2),
                                        yend = rep(1.05, 2)),
                         tibble::tibble(x = c(1.05, 3.95), # CD vs. HC
                                        ystart = rep(1.3, 2),
                                        yend = rep(1.25, 2)),
                         tibble::tibble(x = c(3.05, 3.95), # nonIBD vs. HC
                                        ystart = rep(0.7, 2),
                                        yend = rep(0.65, 2))
) %>% 
  purrr::reduce(rbind)


p <- metadata %>% ggplot(aes(x = Disease, y = `Score (dysbiosis)`)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.25)) +
  # scale_y_continuous(limits = c(-0.73, 1.8)) +
  geom_segment(data = tb_segments, aes(x = xstart, xend = xend,
                                       y = y, yend = y)) +
  geom_segment(data = tb_segments_vert, aes(x = x, xend = x,
                                            y = ystart, yend = yend)) +
  geom_text(data = tb_q, aes(x = x, y = y, label = q_print),
            hjust = 0.5, vjust = 0, size = 5) +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  ylab("Score (dysbiosis)") +
  theme(axis.title.x = element_blank())
ggsave(path,
       p, 
       width = 7.5, height = 6)
