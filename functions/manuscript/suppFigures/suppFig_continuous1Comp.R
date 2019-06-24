load("data/physeq/genera_adj.RData")
tb_loading <- readr::read_tsv("results/6-unsupervised_continuous/genera/cor_cutoff_0.65/avg_loading.tsv")
path <- "supp_materials/suppFigures/"

mat_scores <- physeq_genera_adj %>% 
  to_relativeAbundance() %>% 
  smar::otu_table2() %>% 
  sqrt() %>% asin() %>% t() %>% 
  multiply_by_matrix(as.matrix(tb_loading[, c("Cluster_1", "Cluster_2")])) %>% 
  set_colnames(c("Score (dysbiosis)", "Score (phyla tradeoff)"))
metadata <- smar::sample_data2(physeq_genera_adj) %>% 
  cbind(mat_scores) %>% 
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
  dplyr::filter(q < 0.05) %>% 
  dplyr::mutate(x = c(1.5, 2, 2.5, 2.5, 3),
                y = c(1.02, 1.42, 1.62, 1.02, 1.22),
                q_print = format(q, digits = 2, scientific = TRUE))
tb_segments <- list(tibble::tibble(xstart = 2.05, # UC vs. nonIBD
                                   xend = 2.95,
                                   y = 1),
                    tibble::tibble(xstart = 2.05, # UC vs. HC
                                   xend = 3.95,
                                   y = 1.2),
                    tibble::tibble(xstart = 1.05, # CD vs. UC
                                   xend = 1.95,
                                   y = 1),
                    tibble::tibble(xstart = 1.05, # CD vs. nonIBD
                                   xend = 2.95,
                                   y = 1.4),
                    tibble::tibble(xstart = 1.05, # CD vs. control
                                   xend = 3.95,
                                   y = 1.6)
                    ) %>% 
  purrr::reduce(rbind)
tb_segments_vert <- list(tibble::tibble(x = c(2.05, 2.95), # UC vs. nonIBD
                                        ystart = rep(1, 2),
                                        yend = rep(0.95, 2)),
                         tibble::tibble(x = c(2.05, 3.95), # UC vs. HC
                                        ystart = rep(1.2, 2),
                                        yend = rep(1.15, 2)),
                         tibble::tibble(x = c(1.05, 1.95), # CD vs. UC
                                        ystart = rep(1, 2),
                                        yend = rep(0.95, 2)),
                         tibble::tibble(x = c(1.05, 2.95), # CD vs. nonIBD
                                        ystart = rep(1.4, 2),
                                        yend = rep(1.35, 2)),
                         tibble::tibble(x = c(1.05, 3.95), # CD vs. control
                                        ystart = rep(1.6, 2),
                                        yend = rep(1.55, 2))
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
ggsave(paste0(path, "suppFig_continuous1Comp.pdf"),
       p, 
       width = 7.5, height = 6)
