path <- "supp_materials/suppFigures/"
tb.validation.results <- readr::read_tsv("results/6-unsupervised_continuous/genera/validation.tsv")
tb.validation.summary.sd <- tb.validation.results %>% 
  dplyr::filter(class_data == "Bootstrap") %>% 
  dplyr::group_by(batch, Cluster) %>% 
  dplyr::summarise(sd = sd(InnerProd))
tb.validation.summary <- tb.validation.results %>% 
  dplyr::filter(class_data == "Original") %>% 
  dplyr::left_join(tb.validation.summary.sd, by = c("Cluster", "batch"))
tb_plot <- tb.validation.summary %>% 
  dplyr::left_join(tb_1, by = c("batch" = "study_full_name")) %>% 
  dplyr::mutate(Training = ifelse(batch %in% c("CS-PRISM", "RISK", "Pouchitis", "PROTECT"),
                                  "Training",
                                  "Validation") %>% 
                  factor(levels = c("Training", "Validation")),
                Study_plot = dplyr::case_when(
                  !is.na(Study) ~ Study,
                  batch == "nonIBD" ~ "Control (non-IBD)",
                  batch == "HC" ~ "Control (healthy)",
                  batch == "Neg. control" ~ "Neg. control"
                ),
                control_levels = factor(batch, levels = c("nonIBD", "HC", "Neg. control"))) %>% 
  dplyr::arrange(Training, !is.na(control_levels), control_levels, batch) %>% 
  dplyr::mutate(Study_plot = factor(Study_plot, levels = unique(Study_plot)),
                Cluster = Cluster %>% 
                  dplyr::recode_factor("Cluster 1" = "Score (dysbiosis)",
                                       "Cluster 2" = "Score (phyla tradeoff)"))
tb_anno <- tibble::tibble(x = 0.5,
                              y = 0,
                              label = "Training data",
                              Cluster = "Score (dysbiosis)")
p <- tb_plot %>% 
  ggplot(aes(x = Study_plot, y = InnerProd)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = InnerProd - sd, ymax = InnerProd + sd), width = 0.5) +
  facet_grid(.~Cluster) +
  smar::rotate_xaxis(angle = 30) +
  geom_hline(yintercept = 0.65, linetype = "dashed", color = "grey") +
  ylab("Abs. Cosine Coef.") +
  theme(axis.text.x = element_text(color = c(rep("red", 4), rep("black", 9))),
        axis.title.x = element_blank()) +
  geom_text(data = tb_anno,
            aes(x = x, y = y, label = label),
            color = "red",
            hjust = 0,
            size = 6)
ggsave(paste0(path, "suppFig_continuousValidation.pdf"),
       p, width = 10, height = 5)  
