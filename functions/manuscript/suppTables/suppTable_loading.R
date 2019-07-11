tb_loading <- readr::read_tsv("results/6-unsupervised_continuous/genera/cor_cutoff_0.65/avg_loading.tsv")
suppTable_loading <- function(tb_loading,
                              path = "supp_materials/suppTables/suppTable_loading.csv") {
  tb_loading %>% 
    dplyr::select(Feature = feature,
                  `Loading (dysbiosis)` = Cluster_1,
                  `Loading (phyla trade-off)` = Cluster_2) %>% 
    dplyr::arrange(`Loading (dysbiosis)`) %>% 
    readr::write_csv(path) %>% 
    return()
}
tb_loading <- suppTable_loading(tb_loading)
