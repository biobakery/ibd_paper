#' Calculate feature fold change between conditions, given a phyloseq object and 
#' test parameter specifications
#'
#' @param phylo phyloseq object with feature and metadata table to calculate fold change with
#' @param test_param list of test parameter specifications. Must have components 
#' \code{test_variable} and \code{contrasts}.
#'
#' @return a one feature per row table, specifying the features' fold change between conditions
calculateFoldChange <- function(phylo, test_param) {
  tb_phylo <- phylo %>%
    phyloseq::transform_sample_counts(MMUPHin:::tss) %>% 
    phyloseq_to_tb %>% 
    dplyr::filter(!!test_param$exprs_filter) %>% 
    dplyr::mutate(test_variable = !!sym(test_param$test_variable) %>% 
                    dplyr::recode(!!!test_param$contrasts,
                                  .missing = NA_character_))
  tb_foldChange <- tb_phylo %>% 
    dplyr::group_by_at(c("feature", "test_variable", colnames(tax_table2(phylo)))) %>% 
    dplyr::summarise(mean_abd = mean(abundance)) %>% 
    tidyr::spread(key = test_variable,
                  value = mean_abd) %>% 
    dplyr::mutate(fold_change = `1` / `0`)
  return(tb_foldChange)
}
