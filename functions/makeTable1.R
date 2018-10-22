makeTable1 <- function(df_metadata) {
  # if(anyDuplicated(df_metadata$dataset_name))
  #   stop("Study must be unique!")
  # 
  tibble::tibble(
    n_samples = dplyr::n_distinct(df_metadata$sample_accession_16S),
    n_subject = dplyr::n_distinct(df_metadata$subject_accession),
    repeated = df_metadata %>% 
      dplyr::group_by(subject_accession, body_site) %>% 
      dplyr::summarise(n = n()) %>% extract("n") %>% 
      is_greater_than(1) %>% any()
  )
}
