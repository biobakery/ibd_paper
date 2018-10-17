# This function is used to check which taxa are uniquely present/absent for a study
unique_taxa <- 
  function(phylo, # an aggregated phyloseq objects of multiple studies
           name_study # name of the study variable in sample_data(phylo)
  ) {
    if(!(name_study %in% colnames(sample_data2(phylo))))
      stop("name_study is not one of the metadata variables available in phylo!")
    tb_study <- tibble::tibble(
      rownames = rownames(sample_data2(phylo)),
      study = sample_data2(phylo)[, name_study]
    )
    if(any(is.na(tb_study$rownames) & is.na(tb_study$study)))
      stop("Sample/study indicator cannot be NA!")
    
    phylo_ra <- phylo %>% 
      phyloseq::transform_sample_counts(tss)
    tb_long <- otu_table2(phylo_ra) %>% 
      t() %>% 
      data.frame(check.names = FALSE) %>% 
      tibble::as_tibble(rownames = "rownames") %>% 
      tidyr::gather(key = feature, 
                    value = abundance,
                    dplyr::one_of(phyloseq::taxa_names(phylo))) %>% 
      dplyr::left_join(tb_study, by = "rownames")
    tb_presence <- tb_long %>% 
      dplyr::group_by(feature, study) %>% 
      dplyr::summarise(present = any(abundance > 0)) %>% 
      dplyr::group_by(feature) %>% 
      dplyr::mutate(n_present = sum(present),
                    n_absent = sum(!present))
    tb_long <- tb_long %>% 
      dplyr::left_join(tb_presence,
                       by = c("feature", "study"))
    
    tb_uniqPres <- tb_long %>% 
      dplyr::filter(n_present == 1, present) %>% 
      dplyr::group_by(feature, study) %>% 
      dplyr::summarise(mean_abundance = mean(abundance)) %>% 
      dplyr::ungroup()
    tb_uniqAbs <-tb_long %>% 
      dplyr::filter(n_absent == 1) %>% 
      dplyr::group_by(feature) %>% 
      dplyr::summarise(mean_abundance = mean(abundance)) %>% 
      dplyr::ungroup() %>% 
      dplyr::left_join(
        tb_presence %>% 
          dplyr::filter(n_absent == 1,
                        !present),
        by = "feature"
      ) %>% 
      dplyr::select(feature, study, mean_abundance)
    
    return(list(tb_presence = tb_presence,
                tb_uniqPres = tb_uniqPres,
                tb_uniqAbs = tb_uniqAbs))
  }