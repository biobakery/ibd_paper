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
    
    tb_long <- otu_table2(phylo) %>% 
      t() %>% 
      data.frame(check.names = FALSE) %>% 
      tibble::as_tibble(rownames = "rownames") %>% 
      tidyr::gather(key = feature, 
                    value = abundance,
                    dplyr::one_of(phyloseq::taxa_names(phylo))) %>% 
      dplyr::left_join(tb_study, by = "rownames")
    tb_presence <- tb_long %>% 
      dplyr::group_by(feature, study) %>% 
      dplyr::summarise(present = any(abundance > 0))
    tb_presence <- tb_presence %>% 
      group_by(feature, )
      
    tmp2 <- df_all %>% 
      group_by(genus, dataset_name) %>% 
      summarise(present = any(abundance > 0)) %>% 
      group_by(genus) %>% 
      summarise(npresent = sum(present))
    df_all <- df_all %>% left_join(tmp1) %>% left_join(tmp2)
    tmp1 <- df_all %>% 
      group_by(genus, dataset_name) %>% 
      summarise(notpresent = all(abundance == 0))
    tmp2 <- df_all %>% 
      group_by(genus, dataset_name) %>% 
      summarise(notpresent = all(abundance == 0)) %>% 
      group_by(genus) %>% 
      summarise(nnotpresent = sum(notpresent),
                datasetnotpresent = dataset_name[notpresent][1])
    df_all <- df_all %>% 
      left_join(tmp1) %>% 
      left_join(tmp2)
    gplot <- df_all %>% 
      filter(nnotpresent == 1) %>% 
      group_by(genus) %>% 
      summarise(abundance = mean(abundance),
                datasetnotpresent = datasetnotpresent[1]) %>% 
      ungroup() %>% 
      arrange(datasetnotpresent, abundance) %>% 
      mutate(genus = factor(genus, levels = genus)) %>% 
      ggplot(aes(x = genus, y = abundance, color = datasetnotpresent)) +
      geom_point() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      coord_flip()
    ggsave(gplot, file = "data/phyloseq/fig3.pdf", height = 10, width = 15)
    gplot <- df_all %>% 
      filter(npresent == 1, present == TRUE) %>% 
      group_by(genus, dataset_name) %>% 
      summarise(abundance = mean(abundance)) %>% 
      ungroup() %>% 
      filter(abundance > 2.5e-5) %>% 
      arrange(dataset_name, abundance) %>% 
      mutate(genus = factor(genus, levels = genus)) %>% 
      ggplot(aes(x = genus, y = abundance, color = dataset_name)) +
      geom_point() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      coord_flip()
    ggsave(gplot, file = "data/phyloseq/fig2.pdf", height = 10, width = 15)
    tmp2 %>% filter(npresent == 1) %>% 
      `$`(genus)
    
    
    
    for(study in studies) {
      phylo_tmp <- l_phylo[[study]]
      otu_tmp <- otu_table(phylo_tmp)@.Data
      meta_tmp <- 
        tax_tmp <- tax_table(l_phylo[[study]])@.Data
      otu_new <- matrix(NA, nrow = nrow(mat_tax_all), ncol = nsamples())
    }
  }