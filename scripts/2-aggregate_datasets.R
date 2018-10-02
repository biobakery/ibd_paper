l_phylo <- list()
for(study in studies) {
  print(study)
  load(paste0("data/phyloseq/", study, "/", "genus.RData"))
  l_phylo[[study]] <- phylo_genus
  tax_table(l_phylo[[study]])@.Data[1:5, 6] %>% cat("\n")
}
tax_all <- l_phylo %>% 
  map(function(phylo) tax_table(phylo)@.Data[, 1:6] %>% 
        apply(1, paste, collapse = "|")) %>% 
  unlist() %>% 
  unique()
mat_tax_all <- tax_all %>% 
  sapply(strsplit, split = "|", fixed = TRUE) %>% 
  data.frame(check.names = FALSE) %>% 
  as.matrix() %>% 
  t()
colnames(mat_tax_all) <- paste0("Rank", 1:6)
mat_otu_all <- studies %>% 
  map_dfc(function(study) {
    otu_tmp <- otu_table(l_phylo[[study]])@.Data
    tax_tmp <- tax_table(l_phylo[[study]])@.Data
    rownames(otu_tmp) <- tax_tmp[, 1:6] %>% apply(1, paste, collapse = "|")
    otu_new <- matrix(0, nrow = nrow(mat_tax_all), ncol = ncol(otu_tmp))
    dimnames(otu_new) <- list(rownames(mat_tax_all), colnames(otu_tmp))
    otu_new[rownames(otu_tmp), ] <- otu_tmp
    return(as.data.frame(otu_new, check.names = FALSE))
  }) %>% 
  as.matrix
rownames(mat_otu_all) <- rownames(mat_tax_all)
df_meta_all <- studies %>% 
  map_dfr(function(study) {
    df <- data.frame(sample_data(l_phylo[[study]]),
                     check.names = FALSE) %>% 
      mutate(PMID = as.character(PMID),
             BMI = as.character(BMI),
             calprotectin = as.character(calprotectin),
             time_point = as.character(time_point),
             subject_accession = as.character(subject_accession),
             family = as.character(family),
             PCDAI = as.character(PCDAI))
  }) ## FIXME
rownames(df_meta_all) <- colnames(mat_otu_all)
phylo_all <- phyloseq(
  otu_table(mat_otu_all, taxa_are_rows = TRUE),
  sample_data(df_meta_all),
  tax_table(mat_tax_all)
)