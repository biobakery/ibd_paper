library(MMUPHin)
mat_otu_adj <- MMUPHin::adjust.batch(mat_otu_all, batch = "dataset_name", 
                                     covariates = c("disease", "body_site"),
                                     data = df_meta_all)
phylo_all_adj <- phylo_all
