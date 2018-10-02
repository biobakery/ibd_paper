.rs.restartR()
library(magrittr)
source("scripts/misc/constants.R")

# load the biom files, some initial QC
l_biom <- list()
for(study in studies) {
  print(study)
  biom <- paste0(dir_processed, 
                 "processed/", 
                 study, 
                 "/16s/all_samples_taxonomy_closed_reference.biom") %>% 
    phyloseq::import_biom()
  
  # all otu clusters should have at least some presence in the data,
  # based on the pipeline
  if(any(phyloseq::taxa_sums(biom) == 0)) 
    stop("There are zero sum taxa in the data!")
  l_biom[[study]] <- biom
  
  # 
}

template <- paste0(dir_processed, "data/template.csv") %>% 
  readr::read_csv()
column_setup <- template$var.class %>% 
  dplyr::recode("character" = "c",
                "numeric" = "d",
                "integer" = "i") %>% 
  paste(collapse = "")

for(study in studies) {
  print(study)
  biom <- paste0(dir_processed, 
                 "processed/", 
                 study, 
                 "/16s/all_samples_taxonomy_closed_reference.biom") %>% 
    phyloseq::import_biom()
  # some quality checks on the processed data
  if(any(phyloseq::taxa_sums(biom) == 0)) 
    stop("There are zero sum taxa in the data!")
  # metadata <- paste0(dir_processed, 
  #                    "processed/", 
  #                    study, 
  #                    "/metadata/metadata.txt") %>% 
  #   readr::read_tsv(col_types = column_setup)
  # biom <- paste0(dir_processed, 
  #                "processed/", 
  #                study, 
  #                "/16s/all_samples_taxonomy_closed_reference.biom") %>% 
  #   phyloseq::import_biom()
  # biom %>% 
  #   phyloseq::otu_table() %>% 
  #   apply()
  # sample_common <- intersect(metadata$`16S_sample_accession`,
  #                            sample_names(biom)[apply(otu_table(biom)@.Data > 0, 2, any)]) ## FIXME
  # cat(nrow(metadata),
  #     nsamples(biom),
  #     length(sample_common), "\n")
  # 
  # metadata <- metadata %>% filter(`16S_sample_accession` %in% sample_common) %>% 
  #   as.data.frame()
  # rownames(metadata) <- metadata$`16S_sample_accession`
  # phylo <- phyloseq(
  #   otu_table(biom)[, sample_common],
  #   sample_data(metadata[sample_common, ]),
  #   tax_table(biom)
  # )
  # phylo_genus <- tax_glom(phylo, taxrank = colnames(tax_table(phylo))[6])
  # dir.create(paste0("data/phyloseq/", study, "/"))
  # save(phylo, file = paste0("data/phyloseq/", study, "/", "species.RData"))
  # save(phylo_genus, file = paste0("data/phyloseq/", study, "/", "genus.RData"))
}
