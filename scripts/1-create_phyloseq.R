rm(list = ls())
library(tidyverse)
library(phyloseq)
studies <- c("BIDMC-FMT",
             "CS-PRISM",
             "Jansson_Lamendella_Crohns",
             "jansson_twins_ibd",
             "LSS-PRISM",
             "MucosalIBD",
             "Pouchitis",
             "RISK",
             "Herfarth_CCFA_Microbiome_3B_combined")

for(study in studies) {
  print(study)
  metadata <- paste0("processed/", study, "/metadata/metadata.txt") %>% 
    read_tsv
  biom <- paste0("processed/", study, "/16s/all_samples_taxonomy_closed_reference.biom") %>% 
    import_biom()
  sample_common <- intersect(metadata$`16S_sample_accession`,
                             sample_names(biom)[apply(otu_table(biom)@.Data > 0, 2, any)]) ## FIXME
  cat(nrow(metadata),
      nsamples(biom),
      length(sample_common), "\n")
  
  metadata <- metadata %>% filter(`16S_sample_accession` %in% sample_common) %>% 
    as.data.frame()
  rownames(metadata) <- metadata$`16S_sample_accession`
  phylo <- phyloseq(
    otu_table(biom)[, sample_common],
    sample_data(metadata[sample_common, ]),
    tax_table(biom)
  )
  phylo_genus <- tax_glom(phylo, taxrank = colnames(tax_table(phylo))[6])
  dir.create(paste0("data/phyloseq/", study, "/"))
  save(phylo, file = paste0("data/phyloseq/", study, "/", "species.RData"))
  save(phylo_genus, file = paste0("data/phyloseq/", study, "/", "genus.RData"))
}
