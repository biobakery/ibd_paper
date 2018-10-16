.rs.restartR()
dir.create("results/2-QC/", recursive = TRUE, showWarnings = FALSE)
source("scripts/misc/constants.R")
source("scripts/misc/helpers.R")

l_biom_species <- list()
l_biom_genera <- list()
for(study in studies) {
  load(paste0("data/phyloseq/", study, "/genus.RData"))
  load(paste0("data/phyloseq/", study, "/species.RData"))
  l_biom_species[[study]] <- phylo
  l_biom_genera[[study]] <- phylo_genus
}

l_biom_genera <- 
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

metadata_all <- c("BIDMC-FMT",
                  "Jansson_Lamendella_Crohns",
                  "jansson_twins_ibd",
                  "LSS-PRISM",
                  "MucosalIBD",
                  "Pouchitis",
                  "RISK") %>% 
  purrr::map_dfr(function(study) {
    paste0("processed/", study, "/metadata/metadata.txt") %>% 
      read_tsv %>% 
      mutate(PMID = as.character(PMID),
             BMI = as.numeric(BMI),
             calprotectin = as.numeric(calprotectin),
             age = as.numeric(calprotectin),
             time_point = as.character(time_point),
             subject_accession = as.character(subject_accession),
             family = as.character(family),
             PCDAI = as.numeric(PCDAI)
      )
  })
