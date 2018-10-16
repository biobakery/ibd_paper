.rs.restartR()
dir.create("results/1-create_phyloseq/", recursive = TRUE, showWarnings = FALSE)
source("scripts/misc/constants.R")
source("scripts/misc/helpers.R")

# Load the biom files, some initial QC ------------------------------------

l_biom <- list()
for(study in studies) {
  print(study)
  biom <- paste0(dir_processed, 
                 "processed/", 
                 study, 
                 "/16s/all_samples_taxonomy_closed_reference.biom") %>% 
    phyloseq::import_biom()
  if(any(phyloseq::taxa_sums(biom) == 0)) 
    stop("There are zero sum taxa in the data!") # all otu clusters should have at least some presence 
                                                 # in the data, based on the pipeline
  l_biom[[study]] <- biom
}

# Some lib size sanity check
lib_size <- studies %>% 
  purrr::map_dfr(~ tibble::tibble(
    study = .x,
    sample = l_biom %>% 
      extract2(.x) %>% 
      phyloseq::sample_names(),
    lib.size = l_biom %>% 
      extract2(.x) %>% 
      phyloseq::sample_sums()
  ))  
read_counts <- studies %>% 
  purrr::map_dfr(~ paste0(dir_processed, 
                          "processed/", 
                          .x, 
                          "/16s/all_samples_read_counts.tsv") %>% 
                   readr::read_tsv(col_types = "cddd") %>% 
                   dplyr::mutate(study = .x)) 
lib_size %>% 
  dplyr::filter(sample %>% 
                  is_in(read_counts$`# sample`) %>% 
                  not)
read_counts %>%
  dplyr::filter(`# sample` %>% 
                  is_in(lib_size$sample) %>% 
                  not) # samples that aren't present in the OTU table 
                       # are because of low read counts
lib_size <- lib_size %>% 
  dplyr::left_join(read_counts,
                   c("study" = "study",
                     "sample" = "# sample"))
lib_size %>% 
  ggplot(aes(x = lib.size, y = `reads mapping to OTU with taxonomy`)) +
  geom_point() +
  facet_wrap(~study) # the two are equal
lib_size %>% 
  tidyr::gather(key = "Type of read count",
                value = "N",
                "original read count",
                "lib.size",
                `reads mapping to OTU with taxonomy`,
                `reads mapping to unclassifed OTU`,
                factor_key = TRUE) %>% 
  ggplot(aes(x = study, y = log10(N + 1), color = `Type of read count`)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = "results/1-create_phyloseq/fig_libSize.pdf",
       width = 10,
       height = 8)

# Create phyloseq objects -------------------------------------------------

template <- paste0(dir_processed, "data/template.csv") %>% 
  readr::read_csv()
column_setup <- template$var.class %>% 
  dplyr::recode("character" = "c",
                "numeric" = "d",
                "integer" = "i") %>% 
  paste(collapse = "")

# for(study in studies) {
#   print(study)
#   biom <- paste0(dir_processed, 
#                  "processed/", 
#                  study, 
#                  "/16s/all_samples_taxonomy_closed_reference.biom") %>% 
#     phyloseq::import_biom()
#   # some quality checks on the processed data
#   if(any(phyloseq::taxa_sums(biom) == 0)) 
#     stop("There are zero sum taxa in the data!")
#   # metadata <- paste0(dir_processed, 
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