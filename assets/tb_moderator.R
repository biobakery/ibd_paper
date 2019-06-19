load("data/physeq/genera.RData")
metadata <- physeq_genera %>% smar::sample_data2()
metadata_test <- readxl::read_xlsx("assets/table1_metadata.xlsx") %>% 
  dplyr::right_join(metadata, by = c("study_full_name" = "dataset_name")) %>% 
  dplyr::mutate(sample_type2 = 
                  sample_type %>% 
                  dplyr::recode("biopsy" = "Biopsy",
                                "stool" = "Stool"),
                disease2 = 
                  disease %>% 
                  dplyr::recode("control" = "Control")) %>% 
  dplyr::mutate(study_site = paste(Study,
                                   sample_type2, 
                                   sep = ", "),
                dataset_site = paste(study_full_name,
                                     sample_type, 
                                     sep = "_"),
                study_site_disease = paste(Study,
                                           sample_type2, 
                                           disease2,
                                           sep = ", "),
                dataset_site_disease = paste(study_full_name,
                                             sample_type,
                                             disease,
                                             sep = "_"))
df_moderator_site <- metadata_test %>% 
  dplyr::select(study_site, dataset_site, sample_type, 
                Study, study_full_name) %>% 
  dplyr::group_by(study_site) %>% 
  dplyr::summarise_all(function(x) unique(x)) %>% 
  dplyr::mutate(sample_type = factor(sample_type)) %>% 
  as.data.frame() %>% tibble::column_to_rownames("dataset_site") 
df_moderator_siteDisease <- metadata_test %>% 
  dplyr::select(study_site_disease, dataset_site_disease, sample_type, disease,
                Study, study_full_name) %>% 
  dplyr::filter(disease %in% c("CD", "UC")) %>% 
  dplyr::group_by(study_site_disease) %>% 
  dplyr::summarise_all(function(x) unique(x)) %>% 
  dplyr::mutate(sample_type = factor(sample_type),
                disease = factor(disease)) %>% 
  as.data.frame() %>% tibble::column_to_rownames("dataset_site_disease")