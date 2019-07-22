source("assets/tb_1.R")
load("data/physeq/genera_prefilter.RData")
df_metadata <- smar::sample_data2(physeq_genera_prefilter)


tb_table1 <- df_metadata %>% 
  dplyr::group_by(dataset_name) %>% 
  dplyr::summarise(
    n_subject = dplyr::n_distinct(subject_accession, na.rm = TRUE),
    n_subject_baseline = dplyr::n_distinct(subject_accession[is_baseline], na.rm = TRUE),
    n_sample = dplyr::n_distinct(sample_accession_16S, na.rm = TRUE),
    n_sample_baseline = dplyr::n_distinct(sample_accession_16S[is_baseline], na.rm = TRUE),
    n_CD = dplyr::n_distinct(subject_accession[disease == "CD"], na.rm = TRUE),
    n_UC = dplyr::n_distinct(subject_accession[disease == "UC"], na.rm = TRUE),
    n_control = dplyr::n_distinct(subject_accession[disease == "control"], na.rm = TRUE),
    n_male = dplyr::n_distinct(subject_accession[gender == "m"], na.rm = TRUE),
    n_female = dplyr::n_distinct(subject_accession[gender == "f"], na.rm = TRUE),
    n_missing = dplyr::n_distinct(subject_accession[is.na(gender)], na.rm = TRUE),
    n_biopsy = sum(sample_type == "biopsy", na.rm = TRUE),
    n_stool = sum(sample_type == "stool", na.rm = TRUE)
  )
tb_table1_age <- tb_1 %>% 
  dplyr::right_join(df_metadata, by = c("study_full_name" = "dataset_name")) %>% 
  dplyr::group_by(study_full_name) %>% 
  dplyr::mutate(is_baseline_or_unique = ifelse(`Longitudinal sampling?` == "No", 
                                               !duplicated(subject_accession), 
                                               is_baseline)) %>% 
  dplyr::filter(is_baseline_or_unique) %>% 
  dplyr::summarise(age_mean = mean(age, na.rm = TRUE),
                   age_sd = sd(age, na.rm = TRUE))
tb_table1 <- tb_table1 %>% dplyr::left_join(tb_table1_age, by = c("dataset_name" = "study_full_name"))
tb_table1 <- tb_1 %>% 
  dplyr::right_join(tb_table1, by = c("study_full_name" = "dataset_name"))

tb_table1_toshow <- tb_table1 %>% 
  dplyr::group_by(1:dplyr::n()) %>% 
  dplyr::transmute(Study = paste0(Study, "\n", "{", PMID, "}"),
                   `Brief description` = `Description of study design`,
                   `N subject` = n_subject,
                   `N sample` = ifelse(`Longitudinal sampling?` == "Yes",
                                       paste0(" (", n_sample_baseline, ")"),
                                       "") %>% paste0(n_sample, .),
                   `Phenotype(s)` = c("CD", "UC", "Control") %>% 
                     purrr::map2_chr(c(n_CD, n_UC, n_control), function(disease, n) {
                       if(n == 0) return("none")
                       return(paste0(disease, " ", n))
                     }) %>% setdiff("none") %>% paste0(collapse = "/\n"),
                   Age = ifelse(is.na(age_mean),
                                "",
                                paste0(round(age_mean, 2), "\n(", round(age_sd, 2), ")")),
                   Gender = c("Male", "Female", "Missing") %>% 
                     purrr::map2_chr(c(n_male, n_female, n_missing), function(gender, n) {
                       if(n == 0) return("none")
                       return(paste0(gender, " ", round(n/n_subject*100, 0), "%"))
                     }) %>% setdiff("none") %>% paste0(collapse = "/\n"),
                   `Sample type(s)`= c("Biopsy", "Stool") %>% 
                     purrr::map2_chr(c(n_biopsy, n_stool), function(sample_type, n) {
                       if(n == 0) return("none")
                       if(n == n_sample) return(sample_type)
                       return(paste0(sample_type, " ", round(n/n_sample*100, 0), "%"))
                     }) %>% setdiff("none") %>% paste0(collapse = "/\n")
                   ) %>% 
  dplyr::ungroup() %>% dplyr::select(-`1:dplyr::n()`)
writexl::write_xlsx(tb_table1_toshow, "tables/table1/tabe1.xlsx")
