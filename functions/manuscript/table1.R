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
    n_HC = dplyr::n_distinct(subject_accession[control == "HC"], na.rm = TRUE),
    n_nonIBD = dplyr::n_distinct(subject_accession[control == "nonIBD"], na.rm = TRUE),
    n_L.cat = dplyr::n_distinct(subject_accession[!is.na(L.cat)], na.rm = TRUE),
    n_B.cat = dplyr::n_distinct(subject_accession[!is.na(B.cat)], na.rm = TRUE),
    n_E.cat = dplyr::n_distinct(subject_accession[!is.na(E.cat)], na.rm = TRUE),
    n_male = dplyr::n_distinct(subject_accession[gender == "m"], na.rm = TRUE),
    n_female = dplyr::n_distinct(subject_accession[gender == "f"], na.rm = TRUE),
    n_missing = dplyr::n_distinct(subject_accession[is.na(gender)], na.rm = TRUE),
    n_biopsy = sum(sample_type == "biopsy", na.rm = TRUE),
    n_stool = sum(sample_type == "stool", na.rm = TRUE),
    n_antibiotics = sum(antibiotics == "y", na.rm = TRUE),
    n_antibiotics_n = sum(antibiotics == "n", na.rm = TRUE),
    n_immunosuppressants = sum(immunosuppressants == "y", na.rm = TRUE),
    n_immunosuppressants_n = sum(immunosuppressants == "n", na.rm = TRUE),
    n_steroids = sum(steroids == "y", na.rm = TRUE),
    n_steroids_n = sum(steroids == "n", na.rm = TRUE),
    n_mesalamine_5ASA = sum(mesalamine_5ASA == "y", na.rm = TRUE),
    n_mesalamine_5ASA_n = sum(mesalamine_5ASA == "n", na.rm = TRUE),
    lib_size_mean = mean(read_depth),
    lib_size_sd = sd(read_depth)
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
  dplyr::group_by(1:n()) %>% 
  dplyr::transmute(Study = ifelse(is.na(PMID),
                                  "",
                                  paste0("{", PMID, "}")) %>% 
                     paste0(Study, .),
                   `Brief description` = `Description of study design`,
                   Longitudinal = ifelse(`Longitudinal sampling?` == "Yes",
                                         "x", ""),
                   Pediatric = ifelse(`Pediatric cohort` %in% "Yes",
                                      "x", ""),
                   `N subject` = n_subject,
                   `N sample` = ifelse(`Longitudinal sampling?` == "Yes",
                                       paste0(" (", n_sample_baseline, ")"),
                                       "") %>% paste0(n_sample, .),
                   Disease = c("CD", "UC", "Control") %>% 
                     purrr::map2_chr(c(n_CD, n_UC, n_control), function(disease, n) {
                       if(n == 0) return("none")
                       return(paste0(disease, " ", n))
                     }) %>% setdiff("none") %>% paste0(collapse = "/"),
                   Classification = c("Location", "Behavior", "Extent") %>% 
                     purrr::map2_chr(c(n_L.cat, n_B.cat, n_E.cat), function(classification, n) {
                       if(n == 0) return("none")
                       return(classification)
                     }) %>% setdiff("none") %>% paste0(collapse = "/"),
                   Control = c("Healthy", "non-IBD") %>% 
                     purrr::map2_chr(c(n_HC, n_nonIBD), function(control, n) {
                       if(n == 0) return("none")
                       return(control)
                     }) %>% setdiff("none") %>% paste0(collapse = "/"),
                   Age = ifelse(is.na(age_mean),
                                "",
                                paste0(round(age_mean, 2), " (", round(age_sd, 2), ")")),
                   Gender = c("Male", "Female", "Missing") %>% 
                     purrr::map2_chr(c(n_male, n_female, n_missing), function(gender, n) {
                       if(n == 0) return("none")
                       return(paste0(gender, " ", round(n/n_subject*100, 0), "%"))
                     }) %>% setdiff("none") %>% paste0(collapse = "/"),
                   `Sample type`= c("Biopsy", "Stool") %>% 
                     purrr::map2_chr(c(n_biopsy, n_stool), function(sample_type, n) {
                       if(n == 0) return("none")
                       if(n == n_sample) return(sample_type)
                       return(paste0(sample_type, " ", round(n/n_sample*100, 0), "%"))
                     }) %>% setdiff("none") %>% paste0(collapse = "/"),
                   `Treatment` = c("Antibiotics", "5ASA", "Immunosuppressants", "Steroids") %>% 
                     purrr::map2_chr(list(c(n_antibiotics, n_antibiotics_n),
                                          c(n_mesalamine_5ASA, n_mesalamine_5ASA_n),
                                          c(n_immunosuppressants, n_immunosuppressants_n),
                                          c(n_steroids, n_steroids_n)), function(treatment, n) {
                                            if(any(n == 0)) return("none")
                                            return(paste0(treatment, " (", sum(n), ")"))
                                          }) %>% setdiff("none") %>% paste0(collapse = "/"),
                   `Library size` = paste0(round(lib_size_mean/1000, 0)*1000, 
                                           " (",
                                           round(lib_size_sd/1000, 0)*1000, ")")
                   )
readr::write_tsv(tb_table1_toshow, paste0("tables/table1/table1.tsv"))
