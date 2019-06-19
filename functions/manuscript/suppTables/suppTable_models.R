load("results/4.1-permanova/permanova_config.RData")
source("assets/ma_tests.R")
suppTable_models <- function(tb_permanova, l_tests, l_tests_subtype,
                             path = "supp_materials/suppTables/") {
  
  # PERMANOVA models
  vars_PERMANOVA_map <- c("dataset_name" = "Study", 
                          "batch" = "Batch",
                          "sample_type" = "Biopsy vs. stool", 
                          "body_site" = "Biopsy location",
                          "IBD" = "IBD vs. contorl",
                          "disease" = "CD vs. UC",
                          "L.cat" = "CD Location", 
                          "B.cat" = "CD Behavior", 
                          "E.cat" = "UC Extent",
                          "antibiotics" = "Antibiotics", 
                          "immunosuppressants" = "Immunosuppressants", 
                          "steroids" = "Steroids", 
                          "mesalamine_5ASA" = "5-ASA", 
                          "age.cat" = "Age (< 18)", 
                          "age_at_diagnosis.cat" = "Age at diagnosis", 
                          "gender" = "Gender", 
                          "race" = "Race",
                          "subject_accession" = "Subject")
  
  tb_write_permanova <- tb_permanova %>% 
    dplyr::mutate(Variable = variable %>% 
                    dplyr::recode_factor(!!!vars_PERMANOVA_map)) %>% 
    dplyr::group_by(Variable) %>% 
    dplyr::summarise(variable_adjust = ifelse(all(is.na(variable_adjust)),
                                              NA_character_,
                                              variable_adjust %>% unique %>% 
                                                setdiff(NA_character_)),
                     variable_type = unique(variable_type),
                     variable = unique(variable)) %>% 
    dplyr::mutate(`Conditional on` = 
                    dplyr::case_when(
                      variable_adjust == "sample_type" ~ "Sample is biopsy",
                      variable_adjust == "IBD" ~ "Subject is IBD",
                      variable_adjust == "disease" & variable %in% c("L.cat", "B.cat") ~ "Subject is CD",
                      variable_adjust == "disease" & variable %in% c("E.cat") ~ "Subject is UC",
                      variable_adjust == "disease" ~ "Disease (CD/UC/control)",
                      TRUE ~ NA_character_
                    ),
                  `Permutation for all cohorts test` = dplyr::case_when(
                    Variable == "Study" ~ "Freely across cohorts but along with subjects",
                    Variable == "Batch" ~ NA_character_,
                    Variable == "Subject" ~ "Within cohorts but otherwise freely",
                    variable_type == "subject" ~ "Within cohorts; also along with subjects for longitudinal cohorts",
                    variable_type == "sample" ~ "Within cohorts; also within subjects for longitudinal cohorts"
                  ),
                  `Permutation for individual longitudinal cohort` = dplyr::case_when(
                    Variable == "Study" ~ NA_character_,
                    Variable %in% c("Batch", "Subject") ~ "Permute freely",
                    variable_type == "subject" ~ "Along with subjects",
                    variable_type == "sample" ~ "Within subjects"
                  ),
                  `Permutation for individual cross-sectional cohort` = dplyr::case_when(
                    Variable %in% c("Study", "Subject") ~ NA_character_,
                    TRUE ~ "Permute freely"
                  )) %>% 
    dplyr::select(-variable, -variable_adjust, -variable_type)
  
  tb_MA <- c(l_tests[c(tests_disease, tests_treatment)], 
             l_tests["B2_vs_B1"],
             l_tests["B3_vs_B1"],
             l_tests_subtype["CD Behavior"],
             l_tests["E2_vs_E1"],
             l_tests["E3_vs_E1"],
             l_tests_subtype["UC Extent"]) %>% 
    purrr::imap_dfr(function(parameters, i_test) {
      tibble::tibble(test = i_test,
                     batch_variable = parameters$batch_variable,
                     covariates = list(parameters$covariates))
    })
  
  var_MA_map <- c(vars_PERMANOVA_map,
                  "CD Behavior" = "(B3 vs. B1) - (B23 vs. B1)",
                  "UC Extent" = "(E3 vs. E1) - (E23 vs. E1)")
  tb_write_MA <- tb_MA %>% 
    dplyr::mutate(Test = test %>% dplyr::recode(!!!var_MA_map) %>% 
                    stringr::str_replace_all(stringr::fixed("_vs_"), " vs. "),
                  Subset = test %>% 
                    dplyr::recode("IBD_vs_control" = NA_character_,
                                  "CD_vs_control" = "CD/control subjects",
                                  "UC_vs_control" = "UC/control subjects",
                                  "CD_vs_UC" = "CD/UC subjects",
                                  "antibiotics" = "CD/UC samples with available antibiotics information",
                                  "immunosuppressants" = "CD/UC samples with available immunosuppressants information",
                                  "steroids" = "CD/UC samples with available steroids information",
                                  "mesalamine_5ASA" = "CD/UC samples with available 5-ASA information",
                                  "B2_vs_B1" = "B1/B2 CDs",
                                  "B3_vs_B1" = "B1/B3 CDs",
                                  "CD Behavior" = "CDs with available Behavior classification",
                                  "E2_vs_E1" = "E1/E2 UCs",
                                  "E3_vs_E1" = "E1/E3 UCs",
                                  "UC Extent" = "UCs with available Extent classification"),
                  `Cohort stratification` = batch_variable %>% 
                    dplyr::recode("study_site" = "Cohort, biopsy/stool",
                                  "study_site_disease" = "Cohort, biopsy/stool, CD/UC"),
                  Covariates = covariates %>% 
                    purrr::map_chr(~.x %>% 
                                     stringr::str_replace_all(stringr::fixed("_fill"), "") %>% 
                                     dplyr::recode(!!!var_MA_map) %>% 
                                     paste(collapse = ",")) %>% 
                    paste0(ifelse(test %in% c("CD Behavior", "UC Extent"),
                                  "; see Methods for model details",
                                  ""))) %>% 
    dplyr::select(Test, Subset, `Cohort stratification`, Covariates)
  
  l_sheets <- list("Omnibus testing models" = tb_write_permanova,
                   "Per-feature testing models" = tb_write_MA)
  writexl::write_xlsx(l_sheets, path = paste0(path, "suppTable_models.xlsx"), format_headers = FALSE)
  return(l_sheets)
}
