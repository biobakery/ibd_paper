suppTable_cohortSummary <- function(phylo,
                                    path = "supp_materials/suppTables/suppTable_cohortSummary.csv") {
  df_metadata <- sample_data2(phylo) %>% 
    dplyr::left_join(tb_1,
                     by = c("dataset_name" = "study_full_name")) %>% 
    dplyr::group_by(dataset_name) %>% 
    dplyr::mutate(is_baseline_or_unique = ifelse(`Longitudinal sampling?` == "No", 
                                                 !duplicated(subject_accession), 
                                                 is_baseline)) %>% 
    dplyr::ungroup()
  
  template <- readr::read_csv("../ibd_meta_analysis/data/template.csv")
  variables_toSummarise <- c("control",
                             "L.cat", "E.cat", "B.cat", "perianal",
                             "age_at_diagnosis", "age_at_diagnosis.cat",
                             "race", "gender", "BMI", "alcohol", "smoke",
                             "calprotectin", "PCDAI",
                             "antibiotics", "immunosuppressants", "steroids", 
                             "mesalamine_5ASA")
  template_subset <- template %>% 
    dplyr::filter(col.name %in% variables_toSummarise) %>% 
    dplyr::mutate(subject = factor(`subject specific?`, levels = c("y", "n"))) %>% 
    dplyr::arrange(subject) %>% 
    dplyr::select(col.name, var.class, `subject specific?`, allowedvalues)
  
  df_summary <- (1:nrow(template_subset)) %>% 
    purrr::map_dfc(function(i_template) {
      # for each variable, summarize per-cohort
      variable_name <- template_subset[i_template, ]$col.name
      is_subject <- template_subset[i_template, ]$`subject specific?` %in% "y"
      variable_type <- template_subset[i_template, ]$var.class
      variable_categories <- NULL
      if(variable_type == "character") 
        variable_categories <- template_subset[i_template, ]$allowedvalues %>% 
        strsplit(split = "|", fixed = TRUE) %>% 
        magrittr::extract2(1)
      
      df_result <- df_metadata %>% 
        dplyr::group_by(dataset_name) %>% 
        dplyr::summarise(summary = 
                           ifelse(is_subject,
                                  (!!rlang::sym(variable_name))[is_baseline_or_unique] %>% 
                                    summarise_variable(variable_type = variable_type,
                                                       categories = variable_categories),
                                  !!rlang::sym(variable_name) %>% 
                                    summarise_variable(variable_type = variable_type,
                                                       categories = variable_categories)))
      df_result <- df_result[, "summary", drop = FALSE]
      colnames(df_result) <- c(variable_name)
      return(df_result)
    })
  
  df_summary$dataset_name <- unique(df_metadata$dataset_name)
  
  df_summary_3000 <- df_metadata %>% 
    dplyr::group_by(dataset_name) %>% 
    dplyr::summarise(perc_3000 = (mean(read_depth > 3000)*100) %>% 
                       round(digits = 2) %>%
                       paste0("%"))
  
  tb_cohortSummary <- df_summary %>% 
    dplyr::left_join(df_summary_3000, by = "dataset_name") %>% 
    dplyr::left_join(tb_1,
                     by = c("dataset_name" = "study_full_name")) %>% 
    dplyr::select(Study,
                  Control = control,
                  `CD Montreal Location` = L.cat,
                  `CD Montreal Behavior` = B.cat,
                  `CD perianal` = perianal,
                  `Age at diagnosis` = age_at_diagnosis,
                  `Age at diagnosis Montreal classification` = age_at_diagnosis,
                  Race = race,
                  Gender = gender,
                  BMI = BMI,
                  Alcohol = alcohol,
                  Smoke = smoke,
                  `Fecal calprotectin` = calprotectin,
                  PDAI = PCDAI,
                  Antibiotics = antibiotics,
                  Immunosuppressants = immunosuppressants,
                  Steroids = steroids,
                  `Mesalamine/5ASA` = mesalamine_5ASA,
                  `% samples > 3000 read depth` = perc_3000)
  if(!is.null(path))
    readr::write_csv(tb_cohortSummary, 
                     path)
  
  return(tb_cohortSummary)
}
