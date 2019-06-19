source("assets/ma_tests.R")
source("assets/tb_moderator.R")
l_results <- list()
for(test in c(tests_disease, tests_treatment)) {
  load(paste0("results/4.2-maaslin2/genera/", test, "/result.RData"))
  l_results[[test]] <- result
}
load("data/physeq/genera_adj.RData")
tax_table <- smar::tax_table2(physeq_genera_adj)

suppTable_MA <- function(l_results,
                         tests_disease,
                         tests_treatment,
                         tax_table,
                         path = "supp_materials/suppTables/") {
  ll_tb_write <- list(disease = tests_disease,
                      treatment = tests_treatment) %>% 
    purrr::imap(function(tests, name_tests) {
      # create table with all results to print
      result_all <- l_results[tests] %>% 
        purrr::imap_dfr(function(result, test) {
          tb_ma <- result$tb_summarise %>%
            dplyr::filter(k > 1) %>% 
            dplyr::mutate(qval.fdr_meta = p.adjust(pval, method = "fdr")) %>% 
            dplyr::select(feature, coef, stderr, pval, qval.fdr_meta, I2) %>% 
            dplyr::mutate(batch = "MA")
          tb_ind <- result$result$ind.results %>% 
            purrr::map_dfr(~ .x) %>% 
            dplyr::select(feature, coef, stderr, batch) %>% 
            dplyr::filter(feature %in% tb_ma$feature)
          dplyr::bind_rows(tb_ma, tb_ind) %>% 
            dplyr::filter(!is.na(coef)) %>% 
            dplyr::mutate(test = test) %>% 
            return()
        }) %>% 
        dplyr::group_by(feature) %>% 
        dplyr::filter(min(qval.fdr_meta, na.rm = TRUE) < 0.05) %>% 
        dplyr::ungroup() %>% 
        dplyr::left_join(tax_table %>% 
                           as.data.frame() %>% 
                           tibble::rownames_to_column("feature"),
                         by = "feature") %>% 
        dplyr::filter(Rank5 != "f__") %>% 
        dplyr::mutate(feature = feature %>% 
                        stringr::str_replace_all(stringr::fixed("|NA|NA"),
                                                 ""))
      
      # add batch annotation
      if(name_tests == "disease")
        result_all <- df_moderator_site %>% 
          tibble::rownames_to_column("batch") %>% 
          dplyr::right_join(result_all, by = "batch") %>% 
          dplyr::mutate(strata = study_site %>% 
                          {ifelse(is.na(.),
                                  "MA effect",
                                  .)} %>% 
                          factor(levels = c(df_moderator_site$study_site,
                                            "MA effect")))
      if(name_tests == "treatment")
        result_all <- df_moderator_siteDisease %>% 
          tibble::rownames_to_column("batch") %>% 
          dplyr::right_join(result_all, by = "batch") %>% 
          dplyr::mutate(strata = study_site_disease %>% 
                          {ifelse(is.na(.),
                                  "MA effect",
                                  .)} %>% 
                          factor(levels = c(df_moderator_siteDisease$study_site_disease,
                                            "MA effect")))
      
      # make features ordered the same way as in Fig. 3
      features_ordered <- result_all %>% 
        dplyr::filter(batch == "MA" & test == tests[1]) %>% 
        dplyr::group_by(sign(coef), Rank5) %>%
        dplyr::mutate(max_coef = sign(coef)*max(abs(coef), na.rm = TRUE)) %>% 
        dplyr::ungroup() %>%
        dplyr::arrange(sign(coef), max_coef, coef) %>%
        dplyr::pull(feature) %>% 
        unique()
      result_all <- result_all %>% 
        dplyr::mutate(feature = factor(feature, levels = features_ordered))
      
      # format tb
      tb_write <- result_all %>% 
        dplyr::arrange(test, feature, strata) %>% 
        dplyr::group_by(test, feature) %>% 
        dplyr::mutate(Feature = ifelse((1:dplyr::n()) == 1,
                                       feature %>% as.character(),
                                       NA_character_),
                      `Cohort strata` = strata,
                      `Effect size` = coef,
                      `Standard error` = stderr,
                      `P value` = pval,
                      `Q value` = qval.fdr_meta,
                      `I2` = I2) %>% 
        dplyr::ungroup()
      
      l_tb_write <- tests %>% 
        purrr::map(function(i_test) {
          tb_write %>% 
            dplyr::filter(test == i_test) %>% 
            dplyr::select(Feature, `Cohort strata`, `Standard error`, `P value`, `Q value`, I2)
        })
      
      names(l_tb_write) <- tests %>% 
        dplyr::recode("CD_vs_UC" = "CD vs. UC",
                      "UC_vs_control" = "UC vs. control",
                      "CD_vs_control" = "CD vs. control",
                      "IBD_vs_control" = "IBD vs. control",
                      "mesalamine_5ASA" = "5-ASA",
                      "steroids" = "Steroids",
                      "immunosuppressants" = "Immunosuppressants",
                      "antibiotics" = "Antibiotics")
      
      return(l_tb_write)
    })
  l_sheets <- c(ll_tb_write[[1]], ll_tb_write[[2]])
  writexl::write_xlsx(l_sheets, path = paste0(path, "suppTable_MA.xlsx"), format_headers = FALSE)
  return(l_sheets)
}
