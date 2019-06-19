source("assets/ma_tests.R")
source("assets/tb_moderator.R")
l_results <- list()
for(test in c(tests_disease, tests_treatment)) {
  load(paste0("results/4.2-maaslin2/genera/", test, "/result.RData"))
  l_results[[test]] <- result
}
load("data/physeq/genera_adj.RData")
tax_table <- smar::tax_table2(physeq_genera_adj)

suppTable_interaction <- function(l_results,
                         tests_disease,
                         tests_treatment,
                         tax_table,
                         path = "supp_materials/suppTables/") {
  ll_tb_write <- list(disease = tests_disease,
                      treatment = tests_treatment) %>% 
    purrr::imap(function(tests, name_tests) {
      # create table with all results to print
      purrr::map(tests, function(i_test) {
        result_all <- paste0("results/4.4-interaction/genera/",
                             i_test,
                             "_moderator_results.tsv") %>% 
          readr::read_tsv() %>% 
          dplyr::filter(moderator_level == "difference",
                        !is.na(coef)) %>% 
          dplyr::group_by(feature) %>% 
          dplyr::filter(max(k) >= 2) %>% 
          dplyr::ungroup() %>% 
          dplyr::group_by(moderator_level) %>% 
          dplyr::mutate(qval.fdr_meta = pval %>% p.adjust("fdr")) %>% 
          dplyr::ungroup() %>% 
          dplyr::mutate(batch = moderator %>% 
                          dplyr::recode(sample_type = "stool vs. biopsy",
                                        disease = "UC vs. CD") %>% 
                          paste0("MA effect (", ., ")")) %>% 
          dplyr::group_by(feature) %>% 
          dplyr::filter(min(qval.fdr_meta) < 0.05) %>% 
          dplyr::ungroup() %>% 
          dplyr::select(feature, coef, stderr, pval, qval.fdr_meta, batch)
        if(nrow(result_all) == 0) return(data.frame(Feature = NULL,
                                                        `Cohort strata` = NULL,
                                                        `Effect size` = NULL,
                                                        `Standard error` = NULL,
                                                        `P value` = NULL,
                                                        `Q value` = NULL))
        # individual study effects
        result_all <- l_results[[i_test]]$result$ind.results %>% 
          purrr::map_dfr(~ .x) %>% 
          dplyr::select(feature, coef, stderr, batch) %>% 
          dplyr::filter(feature %in% result_all$feature) %>% 
          dplyr::bind_rows(result_all) %>% 
          dplyr::mutate(feature = feature %>% 
                          stringr::str_replace_all(stringr::fixed("|NA|NA"), ""))
        
        if(name_tests == "disease")
          result_all <- df_moderator_site %>% 
          tibble::rownames_to_column("batch") %>% 
          dplyr::right_join(result_all, by = "batch") %>% 
          dplyr::mutate(strata = study_site %>% 
                          {ifelse(is.na(.),
                                  batch,
                                  .)} %>% 
                          factor(levels = c(df_moderator_site$study_site,
                                            "MA effect (stool vs. biopsy)",
                                            "MA effect (UC vs. CD)")))
        if(name_tests == "treatment")
          result_all <- df_moderator_siteDisease %>% 
          tibble::rownames_to_column("batch") %>% 
          dplyr::right_join(result_all, by = "batch") %>% 
          dplyr::mutate(strata = study_site_disease %>% 
                          {ifelse(is.na(.),
                                  batch,
                                  .)} %>% 
                          factor(levels = c(df_moderator_siteDisease$study_site_disease,
                                            "MA effect (stool vs. biopsy)",
                                            "MA effect (UC vs. CD)")))
        
        # order
        tb_write <- result_all %>% 
          dplyr::arrange(feature, strata) %>% 
          dplyr::group_by(feature) %>% 
          dplyr::transmute(Feature = ifelse((1:dplyr::n()) == 1,
                                            feature,
                                            NA_character_),
                           `Cohort strata` = strata,
                           `Effect size` = coef,
                           `Standard error` = stderr,
                           `P value` = pval,
                           `Q value` = qval.fdr_meta) %>% 
          dplyr::ungroup() %>% 
          dplyr::select(-feature)
        
        return(tb_write)
      })
    })
  l_sheets <- c(ll_tb_write[[1]], ll_tb_write[[2]])
  names(l_sheets) <- c(tests_disease, tests_treatment) %>% 
    dplyr::recode("CD_vs_UC" = "CD vs. UC",
                  "UC_vs_control" = "UC vs. control",
                  "CD_vs_control" = "CD vs. control",
                  "IBD_vs_control" = "IBD vs. control",
                  "mesalamine_5ASA" = "5-ASA",
                  "steroids" = "Steroids",
                  "immunosuppressants" = "Immunosuppressants",
                  "antibiotics" = "Antibiotics")
  l_sheets <- l_sheets[l_sheets %>% purrr::map_lgl(~nrow(.x) > 0)]
  writexl::write_xlsx(l_sheets, path = paste0(path, "suppTable_interaction.xlsx"), format_headers = FALSE)
  return(l_sheets)
}