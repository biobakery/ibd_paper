source("assets/ma_tests.R")
source("assets/tb_moderator.R")
l_results <- list()
for(test in c(tests_CD, tests_UC)) {
  load(paste0("results/4.2-maaslin2/genera/", test, "/result.RData"))
  l_results[[test]] <- result
}
load("data/physeq/genera_adj.RData")
tax_table <- smar::tax_table2(physeq_genera_adj)
suppTable_montreal <- function(l_results,
                               tests_CD,
                               tests_UC,
                               tax_table,
                               path = "supp_materials/suppTables/") {
  l_sheets <- list(`CD Behavior` = tests_CD,
                   `UC Extent` = tests_UC) %>% 
    purrr::imap(function(tests, name_tests) {
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
        })
      
      # load contrasts results
      if(name_tests == "CD Behavior")
        test_comparison <- "(B3 vs. B1) - (B23 vs. B1)"
      if(name_tests == "UC Extent")
        test_comparison <- "(E3 vs. E1) - (E23 vs. B1)"
      load(paste0("results/4.3-montreal/genera/",
                  name_tests, "/result.RData"))
      tb_ma <- result$tb_summarise %>%
        dplyr::mutate(batch = "MA") %>% 
        dplyr::filter(k > 1) %>% 
        dplyr::select(feature, coef, stderr, pval, I2) %>% 
        dplyr::filter(!is.na(coef)) %>% 
        dplyr::mutate(test = test_comparison) %>% 
        dplyr::filter(feature %in% result_all$feature)
      tb_ind <- result$result$ind.results %>% 
        purrr::map_dfr(~ .x) %>% 
        dplyr::filter(feature %in% tb_ma$feature) %>% 
        dplyr::select(feature, coef, stderr, batch) %>% 
        dplyr::filter(!is.na(coef)) %>% 
        dplyr::mutate(test = test_comparison) %>% 
        dplyr::filter(feature %in% result_all$feature)
      
      result_all <- result_all %>% 
        dplyr::bind_rows(tb_ma, tb_ind) %>% 
        dplyr::group_by(feature) %>% 
        dplyr::filter(min(qval.fdr_meta, na.rm = TRUE) < 0.1) %>% 
        dplyr::filter(sign(coef[test == tests[1] & batch == "MA"]) == sign(coef[test == tests[2] & batch == "MA"]) &
                        sign(coef[test == tests[1] & batch == "MA"] - coef[test == tests[2] & batch == "MA"]) == 
                        sign(coef[test == tests[1] & batch == "MA"])) %>% 
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
      result_all <- df_moderator_site %>% 
        tibble::rownames_to_column("batch") %>% 
        dplyr::right_join(result_all, by = "batch") %>% 
        dplyr::mutate(strata = study_site %>% 
                        {ifelse(is.na(.),
                                "MA effect",
                                .)} %>% 
                        factor(levels = c(df_moderator_site$study_site,
                                          "MA effect")))
      
      # format order test
      result_all <- result_all %>% 
        dplyr::mutate(test = test %>% 
                        stringr::str_replace_all(stringr::fixed("_vs_"),
                                                 " vs. ") %>% 
                        factor(levels = c(rev(tests) %>% 
                                            stringr::str_replace_all(stringr::fixed("_vs_"),
                                                                     " vs. "),
                                          test_comparison)))
      
      # make features ordered the same way as in Fig. 4
      features_ordered <- result_all %>% 
        dplyr::filter(batch == "MA", test == levels(test)[2]) %>% 
        dplyr::group_by(sign(coef), Rank5) %>%
        dplyr::mutate(max_coef = sign(coef)*max(abs(coef), na.rm = TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(sign(coef), max_coef, coef) %>%
        dplyr::pull(feature) %>% 
        rev()
      result_all <- result_all %>% 
        dplyr::mutate(feature = factor(feature, levels = features_ordered))
      
      # format tb
      tb_write <- result_all %>% 
        dplyr::arrange(feature, test, strata) %>% 
        dplyr::group_by(feature, test) %>% 
        dplyr::mutate(Test = ifelse((1:dplyr::n()) == 1,
                                    test %>% as.character(),
                                    NA_character_)) %>% 
        dplyr::group_by(feature) %>% 
        dplyr::mutate(Feature = ifelse((1:dplyr::n()) == 1,
                                       feature %>% as.character(),
                                       NA_character_),
                      `Cohort strata` = strata,
                      `Effect size` = coef,
                      `Standard error` = stderr,
                      `P value` = pval,
                      `Q value` = qval.fdr_meta,
                      `I2` = I2) %>% 
        dplyr::ungroup() %>% 
        dplyr::select(Feature, Test, `Cohort strata`,
                      `Effect size`, `Standard error`, `P value`, `Q value`, I2)
    })
  writexl::write_xlsx(l_sheets, path = paste0(path, "suppTable_montreal.xlsx"), format_headers = FALSE)
  return(l_sheets)
}