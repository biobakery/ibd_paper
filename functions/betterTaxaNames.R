#' More concise names for (greengenes format) microbial species
#'
#' @param order (vector of) order names
#' @param family (vector of) family names
#' @param genus (vector of) genus names
#' @param species (vector of) species names
#'
#' @return simplified species names
#' @export
betterSpeciesNames <- function(order, family, genus, species) {
  dplyr::case_when(
    species != "s__" ~ paste(genus %>% 
                               stringr::str_replace_all(stringr::fixed("g__"), ""),
                             species %>% 
                               stringr::str_replace_all(stringr::fixed("s__"), "")),
    genus != "g__" ~ paste(genus %>% 
                             stringr::str_replace_all(stringr::fixed("g__"), ""),
                           "unclassified"),
    family != "f__" ~ paste0(family %>% 
                               stringr::str_replace_all(stringr::fixed("f__"), ""),
                             "(f) unclassified"),
    order != "o__" ~ paste0(order %>% 
                              stringr::str_replace_all(stringr::fixed("o__"), ""),
                            "(o) unclassified"),
    TRUE ~ "unclassified at order"
  )
}

#' More concise names for (greengenes format) microbial genera
#'
#' @param order (vector of) order names
#' @param family (vector of) family names
#' @param genus (vector of) genus names
#'
#' @return simplified genus names
#' @export
betterGeneraNames <- function(order, family, genus) {
  dplyr::case_when(
    genus != "g__" ~ genus %>% stringr::str_replace_all(stringr::fixed("g__"), ""),
    family != "f__" ~ paste0(family %>% 
                               stringr::str_replace_all(stringr::fixed("f__"), ""),
                             "(f) unclassified"),
    order != "o__" ~ paste0(order %>% 
                              stringr::str_replace_all(stringr::fixed("o__"), ""),
                            "(o) unclassified"),
    TRUE ~ "unclassified at order"
  )
}