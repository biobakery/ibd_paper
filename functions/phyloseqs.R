#' Transform count phyloseq object to relative abundance space. Very specific to this project,
#' as it requires a pre-recorded \code{read_depth} metadata column in the phyloseq object. This
#' is necessary because the sample-wise sums do not necessarily equal read depth, after feature
#' filtering
#'
#' @param physeq 
#'
#' @return transformed phyloseq object
#' @export
#'
#' @examples
to_relativeAbundance <- function(physeq) {
  mat_otu <- smar::otu_table2(physeq)
  if(all(mat_otu < 1))
    stop("Warning! It seems the OTU table is already on relative scale!")
  read_depth <- smar::sample_data2(physeq)$read_depth
  mat_otu_ra <- t(t(mat_otu) / read_depth)
  # if read depth is zero then all ras are zero as well
  mat_otu_ra[, read_depth == 0] <- 0
  dimnames(mat_otu_ra) <- dimnames(mat_otu)
  phyloseq::otu_table(physeq) <- otu_table(mat_otu_ra, taxa_are_rows = TRUE)
  return(physeq)
}

#' Combine list of meta-analysis phyloseq objects into one
#'
#' @param l_physeq list of (named!) phyloseq objects, one for each study
#'
#' @return combined phyloseq object
#' @importFrom magrittr "%>%"
#' @export
#'
#' @examples
combine_phyloseq <- function(l_physeq) {
  studies <- names(l_physeq)
  if(is.null(studies)) stop("Names of l_physeq needed to provide study names!")
  
  # Check taxonamy tables are consistent; create overall taxonamy table
  l_tax <- l_physeq %>% 
    purrr::map(smar::tax_table2)
  ranks1 <- l_tax %>% 
    purrr::map_dbl(~ is.na(.x) %>% 
                     not %>% 
                     apply(2, all) %>% 
                     sum)
  ranks2 <- l_tax %>% 
    purrr::map_dbl(~ is.na(.x) %>% 
                     not %>% 
                     apply(2, any) %>% 
                     sum)
  if(!all(ranks1 == ranks2)) 
    stop("Check the tax tables - some columns have sporadic missing values!")
  if(dplyr::n_distinct(ranks1) > 1)
    stop("Not all phyloseq objects have the same taxonomy level!")

  l_tax_all <- l_tax %>% 
    purrr::map(apply, MARGIN=1, FUN=paste, collapse = "|") 
  if(l_tax_all %>% 
     purrr::map_lgl(~anyDuplicated(.x) > 0) %>% 
     any) stop("Taxonamy needs to be uniquely mappable to features to aggregate!")
  tax_all <- l_tax_all %>% 
    unlist() %>% 
    unique()
  mat_tax_all <- tax_all %>% 
    sapply(strsplit, split = "|", fixed = TRUE) %>% 
    data.frame(check.names = FALSE) %>% 
    as.matrix() %>% 
    t()
  colnames(mat_tax_all) <- colnames(l_tax[[1]])
  
  # Create overall OTU table
  mat_otu_all <- studies %>% 
    purrr::map(function(study) {
      otu_tmp <- smar::otu_table2(l_physeq[[study]])
      rownames(otu_tmp) <- l_tax[[study]] %>% apply(1, paste, collapse = "|")
      otu_new <- matrix(0, nrow = nrow(mat_tax_all), ncol = ncol(otu_tmp))
      dimnames(otu_new) <- list(rownames(mat_tax_all), colnames(otu_tmp))
      otu_new[rownames(otu_tmp), colnames(otu_tmp)] <- otu_tmp
      return(otu_new)
    }) %>% 
    purrr::reduce(cbind)
  
  # Check metadata tables are consistent; create overall metadata table
  l_metadata <- l_physeq %>% 
    purrr::map(smar::sample_data2)
  for(i in 1:(length(l_tax) - 1)) {
    if(!all(colnames(l_metadata[[1]]) == colnames(l_metadata[[i + 1]])))
      stop("Not all phyloseq objects have the same metadata names!")
    if(!all(sapply(l_metadata[[1]], class) == sapply(l_metadata[[i + 1]], class)))
      stop("Not all phyloseq objects have the same metadata types!")
  }
  df_metadata_all <- l_metadata %>% 
    purrr::reduce(rbind)
  
  # Dimension names of the aggregated tables need to agree
  if(!all(colnames(mat_otu_all) == rownames(df_metadata_all))) {
    stop("Sample names in the aggregated OTU table and metadata table do not agree!")
  }
  if(!all(rownames(mat_otu_all) == rownames(mat_tax_all))) {
    stop("Feature names in the aggregated OTU table and metadata table do not agree!")
  }
  
  phyloseq(
    phyloseq::otu_table(mat_otu_all, taxa_are_rows = TRUE),
    phyloseq::sample_data(df_metadata_all),
    phyloseq::tax_table(mat_tax_all)
  )
}