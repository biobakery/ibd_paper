# These functions are to replace the otu_table, sample_data, and tax_table
# functions in phyloseq, so that the returned objects are matrix, data.frame, and
# matrix, respectively
otu_table2 <- function(phylo) {
  if(!(class(phylo) == "phyloseq"))
    stop("Must have phyloseq object as input!")
  phyloseq::otu_table(phylo)@.Data
}
sample_data2 <- function(phylo) {
  if(!(class(phylo) == "phyloseq"))
    stop("Must have phyloseq object as input!")
  data.frame(phyloseq::sample_data(phylo),
             check.names = FALSE,
             stringsAsFactors = FALSE)
}
tax_table2 <- function(phylo) {
  phyloseq::tax_table(phylo)@.Data
}

# This funciton is to aggregate a list of (uniformly prepared) phyloseq objects into
# one single object
combine_phyloseq <- function(l_phylo) {
  studies <- names(l_phylo)
  if(is.null(studies)) stop("Names of l_phylo needed to provide study names!")
  
  # Check taxonamy tables are consistent; create overall taxonamy table
  l_tax <- l_phylo %>% 
    purrr::map(tax_table2)
  for(i in 1:(length(l_tax) - 1)) {
    if(all(colnames(l_tax[[1]]) == colnames(l_tax[[i + 1]])) %>% 
       not())
      stop("Not all phyloseq objects have the same taxonomic levels!")
  }
  tax_all <- l_tax %>% 
    purrr::map(.f = apply, MARGIN=1, FUN=paste, collapse = "|") %>% 
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
      otu_tmp <- otu_table2(l_phylo[[study]])
      rownames(otu_tmp) <- l_tax[[study]] %>% apply(1, paste, collapse = "|")
      otu_new <- matrix(0, nrow = nrow(mat_tax_all), ncol = ncol(otu_tmp))
      dimnames(otu_new) <- list(rownames(mat_tax_all), colnames(otu_tmp))
      otu_new[rownames(otu_tmp), colnames(otu_tmp)] <- otu_tmp
      return(otu_new)
    }) %>% 
    purrr::reduce(cbind)
  
  # Check metadata tables are consistent; create overall metadata table
  l_metadata <- l_phylo %>% 
    purrr::map(sample_data2)
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