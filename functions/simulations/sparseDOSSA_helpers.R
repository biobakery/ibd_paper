create_metadataMatrix <- function(df_metadata,
                                  metadata_type,
                                  fDummyData = TRUE,
                                  scale = FALSE) {
  if(any(colnames(df_metadata) != names(metadata_type)))
    stop("Variable names in df_metadata and metadata_type do not agree!")
  if(any(sapply(df_metadata, class) != metadata_type))
    stop("Variable classes in df_metadata and metadata_type do not agree!")
  l_mat_metadata <- lapply(names(metadata_type), function(variable) {
    if(metadata_type[variable] != "factor" | !fDummyData) {
      mat_tmp <- rbind(df_metadata[, variable])
      rownames(mat_tmp) <- variable
      return(mat_tmp)
    } else {
      lvls <- levels(df_metadata[, variable])
      nlvls <- nlevels(df_metadata[, variable])
      mat_tmp <- t(sapply(1:nlvls, function(ilvl) {
        (df_metadata[, variable] == lvls[ilvl]) * 1
      }))
      rownames(mat_tmp) <- paste0(variable, "_", 1:nlvls)
      return(mat_tmp)
    }
  })
  mat_metadata <- Reduce("rbind", l_mat_metadata)
  if(scale) {
    mean_metadata <- apply(mat_metadata, 1, mean)
    sd_metadata <- apply(mat_metadata, 1, sd)
    sd_metadata[sd_metadata == 0] <- 1
    mat_metadata <- (mat_metadata - mean_metadata) / sd_metadata
  }
  return(mat_metadata)
}

create_effectSize <- function(effectSize,
                              df_metadata,
                              metadata_type,
                              fDummyData = TRUE) {
  if(any(names(effectSize) != names(metadata_type)))
    stop("Variable names in effectSize and metadata_type do not agree!")
  if(any(colnames(df_metadata) != names(metadata_type)))
    stop("Variable names in df_metadata and metadata_type do not agree!")
  if(any(sapply(df_metadata, class) != metadata_type))
    stop("Variable classes in df_metadata and metadata_type do not agree!")
  l_effectSize <- lapply(names(metadata_type), function(variable) {
    if(!fDummyData | metadata_type[variable] != "factor") {
      effectSize_tmp <- effectSize[variable]
      names(effectSize_tmp) <- variable
      return(effectSize_tmp)
    } else {
      nlvls <- nlevels(df_metadata[, variable])
      effectSize_tmp <- rep(effectSize[variable], nlvls)
      names(effectSize_tmp) <- paste0(variable, "_", 1:nlvls)
      return(effectSize_tmp)
    }
  })
  return(unlist(l_effectSize))
}

create_spikein.mt <- function(number_features,
                              percent_spiked,
                              effectSize,
                              direction = TRUE,
                              same_features = FALSE,
                              seed) {
  nFeatureSpiked <- floor(number_features * percent_spiked)
  if(nFeatureSpiked == 0)
    stop("No features are spiked in with current configuration!")
  
  if(same_features) features_spike <- sample.int(n = number_features, size = nFeatureSpiked)
  set.seed(seed)
  l_spikein.mt <- lapply(1:length(effectSize), function(i) {
    if(!same_features) features_spike <- sample.int(n = number_features, size = nFeatureSpiked)
    strength <- rep(effectSize[i], nFeatureSpiked)
    if(direction) strength[sample.int(nFeatureSpiked, size = floor(nFeatureSpiked / 2))] <- -effectSize[i]
    data.frame(feature = features_spike,
               metadata = i,
               strength = strength,
               stringsAsFactors = FALSE, row.names = NULL)
  })
  return(Reduce("rbind", l_spikein.mt))
}

format_spikein.mt <- function(df_spikein.mt) {
  df_spikein.mt %>% 
    dplyr::filter(strength != 0) %>% 
    dplyr::group_by(feature) %>% 
    dplyr::summarise(metadata = paste(metadata, collapse = ";"), 
                     strength = paste(strength, collapse = ";")) %>% 
    dplyr::ungroup()
}

extract_sparseDOSSA <- function(sparseDOSSA_fit) {

  # metadata + feature data
  sparsedossa_results <- as.data.frame(sparseDOSSA_fit$OTU_count)
  rownames(sparsedossa_results) <- sparsedossa_results$X1
  nMetadata <- sum(grepl("Metadata", sparsedossa_results$X1, fixed = TRUE))
  nMicrobes <- sum(grepl("Feature_spike", sparsedossa_results$X1, fixed = TRUE))
  nSamples <- ncol(sparsedossa_results) - 1
  sparsedossa_results <- sparsedossa_results[-1, -1]
  colnames(sparsedossa_results) <- paste('Sample', 1:nSamples, sep='')
  data <- as.matrix(sparsedossa_results[-c((nMetadata+1):(2*nMicrobes+nMetadata)), ])
  data <- data.matrix(data)
  class(data) <- "numeric"

  # Spiked-in features and metadata
  # truth <- c(unlist(sparseDOSSA_fit$truth))
  # truth <- truth[!stringi::stri_detect_fixed(truth,":")]
  # significant_features <- truth[grepl("Feature", truth, fixed = TRUE)]
  # significant_metadata <- truth[-(1+nMetadata)] %>%
  #   stringr::str_subset("Metadata") %>%
  #   stringr::str_replace("_Level_.+", "") %>%
  #   unique

  # Extract Metadata
  # metadata <- as.data.frame(t(data[(1:nMetadata), ]))

  # Rename True Positive Metadata - Same Format at Mcmurdie and Holmes (2014)
  # which.TP <- colnames(metadata) %in% significant_metadata
  # meta_newname <- paste0(colnames(metadata)[which.TP], "_TP")
  # colnames(metadata)[which.TP] <- meta_newname

  # Extract Features
  features <- as.data.frame(t(data[-c(1:nMetadata),]))

  # Rename Features and True Positive Features - Same Format at Mcmurdie and Holmes (2014)
  # wh.TP <- colnames(features) %in% significant_features
  colnames(features) <- paste("Feature", 1:nMicrobes, sep = "")
  # newname <- paste0(colnames(features)[wh.TP], "_TP")
  # colnames(features)[wh.TP] <- newname

  # feature table has rows as features
  features <- t(features)

  # Return as list
  return(list(features=features))
}
