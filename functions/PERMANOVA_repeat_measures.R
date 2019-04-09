#' PERMANOVA with repeat measure-aware permutations. Block sizes are allowed to
#' differ.
#'
#' @param D An N-by-N distance matrix (must be a \code{dist} object).
#' @param permute_within Data frame with N rows containing metadata to test per
#' sample
#' @param blocks A length-N vector containing the block structure.
#' @param block_data Data frame with per-block metadata. If \code{blocks} is
#' numeric, \code{block_data}'s rows must match those indices. If \code{blocks}
#' is a \code{factor}, then the row ordering must match the factor levels. If
#' \code{blocks} is a character vector, \code{block_data} must have rownames
#' matching the contents of \code{blocks}.
#' @param permutations Number of permutations to test
#' @param metadata_order Order of the metadata in the model. If not given, this
#' is assumed to be within-block metadata first, followed by block metadata, in
#' the order given in \code{permute_within} and \code{block_data}.
#' @return Same structure as \code{adonis}, with p-values recalculated based on
#' permutations that are aware of the block structure of the data.
#' Metadata in \code{permute_within} are permuted within blocks, whereas
#' metadata in \code{block_data} are first permuted across blocks, and then
#' assigned to samples according to the block structure.
PERMANOVA_repeat_measures <- function(
  D, permute_within, blocks = NULL, block_data = NULL, permutations=999,
  metadata_order = c(names(permute_within), names(block_data)),
  ncores=NULL,
  na.rm=F) {

  # Make sure D is a dist object
  if (class(D) != "dist") {
    stop("D must be a dist object")
  }

  # Default to free permutations if blocks is not given
  if (!is.null(block_data) && is.null(blocks)) {
    stop("blocks must be given if block_data is present")
  } else if (is.null(blocks)) {
    blocks <- rep(1, nrow(permute_within))
    block_data <- as.data.frame(matrix(0, nrow=1, ncol=0))
  } else if (length(unique(blocks)) == 1) {
    warning("blocks only contains one unique value")
  }

  # Ensure no metadata overlap between permute_within and block_data
  if (length(intersect(names(permute_within), names(block_data))) > 0) {
    stop("metadata is repeated across permute_within and block_data")
  }

  # Ensure that metadata_order only contains stuff in permute_within and block_data
  if(length(setdiff(metadata_order, union(names(permute_within), names(block_data)))) > 0) {
    stop("metadata_order contains metadata not in permute_within and block_data")
  }

  # Ensure that the data in permute_within matches that in dist
  ord <- rownames(as.matrix(D))
  if (length(ord) != nrow(permute_within) || length(blocks) != length(ord)) {
    stop("blocks, permute_within, and D are not the same size")
  }
  if (is.null(rownames(permute_within))) {
    warning("permute_within has no rownames - can't verify sample orders")
  } else if (!all(ord == rownames(permute_within))) {
    stop("rownames do not match between permute_within and D")
  }

  # Ensure matching between blocks and block_data
  if (any(is.na(blocks))) {
    stop("NAs are not allowed in blocks")
  }
  if (is.factor(blocks)) {
    if (any(!(levels(blocks) %in% rownames(block_data)))) {
      stop("not all block levels are contained in block_data")
    }
    # Match blocks with block_data and discard level information
    block_data <- block_data[match(levels(blocks), rownames(block_data)), , drop=F]
    blocks <- as.numeric(blocks)
  } else if (is.numeric(blocks)) {
    if (blocks < 1 || max(blocks) > nrow(block_data)) {
      stop("Numeric blocks has indices out of range")
    }
  } else if (is.character(blocks)) {
    if (is.null(rownames(block_data)) || !all(blocks %in% rownames(block_data))) {
      stop("blocks does not match the rownames of block_data")
    }
    # Transform to numeric
    blocks <- match(blocks, rownames(block_data))
  } else {
    stop("blocks must be a numeric, factor, or character vector")
  }

  # Error out on NA metadata rather than allowing adonis to error out with
  # a totally nonsensical error message
  na.removed <- 0
  if (any(is.na(permute_within)) || any(is.na(block_data))) {
    if (na.rm) {
      n_prerm <- length(blocks)

      # Remove NAs in block_data
      hasna <- (rowSums(is.na(block_data)) > 0) | (sapply(split(rowSums(is.na(permute_within)) > 0, blocks), mean) == 1)
      block_data <- block_data[!hasna,, drop=F]
      keep <- !hasna[blocks]
      blocks <- cumsum(!hasna)[blocks]

      blocks <- blocks[keep]
      permute_within <- permute_within[keep,, drop=F]
      D <- as.matrix(D)[keep, keep]
      # block_data is not subset, as the rows with NAs are no longer referenced in blocks

      # Remove NAs in permute_within
      keep <- rowSums(is.na(permute_within)) == 0
      blocks <- blocks[keep]
      permute_within <- permute_within[keep,, drop=F]
      D <- as.dist(D[keep, keep])

      if (length(blocks) < ncol(permute_within) + ncol(block_data)) {
        stop(sprintf("After omitting samples with NAs, the number of samples (%d) is less than the number of metadata (%d)",
                     length(blocks), ncol(permute_within) + ncol(block_data)))
      } else if (length(blocks) < n_prerm * 0.5) {
        warning(sprintf("Removed %d samples with NA metadata", n_prerm - length(blocks)))
      }
      na.removed <- n_prerm - length(blocks)
    } else {
      stop("Some metadata is NA! adonis does not support any NA in the metadata")
    }
  }

  # Warn on some suspicious input
  persample <- apply(permute_within, 1, function(x)is.factor(x) && !any(duplicated(x)))
  if (any(persample)) {
    warning(sprintf("%s in permute_within has one DOF per sample.", colnames(permute_within)[which(persample)[1]]))
  }
  if (length(unique(blocks)) < nrow(block_data)) {
    warning("Not all blocks have a sample associated with them. Block permutations will still be performed over the full set of blocks - if this is not desired, subset block_data to only the blocks which appear in the data.")
  }
  if (!any(duplicated(blocks))) {
    warning("blocks contains no duplicated elements")
  }

  library(vegan)
  library(permute)

  # Test statistic from non-permuted data
  mtdat <- cbind(permute_within, block_data[blocks,,drop=F])
  ad <- adonis(D ~ ., permutations=0, data=mtdat[, metadata_order, drop=F])
  R2 <- ad$aov.tab$R2
  names(R2) <- rownames(ad$aov.tab)

  # Permutations
  if(is.null(ncores)) {
    nullsamples <- matrix(NA, nrow=length(R2), ncol=permutations)
    for (i in seq_len(permutations)) {
      within.i <- shuffle(nrow(permute_within), control=how(blocks=blocks))
      block.i <- sample(seq_len(nrow(block_data)))
      mtdat <- cbind(
        permute_within[within.i,,drop=F],
        block_data[block.i,,drop=F][blocks,,drop=F])
      perm.ad <- adonis(D ~ ., permutations=0, data=mtdat[, metadata_order, drop=F])

      nullsamples[,i] <- perm.ad$aov.tab$R2
    }
  }
  if(!is.null(ncores)){
    library("foreach")
    library("doParallel")

    registerDoParallel(ncores)
    nullsamples <- foreach(i = seq_len(permutations),
                           .combine = cbind) %dopar% {
                             within.i <- shuffle(nrow(permute_within), control=how(blocks=blocks))
                             block.i <- sample(seq_len(nrow(block_data)))
                             mtdat <- cbind(
                               permute_within[within.i,,drop=F],
                               block_data[block.i,,drop=F][blocks,,drop=F])
                             perm.ad <- adonis(D ~ ., permutations=0, data=mtdat[, metadata_order, drop=F])
                             return(perm.ad$aov.tab$R2)
                           }
    stopImplicitCluster()
  }


  # For residuals, test the other direction (i.e. p-value of all covariates)
  n <- length(R2)
  R2[n-1] <- 1 - R2[n-1]
  nullsamples[n-1,] <- 1 - nullsamples[n-1,]

  # P value calculation similar to adonis's
  exceedances <- rowSums(nullsamples > R2)
  P <- (exceedances + 1) / (permutations + 1)

  P[n] <- NA    # No p-values for "Total"
  ad$aov.tab$`Pr(>F)` <- P

  if (na.rm) {
    ad$na.removed <- na.removed
  }

  return (ad)
}

#' PERMANOVA with two blocking structures?
#'
#' @param D An N-by-N distance matrix (must be a \code{dist} object)
#' @param study A per-sample (length-N) character vector containing the study identifiers
#' @param study_longitudinal A per-study logical vector indicate if the studies are
#' longitudinal. Must have names matching the contents of \code{study}
#' @param subject A per-sample (length-N) character vector containing the subject identifiers
#' @param subject_data Data frame with per-subject metadata. These will be permuted within
#' study. Must have rownames matching the contents of \code{blocks}
#' @param sample_data Data frame with per-sample metadata sample. These will be permuted
#' within subject for longitudinal studies, and permuted within study for cross sectionals.
#' @param metadata_order Order of the metadata in the model. If not given, this
#' is assumed to be within-block metadata first, followed by block metadata, in
#' the order given in \code{subject_data} and \code{sample_data}
#' @param permutations Number of permutations to test
#' @param ncores Number of cores (for parallel computing)
#' @return Same structure as \code{adonis}, with p-values recalculated based on
#' permutations that are aware of the block structure of the data.
PERMANOVA_repeat_measures_meta <- function(
  D,
  study, study_longitudinal,
  subject, subject_data = NULL,
  sample_data = NULL,
  metadata_order = c(names(subject_data), names(sample_data)),
  permutations=999, ncores=1)
{

  # Make sure D is a dist object
  if (class(D) != "dist") {
    stop("D must be a dist object")
  }

  # check sample identifiers in D, study, and subject
  if(nrow(as.matrix(D)) != length(study) |
     nrow(as.matrix(D)) != length(subject))
    stop("sample number from D, study, and subject must match!")
  if(any(rownames(as.matrix(D)) != names(study)) |
     any(rownames(as.matrix(D)) != names(subject)))
    stop("sample identifiers from D, study, and subject must match!")

  # check study data
  if(!is.character(study))
    stop("study identifiers must be character!")
  if(!is.logical(study_longitudinal))
    stop("study_longitudinal must be TRUE/FALSE!")
  if(length(unique(study)) != length(study_longitudinal))
    stop("study number from study and study_longitudinal must match!")
  if(!setequal(study, names(study_longitudinal)))
    stop("study identifiers from study and study_longitudinal must match!")

  # subject_data and sample_data cannot both be missing
  if(is.null(subject_data) & is.null(sample_data)) {
    stop("At least one of subject_data and sample_data must be provided!")
  }

  # check subject data, create if missing
  if(!is.character(subject))
    stop("subject must be character!")
  if(is.null(subject_data)) {
    subject_data <- data.frame(place_holder_subject = rep(1, length(unique(subject))),
                               row.names = unique(subject))
  }
  if(length(unique(subject)) != nrow(subject_data))
    stop("subject number from subject and subject_data must match!")
  if(!setequal(subject, rownames(subject_data)))
    stop("subject identifiers from subject and subject_data must match!")

  # study identifiers and subject identifiers cannot overlap for the funciton to work
  if(length(intersect(study, subject)) > 0)
    stop("study identifiers and subject identifiers cannot overlap!")
  # check that subject are study-specific
  if(any(apply(table(study, subject) > 0, 2, sum) > 1))
    stop("subject must be study-specific!")
  # map studies to subjects, and create block variable to permute
  # subject-specific data with
  study_subject <- study[!duplicated(subject)]
  names(study_subject) <- subject[!duplicated(subject)]
  study_permute <- study_subject[rownames(subject_data)]

  # create sample data if missing
  if(is.null(sample_data)) {
    sample_data <- data.frame(place_holder_sample = rep(1, nrow(as.matrix(D))),
                              row.names = rownames(as.matrix(D)))
  }
  if(nrow(as.matrix(D)) != nrow(sample_data))
    stop("sample number from D and sample_data must match!")
  if(any(rownames(as.matrix(D)) != rownames(sample_data)))
    stop("sample identifiers from D and sample_data must match!")
  # create subject identifiers to permute sample-specific data with
  # these are just subject identifiers for longitudinal data
  # and study identifiers for cross-sectional data
  subject_permute <- ifelse(study_longitudinal[study],
                            subject,
                            study)

  # Ensure no metadata overlap between subject_data and sample_data
  if(length(intersect(names(subject_data), names(sample_data))) > 0) {
    stop("metadata is repeated across subject_data and sample_data")
  }

  # Ensure that metadata_order only contains stuff in subject_data and sample_data
  if(length(setdiff(metadata_order, union(names(subject_data), names(sample_data)))) > 0) {
    stop("metadata_order contains metadata not in subject_data and sample_data!")
  }


  # Warn on some suspicious input
  # persample <- apply(sample_data, 1, function(x)is.factor(x) && !any(duplicated(x)))
  # if (any(persample)) {
  #   warning(sprintf("%s in sample_data has one DOF per sample.", colnames(sample_data)[which(persample)[1]]))
  # }
  # if (length(unique(subject)) < nrow(subject_data)) {
  #   warning("Not all subject have a sample associated with them. Block permutations will still be performed over the full set of subject - if this is not desired, subset subject_data to only the subject which appear in the data.")
  # }
  # if (!any(duplicated(subject))) {
  #   warning("subject contains no duplicated elements")
  # }

  library(vegan)
  library(permute)
  library(foreach)
  library(doParallel)

  mtdat <- cbind(subject_data[subject,,drop=F], sample_data)
  # Error out on NA metadata rather than allowing adonis to error out with a totally
  # nonsensical error message
  if (any(is.na(mtdat[, metadata_order, drop=F]))) {
    stop("Some metadata is NA! adonis does not support any NA in the metadata")
  }
  # Test statistic from non-permuted data
  ad <- adonis(D ~ ., permutations=0, data=mtdat[, metadata_order, drop=F])
  R2 <- ad$aov.tab$R2
  names(R2) <- rownames(ad$aov.tab)

  registerDoParallel(ncores)
  nullsamples <- foreach(i = seq_len(permutations),
                         .combine = cbind) %dopar%
    {
      subject.i <- shuffle(nrow(subject_data), control=how(blocks=study_permute))
      sample.i <- shuffle(nrow(sample_data), control=how(blocks=subject_permute))
      mtdat <- cbind(subject_data[subject.i,,drop=F][subject,,drop=F],
                     sample_data[sample.i,,drop=F])
      perm.ad <- adonis(D ~ ., permutations=0, data=mtdat[, metadata_order, drop=F])
      return(perm.ad$aov.tab$R2)
    }
  stopImplicitCluster()


  # For residuals, test the other direction (i.e. p-value of all covariates)
  n <- length(R2)
  R2[n-1] <- 1 - R2[n-1]
  nullsamples[n-1,] <- 1 - nullsamples[n-1,]

  # P value calculation similar to adonis's
  exceedances <- rowSums(nullsamples > R2)
  P <- (exceedances + 1) / (permutations + 1)

  P[n] <- NA    # No p-values for "Total"
  ad$aov.tab$`Pr(>F)` <- P

  return (ad)
}

fit_permanova_overall(D = dist_genera_adj,
                       variable = variable,
                       variable_class = "subject",
                       covariates = NULL,
                       block_covariates = NULL,
                       block_variable = "subject_accession",
                       data = metadata_test,
                       permutations = n_permutations,
                       ncores = ncores)

fit_permanova_per_study(D = dist_genera_adj,
                        variable = variable,
                        variable_class = "subject",
                        covariates = NULL,
                        block_covariates = NULL,
                        block_variable = "subject_accession",
                        data = metadata_test,
                        permutations = n_permutations,
                        ncores = ncores)


fit_permanova_variable <- function(
  D,
  variable, # variable to test for R2
  variable_class, # is this a per-sample or per-subject variable?
  covariates = NULL,
  block_covariates = NULL, # additional covariates (subject specific) to test (come before variable in the model)
  block_variable,
  data,
  permutations = 999,
  ncores = NULL
) {
  if(!all(c(variable, covariates, block_covariates, block_variable) %in% names(data)))
    stop("variable and/or block_covariates and/or  block_variable not present in data!")
  if(!(variable_class %in% c("sample", "subject")))
    stop("variable_class must be either sample or subject!")

  if(any(is.na(data[, c(block_variable, block_covariates, covariates), drop = FALSE])))
    stop("Cannot have missing values in block_covariates or block variable!")

  # Fill in is_na indicator if variable has any missing values
  variable_na <- NULL
  if(any(is.na(data[, variable]))) {
    variable_na <- paste0(variable, "_is_NA")
    data[, variable_na] <- is.na(data[, variable])
    data[, variable] <- fill_na(data[, variable])
  }

  permute_within <- data.frame(rows_sample = 1:nrow(data)) # this one has to be there no matter what
  rownames(permute_within) <- rownames(data)
  permute_within <- cbind(permute_within, data[, covariates, drop = FALSE])
  if(variable_class == "sample")
    permute_within <- cbind(permute_within, data[, c(variable_na, variable), drop = FALSE])

  if(!is.null(block_variable)) {
    blocks <- data[, block_variable, drop = TRUE]
    block_data <- data.frame(rows_subject = blocks)
    block_data <- cbind(block_data, data[, block_covariates, drop = FALSE])
    if(variable_class == "subject")
      block_data <- cbind(block_data, data[, c(variable_na, variable), drop = FALSE])

    # make sure that all variables in block data are indeed block specific
    block_data <-
      block_data %>%
      dplyr::group_by(rows_subject) %>%
      dplyr::distinct()
    test_block_data <- block_data %>%
      dplyr::summarise(n_distinct = n())
    if(!all(test_block_data$n_distinct == 1))
      stop("Block variables aren't unique!")
    block_data <- as.data.frame(dplyr::ungroup(block_data))
    rownames(block_data) <- block_data$rows_subject
  }
  if(is.null(block_variable)) {
    blocks <- NULL
    block_data <- NULL
  }

  metadata_order <- c(block_covariates, covariates, variable_na, variable)
  fit_adonis <- PERMANOVA_repeat_measures(D = D,
                                          permute_within = permute_within,
                                          blocks = blocks,
                                          block_data = block_data,
                                          permutations = permutations,
                                          metadata_order = metadata_order,
                                          ncores = ncores, na.rm = FALSE)

  return(fit_adonis)
}

fit_permanova_variable <- function(
  D,
  variable, # variable to test for R2
  variable_class, # is this a per-sample or per-subject variable?
  covariates = NULL,
  subject_covariates = NULL, # additional covariates (subject specific) to test (come before variable in the model)
  subject,
  study,
  study_longitudinal,
  data,
  permutations = 999,
  ncores = NULL
) {
  if(!all(c(variable, covariates, subject, block_variable) %in% names(data)))
    stop("variable and/or block_covariates and/or  block_variable not present in data!")
  if(!(variable_class %in% c("sample", "subject")))
    stop("variable_class must be either sample or subject!")

  if(any(is.na(data[, c(block_variable, block_covariates, covariates), drop = FALSE])))
    stop("Cannot have missing values in block_covariates or block variable!")

  # Fill in is_na indicator if variable has any missing values
  variable_na <- NULL
  if(any(is.na(data[, variable]))) {
    variable_na <- paste0(variable, "_is_NA")
    data[, variable_na] <- is.na(data[, variable])
    data[, variable] <- fill_na(data[, variable])
  }

  permute_within <- data.frame(rows_sample = 1:nrow(data)) # this one has to be there no matter what
  rownames(permute_within) <- rownames(data)
  permute_within <- cbind(permute_within, data[, covariates, drop = FALSE])
  if(variable_class == "sample")
    permute_within <- cbind(permute_within, data[, c(variable_na, variable), drop = FALSE])

  if(!is.null(block_variable)) {
    blocks <- data[, block_variable, drop = TRUE]
    block_data <- data.frame(rows_subject = blocks)
    block_data <- cbind(block_data, data[, block_covariates, drop = FALSE])
    if(variable_class == "subject")
      block_data <- cbind(block_data, data[, c(variable_na, variable), drop = FALSE])

    # make sure that all variables in block data are indeed block specific
    block_data <-
      block_data %>%
      dplyr::group_by(rows_subject) %>%
      dplyr::distinct()
    test_block_data <- block_data %>%
      dplyr::summarise(n_distinct = n())
    if(!all(test_block_data$n_distinct == 1))
      stop("Block variables aren't unique!")
    block_data <- as.data.frame(dplyr::ungroup(block_data))
    rownames(block_data) <- block_data$rows_subject
  }
  if(is.null(block_variable)) {
    blocks <- NULL
    block_data <- NULL
  }

  metadata_order <- c(block_covariates, covariates, variable_na, variable)
  fit_adonis <- PERMANOVA_repeat_measures(D = D,
                                          permute_within = permute_within,
                                          blocks = blocks,
                                          block_data = block_data,
                                          permutations = permutations,
                                          metadata_order = metadata_order,
                                          ncores = ncores, na.rm = FALSE)

  return(fit_adonis)
}