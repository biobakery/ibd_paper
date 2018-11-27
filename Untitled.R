fit_rma.mod <- function(l.Maaslin.fit, data.moderator,
                            method = "REML"){
  lvl.batch <- names(l.Maaslin.fit)
  if(!all(lvl.batch %in% rownames(data.moderator)))
    stop("data.moderator must have all the batches fitted in Maaslin!")
  data.moderator <- data.moderator[lvl.batch, , drop = FALSE]
  exposure <- unique(l.Maaslin.fit[[1]]$metadata)
  values.exposure <- unique(l.Maaslin.fit[[1]]$value)
  features <- unique(l.Maaslin.fit[[1]]$feature)
  l.results <- list()
  for(value.exposure in values.exposure) {
    i.result <- data.frame(feature = features,
                           exposure = value.exposure,
                           tau2 = NA,
                           se.tau2 = NA,
                           p.tau2 = NA,
                           p.moderator = NA,
                           I2 = NA,
                           H2 = NA,
                           R2 = NA, stringsAsFactors = FALSE)
    rownames(i.result) <- i.result$feature
    # sanity check
    if(any(features != l.Maaslin.fit[[2]][l.Maaslin.fit[[2]]$value == value.exposure, "feature"]))
      stop("Feature names don't match between l.Maaslin.fit components!")
    betas <- sapply(l.Maaslin.fit, function(i.Maaslin.fit) {
      i.Maaslin.fit[i.Maaslin.fit$value == value.exposure, "coef"]
    })
    sds <- sapply(l.Maaslin.fit, function(i.Maaslin.fit) {
      i.Maaslin.fit[i.Maaslin.fit$value == value.exposure, "stderr"]
    })
    rownames(betas) <- rownames(sds) <- features
    ind.feature <- !is.na(betas) & !is.na(sds) & (sds != 0)
    count.feature <- apply(ind.feature, 1, sum)
    for(feature in features) {
      if(count.feature[feature] <= 1) next
      suppressWarnings(tmp.rma.fit <-
                         try(metafor::rma.uni(yi = betas[feature, ind.feature[feature, ]],
                                              sei = sds[feature, ind.feature[feature, ]],
                                              mod = ~.,
                                              data = data.moderator[ind.feature[feature, ], ,
                                                                    drop = FALSE],
                                              method = method,
                                              control = list(threshold = 1e-10,
                                                             maxiter = 1000)),
                             silent = TRUE)) # FIXME
      if("try-error" %in% class(tmp.rma.fit))
        next
      if(is.null(tmp.rma.fit$R2))
        next
      i.result[feature, c("tau2",
                          "se.tau2",
                          "p.tau2",
                          "p.moderator",
                          "I2",
                          "H2",
                          "R2")] <-
        unlist(tmp.rma.fit[c("tau2",
                             "se.tau2",
                             "QEp",
                             "QMp",
                             "I2",
                             "H2",
                             "R2")])
    }
    l.results[[value.exposure]] <- i.result
  }
  results <- Reduce("rbind", l.results)
  results$R2[is.na(results$R2) & !is.na(results$tau2)] <- 0
  return(results)
}