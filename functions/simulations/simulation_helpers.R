create_batch <- function(nSample_perBatch, nBatch) {
  batch.samples <- sample.int(nSample_perBatch * nBatch)
  batch.samples <- cut(batch.samples, 
                       breaks = seq(0, nSample_perBatch * nBatch, by = nSample_perBatch),
                       labels = 1:nBatch)
  return(batch.samples)
}

create_imbalancedBatches <- function(nSample_perBatch, nBatch, imbalance) {
  if(nBatch != 2) 
    stop("Currently only support 2 batches!")
  batch <- create_batch(nSample_perBatch, nBatch)
  exposure <- rep(NA_real_, length = length(batch))
  exposure[batch == 1] <- rbinom(n = nSample_perBatch, size = 1, prob = 0.5 - imbalance)
  exposure[batch == 2] <- rbinom(n = nSample_perBatch, size = 1, prob = 0.5 + imbalance)
  data.frame(batch = as.factor(batch),
             exposure = as.factor(exposure),
             row.names = paste0("Sample", 1:(nSample_perBatch * nBatch)))
}

create_continuousStructure <- function(nSample_perBatch, nBatch) {
  batch <- create_batch(nSample_perBatch, nBatch)
  score <- runif(nSample_perBatch * nBatch, min = -1, max = 1)
  data.frame(batch = as.factor(batch),
             score = score,
             row.names = paste0("Sample", 1:(nSample_perBatch * nBatch)))
}