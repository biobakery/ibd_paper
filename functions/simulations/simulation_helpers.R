create_batch <- function(nSample_perBatch, nBatch) {
  batch.samples <- sample.int(nSample_perBatch * nBatch)
  batch.samples <- cut(batch.samples, 
                       breaks = seq(0, nSample_perBatch * nBatch, by = nSample_perBatch),
                       labels = 1:nBatch)
  return(batch.samples)
}

create_continuousStructure <- function(nSample_perBatch, nBatch) {
  batch <- create_batch(nSample_perBatch, nBatch)
  score <- runif(nSample_perBatch * nBatch, min = -1, max = 1)
  data.frame(batch = as.factor(batch),
             score = score,
             row.names = paste0("Sample", 1:(nSample_perBatch * nBatch)))
}