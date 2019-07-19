library(magrittr)
ii <- commandArgs(trailingOnly = TRUE) %>% as.integer()

load("/n/hutlab11_nobackup/users/syma/ibd_paper/results/simulations/lm.meta/tb_sim.RData")
tb_sim_subset <- tb_sim %>% 
  dplyr::filter(nSample_perBatch == 100,
                nMicrobe == 200,
                spikeMicrobes == 0.05)
i <- tb_sim_subset$i[ii]
mat_otu <- read.table(paste0("/n/hutlab11_nobackup/users/syma/ibd_paper/results/simulations/lm.meta/sparseDOSSA_sets/",
                      i, ".tsv"),
                      header = TRUE,
                      sep = "\t",
                      row.names = 1) %>%
  as.matrix()

print(i)

# i_simSetup <- tb_sim[i, ]
# df_metadata <- i_simSetup$df_metadata[[1]] %>% 
#   dplyr::transmute(main = as.numeric(as.character(exposure)),
#                    confounder = NA_real_,
#                    batch = as.numeric(as.character(batch)))
# sSet <- SummarizedExperiment::SummarizedExperiment(
#   assays = list(mat_otu),
#   colData = df_metadata
# )
# 
# fit.BDMMA <- BDMMAcorrect::BDMMA(Microbiome_dat = sSet)
# save(sSet, file = paste0("/n/hutlab11_nobackup/users/syma/ibd_paper/results/simulations/lm.meta/BDMMA/",
#                          i, ".RData"))
