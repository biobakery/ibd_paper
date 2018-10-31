lib.toLoad <- 
  c("magrittr", # for easy operators
    "ggplot2", "cowplot", # for plotting
    "phyloseq" # for phyloseq objects
  )
lib.required <- c(lib.toLoad)
studies <- c("BIDMC-FMT",
             "CS-PRISM",
             "Herfarth_CCFA_Microbiome_3B_combined",
             "HMP2",
             "Jansson_Lamendella_Crohns",
             "jansson_twins_ibd",
             "LSS-PRISM",
             "MucosalIBD",
             "Pouchitis",
             "PROTECT",
             "RISK")
dir_processed <- "../ibd_meta_analysis/"
