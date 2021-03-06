# required packages -------------------------------------------------------

lib.toLoad <- 
  c("magrittr", # for easy operators
    "ggplot2", 
    "phyloseq" # for phyloseq objects
  )
lib.required <- c(lib.toLoad,
                  "tidyverse",
                  "cowplot", 
                  "ggrepel", # for plotting
                  "DT", # for table formatting
                  "hash", "optparse" # for Maaslin2
                  )

# study information -------------------------------------------------------

studies <- c("BIDMC-FMT",
             "CS-PRISM",
             "Herfarth_CCFA_Microbiome_3B_combined",
             "HMP2",
             "Jansson_Lamendella_Crohns",
             "LSS-PRISM",
             "MucosalIBD",
             "Pouchitis",
             "PROTECT",
             "RISK")
studies_longitudinal <- c("Herfarth_CCFA_Microbiome_3B_combined",
                          "HMP2",
                          "Jansson_Lamendella_Crohns",
                          "LSS-PRISM",
                          "PROTECT")

# directory for processed samples -----------------------------------------

dir_processed <- "../ibd_meta_analysis/"
