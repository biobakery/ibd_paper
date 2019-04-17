setup <- function(dir) {
  # Load constants
  source(paste0(dir, "/assets/constants.R"))
  
  # Load study information
  source(paste0(dir, "/assets/tb_1.R"))
  
  # Load global packages
  for(i.lib in lib.toLoad) library(i.lib, character.only = TRUE)
  
  # Plotting options
  theme_set(theme_bw() +
              theme(panel.grid.major = element_blank(), 
                    panel.grid.minor = element_blank(),
                    title = element_text(size = 16),
                    axis.title = element_text(size = 16),
                    axis.text = element_text(size = 14),
                    legend.text = element_text(size = 14),
                    strip.text = element_text(size = 14)))
}