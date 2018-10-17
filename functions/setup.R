setup <- function(dir) {
  # Load constants
  source(paste0(dir, "/assets/constants.R"))
  
  # Load global packages
  for(i.lib in lib.toLoad) library(i.lib, character.only = TRUE)
  
  # Plotting options
  theme_set(theme_bw())
}