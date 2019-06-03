setup <- function(dir) {
  # Load constants
  source(paste0(dir, "/assets/constants.R"))
  
  # Load study information
  source(paste0(dir, "/assets/tb_1.R"))
  
  # Load global packages
  for(i.lib in lib.toLoad) library(i.lib, character.only = TRUE)
  
  # set ggplot theme
  smar::set_ggplot()
}