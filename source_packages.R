# Source file: packages used in wolf pva project

####--------------------------------------------####
## PACKAGES
####--------------------------------------------####

mylibraries <- c("dplyr", "ggeffects", "vegan", "tibble",
                 "tidyr", "forcats", "ggplot2", "cowplot", 
                 "janitor", "scico", "sf", "raster", 
                 "png", "grid", "ggspatial", "ggrepel")

for (i in 1:length(mylibraries)) {
  if(mylibraries[i] %in% rownames(installed.packages()) == FALSE) {install.packages(mylibraries[i])}
}
lapply(mylibraries, require, character.only = TRUE)
