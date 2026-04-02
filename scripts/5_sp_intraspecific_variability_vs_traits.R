##Relating species level intraspecific variability to species traits 


library(tidyverse)
library(ggplot2)
####Load intraspecific variability at site level data 
load("data/SpeciesIntraspecific_variance_perSite_Env.RData")


##Calculate the mean instraspecific variability of each species across all food webs 
sp_var_mean <- SpVar_env %>%
  group_by(fish_species, waterbody_type) %>% ##separated by ecosystem type- but this may not be necessary 
  filter(!is.na(sp_site_var_N)) %>% 
  summarise_at(vars(sp_site_var_C, sp_site_var_N), list(mean = mean))

##then you would want to join this to trait data 

