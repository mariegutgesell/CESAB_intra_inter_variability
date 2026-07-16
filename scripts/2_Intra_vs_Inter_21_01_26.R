################################################################################
# Title:        Add the environmental drivers dataset as to intraspecific dataset 
# Description:  From the environmental dataset, computes:
#               - a more elaborated category for Ecosystem_type
#               - xxxxxxx
#               - xxxxxxx
# Date:         2026-01-21
# Version:      1.0
# Notes:        Any additional information or context
# Authors:       @Arnaud Sentis
################################################################################

######################################################
# 1. SET UP -------
######################################################

## Clear workspace (remove all objects)
rm(list = ls())


## libraries---------
library(dplyr)
library(rfishbase)
library(openxlsx)
library(readxl)
library(tidyr)
library(DHARMa)
library(ggpubr)
library(rprojroot)
library(forcats)
library(FactoMineR)
library(factoextra)
library(cluster)
library(ggrepel)
library(sf)
library(viridis)
library(ggplot2)
library(ggpubr)
library(lubridate)
library(mapview)
library(stringr)
##  Set the working directory to the root of the project ------
#root.dir = find_rstudio_root_file()
#data.dir = paste0(root.dir,'/data')
#script.dir = paste0(root.dir,'/src')
#figures.dir = paste0(root.dir,'/figures')

#setwd(script.dir)

### choice of color palette---------
col_pal<-c("darkgrey","deepskyblue1","deepskyblue2","deepskyblue3","darkolivegreen","darkolivegreen2","darkolivegreen3",
           "coral","coral2","brown1","brown2","brown3","goldenrod1","goldenrod2","goldenrod3")

## import/load datasets and prep them---------
#load("../data/Env.RData") 
load("data/Env_27Mar26.RData") 

##plot out on map
Env_sf <-  st_as_sf(Env_27Mar26,
  coords = c("Latitude",
             "Longitude"),
  crs = 4326,
  remove = FALSE
)
mapview(Env_sf)



#Load data and "normalize" d13C and d15N per site for same "baseline"
##per site: d13C_i=(d13C_i - d13C_min)/(d13C_max-d13C_min)
##per site: d15N_i=(d15N_i - d15N_min)/(d15N_max-d15N_min)

#ALLindiv_June2025 <- read_excel("../data/ALLindiv_June2025.xlsx")
#ALLindiv_June2025 <- read_excel("data/ALLindiv_June2025.xlsx")
##ALLindiv_June2025 <- read_excel("ALLindiv_January2026.xlsx")

df_all <- read_excel("data/FINAL_ALLindiv_February2026.xlsx") 


##clean data
##select only fish, where have C/N data and associated scientific name
#DataFish<-subset(df_all,!is.na(d15N) & !is.na(d13C) & organism_type=="fish" & fish_species != "NA")
##note: when we remove fish_species is NA we are removing any fish that are identified to only to genus or higher level - this can create false diversity and intravar. because we are removing individuals from different speceis from a community .. 

##alternative way to retain taxonomic resolution variation -- new species column to use going forward then is fish_scientific_name
DataFish <- df_all %>%
  filter(!is.na(d15N) & !is.na(d13C) & organism_type=="fish") %>%
  mutate(across(c(fish_species, fish_genus, fish_family, fish_order), ~ na_if(.x, "NA"))) %>%
  mutate(fish_name_level = case_when(
    !is.na(fish_species) ~ "species",
    !is.na(fish_genus) ~ "genus",
    !is.na(fish_family) ~ "family",
    !is.na(fish_order) ~ "order"),
    fish_scientific_name = coalesce(fish_species, fish_genus, fish_family, fish_order)) %>%
  filter(!is.na(fish_scientific_name))



##add species-site identifier column
DataFish$sp_site<-paste(DataFish$fish_scientific_name,DataFish$collection_site_id,sep="_")
##Assign 1 - number of fish
DataFish$num<-1


##normalise data
DataFish<-DataFish %>% 
  mutate(year = str_extract(collection_date, "^\\d{4}") %>% as.numeric()) %>% ##create column of year
  mutate(site_year_code = paste(FWB_id, year, sep = "_")) %>%
  group_by(site_year_code) %>% 
  mutate(d15N_norm = (d15N - min(d15N, na.rm = TRUE))/(max(d15N, na.rm = TRUE)-min(d15N, na.rm = TRUE)),
          d13C_norm = (d13C - min(d13C, na.rm = TRUE))/(max(d13C, na.rm = TRUE)-min(d13C, na.rm = TRUE)))


##Read in cutoff site list, and then select the cutoff/criteria 
cutoff_site_list_all <- read.csv("data/cuttoff_site_lists.csv")  


cutoff_site_list <- cutoff_site_list_all %>%
  mutate(site_year_code = paste(FWB_id, year, sep = "_")) %>%
  filter(min_sp_num == 2) %>% ##select minimum species number
  filter(cutoff == "75pct") ##select cutoff 
##453 sites 

##select site/years from chosen cutoff 
DataFish_final <- DataFish %>%
  filter(site_year_code %in% cutoff_site_list$site_year_code)


##what, why are there 3 food webs without environmental data? 
#I. Intraspecific variation per species------------
#some species have var = NA when only one individual sampled
#we ignore such species in average intra var but still include it in inter var
colnames(DataFish_final)
str(DataFish_final)
#levels(DataFish$Diet)
DataFish_final$collected_sample_length_mm <- as.numeric(DataFish_final$collected_sample_length_mm)



SpVar<-DataFish_final %>% 
  group_by(FWB_id, site_year_code, sp_site,collection_site_id,fish_scientific_name, fish_name_level,
           waterbody_type, 
           #ecosystem_area_km2, ecosystem_width_m,  ##for 3 FWB_id there are two sizes, so i dont think we want to group by this .. not consistent across site
           collection_decimal_longitude, collection_decimal_latitude) %>% 
  summarise(sp_site_mean_N = mean(d15N_norm, na.rm = TRUE),
            sp_site_var_N = var(d15N_norm, na.rm = TRUE),
            sp_site_mean_C = mean(d13C_norm, na.rm = TRUE),
            sp_site_var_C = var(d13C_norm, na.rm = TRUE),
            sp_site_mean_length = mean(collected_sample_length_mm, na.rm = TRUE), ##do we want to only use fork or total length? or group by this? 
            sp_site_var_length = var(collected_sample_length_mm, na.rm = TRUE),
            collection_decimal_longitude = mean(collection_decimal_longitude, na.rm = TRUE),
            collection_decimal_latitude = mean(collection_decimal_latitude, na.rm = TRUE),
            sp_site_num_ind = sum(num), 
            sp_site_body_size_range = max(collected_sample_length_mm) - min(collected_sample_length_mm)) ##number of samples per species

SpVar$VarTot<-SpVar$sp_site_var_N+SpVar$sp_site_var_C
colnames(SpVar)
SpVar$collection_decimal_longitude
library(sf)
library(dplyr)

na_fish <- SpVar %>%
  filter(fish_scientific_name == "NA")
na_fish <- SpVar %>%
  filter(is.na(fish_scientific_name))
##should be 0

test <- SpVar %>%
  select(collection_site_id) %>%
  group_by(collection_site_id) %>%
  count()
#453 sites - so i think okay .. 


############################################################
# Spatial join between fish site data and environmental data
# Purpose: Merge fish variables with environmental variables
#          using spatial proximity (robust to coordinate mismatch)
############################################################
#-----------------------------------------------------------
# 0. Keep a copy of the original data
#-----------------------------------------------------------
SpVar_original <- SpVar

#-----------------------------------------------------------
# 1. Keep only rows with valid coordinates
#-----------------------------------------------------------
SpVar_clean <- SpVar %>%
  filter(
    !is.na(collection_decimal_longitude),
    !is.na(collection_decimal_latitude)
  )

#-----------------------------------------------------------
# 2. Convert to sf
#-----------------------------------------------------------
SpVar_sf <- SpVar_clean %>%
  st_as_sf(
    coords = c("collection_decimal_longitude",
               "collection_decimal_latitude"),
    crs = 4326,
    remove = FALSE
  )
mapview(SpVar_sf)
#-----------------------------------------------------------
# 3. Ensure Env uses same CRS
#-----------------------------------------------------------
Env_sf <- st_transform(Env_27Mar26, 4326)
mapview(Env_sf)
#-----------------------------------------------------------
# 4. Spatial join (nearest environmental site)
#-----------------------------------------------------------
SpVar_env_sf <- st_join(
  SpVar_sf,
  Env_sf,
  join = st_nearest_feature
)

test <- SpVar_env_sf %>%
  select(FWB_id) %>%
  distinct()
#-----------------------------------------------------------
# 5. (Optional) Distance check
#-----------------------------------------------------------
SpVar_m <- st_transform(SpVar_sf, 3857)
Env_m   <- st_transform(Env_sf, 3857)

nearest_id <- st_nearest_feature(SpVar_m, Env_m)

SpVar_env_sf$distance_to_env_m <- as.numeric(
  st_distance(SpVar_m, Env_m[nearest_id, ], by_element = TRUE)
)


############################################################
#  Optional check: how many matches are far away?
#-----------------------------------------------------------
summary(SpVar_env_sf$distance_to_env_m)
hist(SpVar_env_sf$distance_to_env_m)
# Example: flag suspicious matches (> 5 km)
SpVar_env_sf <- SpVar_env_sf %>%
  mutate(
    match_flag = ifelse(distance_to_env_m > 5000, "check", "ok")
  )
sum(SpVar_env_sf$match_flag == "check", na.rm = TRUE)

verify<-SpVar_env_sf %>%
  filter(match_flag == "check") %>%
  select(sp_site,
         collection_site_id,
         distance_to_env_m,
         collection_decimal_latitude,
         collection_decimal_longitude)

#-----------------------------------------------------------
# 6. Drop geometry if you want a regular data.frame
#-----------------------------------------------------------
SpVar_env <- SpVar_env_sf %>%
  st_drop_geometry()

# Final object: SpVar_env


colnames(SpVar_env)
p1bis<-ggplot(SpVar_env,aes(x=fish_scientific_name,y=sp_site_var_N,col=waterbody_type))+geom_point(alpha=0.5) +
  theme_bw()+ 
  geom_smooth()+xlab("abs latitude")+ylab("specioes intraspecific variance")
p1bis


#II. Intraspecific vs. interspecific variation per site with environemental data
colnames(SpVar)
colnames(SpVar_env)
SpVar_env$num<-1

SiteVar <- SpVar_env %>%
  group_by(
    FWB_id,
    site_year_code,
    collection_site_id,
    waterbody_type,
   collection_decimal_latitude,
  collection_decimal_longitude,
#  ecosystem_area_km2, ##for 3 FWB there are two different ecosystem sizes (see below), but since not using this i dont think then dont need to group by this here ... 
#  ecosystem_width_m
  ) %>%
  summarise(
    
    # ----------------------------
    # Fish community statistics
    # ----------------------------
    site_interspe_var_N = var(sp_site_mean_N, na.rm = TRUE), ##site level interspecific variation: variation in mean dN between species
    site_intraspe_var_N = mean(sp_site_var_N, na.rm = TRUE), ##site level intraspecific variation: mean within species variation across all species at the site
    
    site_interspe_var_C = var(sp_site_mean_C, na.rm = TRUE),
    site_intraspe_var_C = mean(sp_site_var_C, na.rm = TRUE),
    
    site_nbspe = sum(num, na.rm = TRUE), ##species richness
    
   collection_decimal_longitude = mean(collection_decimal_longitude, na.rm = TRUE),
    collection_decimal_latitude  = mean(collection_decimal_latitude,  na.rm = TRUE),
    
    site_mean_sample_id = mean(sp_site_num_ind, na.rm = TRUE), ##mean number of individuals in sample
    site_min_sample_id  = min(sp_site_num_ind,  na.rm = TRUE), ##min number of individuals in sample 
    
    # ----------------------------
    # Environmental variables
    # (identical within site)
    # ----------------------------
   # Food_web_ID        = first(Food web_ID),
    Climate_zone       = first(Climate_zone),
    Type               = first(Type),
    Latitude           = first(Latitude),
    Longitude          = first(Longitude),
    Size               = first(Size),
    dis_r_sv           = first(dis_r_sv),
    #npp_mean           = first(npp_mean),
    TP                 = first(TP),
    TN                 = first(TN),
    Distance_km        = first(Distance_km),
    temp               = first(temp),
    hft                = first(hft),
    prec               = first(prec),
    pop                = first(pop),
    crp                = first(crp),
    urb                = first(urb),
    upstr              = first(upstr),
    regul              = first(regul),
    pop_den            = first(pop_den),
  #  npp_class          = first(npp_class),
    size_class         = first(size_class),
    size_z_scored      = first(size_z_scored),
    hydro_dis_z_scored = first(hydro_dis_z_scored),
    TP_class           = first(TP_class),
    Climate_zone_e     = first(Climate_zone_e),
    Climate_zone_e2    = first(Climate_zone_e2),
    distance_to_env_m  = first(distance_to_env_m),
    match_flag         = first(match_flag),
    
    .groups = "drop"
  )

##check for any duplicates
#test <- SiteVar %>%
#  select(FWB_id, collection_site_id, waterbody_type) %>%
#  group_by(FWB_id, collection_site_id, waterbody_type) %>%
#  count()
##FWB_0175, FWB_0215, FWB_0541 - do we know why in DataFish there are two different ecosystem_area_m or ecosystem_width_m?
##NOTE: these issues have been resolved by not grouping by ecosystem type - now only have 1 value per site



##Calculate proportion of intraspecific variation from C and N variance
SiteVar$propintraspecific_N<-SiteVar$site_intraspe_var_N/(SiteVar$site_interspe_var_N+SiteVar$site_intraspe_var_N)
SiteVar$propintraspecific_C<-SiteVar$site_intraspe_var_C/(SiteVar$site_interspe_var_C+SiteVar$site_intraspe_var_C)
SiteVar$propintraspecific_Total<-(SiteVar$site_intraspe_var_N+SiteVar$site_intraspe_var_C)/(SiteVar$site_interspe_var_N+SiteVar$site_intraspe_var_N+SiteVar$site_interspe_var_C+SiteVar$site_intraspe_var_C)


##save - NOTE: remember, this is filtered by cutoff
#save(SiteVar_final, file = "data/Intraspecific_contribution_perSite_Env.RData")
#save(SpVar_final, file = "data/SpeciesIntraspecific_variance_perSite_Env.RData")
write.csv(SiteVar, "data/SiteVar_2sp_75cutoff.csv", row.names = FALSE) ##update file name if using different cutoffs
write.csv(SpVar, "data/SpVar_2sp_75cutoff.csv", row.names = FALSE) 


##how many species have at least 10 species?
sites_10sp <- SiteVar %>%
  filter(site_nbspe >= 8)
##potentially focus 

##Calculate the proporton of total individuals in a food web that the species with less than 3 replicates make up (i.e., calculating rare species)
##Okay so -- current cutoff, 2 species, 75% -- for those 25% w/ less than 3 replicates, what proportion of total number of individuals in that food web are they ? 
cutoff_site_list_25 <- cutoff_site_list %>%
  filter(min_num_ind_per_sample < 3)

datafish_25 <- DataFish %>%
  filter(site_year_code %in% cutoff_site_list_25$site_year_code)


datafish_25_total <- datafish_25 %>%
  group_by(site_year_code) %>%
  count() %>%
  rename(num_ind_total = "n")


datafish_25_sp_total <- datafish_25 %>%
  group_by(site_year_code, fish_scientific_name) %>%
  count() %>%
  rename(num_ind_per_species = "n") %>%
  left_join(datafish_25_total, by = "site_year_code") %>%
  mutate(percent_total_ind = num_ind_per_species/num_ind_total*100)



rare_sp <- datafish_25_sp_total %>%
  filter(num_ind_per_species < 3)

rare_sp_mean <- rare_sp %>%
  ungroup() %>%
  summarise(mean_percent_total_ind = mean(percent_total_ind),
            sd_percent_total_ind = sd(percent_total_ind))

ggplot(rare_sp, aes(x = percent_total_ind)) +
  geom_histogram()

rare_sp %>%
  group_by(site_year_code) %>%
  summarise(percent_total_ind_sum = sum(percent_total_ind)) %>%
  ggplot(aes(x = percent_total_ind_sum)) +
  geom_histogram()


