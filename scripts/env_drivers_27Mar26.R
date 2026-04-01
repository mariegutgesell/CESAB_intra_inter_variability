################################################################################
# Title:        Generates the environmental drivers dataset as predictors for PPMR 
# Description:  From the environmental dataset, computes:
#               - a more elaborated category for Ecosystem_type
#               - xxxxxxx
#               - xxxxxxx
# Date:         2026-01-20
# Version:      1.0
# Notes:        Any additional information or context
# Authors:       @Marie-Elodie Perga  @Arnaud Sentis
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
library(dplyr)
library(DHARMa)
library(ggpubr)
library(rprojroot)
library(forcats)
library(FactoMineR)
library(factoextra)
library(cluster)
library(ggrepel)
library(sf)
library(dplyr)
library(viridis)



##  Set the working directory to the root of the project ------
root.dir = find_rstudio_root_file()
#root.dir =("C:/Users/asentis/INRAE/AAAprojects/Food_webs_CESAB/Intraspecific")
data.dir = paste0(root.dir,'/data')
script.dir = paste0(root.dir,'/scripts')
figures.dir = paste0(root.dir,'/figures')

setwd(script.dir)

### choice of color palette---------
col_pal<-c("darkgrey","deepskyblue1","deepskyblue2","deepskyblue3","darkolivegreen","darkolivegreen2","darkolivegreen3",
           "coral","coral2","brown1","brown2","brown3","goldenrod1","goldenrod2","goldenrod3")

## import/load datasets and prep them---------
#EnvData_len <- read_excel("../data/FOODWEBS_Individual_EnvironmentalData_28Aug25.xlsx",sheet = "Lentic")
EnvData_len <- read_excel("../data/FOODWEBS_Individual_EnvironmentalData_23Mar26.xlsx",sheet = "Lentic environment")
##nutrient data
df_all <- read_excel("../data/FINAL_ALLindiv_February2026.xlsx")  
str(df_all) ##issues with a lot of variables as character instead of numeric.
colnames(df_all)

##keep only Fodd web id and longitude and latitude
df_clean <- df_all %>%
  group_by(FWB_id) %>%
  summarise(
    latitude = mean(collection_decimal_latitude, na.rm = TRUE),
    longitude = mean(collection_decimal_longitude, na.rm = TRUE)
  )

##add longitude and latitude to EnvData
EnvData_len <- EnvData_len %>%
  left_join(df_clean, by = "FWB_id")


#load nutrient data
NutData <- read.csv("../data/Global_Nutrient_Data.csv")
str(NutData)

##to add the nutrient data from the nearest location
# sites
sites_sf <- st_as_sf(
  EnvData_len,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

# nutriments
nut_sf <- st_as_sf(
  NutData,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)

##join nutrient data
EnvData_len_nut <- st_join(
  sites_sf,
  nut_sf,
  join = st_nearest_feature
)

##add distance in degrees
EnvData_len_nut$Distance_deg <- st_distance(
  sites_sf,
  nut_sf[st_nearest_feature(sites_sf, nut_sf), ],
  by_element = TRUE
) %>% as.numeric()


##add distance in kilometers
sites_sf_m <- st_transform(sites_sf, 3857)
nut_sf_m   <- st_transform(nut_sf, 3857)

EnvData_len_nut$Distance_km <- as.numeric(
  st_distance(
    sites_sf_m,
    nut_sf_m[st_nearest_feature(sites_sf_m, nut_sf_m), ],
    by_element = TRUE
  )
) / 1000




##lotic systems
#EnvData_lot <- read_excel("../data/FOODWEBS_Individual_EnvironmentalData_28Aug25.xlsx",sheet = "Lotic")
EnvData_lot <- read_excel("../data/FOODWEBS_Individual_EnvironmentalData_23Mar26.xlsx",sheet = "Lotic environment")
colnames(EnvData_lot)

##add longitude and latitude to EnvData
EnvData_lot <- EnvData_lot %>%
  left_join(df_clean, by = "FWB_id")


# sites
sites_sf_lot <- st_as_sf(
  EnvData_lot,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)
##join nutrient data
EnvData_lot_nut <- st_join(
  sites_sf_lot,
  nut_sf,
  join = st_nearest_feature
)

##add distance in degrees
EnvData_lot_nut$Distance_deg <- st_distance(
  sites_sf_lot,
  nut_sf[st_nearest_feature(sites_sf_lot, nut_sf), ],
  by_element = TRUE
) %>% as.numeric()
colnames(EnvData_lot_nut)

##add distance in kilometers
sites_sf_m_lot <- st_transform(sites_sf_lot, 3857)

EnvData_lot_nut$Distance_km <- as.numeric(
  st_distance(
    sites_sf_m_lot,
    nut_sf_m[st_nearest_feature(sites_sf_m_lot, nut_sf_m), ],
    by_element = TRUE
  )
) / 1000

save(EnvData_len_nut, file = "../data/EnvData_len_nut_March26.RData")
save(EnvData_lot_nut, file = "../data/EnvData_lot_nut_March26.RData")



####load the data
# load("../data/EnvData_len_nut_March26.RData") 
# load("../data/EnvData_lot_nut_March26.RData") 
EnvData_lot_nut$ria_ha_csu

EnvData_lot_nut$dis_r_sv=(EnvData_lot_nut$dis_m3_pmx - EnvData_lot_nut$dis_m3_pmn)/(0.001*EnvData_lot_nut$riv_tc_csu)

colnames(EnvData_lot_nut)
Env_river_b<-EnvData_lot_nut %>% 
  dplyr::select ('Food web_ID'=FWB_id,
                 #Region,
                 Climate_zone=clz_cl_cmj,
                 Type=waterbody_type,
                 #Country=Country,
                 Latitude=latitude,
                 Longitude=longitude,
                 Size=ria_ha_csu,
                 #Size2=ria_ha_csu,
                 dis_r_sv=dis_r_sv,
                 #npp_mean=NPP,
                 #DRP=DRP,
                 TP=TP,#productivity as P yield
                 TN=TN,
                 Distance_km =Distance_km,
                 upstr=ria_ha_usu,
                 pop=pop_ct_csu,
                 pop_den=ppd_pk_cav,
                 temp=tmp_dc_cyr,
                 hft=hft_ix_c09, #human footprint
                 prec=pre_mm_cyr,
                 crp=crp_pc_cse,
                 urb=urb_pc_cse,
                 regul=dor_pc_pva,
                ) %>%
  mutate(Size=log(Size),
         dis_r_sv=log(dis_r_sv+0.00001),#avoid zeros
         pop=log(pop+0.001),
         upstr=log(upstr+0.00001),
         pop_den=log(pop_den+0.00001),
         TP=log(TP))%>%
  # mutate(npp_class=case_when(
  #   npp_mean<median(npp_mean,na.rm=TRUE) ~ "low",
  #   npp_mean>median(npp_mean,na.rm=TRUE) ~ "high"
  # ))%>%
  mutate(size_class=case_when(
    Size<median(Size,na.rm=TRUE) ~ "small",
    Size>median(Size,na.rm=TRUE) ~ "large"
  ))%>%
  filter(Size!= "-Inf") %>%
  mutate(size_z_scored=scale(Size, center = TRUE, scale = TRUE))%>%#scale size
  mutate(hydro_dis_z_scored=scale(dis_r_sv, center = TRUE, scale = TRUE)) %>%#scale flow variability
  mutate(TP_class=case_when(
    TP<median(TP,na.rm=TRUE) ~ "low_prod",
    TP>median(TP,na.rm=TRUE) ~ "high_prod"))
mean(Env_river_b$Size)
Env_river_b$size_z_scored
##lentic system
#flow variability
#compute flow variability
EnvData_len_nut$dis_r_sv=(EnvData_len_nut$dis_m3_pmx - EnvData_len_nut$dis_m3_pmn)/(EnvData_len_nut$Vol_total)

Env_lake_b<-EnvData_len_nut %>% dplyr::select (
  'Food web_ID'=FWB_id,
  #Region,
  Climate_zone=clz_cl_lmj,
  Type=waterbody_type,
  #Country=Country,
  Latitude=latitude,
  Longitude=longitude,
  Size=Lake_area,#lake size in km2
  dis_r_sv=dis_r_sv,#flow variability
  #npp_mean=NPP,
  #DRP=DRP,
  TP=TP,#productivity as P yield
  TN=TN,
  Distance_km =Distance_km,
  temp=tmp_dc_lyr,
  hft=hft_ix_u09,#human footprint
  prec=pre_mm_lyr,
  pop=pop_ct_vsu,
  crp=crp_pc_vse,
  urb=urb_pc_vse,
  upstr=ria_ha_usu,
  regul=dor_pc_pva,  pop_den=ppd_pk_vav)%>%
  mutate(Size=log(Size),
         dis_r_sv=log(dis_r_sv+0.00001),
         pop=log(pop+0.001),
         upstr=log(upstr+0.00001),
         pop_den=log(pop_den+0.00001),
         TP=log(TP))%>%
  # mutate(npp_class=case_when(
  #   npp_mean<median(npp_mean,na.rm=TRUE) ~ "low",
  #   npp_mean>median(npp_mean,na.rm=TRUE) ~ "high"
  # ))%>%
  mutate(size_class=case_when(
    Size<median(Size,na.rm=TRUE) ~ "small",
    Size>median(Size,na.rm=TRUE) ~ "large"
  ))%>%
  mutate(size_z_scored=scale(Size, center = TRUE, scale = TRUE)) %>%#scale size
  mutate(hydro_dis_z_scored=scale(dis_r_sv, center = TRUE, scale = TRUE))%>%#scale flow variability
  mutate(TP_class=case_when(
    TP<median(TP,na.rm=TRUE) ~ "low_prod",
    TP>median(TP,na.rm=TRUE) ~ "high_prod"))

Env<-rbind(Env_lake_b,Env_river_b)
Env_27Mar26<-Env %>% 
  mutate(Climate_zone_e=as.character(as.factor(Climate_zone))) %>%
  mutate(Climate_zone_e = fct_recode(Climate_zone_e,
                                     "Cold and wet"="5",
                                     "Extremely cold and mesic"="6",
                                     "Cold and mesic" = "7",
                                     "Cool temperate and dry"="8",
                                     "Cool temperate and xeric"="9",
                                     "Cool temperate and moist"="10",
                                     "Warm temperate and mesic"="11",
                                     "Warm temperate and xeric"="12",
                                     "Hot and mesic"="13",
                                     "Hot and dry"="14",
                                     "Hot and arid"="15",
                                     "Extremely hot and arid"="16",
                                     "Extremely hot and xeric"="17",
                                     "Extremely hot and moist"="18"))%>%
  mutate(Climate_zone_e2 = fct_recode(Climate_zone_e,#aggregate some categories of cliamte
                                      "Cold and wet/mesic"= "Cold and wet",
                                      "Cold and wet/mesic"="Extremely cold and mesic",
                                      "Cold and wet/mesic"="Cold and mesic" ,
                                      "Cool temperate and dry/xeric"="Cool temperate and dry",
                                      "Cool temperate and dry/xeric"="Cool temperate and xeric",
                                      "Cool and moist"="Cool temperate and moist",
                                      "Warm temperate"="Warm temperate and mesic",
                                      "Warm temperate"="Warm temperate and xeric",
                                      "Hot and moist"="Hot and mesic",
                                      "Hot and dry"="Hot and dry",
                                      "Hot and dry"="Hot and arid",
                                      "Hot and dry"="Extremely hot and arid",
                                      "Hot and dry"="Extremely hot and xeric",
                                      "Hot and moist"="Extremely hot and moist"))  






#save(Env_river_b,file = "Env_river_b.RData")
#save(Env_lake_b,file = "Env_lake_b.RData")
save(Env_27Mar26,file = "../data/Env_27Mar26.RData")

load("../data/Env_27Mar26.RData")



ggplot(Env, aes(x = TP, y = Distance_km)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ Type) +
  geom_smooth(method='lm',formula= y~x)+
  theme_bw() +
  labs(
    x = "Total phosphorus (TP)",
    y = "Distance")

Env_far <- Env %>% filter(Distance_km > 1000)


ggplot(Env, aes(y = (TP), x = log(Distance_km))) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE) +
  facet_wrap(~ Type) +
  scale_color_viridis_c(option = "magma", name = "Distance (km)") +
  theme_bw() +
  labs(
    y = "Total phosphorus (TP)",
    x = "log(Distance_km)")


Env <- Env %>%
  mutate(
    Dist_class = cut(
      Distance_km,
      breaks = quantile(Distance_km, probs = seq(0, 1, 0.25), na.rm = TRUE),
      include.lowest = TRUE,
      labels = c("Q1 (close)", "Q2", "Q3", "Q4 (far)")
    )
  )
quantile(Env$Distance_km, probs = seq(0, 1, 0.25), na.rm = TRUE)

ggplot(Env, aes(y = TP, x = Dist_class)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  facet_wrap(Dist_class~ Type) +
  #scale_color_viridis_c(option = "magma", name = "Distance (km)") +
  theme_bw() +
  labs(
    x = "Total phosphorus (TP)",
    y = "Mean NPP",
    title = "NPP vs TP across ecosystems"
  )
