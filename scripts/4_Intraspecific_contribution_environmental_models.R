##Relating intrasp. variability contribution to environmental drivers

library(tidyverse)
library(cowplot)
library(corrplot)
library(tidyverse)
library(readxl)
library(mapview)
library(sf)
library(leaflet)
library(scales)
library(ggpubr)
library(ggplot2)
library(rnaturalearth)
library(ggspatial)
library(lme4)
library(lmerTest)
library(MuMIn)
library(visreg)
library(car)


##Load intraspecific variability at site level data 
load("data/Intraspecific_contribution_perSite_Env.RData")

##checking duplicate site names/lat longs -- don't worry about for now, come back to once using final data
##sometimes have same lat/long with 2 names, and sometimes have same site name with different lat/longs 
latlong1 <- SiteVar %>%
  select(collection_decimal_latitude, collection_decimal_longitude) %>%
  unique()

sites1 <- SiteVar %>%
  select(collection_site_id) %>%
  group_by(collection_site_id) %>%
  count()

##productivity = TP
##type = ecosystem type (lake vs. stream)
##human ftp = hft
##size = size_z_scored
##flow variability = hydro_dis_z_scored
##climate = temp , climate zone 

##potential random effects: species richness, C:N, tissue type 
##Make correlation plot of numerical explanatory variables
exp_var <- SiteVar%>%
  filter(Size != "-Inf") %>%
  ungroup() %>%
  select(TP, size_z_scored, hydro_dis_z_scored, hft)# %>%
#  unique() #
str(exp_var)

duplicated_exp_variables <- SiteVar %>%
  select(collection_site_id,collection_decimal_longitude, collection_decimal_latitude,
         waterbody_type, Climate_zone_e2, TP,
         size_z_scored, hydro_dis_z_scored, hft) %>%
  filter(
    duplicated(select(., waterbody_type, Climate_zone_e2, TP,
                      size_z_scored, hydro_dis_z_scored, hft)) |
      duplicated(select(., waterbody_type, Climate_zone_e2, TP,
                        size_z_scored, hydro_dis_z_scored, hft), fromLast = TRUE)
  )

##so there are a fair amount of duplicated variables, likely because sites are close togehter  - need to check this carefully when have final data
##also, some points are in the ocean, again check if this issue remains with final data 
sites_coord <- st_as_sf(duplicated_exp_variables, coords = c("collection_decimal_longitude", "collection_decimal_latitude"), crs = 4326)
site_map <- mapview(sites_coord, map.types = "Esri.NatGeoWorldMap", legend = FALSE, cex = 2)
site_map

corr_mat <- cor(exp_var, use = "pairwise.complete.obs")

# visualize
corrplot(corr_mat, method = "color", type = "upper",  addCoef.col = "black",
         tl.col = "black", tl.srt = 45)


##Relationships between C/N prop variatibility and species richness
SiteVar <- SiteVar %>%
  mutate(site_nbspe_log = log(site_nbspe))
C_sp_beta <- betareg(propintraspecific_C ~ site_nbspe_log, data = SiteVar)
summary(C_sp_beta)

# create sequence for TP
newdat <- data.frame(
  site_nbspe_log = seq(min(SiteVar$site_nbspe_log, na.rm = TRUE),
           max(SiteVar$site_nbspe_log, na.rm = TRUE),
           length.out = 100)
)

# predictions
newdat$pred <- predict(C_sp_beta, newdata = newdat, type = "response")

ggplot(SiteVar, aes(x =site_nbspe_log , y = propintraspecific_C)) +
  geom_point() +
  geom_line(data = newdat, aes(x = site_nbspe_log, y = pred), linewidth = 1)

N_sp_beta <- betareg(propintraspecific_N ~ site_nbspe_log, data = SiteVar)
summary(N_sp_beta)

# create sequence for TP
newdat <- data.frame(
  site_nbspe_log = seq(min(SiteVar$site_nbspe_log, na.rm = TRUE),
                       max(SiteVar$site_nbspe_log, na.rm = TRUE),
                       length.out = 100)
)

# predictions
newdat$pred <- predict(N_sp_beta, newdata = newdat, type = "response")

ggplot(SiteVar, aes(x =site_nbspe_log , y = propintraspecific_N)) +
  geom_point() +
  geom_line(data = newdat, aes(x = site_nbspe_log, y = pred), linewidth = 1)




##Simple Univariate Models 
##total phosphorous
library(betareg)

C_TP_beta <- betareg(propintraspecific_C ~ TP + site_nbspe_log, data = SiteVar)
summary(C_TP_beta)

# create sequence for TP
newdat <- data.frame(
  TP = seq(min(SiteVar$TP, na.rm = TRUE),
           max(SiteVar$TP, na.rm = TRUE),
           length.out = 100),
  site_nbspe_log = mean(SiteVar$site_nbspe_log, na.rm = TRUE)  # hold constant
)

# predictions
newdat$pred <- predict(C_TP_beta, newdata = newdat, type = "response")

ggplot(SiteVar, aes(x = TP, y = propintraspecific_C)) +
  geom_point() +
  geom_line(data = newdat, aes(x = TP, y = pred), linewidth = 1)


##regular linear regressions
C_TP <- lm(propintraspecific_C ~ TP, data = SiteVar)
summary(C_TP)
c_tp_1 <- ggplot(SiteVar, aes(x = TP, y = propintraspecific_C)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()
c_tp_1
c_tp_2 <- ggplot(SiteVar, aes(x = TP, y = propintraspecific_C, group = Climate_zone_e2, color = Climate_zone_e2)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()
c_tp_2

c_tp_3 <- ggplot(SiteVar, aes(x = TP, y = propintraspecific_C, group = waterbody_type, color = waterbody_type)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()
c_tp_3

N_TP_beta <- betareg(propintraspecific_N ~ TP + site_nbspe_log, data = SiteVar)
summary(N_TP_beta)

# create sequence for TP
newdat <- data.frame(
  TP = seq(min(SiteVar$TP, na.rm = TRUE),
           max(SiteVar$TP, na.rm = TRUE),
           length.out = 100),
  site_nbspe = mean(SiteVar$site_nbspe, na.rm = TRUE)  # hold constant
)

# predictions
newdat$pred <- predict(N_TP_beta, newdata = newdat, type = "response")

ggplot(SiteVar, aes(x = TP, y = propintraspecific_N)) +
  geom_point() +
  geom_line(data = newdat, aes(x = TP, y = pred), linewidth = 1)



N_TP <- lm(propintraspecific_N ~ TP, data = SiteVar)
summary(N_TP)
n_tp_1 <- ggplot(SiteVar, aes(x = TP, y = propintraspecific_N)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()
n_tp_1

n_tp_2 <- ggplot(SiteVar, aes(x = TP, y = propintraspecific_N, group = Climate_zone_e2, color = Climate_zone_e2)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()
n_tp_2

n_tp_3 <- ggplot(SiteVar, aes(x = TP, y = propintraspecific_N, group = waterbody_type, color = waterbody_type)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()
n_tp_3

ggarrange(c_tp_1, c_tp_2, c_tp_3, nrow = 3, ncol = 1)
ggarrange(n_tp_1, n_tp_2, n_tp_3, nrow = 3, ncol = 1)

##include interaction for TP and climate 
##so would want a climate zone * TP interaction 
##human footprint
C_hft_beta <- betareg(propintraspecific_C ~ hft + site_nbspe_log, data = SiteVar)
summary(C_hft_beta)

# create sequence for TP
newdat <- data.frame(
  hft = seq(min(SiteVar$hft, na.rm = TRUE),
           max(SiteVar$hft, na.rm = TRUE),
           length.out = 100),
  site_nbspe = mean(SiteVar$site_nbspe, na.rm = TRUE)  # hold constant
)

# predictions
newdat$pred <- predict(C_hft_beta, newdata = newdat, type = "response")

ggplot(SiteVar, aes(x = hft, y = propintraspecific_C)) +
  geom_point() +
  geom_line(data = newdat, aes(x = hft, y = pred), linewidth = 1)



N_hft_beta <- betareg(propintraspecific_N ~ hft + site_nbspe_log, data = SiteVar)
summary(N_hft_beta)

# create sequence for TP
newdat <- data.frame(
  hft = seq(min(SiteVar$hft, na.rm = TRUE),
            max(SiteVar$hft, na.rm = TRUE),
            length.out = 100),
  site_nbspe = mean(SiteVar$site_nbspe, na.rm = TRUE)  # hold constant
)

# predictions
newdat$pred <- predict(N_hft_beta, newdata = newdat, type = "response")

ggplot(SiteVar, aes(x = hft, y = propintraspecific_N)) +
  geom_point() +
  geom_line(data = newdat, aes(x = hft, y = pred), linewidth = 1)






C_hft <- lm(propintraspecific_C ~ hft, data = SiteVar)
summary(C_hft)
c_hft_1 <- ggplot(SiteVar, aes(x = hft, y = propintraspecific_C)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()
c_hft_1
c_hft_2 <- ggplot(SiteVar, aes(x = hft, y = propintraspecific_C, group = Climate_zone_e2, color = Climate_zone_e2)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()

c_hft_3 <- ggplot(SiteVar, aes(x = hft, y = propintraspecific_C, group = waterbody_type, color = waterbody_type)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()

N_hft <- lm(propintraspecific_N ~ hft, data = SiteVar)
summary(N_hft)
n_hft_1 <- ggplot(SiteVar, aes(x = hft, y = propintraspecific_N)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()
n_hft_1
n_hft_2 <- ggplot(SiteVar, aes(x = hft, y = propintraspecific_N,  group = Climate_zone_e2, color = Climate_zone_e2)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()

n_hft_3 <- ggplot(SiteVar, aes(x = hft, y = propintraspecific_N,  group = waterbody_type, color = waterbody_type)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()

ggarrange(n_hft_1, n_hft_2, n_hft_3, nrow = 3, ncol = 1)
ggarrange(c_hft_1, c_hft_2, c_hft_3, nrow = 3, ncol = 1)

##flow variability
C_fv_beta <- betareg(propintraspecific_C ~ hydro_dis_z_scored + site_nbspe_log, data = SiteVar)
summary(C_fv_beta)

# create sequence for TP
newdat <- data.frame(
  hydro_dis_z_scored = seq(min(SiteVar$hydro_dis_z_scored, na.rm = TRUE),
            max(SiteVar$hydro_dis_z_scored, na.rm = TRUE),
            length.out = 100),
  site_nbspe = mean(SiteVar$site_nbspe, na.rm = TRUE)  # hold constant
)

# predictions
newdat$pred <- predict(C_fv_beta, newdata = newdat, type = "response")

ggplot(SiteVar, aes(x = hydro_dis_z_scored, y = propintraspecific_C)) +
  geom_point() +
  geom_line(data = newdat, aes(x = hydro_dis_z_scored, y = pred), linewidth = 1)

N_fv_beta <- betareg(propintraspecific_N ~ hydro_dis_z_scored + site_nbspe, data = SiteVar)
summary(N_fv_beta)

# create sequence for TP
newdat <- data.frame(
  hydro_dis_z_scored = seq(min(SiteVar$hydro_dis_z_scored, na.rm = TRUE),
                           max(SiteVar$hydro_dis_z_scored, na.rm = TRUE),
                           length.out = 100),
  site_nbspe = mean(SiteVar$site_nbspe, na.rm = TRUE)  # hold constant
)

# predictions
newdat$pred <- predict(N_fv_beta, newdata = newdat, type = "response")

ggplot(SiteVar, aes(x = hydro_dis_z_scored, y = propintraspecific_N)) +
  geom_point() +
  geom_line(data = newdat, aes(x = hydro_dis_z_scored, y = pred), linewidth = 1)



C_fv <- lm(propintraspecific_C ~ hydro_dis_z_scored, data = SiteVar)
summary(C_fv)
c_fv_1 <- ggplot(SiteVar, aes(x = hydro_dis_z_scored, y = propintraspecific_C)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()
c_fv_1

c_fv_2 <- ggplot(SiteVar, aes(x = hydro_dis_z_scored, y = propintraspecific_C,  group = Climate_zone_e2, color = Climate_zone_e2)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()

c_fv_3 <- ggplot(SiteVar, aes(x = hydro_dis_z_scored, y = propintraspecific_C,  group = waterbody_type, color = waterbody_type)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()



N_fv <- lm(propintraspecific_N ~ hydro_dis_z_scored, data = SiteVar)
summary(N_fv)
n_fv_1 <- ggplot(SiteVar, aes(x = hydro_dis_z_scored, y = propintraspecific_N)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()
n_fv_1

n_fv_2 <- ggplot(SiteVar, aes(x = hydro_dis_z_scored, y = propintraspecific_N,  group = Climate_zone_e2, color = Climate_zone_e2)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()

n_fv_3 <- ggplot(SiteVar, aes(x = hydro_dis_z_scored, y = propintraspecific_N,  group = waterbody_type, color = waterbody_type)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()


ggarrange(n_fv_1, n_fv_2, n_fv_3, nrow = 3, ncol = 1)
ggarrange(c_fv_1, c_fv_2, c_fv_3, nrow = 3, ncol = 1)

##ecosystem size
C_es_beta <- betareg(propintraspecific_C ~ size_z_scored + site_nbspe, data = SiteVar)
summary(C_es_beta)

# create sequence for TP
newdat <- data.frame(
  size_z_scored = seq(min(SiteVar$size_z_scored, na.rm = TRUE),
                           max(SiteVar$size_z_scored, na.rm = TRUE),
                           length.out = 100),
  site_nbspe = mean(SiteVar$site_nbspe, na.rm = TRUE)  # hold constant
)

# predictions
newdat$pred <- predict(C_es_beta, newdata = newdat, type = "response")

ggplot(SiteVar, aes(x = size_z_scored, y = propintraspecific_C)) +
  geom_point() +
  geom_line(data = newdat, aes(x = size_z_scored, y = pred), linewidth = 1)



N_es_beta <- betareg(propintraspecific_N ~ size_z_scored + site_nbspe, data = SiteVar)
summary(N_es_beta)

# create sequence for TP
newdat <- data.frame(
  size_z_scored = seq(min(SiteVar$size_z_scored, na.rm = TRUE),
                      max(SiteVar$size_z_scored, na.rm = TRUE),
                      length.out = 100),
  site_nbspe = mean(SiteVar$site_nbspe, na.rm = TRUE)  # hold constant
)

# predictions
newdat$pred <- predict(N_es_beta, newdata = newdat, type = "response")

ggplot(SiteVar, aes(x = size_z_scored, y = propintraspecific_N)) +
  geom_point() +
  geom_line(data = newdat, aes(x = size_z_scored, y = pred), linewidth = 1)



##-Inf for 1 site: Swartspruit - would mean size of 0, some error 
C_es <- lm(propintraspecific_C ~ size_z_scored, data = SiteVar)
summary(C_es)
c_es_1 <- ggplot(SiteVar, aes(x = size_z_scored, y = propintraspecific_C)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()

c_es_2 <- ggplot(SiteVar, aes(x = size_z_scored, y = propintraspecific_C,  group = Climate_zone_e2, color = Climate_zone_e2)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()
c_es_2

c_es_3 <- ggplot(SiteVar, aes(x = size_z_scored, y = propintraspecific_C,  group = waterbody_type, color = waterbody_type)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()

N_es <- lm(propintraspecific_N ~ size_z_scored, data = SiteVar)
summary(N_es)
n_es_1 <- ggplot(SiteVar, aes(x = size_z_scored, y = propintraspecific_N)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()

n_es_2 <- ggplot(SiteVar, aes(x = size_z_scored, y = propintraspecific_N,  group = Climate_zone_e2, color = Climate_zone_e2)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()

n_es_3 <- ggplot(SiteVar, aes(x = size_z_scored, y = propintraspecific_N,  group = waterbody_type, color = waterbody_type)) +
  geom_point() +
  geom_smooth(method= "lm") +
  theme_classic()


ggarrange(n_es_1, n_es_2, n_es_3, nrow = 3, ncol = 1)
ggarrange(c_es_1, c_es_2, c_es_3, nrow = 3, ncol = 1)

##climate class 
C_climate <- aov(propintraspecific_C ~ Climate_zone_e2, data = SiteVar)
summary(C_climate)
TukeyHSD(C_climate)

SiteVar$Climate_zone_e2 <- ordered(SiteVar$Climate_zone_e2, levels = 
                                     c("Cold and wet/mesic", "Cool and moist", "Cool temperate and dry/xeric", "Warm temperate", "Hot and moist", "Hot and dry"))
ggplot(SiteVar, aes(x = Climate_zone_e2, y = propintraspecific_C)) +
  geom_boxplot()+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

N_climate <- aov(propintraspecific_N ~ Climate_zone_e2, data = SiteVar)
summary(N_climate)
TukeyHSD(N_climate)
ggplot(SiteVar, aes(x = Climate_zone_e2, y = propintraspecific_N)) +
  geom_boxplot()+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##ecosystem type - lentic / lotic
C_type <- aov(propintraspecific_C ~ waterbody_type, data = SiteVar)
summary(C_type)
TukeyHSD(C_type)
ggplot(SiteVar, aes(x = waterbody_type, y = propintraspecific_C)) +
  geom_boxplot()+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

N_type <- aov(propintraspecific_N ~ waterbody_type, data = SiteVar)
summary(N_type)
TukeyHSD(N_type)
ggplot(SiteVar, aes(x = waterbody_type, y = propintraspecific_N)) +
  geom_boxplot()+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



###Statistical Analysis 
options(na.action = "na.fail")


##CARBON -----------------
##Multiple Linear Model 
##full model
lm_full_C <- lm(propintraspecific_C ~ TP + waterbody_type + Size + Climate_zone_e2 + hydro_dis_z_scored + hft, data = SiteVar)
summary(lm_full_C)

lm_function_C <- function(x) {
  df <- x 
  
  ##full model
  lm_full_C_int <- lm(
    propintraspecific_C ~ (TP + waterbody_type + size_z_scored + Climate_zone_e2 + hydro_dis_z_scored + hft)^2,
    data = df
  )
  
  ##run AIC on all models using dredging approach -- compares all possible nest models constrained to those w/ CTmax
  model_set <- dredge(lm_full_C_int)
  model_set
  
  ##potential averagin approach: 
   avg_mod <- model.avg(model_set, subset = delta < 2)
    summary(avg_mod)
  
  ##select the best model - don't just do this if more than 1 good model 
  top_models <- subset(model_set, delta < 2)
  top_models
  # index of the most parsimonious model (smallest df)
  best_idx <- which.min(top_models$df)
  
  # get all models in the ΔAIC < 2 set:
  top_model_list <- get.models(model_set, subset = delta < 2)
  
  # pick the most parsimonious one:
  best_model <- top_model_list[[best_idx]]
  best_model_summary <- summary(best_model)
  best_model
  
  

  return(list(
    data = df,
    model_set = model_set,
    top_models = top_models,
    best_model = best_model,
    best_model_summary = best_model_summary,
    avg_mod = avg_mod
  ))
  
}


C_iv = lm_function_C(SiteVar)
C_iv$best_model_summary


##Linear Mixed Effects Model w/ Random Effect for Species 
lmer_function_C <- function(x) {
  df <- x 
  
  ##full model
  lmer_full_C_int <- lmer(
    propintraspecific_C ~
      (TP + waterbody_type + size_z_scored + Climate_zone_e2 + hydro_dis_z_scored + hft)^2 +
      (1 | site_nbspe),
    data = df,
    REML = FALSE
  )
  
  ##run AIC on all models using dredging approach -- compares all possible nest models constrained to those w/ CTmax
  model_set <- dredge(lmer_full_C_int)
  model_set
  
  ##potential averagin approach: 
   avg_mod <- model.avg(model_set, subset = delta < 2)
    summary(avg_mod)
  
  ##select the best model - don't just do this if more than 1 good model 
  top_models <- subset(model_set, delta < 2)
  top_models
  # index of the most parsimonious model (smallest df)
  best_idx <- which.min(top_models$df)
  
  # get all models in the ΔAIC < 2 set:
  top_model_list <- get.models(model_set, subset = delta < 2)
  
  # pick the most parsimonious one:
  best_model <- top_model_list[[best_idx]]
  best_model_summary <- summary(best_model)
  best_model
  
  
  
  return(list(
    data = df,
    model_set = model_set,
    top_models = top_models,
    best_model = best_model,
    best_model_summary = best_model_summary,
    avg_mod = avg_mod
  ))
  
}


C_iv_lmer = lmer_function_C(SiteVar) 
C_iv_lmer$best_model_summary
C_iv_lmer$model_set

saveRDS(C_iv_lmer, file = "outputs/C_iv_lmer.rds")

##NITROGEN -----------------

lm_full_N <- lm(propintraspecific_N ~ TP + waterbody_type + Size + Climate_zone_e2 + hydro_dis_z_scored + hft, data = SiteVar_2)
summary(lm_full_N)

lm_function_N <- function(x) {
  df <- x 
  
  ##full model
  lm_full_N_int <- lm(
    propintraspecific_N ~ (TP + waterbody_type + size_z_scored + Climate_zone_e2 + hydro_dis_z_scored + hft)^2,
    data = df
  )
  
  ##run AIC on all models using dredging approach -- compares all possible nest models constrained to those w/ CTmax
  model_set <- dredge(lm_full_N_int)
  model_set
  
  ##potential averagin approach: 
   avg_mod <- model.avg(model_set, subset = delta < 2)
    summary(avg_mod)
  
  ##select the best model - don't just do this if more than 1 good model 
  top_models <- subset(model_set, delta < 2)
  top_models
  # index of the most parsimonious model (smallest df)
  best_idx <- which.min(top_models$df)
  
  # get all models in the ΔAIC < 2 set:
  top_model_list <- get.models(model_set, subset = delta < 2)
  
  # pick the most parsimonious one:
  best_model <- top_model_list[[best_idx]]
  best_model_summary <- summary(best_model)
  best_model
  
  
  
  return(list(
    data = df,
    model_set = model_set,
    top_models = top_models,
    best_model = best_model,
    best_model_summary = best_model_summary,
    avg_mod = avg_mod
  ))
  
}


N_iv = lm_function_N(SiteVar)
N_iv$best_model_summary


##Linear Mixed Effects Model w/ Random Effect for Species 
lmer_function_N <- function(x) {
  df <- x 
  
  lmer_full_N_int <- lmer(
    propintraspecific_N ~
      (TP + waterbody_type + size_z_scored + Climate_zone_e2 + hydro_dis_z_scored + hft)^2 +
      (1 | site_nbspe),
    data = df,
    REML = FALSE
  )
  
 # full_model <- summary(lmer_full_N_int)
  ##run AIC on all models using dredging approach -- compares all possible nest models constrained to those w/ CTmax
  model_set <- dredge(lmer_full_N_int)
  model_set
  
  ##potential averagin approach: 
   avg_mod <- model.avg(model_set, subset = delta < 2)
    summary(avg_mod)
  
  ##select the best model - don't just do this if more than 1 good model 
  top_models <- subset(model_set, delta < 2)
  top_models
  # index of the most parsimonious model (smallest df)
  best_idx <- which.min(top_models$df)
  
  # get all models in the ΔAIC < 2 set:
  top_model_list <- get.models(model_set, subset = delta < 2)
  
  # pick the most parsimonious one:
  best_model <- top_model_list[[best_idx]]
  best_model_summary <- summary(best_model)
  best_model
  
  
  
  return(list(
    data = df,
    full_model = lmer_full_N_int,
    model_set = model_set,
    top_models = top_models,
    best_model = best_model,
    best_model_summary = best_model_summary, 
    avg_mod = avg_mod
  ))
  
}


N_iv_lmer = lmer_function_N(SiteVar) 
N_iv_lmer$best_model_summary
N_iv_lmer$model_set


##Look at distribution of residuals and error
res  <- resid(N_iv_lmer)
fit  <- fitted(lmer_full_N_int)







v <- visreg(best_model, "propintraspecific_C", partial = TRUE, plot = FALSE)
plot_sig <- ggplot() +
  geom_line(data = v$fit, aes(x = propintraspecific_C, y = visregFit), color = "red", linewidth = 1) +
  geom_ribbon(data = v$fit, aes(x = propintraspecific_C, ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  geom_point(data = v$res, aes(x = CTmax, y = visregRes), size = 2) +
  theme_classic() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
  labs(x = "Mean Estimated CTmax (˚C)", y = paste("Log", nutrient))


plot_marginal_sig <- ggplot() +
  geom_line(data = v$fit, aes(x = CTmax, y = visregFit), color = "red", linewidth = 1, linetype = "dashed") +
  geom_ribbon(data = v$fit, aes(x = CTmax, ymin = visregLwr, ymax = visregUpr), alpha = 0.2) +
  geom_point(data = v$res, aes(x = CTmax, y = visregRes), size = 2) +
  theme_classic() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
  labs(x = "Mean Estimated CTmax (˚C)", y = paste("Log", nutrient))

plot_no_sig <- ggplot() +
  geom_point(data = v$res, aes(x = CTmax, y = visregRes), size = 2) +
  theme_classic() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16)) +
  labs(x = "Mean Estimated CTmax (˚C)", y = paste("Log", nutrient))



