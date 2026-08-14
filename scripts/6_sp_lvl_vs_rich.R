###Species-level changes with diversity 

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
library(cowplot)
library(marginaleffects)

SiteVar <- read.csv("data/SiteVar_2sp_75cutoff.csv") %>%
  mutate(Climate_zone_2cat = case_when(
    startsWith(Climate_zone_e2, "Co") ~"Cold/Cool",
    startsWith(Climate_zone_e2, "Ho") ~"Warm/Hot",
    startsWith(Climate_zone_e2, "Wa") ~"Warm/Hot",
  ))


test <- SiteVar %>%
  filter(is.na(Climate_zone_e2))

ggplot(SiteVar, aes(y = site_interspe_var_N, x = log(site_nbspe))) +
  geom_point() +
  geom_smooth(method = "lm")

lm1_N <- lm(site_interspe_var_N ~ log(site_nbspe), data = SiteVar)
summary(lm1_N)

ggplot(SiteVar, aes(y = site_interspe_var_C, x = log(site_nbspe))) +
  geom_point() +
  geom_smooth(method = "lm")

lm1_C <- lm(site_interspe_var_C ~ log(site_nbspe), data = SiteVar)
summary(lm1_C)

site_var_env <- SiteVar %>%
  select(FWB_id, waterbody_type, Climate_zone_2cat, site_nbspe)


##Number of communities -- 453 (# of sites) --455
##Number of species:
SpVar <- read.csv("data/SpVar_2sp_75cutoff.csv") %>%
  left_join(site_var_env, by = "FWB_id")

test <- SpVar %>%
 # filter(fish_name_level == "species") %>%
  select(fish_scientific_name) %>%
  group_by(fish_scientific_name) %>%
  count() %>%
  filter(n >1)
  distinct()


test %>%
  summarise(
    unique_taxa = n_distinct(fish_scientific_name),
    unique_species = n_distinct(
      fish_scientific_name[fish_name_level == "species"]
    )
  )
##can test for species that are in both -- do they have the different slopes when found in different climates ? 
##do we have fish that we select for this from both lentic and lotic sites? 

##habitat
sp_habitat_list <- SpVar %>%
  select(waterbody_type.x, fish_scientific_name) %>%
  distinct() %>%
  group_by(fish_scientific_name) %>%
  count()
###look at the species and number of species across sampling sites and climate 
sp_climate_list <- SpVar %>%
  select(Climate_zone_2cat, fish_scientific_name, FWB_id) %>%
  distinct() %>%
  group_by(Climate_zone_2cat, fish_scientific_name) %>%
  count() %>%
  rename(num_fw_present = "n")

##select only species with a single climate zone
single_sp_climate <- sp_climate_list %>%
  select(Climate_zone_2cat, fish_scientific_name) %>%
  group_by(fish_scientific_name) %>%
  count() %>%
  filter(n == 1) %>%
  left_join(sp_climate_list, by = "fish_scientific_name")



sp_interest_list <- single_sp_climate %>%
  filter(num_fw_present >= 8)


##Looking at relationship between species richness and variance within certain species
##salmo trutta , pike, roach, carp? 
sp_var_sp_interest <- SpVar %>%
  filter(fish_scientific_name %in% sp_interest_list$fish_scientific_name)  %>%
  mutate(site_nbspe_log = log(site_nbspe),
         site_nbspe_log_scale = scale(site_nbspe_log),
         sp_site_var_C_log = log(sp_site_var_C),
         sp_site_var_N_log = log(sp_site_var_N),
         )
sp_habitat_list <- sp_var_sp_interest %>%
  select(waterbody_type.x, Climate_zone_2cat, fish_scientific_name) %>%
  distinct() #%>%

clim_hab_num <- sp_habitat_list %>%
  select(Climate_zone_2cat, waterbody_type.x) %>%
  group_by(Climate_zone_2cat, waterbody_type.x) %>%
  count()

group_by(fish_scientific_name) %>%
  count()


hist(sp_var_sp_interest$sp_site_var_C)
hist(sp_var_sp_interest$sp_site_var_C_log)
hist(sp_var_sp_interest$sp_site_var_N)
hist(sp_var_sp_interest$sp_site_var_N_log)

test <- sp_var_sp_interest %>%
  filter(fish_scientific_name == "Amniataba percoides")

num_sp <- sp_var_sp_interest %>%
  select(fish_scientific_name, Climate_zone_2cat) %>%
  distinct() %>%
  group_by(Climate_zone_2cat) %>%
  count()

ggplot(sp_var_sp_interest, aes(x = log(site_nbspe), y = sp_site_var_C, group = fish_scientific_name, color = Climate_zone_2cat)) +
  geom_point() +
  geom_smooth(method = "lm")

glmm1_C <- glmmTMB(sp_site_var_C_log ~ site_nbspe_log_scale * Climate_zone_2cat + waterbody_type.x + (site_nbspe_log_scale|fish_scientific_name), data = sp_var_sp_interest)
summary(glmm1_C)
Anova(glmm1_C)
diagnose(glmm1_C)
simulationOutput <- simulateResiduals(fittedModel = glmm1_C, plot = FALSE)
plot(simulationOutput)

sp_cold <- sp_var_sp_interest %>%
  select(fish_scientific_name, Climate_zone_2cat) %>%
  distinct() %>%
  filter(Climate_zone_2cat == "Cold/Cool")
sp_warm <- sp_var_sp_interest %>%
  select(fish_scientific_name, Climate_zone_2cat) %>%
  distinct() %>%
  filter(Climate_zone_2cat == "Warm/Hot")


nd_C_fix <- expand.grid(site_nbspe_log_scale = seq(-2, 4, 0.1), Climate_zone_2cat = c("Cold/Cool", "Warm/Hot"), waterbody_type.x = "lotic")
nd_C_fix$fish_scientific_name <- sp_var_sp_interest$fish_scientific_name[1]
# Population-level predictions (exclude random effects)
pred_C <- predict(glmm1_C, newdata = nd_C_fix, re.form = NA,se.fit = TRUE)

nd_C_fix <- nd_C_fix %>%
  mutate(
    prediction = pred_C$fit,
    se = pred_C$se.fit,
    lower = prediction - 1.96 * se,
    upper = prediction + 1.96 * se
  )
#nd_C_fix$prediction1 <- predict(glmm1_N, newdata = nd_N_fix)



nd_C_cold <- datagrid(site_nbspe_log_scale = seq(-2, 4, 0.1), Climate_zone_2cat = c("Cold/Cool"),waterbody_type.x = "lotic", fish_scientific_name = unique(sp_cold$fish_scientific_name), model = glmm1_C)
nd_C_cold$prediction1 <- predict(glmm1_C, newdata = nd_C_cold)

nd_C_warm <- datagrid(site_nbspe_log_scale = seq(-2, 4, 0.1), Climate_zone_2cat = c("Warm/Hot"),waterbody_type.x = "lotic",  fish_scientific_name = unique(sp_warm$fish_scientific_name), model = glmm1_C)
nd_C_warm$prediction1 <- predict(glmm1_C, newdata = nd_C_warm)

#nd_N_fix <- datagrid(site_nbspe_log_scale = seq(-2, 4, 0.1), Climate_zone_2cat = c("Cold/Cool", "Warm/Hot"),  model = glmm1_N)
#nd_N_fix$prediction1 <- predict(glmm1_N, newdata = nd_N_fix)

c_sp <- ggplot() +
  geom_line(data = nd_C_cold,aes( x = site_nbspe_log_scale, y = prediction1, group = fish_scientific_name ), color = "#99CCCC", alpha = 0.7) +
  geom_line(data = nd_C_warm,aes( x = site_nbspe_log_scale, y = prediction1, group = fish_scientific_name ), color = "#993333", alpha = 0.7) +
  geom_ribbon(data = nd_C_fix, aes(x = site_nbspe_log_scale, ymin = lower, ymax = upper, fill = Climate_zone_2cat), alpha = 0.4)+
   geom_line(data = nd_C_fix,aes( x = site_nbspe_log_scale, y = prediction, group = Climate_zone_2cat, color = Climate_zone_2cat), linewidth = 2)+
  #geom_abline(slope = -0.41131, intercept =  -4.05148, color = "#99CCCC", linewidth = 2) +
  #geom_abline(slope = -0.41131 -0.03275 , intercept = -4.05148 + 0.14264, color = "#993333", linewidth = 2) +
  scale_color_manual(name = "Climate zone", values = c("Cold/Cool" = "#99CCCC","Warm/Hot" = "#993333")) +
  scale_fill_manual(name = "Climate zone", values = c("Cold/Cool" = "#99CCCC","Warm/Hot" = "#993333")) +
  theme_classic() +
  xlab("Species Richness (log scaled)") +
  ylab("Intraspecific Variance in C (log)") +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))
c_sp

##need to add the confidence interval to the fixed effect lines
#DHARMa::residuals(simulationOutput)

glmm1_N <- glmmTMB(sp_site_var_N_log ~ site_nbspe_log_scale * Climate_zone_2cat +  waterbody_type.x + (site_nbspe_log_scale|fish_scientific_name), data = sp_var_sp_interest)
summary(glmm1_N)
Anova(glmm1_N)
ranef(glmm1_N)

simulationOutput <- simulateResiduals(fittedModel = glmm1_N, plot = FALSE)
plot(simulationOutput)
#DHARMa::residuals(simulationOutput)

min(sp_var_sp_interest$site_nbspe_log_scale)

max(sp_var_sp_interest$site_nbspe_log_scale)



nd_N_cold <- datagrid(site_nbspe_log_scale = seq(-2, 4, 0.1), Climate_zone_2cat = c("Cold/Cool"),waterbody_type.x = "lotic",  fish_scientific_name = unique(sp_cold$fish_scientific_name), model = glmm1_N)
nd_N_cold$prediction1 <- predict(glmm1_N, newdata = nd_N_cold)

nd_N_warm <- datagrid(site_nbspe_log_scale = seq(-2, 4, 0.1), Climate_zone_2cat = c("Warm/Hot"),waterbody_type.x = "lotic",  fish_scientific_name = unique(sp_warm$fish_scientific_name), model = glmm1_N)
nd_N_warm$prediction1 <- predict(glmm1_N, newdata = nd_N_warm)


nd_N_fix <- expand.grid(site_nbspe_log_scale = seq(-2, 4, 0.1), Climate_zone_2cat = c("Cold/Cool", "Warm/Hot"), waterbody_type.x = "lotic")
nd_N_fix$fish_scientific_name <- sp_var_sp_interest$fish_scientific_name[1]
# Population-level predictions (exclude random effects)
pred_N <- predict(glmm1_N, newdata = nd_N_fix, re.form = NA,se.fit = TRUE)

nd_N_fix <- nd_N_fix %>%
  mutate(
    prediction = pred_N$fit,
    se = pred_N$se.fit,
    lower = prediction - 1.96 * se,
    upper = prediction + 1.96 * se
  )
#nd_N_fix <- datagrid(site_nbspe_log_scale = seq(-2, 4, 0.1), Climate_zone_2cat = c("Cold/Cool", "Warm/Hot"),  model = glmm1_N)
#nd_N_fix$prediction1 <- predict(glmm1_N, newdata = nd_N_fix)

n_sp <- ggplot() +
  geom_line(data = nd_N_cold,aes( x = site_nbspe_log_scale, y = prediction1, group = fish_scientific_name ), color = "#99CCCC", alpha = 0.7) +
  geom_line(data = nd_N_warm,aes( x = site_nbspe_log_scale, y = prediction1, group = fish_scientific_name ), color = "#993333", alpha = 0.7) +
  geom_ribbon(data = nd_N_fix, aes(x = site_nbspe_log_scale, ymin = lower, ymax = upper, fill = Climate_zone_2cat), alpha = 0.4)+
  geom_line(data = nd_N_fix,aes( x = site_nbspe_log_scale, y = prediction, group = Climate_zone_2cat, color = Climate_zone_2cat), linewidth = 2)+
  
 # geom_line(data = nd_N_fix,aes( x = site_nbspe_log_scale, y = prediction1, group = Climate_zone_2cat, color = Climate_zone_2cat), linewidth = 2)+
 # geom_abline(slope = -0.7438, intercept = -4.4769, color = "#99CCCC", linewidth = 2) +
#  geom_abline(slope = -0.7438 +  0.4380, intercept = -4.4769 + 0.3614, color = "#993333", linewidth = 2) +
#  scale_color_manual(values = c("#99CCCC", "#993333")) +
  scale_color_manual(name = "Climate zone", values = c("Cold/Cool" = "#99CCCC","Warm/Hot" = "#993333")) +
  scale_fill_manual(name = "Climate zone", values = c("Cold/Cool" = "#99CCCC","Warm/Hot" = "#993333")) +
  theme_classic() +
  xlab("Species Richness (log scaled)") +
  ylab("Intraspecific Variance in N (log)") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))

n_sp





##see if model runs if species richness is not scaled 
ggarrange(c_sp, n_sp, labels = c("a", "b"), font.label = list(size = 18, face = "bold"),nrow =1, ncol = 2, legend = c("none"))

##account for body size in model, just to show that this isn't driven by body size, add into model , then remove 

##OLD --- 
ggplot(sp_var_sp_interest, aes(x = log(site_nbspe), y = log(sp_site_var_N), group = fish_species, color = Climate_zone_2cat)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  theme(legend.position = "none") +
  facet_wrap(~fish_species)



ggplot(sp_var_sp_interest, aes(x = log(site_nbspe), y = log(sp_site_var_C), group = fish_species, color = Climate_zone_2cat)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  theme(legend.position = "none") +
  facet_wrap(~fish_species, scale = "free")



glmm1_N <- glmmTMB(sp_site_var_N ~ log(site_nbspe) + Climate_zone_2cat + (1|fish_species), data = sp_var_sp_interest)
summary(glmm1_N)
Anova(glmm1_N)




##set up data for running model
glmm1_N <- glmmTMB(sp_site_var_N ~ log(site_nbspe)*Climate_zone_2cat + (1|fish_species), data = sp_var_sp_interest)
summary(glmm1_N)
Anova(glmm1_N)

glmm1_C <- glmmTMB(sp_site_var_C ~ log(site_nbspe)*Climate_zone_2cat + (1|fish_species), data = sp_var_sp_interest)
summary(glmm1_C)
Anova(glmm1_C)

lm_func_C <- function(df){
  lm1 <- lm(sp_site_var_C ~ log(site_nbspe), data = df)
 # summary(lm1)
  coef(lm1)[["log(site_nbspe)"]]
  
}

lm_C_all <- split(sp_var_sp_interest, paste0(sp_var_sp_interest$fish_species)) %>%
  map(lm_func_C) %>%
  bind_rows() %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(isotope = "C") %>%
  rename(fish_species = "variable") %>%
  left_join(sp_var_sp_interest %>% select(fish_species, site_nbspe, Climate_zone_2cat) %>% distinct, by = "fish_species")

ggplot(lm_C_all, aes(x = value, fill = Climate_zone_2cat)) +
  geom_histogram()

ggplot(lm_C_all, aes(x = Climate_zone_2cat, y = value)) +
  geom_boxplot() +
  ylab("Slope")

hist(lm_C_all$value)
  
lm_func_N <- function(df){
  lm1 <- lm(sp_site_var_N ~ log(site_nbspe), data = df)
  # summary(lm1)
  coef(lm1)[["log(site_nbspe)"]]
  
}

lm_N_all <- split(sp_var_sp_interest, paste0(sp_var_sp_interest$fish_species)) %>%
  map(lm_func_N) %>%
  bind_rows() %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(isotope = "N") %>%
  rename(fish_species = "variable") %>%
  left_join(sp_var_sp_interest %>% select(fish_species, site_nbspe, Climate_zone_2cat) %>% distinct, by = "fish_species")

ggplot(lm_N_all, aes(x = value, fill = Climate_zone_2cat)) +
  geom_histogram()

ggplot(lm_N_all, aes(x = Climate_zone_2cat, y = value)) +
  geom_boxplot() +
  ylab("Slope")


hist(lm_N_all$value)


lm_all <- rbind(lm_C_all, lm_N_all)



ggplot(lm_all, aes(x = isotope, y = value, group = )) +
  geom_boxplot() +
  theme_classic() +
  ylab("Slope")


ggplot(sp_var_sp_interest, aes(x = sp_site_body_size_range, y = sp_site_var_C)) +
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~fish_species, scale = "free")

ggplot(sp_var_sp_interest, aes(x = sp_site_body_size_range, y = sp_site_var_N)) +
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~fish_species, scale = "free")

ggplot(sp_var_sp_interest, aes(x = sp_site_var_length, y = sp_site_var_C)) +
  geom_point() +
  facet_wrap(~fish_species, scale = "free")

##filter out the species of interest