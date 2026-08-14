##Generate map and descriptive statistics of intraspecific variability contribution 
##

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


##read in contribution of site level intraspecific contribution RData
#load("data/Intraspecific_contribution_perSite_Env.RData") ##produced by 2_Intra_vs_Inter_21_01_26.R
SiteVar <- read.csv("data/SiteVar_2sp_75cutoff.csv") %>%
  mutate(Climate_zone_2cat = case_when(
    startsWith(Climate_zone_e2, "Co") ~"Cold/Cool",
    startsWith(Climate_zone_e2, "Ho") ~"Warm/Hot",
    startsWith(Climate_zone_e2, "Wa") ~"Warm/Hot",
  ))

#quick look at climate categories:
climate <- SiteVar %>%
  select(Climate_zone, Climate_zone_e, Climate_zone_e2, Climate_zone_2cat) %>%
  distinct()
#fwb_climate <- SiteVar %>%
#  select(FWB_id, Climate_zone_2cat)
#write.csv(fwb_climate, "data/fwb_climate_list.csv")
##C:N ratio -- look for effect especially in muscle, but may not be available 

##also looking for relationships between site environmental variables and community characteristic (e.g., species richness, functional diversity, )

##1) Description of variation -- how much intraspecific. variability contribution in C and N 
## -- also, do we see correlations between contribution of C and N - i.e., when more intrasp. var in C also see intrasp var in N? 

##2) Between site -- Relationship between environmental drivers and proportion of intraspecific variability 

##3) Between species -- drivers of dufferences in intraspecific variability between species -- trait drivers, community level drivers 
##potentially using site as a random effect 

##Looking at decomposition of total community variance - into between-species variance and within species variance
##total community variance = site_intraspe_var + site_intersp_var 
##total within species var = site_intraspe_var (mean within species variation across all species at the site)
##total between species var = site_interspe_var ( variation in mean dN between species at the site)

##then we effectively calculated the proportion of variation that is attributed to intraspecific variation - this is how much of the entire isotopic variability within a food web was associated with intraspecific variability 
#SiteVar$propintraspecific_N<-SiteVar$site_intraspe_var_N/(SiteVar$site_interspe_var_N+SiteVar$site_intraspe_var_N)
#SiteVar$propintraspecific_C<-SiteVar$site_intraspe_var_C/(SiteVar$site_interspe_var_C+SiteVar$site_intraspe_var_C)
#SiteVar$propintraspecific_Total<-(SiteVar$site_intraspe_var_N+SiteVar$site_intraspe_var_C)/(SiteVar$site_interspe_var_N+SiteVar$site_intraspe_var_N+SiteVar$site_interspe_var_C+SiteVar$site_intraspe_var_C)
##these are already calculated, just keeping here for formulas



###1) Create map of sites and generate boxplots that visualize % of total variance that comes from intraspecific variation in both C and N ----
#convert points to sf objects and make base world map
sites_coord <- st_as_sf(SiteVar, coords = c("collection_decimal_longitude", "collection_decimal_latitude"), crs = 4326)  %>%
  st_make_valid() %>%
  st_transform(4326) 

countries <- ne_countries(scale = "medium", returnclass = "sf") %>%
 # st_make_valid() %>%
  st_transform(4326)


sites_geo <- st_join(
  sites_coord,
  countries[, c("admin","iso_a3","continent")],
  join = st_intersects,
  left = TRUE
) %>% rename(country = admin)

na_idx <- which(is.na(sites_geo$country))
if (length(na_idx) > 0) {
  nearest <- st_nearest_feature(sites_geo[na_idx,], countries)
  sites_geo$country[na_idx]   <- countries$admin[nearest]
  sites_geo$iso_a3[na_idx]    <- countries$iso_a3[nearest]
  sites_geo$continent[na_idx] <- countries$continent[nearest]
}

test <- sites_geo %>%
  select(collection_site_id, waterbody_type, country, continent, propintraspecific_N, propintraspecific_C) %>%
  filter(is.na(continent))

sites_geo <- sites_geo %>%
  select(collection_site_id, waterbody_type, Climate_zone_2cat,  country, continent, site_nbspe, propintraspecific_N, propintraspecific_C, propintraspecific_Total) %>%
  pivot_longer(cols= c(propintraspecific_N, propintraspecific_C, propintraspecific_Total), names_to = "prop_intraspecific_var_type", values_to = "prop_intraspecific_var") %>%
  mutate(prop_type = case_when(
    startsWith("propintraspecific_N", prop_intraspecific_var_type) ~ "N",
    startsWith("propintraspecific_C", prop_intraspecific_var_type) ~ "C",
    startsWith("propintraspecific_Total", prop_intraspecific_var_type) ~ "Total",
  ))


world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(admin != "Antarctica")
##Update projection for plotting map -- looks good
# Robinson projection (good-looking global)
crs_robin <- "+proj=robin"

world_r <- st_transform(world, crs_robin)
sites_r <- st_transform(sites_coord, crs_robin)

p_map <- ggplot() +
  geom_sf(data = world_r, fill = "lightgrey", color = "darkgrey", linewidth = 0.1) +
  geom_sf(data = sites_r, aes(fill = Climate_zone_2cat), color = "black", shape = 21, size = 1.5, alpha = 0.9) +
  coord_sf(expand = FALSE) +
  scale_fill_manual(name = "Climate zone", values = c("#99CCCC", "#993333"))+
  # annotation_scale(location = "bl", width_hint = 0.25) +
  #  annotation_north_arrow(location = "bl", which_north = "true",
  #                        style = north_arrow_fancy_orienteering) +
  theme_void() +
  theme(
    # plot.margin = margin(5, 5, 5, 5),
    legend.position = "right"
  ) +
  theme(plot.background = element_rect(fill = "white", color = NA))

p_map

p_map_2 <- ggplot() +
  geom_sf(data = world_r, fill = "lightgrey", color = "darkgrey", linewidth = 0.1) +
  geom_sf(data = sites_r, aes(color = Climate_zone_e2),size = 1, alpha = 0.9) +
  coord_sf(expand = FALSE) +
#  scale_color_manual(values = c("#4C72B0", "#DD8452"))+
  # annotation_scale(location = "bl", width_hint = 0.25) +
  #  annotation_north_arrow(location = "bl", which_north = "true",
  #                        style = north_arrow_fancy_orienteering) +
  theme_void() +
  theme(
    # plot.margin = margin(5, 5, 5, 5),
    legend.position = "right"
  ) +
  theme(plot.background = element_rect(fill = "white", color = NA))

p_map_2




##Brief descriptive stats ----------------

site_var_env <- SiteVar %>%
  select(FWB_id, waterbody_type, Climate_zone_2cat)
##Number of communities -- 453 (# of sites) --455
##Number of species:
SpVar <- read.csv("data/SpVar_2sp_75cutoff.csv") %>%
  left_join(site_var_env, by = "FWB_id")


##climate by category
clim_cat <- SpVar %>%
  select(fish_scientific_name, Climate_zone_2cat) %>%
  distinct() %>%
  group_by(fish_scientific_name) %>%
  count()

##Number of countries:
country_num <- as.data.frame(sites_geo) %>%
  select(country) %>%
  ungroup() %>%
  distinct()

##mean prop intravar
mean_prop_iv <- SiteVar %>%
  summarise_at(vars(propintraspecific_N, propintraspecific_C), list(mean = mean))


##Pairwise t-tests 
ttest_df <- SiteVar %>%
  select(FWB_id, propintraspecific_N, propintraspecific_C) %>%
  pivot_longer(cols = c(propintraspecific_C, propintraspecific_N), names_to = "isotope", values_to = "pro_intrasp")

ttest_result <- t.test(SiteVar$propintraspecific_C,
                       SiteVar$propintraspecific_N, paired= TRUE)
summary(ttest_result)


ttest_result


ggplot(SiteVar, aes(y = propintraspecific_N, propintraspecific_C)) + 
  geom_point()

ttest_df <- SiteVar %>%
  select(FWB_id, propintraspecific_N, propintraspecific_C) %>%
  drop_na() %>%
  mutate(
    difference = propintraspecific_C - propintraspecific_N,
    higher = case_when(
      difference > 0 ~ "C higher",
      difference < 0 ~ "N higher",
      difference == 0 ~ "Equal"
    )
  )

table(ttest_df$higher)
prop.table(table(ttest_df$higher))

ggplot(ttest_df, aes(x = difference)) +
  geom_histogram() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(
    x = "Within-site difference (C - N)",
    y = "Number of sites"
  )

##MODEL:
##glmm w/ ecosystem type, climate 

##Looking at relationship between species richness and variance within certain species
##salmo trutta , pike, roach, carp? 
sp_var_sp_interest <- SpVar %>%
  left_join(SiteVar %>% select(FWB_id, site_nbspe), by = "FWB_id") %>%
  filter(fish_scientific_name %in% c("Salmo trutta", "Esox lucius", "Rutilus rutilus", "Cyprinus carpio"))


ggplot(sp_var_sp_interest, aes(x = site_nbspe, y = sp_site_var_C, color = log(sp_site_num_ind))) +
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~fish_scientific_name, scale = "free")
ggplot(sp_var_sp_interest, aes(x = site_nbspe, y = sp_site_var_N, color = log(sp_site_num_ind))) +
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~fish_scientific_name, scale = "free")

ggplot(sp_var_sp_interest, aes(x = sp_site_body_size_range, y = sp_site_var_C)) +
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~fish_scientific_name, scale = "free")

ggplot(sp_var_sp_interest, aes(x = sp_site_body_size_range, y = sp_site_var_N)) +
  geom_point() +
  geom_smooth(method = "lm")+
  facet_wrap(~fish_scientific_name, scale = "free")

ggplot(sp_var_sp_interest, aes(x = sp_site_var_length, y = sp_site_var_C)) +
  geom_point() +
  facet_wrap(~fish_scientific_name, scale = "free")
##Mean richness:
mean_rich <- SiteVar %>%
  summarise(mean_richness = mean(site_nbspe),
            sd_richness = sd(site_nbspe),
            median_richness = median(site_nbspe))

samp_mean <- SiteVar %>%
  summarise(mean_samplesize = mean(site_mean_sample_id),
            sd_samplesize = sd(site_mean_sample_id),
            median_samplesize = median(site_mean_sample_id))

##mean by continent: 
mean_continent <- as.data.frame(sites_geo) %>%
  group_by(continent, prop_intraspecific_var_type) %>%
  summarise(mean_CIV = mean(prop_intraspecific_var),
            sd_CIV = sd(prop_intraspecific_var),
            median_CIV = median(prop_intraspecific_var))

##run aov for each continent
library(dplyr)
library(purrr)
library(broom)

continent_aov <- sites_geo %>%
  group_by(continent) %>%
  nest() %>%
  mutate(
    model = map(data, ~ aov(prop_intraspecific_var ~ prop_intraspecific_var_type,
                            data = .x)),
    anova = map(model, tidy)
  ) %>%
  select(continent, anova) %>%
  unnest(anova)

continent_aov

continents_aov_2 <- aov(prop_intraspecific_var ~ prop_intraspecific_var_type * continent,
    data = sites_geo)
continents_aov_2
TukeyHSD(continents_aov_2)
##Regression and ANOVA for sp richness, habitat and climate
c_rich_lm <- lm(propintraspecific_C ~ log(site_nbspe), data = SiteVar)
summary(c_rich_lm)

n_rich_lm <- lm(propintraspecific_N ~ log(site_nbspe), data = SiteVar)
summary(n_rich_lm)

c_hab_aov <- aov(propintraspecific_C ~ Type, data = SiteVar)
summary(c_hab_aov)

n_hab_aov <- aov(propintraspecific_N ~ Type, data = SiteVar)
summary(c_hab_aov)

c_clim_aov <- aov(propintraspecific_C ~ Climate_zone_2cat, data = SiteVar)
summary(c_clim_aov)

n_clim_aov <- aov(propintraspecific_N ~ Climate_zone_2cat, data = SiteVar)
summary(n_clim_aov)

##Correlation between % CIV in C and N 
c_n <- lm(propintraspecific_N ~ propintraspecific_C, data = SiteVar)
summary(c_n)
###Some exploratory plots to look at differences in variation contribution ----------------

##plotting distributions
hist_N_intra <- ggplot(SiteVar, aes(x = propintraspecific_N)) +
  geom_histogram() +
  theme_classic() 
hist_N_intra

hist_C_intra <- ggplot(SiteVar, aes(x = propintraspecific_C)) +
  geom_histogram() +
  theme_classic()
hist_C_intra


##comparing variation in proportion 

color1 <- c("darkgrey", "black")

box_c_n_cn <- sites_geo %>%
  ggplot(aes(x = prop_type, y = prop_intraspecific_var)) +
  geom_boxplot() +
#  scale_fill_manual(values = color1)+
  theme_classic() +
  ylab("CIV Proportion") +
  theme(axis.title.x = element_blank())
box_c_n_cn

box_c_n_habitat <- sites_geo %>%
 # filter(prop_type != "Total") %>%
  ggplot(aes(x = prop_type, y = prop_intraspecific_var, fill = waterbody_type)) +
  geom_boxplot() +
  scale_fill_manual(values = color1)+
  theme_classic() +
  ylab("CIV Proportion") +
  theme(axis.title.x = element_blank())

box_c_n_habitat

##show just C and N -- run everything with all 3 for the total, then put in the methods that it is basically averaging between the two 

box_c_n_habitat_sprich_C <- sites_geo %>%
   filter(prop_type == "C") %>%
  ggplot(aes(x = log(site_nbspe), y = prop_intraspecific_var,  color = waterbody_type)) +
  geom_point() +
  geom_smooth(method = "lm")+
  scale_color_manual(values = color1)+
  theme_classic() +
  ylab("CIV Proportion C") +
  theme(axis.title.x = element_blank())
box_c_n_habitat_sprich_C

box_c_n_habitat_sprich_N <- sites_geo %>%
  filter(prop_type == "N") %>%
  ggplot(aes(x = log(site_nbspe), y = prop_intraspecific_var,  color = waterbody_type)) +
  geom_point() +
  geom_smooth(method = "lm")+
  scale_color_manual(values = color1)+
  theme_classic() +
  ylab("CIV Proportion N") +
  theme(axis.title.x = element_blank())
box_c_n_habitat_sprich_N

box_c_n_habitat_sprich_N <- sites_geo %>%
  filter(prop_type == "N") %>%
  ggplot(aes(x = as.factor(site_nbspe), y = prop_intraspecific_var,  fill = waterbody_type)) +
  geom_boxplot() +
  scale_fill_manual(values = color1)+
  theme_classic() +
  ylab("CIV Proportion N") +
  theme(axis.title.x = element_blank())
box_c_n_habitat_sprich_N



color2 <- c("lightblue", "red")
box_climate <- sites_geo %>%
 # filter(prop_type != "Total") %>%
  ggplot(aes(x = prop_type, y = prop_intraspecific_var, fill = Climate_zone_2cat)) +
  geom_boxplot() +
  scale_fill_manual(values = color2)+
  theme_classic() +
  ylab("CIV Proportion") +
  theme(axis.title.x = element_blank())

box_climate

box_c_n_climate_sprich_C <- sites_geo %>%
  filter(prop_type == "C") %>%
  ggplot(aes(x = log(site_nbspe), y = prop_intraspecific_var,  color = Climate_zone_2cat)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_fill_manual(values = color1)+
  theme_classic() +
  ylab("CIV Proportion C") +
  theme(axis.title.x = element_blank())
box_c_n_climate_sprich_C

box_c_n_climate_sprich_N <- sites_geo %>%
  filter(prop_type == "N") %>%
  ggplot(aes(x = log(site_nbspe), y = prop_intraspecific_var,  color = Climate_zone_2cat)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_color_manual(values = color1)+
  theme_classic() +
  ylab("CIV Proportion N") +
  theme(axis.title.x = element_blank())
box_c_n_climate_sprich_N

sp_rich_reg <- sites_geo %>%
 # filter(prop_type != "Total") %>%
  ggplot(aes(x = log(site_nbspe), y = prop_intraspecific_var, color = prop_type)) +
  geom_point() +
  geom_smooth(method = "lm")+
#  scale_color_manual(values = color1)+
  theme_classic() +
  ylab("CIV Proportion") +
  xlab("Species Richness (log)")
 # theme(axis.title.x = element_blank())
sp_rich_reg


plot_sup <- ggarrange(box_c_n_cn, box_c_n_habitat, box_climate, sp_rich_reg, legend = "bottom", nrow = 2, ncol = 2, labels = c("a)", "b)", "c)", "d)"), font.label = list(colour = "black", size = 12))
plot_sup

glmm1 <- g


cn_correlation <- ggplot(SiteVar, aes(x = propintraspecific_C, y = propintraspecific_N)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", color = "darkred") +
  theme_classic() +
  ylab("% CIV in N") +
  xlab("% CIV in C") 

cn_correlation

lm1 <- lm(propintraspecific_N ~ propintraspecific_C, data = SiteVar)
summary(lm1)


plot_1 <- ggarrange(box_c_n_2, cn_correlation, legend = "none", nrow = 1, ncol = 2, labels = c("b)", "c)"), font.label = list(colour = "black", size = 12))
plot_1


##Variation in proportion contribution by species richness
SiteVar_long <- SiteVar %>%
  pivot_longer(cols= c(propintraspecific_N, propintraspecific_C), names_to = "prop_intraspecific_var_type", values_to = "prop_intraspecific_var") %>%
  mutate(prop_type = case_when(
    startsWith("propintraspecific_N", prop_intraspecific_var_type) ~ "N",
    startsWith("propintraspecific_C", prop_intraspecific_var_type) ~ "C",
  ))
cn_sprich_plot <- ggplot(SiteVar_long, aes(x = log(site_nbspe), y = prop_intraspecific_var, group = prop_type, color = prop_type)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  ylab("% CIV") +
  xlab("Species Richness (log)")


plot_2 <- ggarrange(cn_correlation, cn_sprich_plot, legend = "right", nrow = 1, ncol = 2, labels = c("a)", "b)"), font.label = list(colour = "black", size = 12))
plot_2

ggplot(SiteVar, aes(x = log(site_nbspe), y = propintraspecific_N)) +
  geom_point() +
  geom_smooth(method = "lm")


##Variation in proportion contribution by mean and min sample size
SiteVar %>%
  filter(site_min_sample_id >= 2 & site_min_sample_id <= 100) %>%
ggplot(aes(x = log(site_mean_sample_id), y = propintraspecific_C)) +
  geom_point() +
  geom_smooth(method = "lm")




ggplot(SiteVar, aes(x = log(site_mean_sample_id), y = propintraspecific_N)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(SiteVar, aes(x = site_min_sample_id)) +
  geom_histogram(binwidth = 1)

ggplot(SiteVar, aes(x = site_mean_sample_id)) +
  geom_histogram(binwidth = 1)

ggplot(SiteVar, aes(x = site_nbspe)) +
  geom_histogram(binwidth = 1)

##Look at individual species x site level
load("data/SpeciesIntraspecific_variance_perSite_Env.RData")
SpVar_env %>%
  filter(sp_site_num_ind >5) %>%
  ggplot(aes(x = log(sp_site_num_ind), y = sp_site_var_C)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(SpVar_env, aes(x = log(sp_site_num_ind), y = sp_site_var_N)) +
  geom_point() +
  geom_smooth(method = "lm")

##mean d13C/d15N - normalized d15N and d13C


SpVar_env %>%
 # filter(sp_site_num_ind >5) %>%
  ggplot(aes(x = log(sp_site_num_ind), y = sp_site_mean_C)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(SpVar_env, aes(x = log(sp_site_num_ind), y = sp_site_mean_N)) +
  geom_point() +
  geom_smooth(method = "lm")


##Looking at mean / sd of raw isotope values
df_all <- read_excel("data/FINAL_ALLindiv_February2026.xlsx") 



##clean data
##select only fish, where have C/N data and associated scientific name
DataFish<-subset(df_all,!is.na(d15N) & !is.na(d13C) & organism_type=="fish" & !is.na(scientific_name))

##add species-site identifier column
DataFish$sp_site<-paste(DataFish$fish_species,DataFish$collection_site_id,sep="_")
##Assign 1 - number of fish
DataFish$num<-1

##normalise data
DataFish<-DataFish %>% 
  group_by(collection_site_id) %>% 
  mutate(d15N_norm = (d15N - min(d15N, na.rm = TRUE))/(max(d15N, na.rm = TRUE)-min(d15N, na.rm = TRUE)),
         d13C_norm = (d13C - min(d13C, na.rm = TRUE))/(max(d13C, na.rm = TRUE)-min(d13C, na.rm = TRUE)))



##if we remove small sample size, when does slope become 0? 

##need to keep in mind that species richness may bias intrapsecific variability -- could potentially ahndle by taking residuals and then using residuals from this relationships w/ envirnmental variable s 



plot_2 <- ggarrange(p_map, plot_1, legend = "none", nrow = 2, ncol = 1, labels = c("a)", ""), font.label = list(colour = "black", size = 12))
plot_2

