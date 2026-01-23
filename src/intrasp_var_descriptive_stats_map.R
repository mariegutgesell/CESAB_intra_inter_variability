##Relative contribution of intraspecific variability to community variability relative to environmental drivers
##Looking at drivers of variability between sites 

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
##read in individual level si data
#ind_df <- read_excel("data/ALLindiv_January2026.xlsx")

ind_df <- read_excel("data/ALLindiv_June2025.xlsx")

##read in environmental data
#env_df <- read_excel("data/FOODWEBS_Individual_EnvironmentalData_28Aug25.xlsx")
##for individual environmental data from august, need to find list of FW_ID for collection_site_id in the stable isotope data

env_df_lentic <- read_excel("data/FOODWEBS_Individual_EnvironmentalData_13Jun25.xlsx", sheet = "Lenthic") %>%
  rename(Ecosystem = "Ecosystem_Type")
#env_df_lentic$Ecosystem_Type <- as.character(env_df_lentic$Ecosystem_Type)

env_df_lotic <- read_excel("data/FOODWEBS_Individual_EnvironmentalData_13Jun25.xlsx", sheet = "Lotic")
#env_df_lotic$Ecosystem_Type <- as.character(env_df_lotic$Ecosystem_Type)
##for environmental data, follow same way they did in the FCL 

env_df <- full_join(env_df_lentic, env_df_lotic)
  


##read in calculated relative contribution data 

site_is_df <- read_tsv("data/intra_sp_variability/Intraspecific_contribution_perSite.txt") %>%
  rename(longitude = "collection_decimal_longitude", latitude = "collection_decimal_latitude", Ecosystem = "Ecosystem_Type")


##test merge environmental data (still need to transform metrics but lets test)
df <- left_join(site_is_df, env_df, by = c("collection_site_id", "Ecosystem"))




##C:N ratio -- look for effect especially in muscle, but may not be available 

##also looking for relationships between site environmental variables and community characteristic (e.g., species richness, functional diversity, )

##1) Description of variation -- how much intraspecific. variability contribution in C and N 
## -- also, do we see correlations between contribution of C and N - i.e., when more intrasp. var in C also see intrasp var in N? 

##2) Between site -- Relationship between environmental drivers and proportion of intraspecific variability 

##3) Between species -- drivers of dufferences in intraspecific variability between species -- trait drivers, community level drivers 
##potentially using site as a random effect 


##plotting distributions
hist_N_intra <- ggplot(df, aes(x = propintraspecific_N, group = Ecosystem)) +
  geom_histogram() +
  theme_classic() 
hist_N_intra

hist_C_intra <- ggplot(df, aes(x = propintraspecific_C)) +
  geom_histogram(fill = "darkcyan") +
  theme_classic()
hist_C_intra


##getting country names and continents for each site
countries <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_make_valid() %>%
  st_transform(4326)

sites <- st_as_sf(
  df,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

sites_sf <- sites %>%   # whatever your polygon sf object is called
  st_make_valid() %>%
  st_transform(4326)

sites_pt <- st_point_on_surface(sites_sf)


countries2 <- countries %>% st_transform(4326)
sites2 <- st_make_valid(sites_sf) %>% st_transform(4326)
sites_pt <- st_point_on_surface(sites2)

sites_geo <- st_join(
  sites_pt,
  countries2[, c("admin","iso_a3","continent")],
  join = st_intersects,
  left = TRUE
) %>% rename(country = admin)

na_idx <- which(is.na(sites_geo$country))
if (length(na_idx) > 0) {
  nearest <- st_nearest_feature(sites_geo[na_idx,], countries2)
  sites_geo$country[na_idx]   <- countries2$admin[nearest]
  sites_geo$iso_a3[na_idx]    <- countries2$iso_a3[nearest]
  sites_geo$continent[na_idx] <- countries2$continent[nearest]
}

test <- sites_geo %>%
  select(collection_site_id, Ecosystem, country, continent, propintraspecific_N, propintraspecific_C) %>%
  filter(is.na(continent))

sites_geo <- sites_geo %>%
  select(collection_site_id, Ecosystem, country, continent, propintraspecific_N, propintraspecific_C) %>%
  pivot_longer(cols= c(propintraspecific_N, propintraspecific_C), names_to = "prop_intraspecific_var_type", values_to = "prop_intraspecific_var") %>%
  mutate(prop_type = case_when(
    startsWith("propintraspecific_N", prop_intraspecific_var_type) ~ "N",
    startsWith("propintraspecific_C", prop_intraspecific_var_type) ~ "C",
  ))

##comparing variation in proportion 

#color1 <- c("darkred", "darkcyan")
color1 <- c("white", "darkgrey")
violin_c_n <- sites_geo  %>%
  ggplot(aes(x = Ecosystem, y = prop_intraspecific_var,  fill = prop_type)) +
  geom_violin() +
  scale_fill_manual(values = color1)+
  theme_classic() +
  ylab("Relative Contribution to Intraspecific Variability") +
  theme(axis.title.x = element_blank(), legend.position = "right")

violin_c_n

violin_c_n_2 <- sites_geo %>%
  ggplot(aes(x = prop_type, y = prop_intraspecific_var, fill = prop_type)) +
  geom_violin() +
  scale_fill_manual(values = color1)+
  theme_classic() +
  ylab("Relative Contribution to \nIntraspecific Variability") +
  theme(axis.title.x = element_blank(), legend.position = "none")

violin_c_n_2

box_c_n <- sites_geo %>%
  ggplot(aes(x = prop_type, y = prop_intraspecific_var, fill = Ecosystem)) +
  geom_boxplot() +
  scale_fill_manual(values = color1)+
  theme_classic() +
  ylab("Relative Contribution to \nIntraspecific Variability") +
  theme(axis.title.x = element_blank())

box_c_n

box_c_n_2 <- sites_geo %>%
  ggplot(aes(x = prop_type, y = prop_intraspecific_var, fill = prop_type)) +
  geom_boxplot() +
  scale_fill_manual(values = color1)+
  theme_classic() +
  ylab("Relative Contribution to \nIntraspecific Variability") +
  theme(axis.title.x = element_blank(), legend.position = "none")

box_c_n_2

box_c_n_bycontinent <- sites_geo %>%
  ggplot(aes(x = prop_type, y = prop_intraspecific_var, fill = prop_type)) +
  geom_boxplot() +
 scale_fill_manual(values = color1)+
  theme_classic() +
  ylab("Relative Contribution to \nIntraspecific Variability") +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  facet_wrap(~continent)

box_c_n_bycontinent

box_c_n_bycontinent_2 <- sites_geo %>%
  ggplot(aes(x = prop_type, y = prop_intraspecific_var, fill = continent)) +
  geom_boxplot() +
 # scale_fill_manual(values = color1)+
  theme_classic() +
  ylab("Relative Contribution to \nIntraspecific Variability") +
  theme(axis.title.x = element_blank(), legend.position = "right") 

box_c_n_bycontinent_2

##could be cool for each continent, to have a boxplot and then put that on the map 
make_box <- function(dat_ct, title = NULL) {
  ggplot(dat_ct, aes(x = prop_type, y = prop_intraspecific_var, fill = prop_type)) +
    geom_boxplot(width = 0.7, outlier.size = 0.6) +
    labs(title = title, x = NULL, y = NULL) +
   # theme_minimal(base_size = 9) +
    theme_classic()+
    theme(
      legend.position = "none",
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      panel.grid.minor = element_blank(),
      
      ## KEY: Make both backgrounds transparent
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background  = element_rect(fill = "transparent", color = NA)
    ) +
    ylim(0,1)
}

continents <- c("North America","South America","Europe","Africa","Asia","Oceania")

boxplots <- lapply(continents, function(ct) {
  dat_ct <- sites_geo |> st_drop_geometry() |> filter(continent == ct)
  make_box(dat_ct)
})
names(boxplots) <- continents



##save with transparent background
#ggsave("africa.png", plot = africa, bg = "transparent")



cn_correlation <- ggplot(df, aes(x = propintraspecific_C, y = propintraspecific_N)) +
  geom_point(color = "black") +
  geom_smooth(method = "lm", color = "darkred") +
  theme_classic() +
  ylab("Relative Contribution to \nIntraspecific Variability in N") +
  xlab("Relative Contribution to \nIntraspecific Variability in C") 

cn_correlation

lm1 <- lm(propintraspecific_N ~ propintraspecific_C, data = df)
summary(lm1)


plot_1 <- ggarrange(box_c_n_2, cn_correlation, legend = "none", nrow = 1, ncol = 2, labels = c("b)", "c)"), font.label = list(colour = "black", size = 12))
plot_1


ggplot(df, aes(x = log(site_nbspe), y = propintraspecific_C)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(df, aes(x = log(site_nbspe), y = propintraspecific_N)) +
  geom_point() +
  geom_smooth(method = "lm")




##need to keep in mind that species richness may bias intrapsecific variability -- could potentially ahndle by taking residuals and then using residuals from this relationships w/ envirnmental variable s 

##trying 
##Plotting out some maps 

sites_coord <- st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326)
site_map <- mapview(sites_coord, map.types = "Esri.NatGeoWorldMap", legend = FALSE)
site_map




##make map manuscript style

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

# Robinson projection (good-looking global)
crs_robin <- "+proj=robin"

world_r <- st_transform(world, crs_robin)
sites_r <- st_transform(sites_coord, crs_robin)

p_map <- ggplot() +
  geom_sf(data = world_r, fill = "lightgrey", color = "darkgrey", linewidth = 0.1) +
  geom_sf(data = sites_r, size = 1, alpha = 0.9) +
  coord_sf(expand = FALSE) +
 # annotation_scale(location = "bl", width_hint = 0.25) +
#  annotation_north_arrow(location = "bl", which_north = "true",
 #                        style = north_arrow_fancy_orienteering) +
  theme_void() +
  theme(
   # plot.margin = margin(5, 5, 5, 5),
    legend.position = "none"
  ) +
  theme(plot.background = element_rect(fill = "white", color = NA))

p_map

#ggsave("sites_global_map.pdf", p, width = 7.0, height = 4.2, units = "in")
#ggsave("sites_global_map.png", p, width = 7.0, height = 4.2, units = "in", dpi = 600)


##try putting boxplots on the map
positions <- tibble::tribble(
  ~continent,        ~x,   ~y,   ~w,   ~h,
  "North America",   0.05, 0.58, 0.14, 0.18,
  "South America",   0.17, 0.33, 0.14, 0.18,
  "Europe",          0.35, 0.61, 0.14, 0.18,
  "Africa",          0.41, 0.33, 0.14, 0.18,
  "Asia",            0.63, 0.53, 0.14, 0.18,
  "Oceania",         0.75, 0.26, 0.14, 0.18
)

p_final <- ggdraw(p_map)

for (i in seq_len(nrow(positions))) {
  ct <- positions$continent[i]
  p_final <- p_final +
    draw_plot(
      boxplots[[ct]],
      x = positions$x[i], y = positions$y[i],
      width = positions$w[i], height = positions$h[i]
    )
}

p_final






plot_2 <- ggarrange(p_map, plot_1, legend = "none", nrow = 2, ncol = 1, labels = c("a)", ""), font.label = list(colour = "black", size = 12))
plot_2

##getting country names and continents for each site

countries <- ne_countries(scale = "medium", returnclass = "sf")

sites <- st_as_sf(
  df,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

sites_geo <- st_join(
  sites,
  countries[, c("admin", "iso_a3", "continent")],
  join = st_within
) %>%
  rename(country = "admin")



##Do we see correlation between proportion intraspecifc variability in C and N? 
ggplot(df, aes(x = propintraspecific_C, y = propintraspecific_N)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

lm1 <- lm(propintraspecific_N ~ propintraspecific_C, data = df)
summary(lm1)

gamm1.1 <- gam(LC_pred_mean ~ summer_air_temp_mean + log_area_ha + log_lt_cue + s(sampling_year_mean, bs = "re"),
               data = df_summary, method = "REML")
summary(gamm1.1)
visreg(gamm1.1, "summer_air_temp_mean", partial = TRUE)
visreg(gamm1.1, "log_area_ha", partial = TRUE)
visreg(gamm1.1, "log_lt_cue", partial = TRUE)

lm3 <- lm(SEAc ~ summer_air_temp + log_area_ha + log_pred_cue, data = siber_df_final)
summary(lm3)



##add species richness as a random effect

###OLD --------------------------------
#Map with no color indication just sites 
##Plot of proportion intraspecific variability in C
leaflet(sites_coord) |>
  addProviderTiles(providers$CartoDB.PositronNoLabels) |>
  addCircleMarkers(radius = 2,  color = "black", fillColor = "black",  fillOpacity = 0.9,
                   weight = 1,
                   popup = ~paste0("propintraspecific_C: ", signif(vals, 3))
  ) %>%
  addScaleBar(position = c("bottomleft"))


# color palette 
pal <- colorNumeric(
  palette = "viridis",
  domain  = vals,
  na.color = "lightgray"
)

pal_1 <- colorNumeric(
  palette = c("cyan1", "darkcyan", "darkblue"),   # <–– choose any colors
  domain  = vals,
  na.color = "lightgray"
)
##Plot of proportion intraspecific variability in C
leaflet(sites_coord) |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addCircleMarkers(radius = 4,  color = ~pal_1(vals), fillColor = ~pal_1(vals),  fillOpacity = 0.9,
                   weight = 1,
                   popup = ~paste0("propintraspecific_C: ", signif(vals, 3))
  ) |>
  addLegend(
    "bottomright",
    pal = pal_1,
    values = vals,
    title = "Intraspecific\nproportion\nC"
  )


##Plot of proportion intraspecific variability in N
pal_2 <- colorNumeric(
  palette = c("lightpink", "pink4", "darkred"),   # <–– choose any colors
  domain  = vals,
  na.color = "lightgray"
)
leaflet(sites_coord) |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addCircleMarkers(radius = 4,  color = ~pal_2(vals), fillColor = ~pal_2(vals),  fillOpacity = 0.9,
                   weight = 1,
                   popup = ~paste0("propintraspecific_N: ", signif(vals, 3))
  ) |>
  addLegend(
    "bottomright",
    pal = pal_2,
    values = vals,
    title = "Intraspecific\nproportion\nN"
  )

