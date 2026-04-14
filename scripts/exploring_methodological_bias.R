##Exploring potential methodological bias 

library(tidyverse)
library(ggplot2)


#load("data/Env_27Mar26.RData") 

##Potential approach to exploring methodological bias:

##1. Resolve temporal pooling as much as possible.
##2. Rarefy well-sampled species-site groups to examine variance stabilization.
##3. Rarefy whole sites with equal per-species draws and recalculate the partitioning metric.
##4. Run sensitivity analyses across minimum n thresholds.
##5. Add sampling structure covariates to the model: species richness, mean n per species, minimum n per species, proportion of species with n above threshold
##6. Check whether major conclusions persist.

##is there some type of bayesian approach we can use that accounts for small sample size? because otherwise we are going to throw out a lot of data ... 

df_all <- read_excel("data/FINAL_ALLindiv_February2026.xlsx") 

load("data/Intraspecific_contribution_perSite_Env.RData") ##produced by 2_Intra_vs_Inter_21_01_26.R


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


##create a species-site sample sizes from individual-level data
sp_site_counts <- DataFish %>%
  group_by(sp_site, collection_site_id, fish_species) %>%
  summarise(n_ind = n(), .groups = "drop")

# choose the subset for rarefaction
sample_sizes <- c(3, 5, 10, 20)

sp_site_subset <- sp_site_counts %>%
  filter(n_ind >= max(sample_sizes))

# keep only those species-site groups in the individual-level data
DataFish_subset <- DataFish %>%
  semi_join(sp_site_subset, by = c("sp_site", "collection_site_id", "fish_species"))


rarefy_variance <- function(x, sample_sizes = c(3, 5, 10, 20), n_reps = 200) {
  x <- x[!is.na(x)]
  full_n <- length(x)
  full_var <- var(x, na.rm = TRUE)
  
  map_dfr(sample_sizes[sample_sizes <= full_n], function(n_sub) {
    sub_vars <- replicate(n_reps, {
      x_sub <- sample(x, size = n_sub, replace = FALSE)
      var(x_sub, na.rm = TRUE)
    })
    
    tibble(
      n_sub = n_sub,
      full_n = full_n,
      full_var = full_var,
      mean_sub_var = mean(sub_vars, na.rm = TRUE),
      median_sub_var = median(sub_vars, na.rm = TRUE),
      sd_sub_var = sd(sub_vars, na.rm = TRUE),
      q025_sub_var = quantile(sub_vars, 0.025, na.rm = TRUE),
      q975_sub_var = quantile(sub_vars, 0.975, na.rm = TRUE),
      rel_mean_var = mean_sub_var / full_var,
      rel_median_var = median_sub_var / full_var
    )
  })
}


# d13C rarefaction
rare_d13C <- DataFish_subset %>%
  group_by(sp_site, collection_site_id, fish_species) %>%
  summarise(values = list(d13C_norm), .groups = "drop") %>%
  mutate(rarefied = map(values, rarefy_variance, sample_sizes = sample_sizes, n_reps = 200)) %>%
  select(-values) %>%
  unnest(rarefied) %>%
  mutate(isotope = "d13C")

# d15N rarefaction
rare_d15N <- DataFish_subset %>%
  group_by(sp_site, collection_site_id, fish_species) %>%
  summarise(values = list(d15N_norm), .groups = "drop") %>%
  mutate(rarefied = map(values, rarefy_variance, sample_sizes = sample_sizes, n_reps = 200)) %>%
  select(-values) %>%
  unnest(rarefied) %>%
  mutate(isotope = "d15N")

# combine
rare_results <- bind_rows(rare_d13C, rare_d15N)


rare_summary <- rare_results %>%
  group_by(isotope, n_sub) %>%
  summarise(
    n_groups = n(),
    mean_rel_var = mean(rel_mean_var, na.rm = TRUE),
    median_rel_var = median(rel_median_var, na.rm = TRUE),
    q25_rel_var = quantile(rel_median_var, 0.25, na.rm = TRUE),
    q75_rel_var = quantile(rel_median_var, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

rare_summary


ggplot(rare_summary, aes(x = n_sub, y = median_rel_var)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = q25_rel_var, ymax = q75_rel_var), width = 0.3) +
  facet_wrap(~ isotope) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    x = "Rarefied sample size",
    y = "Relative variance estimate (rarefied / full)",
    title = "Effect of sample size on within-species isotope variance"
  ) +
  theme_bw()


##Site-level rarefaction: ------------
#For each site, standardize sampling to exactly k individuals per species, where k = 3, 5, 8, 10, then recalculate the site-level inter- and intraspecific variance metrics and the proportion due to intraspecific variability

##function to rarefy one site at one k
rarefy_site_once <- function(df_site, k, min_species = 3) {
  
  # keep only species with at least k individuals in this site
  df_keep <- df_site %>%
    group_by(fish_species) %>%
    filter(n() >= k) %>%
    ungroup()
  
  # count retained species
  n_species <- df_keep %>%
    distinct(fish_species) %>%
    nrow()
  
  # if too few species remain, return NA
  if (n_species < min_species) {
    return(tibble(
      collection_site_id = unique(df_site$collection_site_id)[1],
      k = k,
      site_nbspe = NA_integer_,
      site_interspe_var_N = NA_real_,
      site_intraspe_var_N = NA_real_,
      site_interspe_var_C = NA_real_,
      site_intraspe_var_C = NA_real_,
      propintraspecific_N = NA_real_,
      propintraspecific_C = NA_real_,
      propintraspecific_Total = NA_real_
    ))
  }
  
  # subsample exactly k individuals per species
  df_sub <- df_keep %>%
    group_by(fish_species) %>%
    slice_sample(n = k, replace = FALSE) %>%
    ungroup()
  
  # species-level summaries within the rarefied site
  sp_summary <- df_sub %>%
    group_by(collection_site_id, fish_species) %>%
    summarise(
      sp_site_mean_N = mean(d15N_norm, na.rm = TRUE),
      sp_site_var_N  = var(d15N_norm, na.rm = TRUE),
      sp_site_mean_C = mean(d13C_norm, na.rm = TRUE),
      sp_site_var_C  = var(d13C_norm, na.rm = TRUE),
      sp_site_num_ind = n(),
      .groups = "drop"
    )
  
  # site-level summaries
  site_summary <- sp_summary %>%
    summarise(
      collection_site_id = first(collection_site_id),
      k = k,
      site_nbspe = n(),
      site_interspe_var_N = var(sp_site_mean_N, na.rm = TRUE),
      site_intraspe_var_N = mean(sp_site_var_N, na.rm = TRUE),
      site_interspe_var_C = var(sp_site_mean_C, na.rm = TRUE),
      site_intraspe_var_C = mean(sp_site_var_C, na.rm = TRUE)
    ) %>%
    mutate(
      propintraspecific_N = site_intraspe_var_N / (site_interspe_var_N + site_intraspe_var_N),
      propintraspecific_C = site_intraspe_var_C / (site_interspe_var_C + site_intraspe_var_C),
      propintraspecific_Total = (site_intraspe_var_N + site_intraspe_var_C) /
        (site_interspe_var_N + site_intraspe_var_N + site_interspe_var_C + site_intraspe_var_C)
    )
  
  return(site_summary)
}

##function to repeat many times for one site and one k
rarefy_site_reps <- function(df_site, k, n_reps = 100, min_species = 3) {
  map_dfr(seq_len(n_reps), function(i) {
    rarefy_site_once(df_site, k = k, min_species = min_species) %>%
      mutate(rep = i)
  })
}

##apply across al sites and all rarefaction levels

# split individual-level data by site
site_list <- DataFish %>%
  group_by(collection_site_id) %>%
  group_split()

# rarefaction levels to test
k_values <- c(3, 5, 8, 10)

# run rarefaction across all sites and k values
site_rarefaction_results <- map_dfr(k_values, function(k_val) {
  map_dfr(site_list, function(df_site) {
    rarefy_site_reps(df_site, k = k_val, n_reps = 10, min_species = 3)
  })
})

##summarize across replicates for each site and k
site_rarefaction_summary <- site_rarefaction_results %>%
  group_by(collection_site_id, k) %>%
  summarise(
    site_nbspe = median(site_nbspe, na.rm = TRUE),
    
    propintraspecific_N_med = median(propintraspecific_N, na.rm = TRUE),
    propintraspecific_N_q25 = quantile(propintraspecific_N, 0.25, na.rm = TRUE),
    propintraspecific_N_q75 = quantile(propintraspecific_N, 0.75, na.rm = TRUE),
    
    propintraspecific_C_med = median(propintraspecific_C, na.rm = TRUE),
    propintraspecific_C_q25 = quantile(propintraspecific_C, 0.25, na.rm = TRUE),
    propintraspecific_C_q75 = quantile(propintraspecific_C, 0.75, na.rm = TRUE),
    
    propintraspecific_Total_med = median(propintraspecific_Total, na.rm = TRUE),
    propintraspecific_Total_q25 = quantile(propintraspecific_Total, 0.25, na.rm = TRUE),
    propintraspecific_Total_q75 = quantile(propintraspecific_Total, 0.75, na.rm = TRUE),
    
    .groups = "drop"
  )


##summarize across all sites for each k
overall_rarefaction_summary <- site_rarefaction_summary %>%
  group_by(k) %>%
  summarise(
    n_sites = n(),
    
    median_prop_N = median(propintraspecific_N_med, na.rm = TRUE),
    q25_prop_N = quantile(propintraspecific_N_med, 0.25, na.rm = TRUE),
    q75_prop_N = quantile(propintraspecific_N_med, 0.75, na.rm = TRUE),
    
    median_prop_C = median(propintraspecific_C_med, na.rm = TRUE),
    q25_prop_C = quantile(propintraspecific_C_med, 0.25, na.rm = TRUE),
    q75_prop_C = quantile(propintraspecific_C_med, 0.75, na.rm = TRUE),
    
    median_prop_Total = median(propintraspecific_Total_med, na.rm = TRUE),
    q25_prop_Total = quantile(propintraspecific_Total_med, 0.25, na.rm = TRUE),
    q75_prop_Total = quantile(propintraspecific_Total_med, 0.75, na.rm = TRUE),
    
    .groups = "drop"
  )

overall_rarefaction_summary

site_rarefaction_summary %>%
  count(k)


##plot site-level rarefaction curves
library(ggplot2)

ggplot(overall_rarefaction_summary, aes(x = k, y = median_prop_Total)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = q25_prop_Total, ymax = q75_prop_Total), width = 0.3) +
  labs(
    x = "Individuals sampled per species",
    y = "Median proportion of intraspecific variability",
    title = "Site-level rarefaction of intraspecific contribution"
  ) +
  theme_bw()


ggplot(overall_rarefaction_summary, aes(x = k, y = median_prop_C)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = q25_prop_C, ymax = q75_prop_C), width = 0.3) +
  labs(
    x = "Individuals sampled per species",
    y = "Median proportion of intraspecific variability of d13C",
    title = "Site-level rarefaction of intraspecific contribution"
  ) +
  theme_bw()

ggplot(overall_rarefaction_summary, aes(x = k, y = median_prop_N)) +
  geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin = q25_prop_N, ymax = q75_prop_N), width = 0.3) +
  labs(
    x = "Individuals sampled per species",
    y = "Median proportion of intraspecific variability of d15N",
    title = "Site-level rarefaction of intraspecific contribution"
  ) +
  theme_bw()


##can also make for median prop C and prop N 

##compare rarefied vs. original site metrics
compare_full_vs_rare <- site_rarefaction_summary %>%
  left_join(
    SiteVar %>%
      select(collection_site_id,
             propintraspecific_N,
             propintraspecific_C,
             propintraspecific_Total,
             site_nbspe,
             site_mean_sample_id,
             site_min_sample_id),
    by = "collection_site_id"
  )

ggplot(compare_full_vs_rare %>% filter(k == 5),
       aes(x = propintraspecific_Total, y = propintraspecific_Total_med)) +
  geom_point(alpha = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Full-data proportion",
    y = "Rarefied proportion (k = 5)",
    title = "Full vs rarefied site-level intraspecific contribution"
  ) +
  theme_bw()


retention_summary <- site_rarefaction_results %>%
  group_by(k, collection_site_id) %>%
  summarise(site_nbspe = median(site_nbspe, na.rm = TRUE), .groups = "drop") %>%
  group_by(k) %>%
  summarise(
    n_sites = sum(!is.na(site_nbspe)),
    median_species_retained = median(site_nbspe, na.rm = TRUE),
    q25_species_retained = quantile(site_nbspe, 0.25, na.rm = TRUE),
    q75_species_retained = quantile(site_nbspe, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

retention_summary

#I. Intraspecific variation per species------------
#some species have var = NA when only one individual sampled
#we ignore such species in average intra var but still include it in inter var
colnames(DataFish)
str(DataFish)
#levels(DataFish$Diet)
DataFish$collected_sample_length_mm <- as.numeric(DataFish$collected_sample_length_mm)

SpVar<-DataFish %>% 
  group_by(sp_site,collection_site_id,fish_species,fish_family,
           waterbody_type, ecosystem_area_km2, ecosystem_width_m, 
           collection_decimal_longitude, collection_decimal_latitude) %>% 
  summarise(sp_site_mean_N = mean(d15N_norm, na.rm = TRUE),
            sp_site_var_N = var(d15N_norm, na.rm = TRUE),
            sp_site_mean_C = mean(d13C_norm, na.rm = TRUE),
            sp_site_var_C = var(d13C_norm, na.rm = TRUE),
            sp_site_mean_length = mean(collected_sample_length_mm, na.rm = TRUE), ##do we want to only use fork or total length? or group by this? 
            sp_site_var_length = var(collected_sample_length_mm, na.rm = TRUE),
            collection_decimal_longitude = mean(collection_decimal_longitude, na.rm = TRUE),
            collection_decimal_latitude = mean(collection_decimal_latitude, na.rm = TRUE),
            sp_site_num_ind = sum(num)) ##number of samples per species

SpVar$VarTot<-SpVar$sp_site_var_N+SpVar$sp_site_var_C
colnames(SpVar)
SpVar$collection_decimal_longitude


##How many sites do we have where the minimum sample size is less than 3 individuals? 
SpVar_min3 <- SpVar %>%
  group_by(collection_site_id) %>%
  filter(all(sp_site_num_ind >= 3)) ##change this to only filter sites where ALL species have at least 3 individuals ... 


test <- SpVar_min3 %>%
  ungroup() %>%
  select(waterbody_type, collection_site_id) %>%
  unique() %>%
  group_by(waterbody_type) %>%
  count()
##this does still give us 382 food webs 

SpVar_min5 <- SpVar %>%
  group_by(collection_site_id) %>%
  filter(all(sp_site_num_ind >= 5))

test <- SpVar_min5 %>%
  ungroup() %>%
  select(waterbody_type, collection_site_id) %>%
  unique() %>%
  group_by(waterbody_type) %>%
  count()

SpVar_min8 <- SpVar %>%
  group_by(collection_site_id) %>%
  filter(all(sp_site_num_ind >= 8))

test <- SpVar_min8 %>%
  select(waterbody_type, collection_site_id) %>%
  unique()

##but then could we only use sites where ALL the species have at least 3 individuals? 


##Try rarefraction within species approach -- 
##1) randomly sample n = 3, n = 5, n = 10, n = 20 
##2) calculate mean/var within site, species
##3) calculate proportion contribution from 


#II. Intraspecific vs. interspecific variation per site with environemental data
colnames(SpVar)
colnames(SpVar_env)
SpVar_env$num<-1

SiteVar <- SpVar_env %>%
  group_by(
    collection_site_id,
    waterbody_type,
    collection_decimal_latitude,
    collection_decimal_longitude,
    ecosystem_area_km2, 
    ecosystem_width_m
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
  )


##Calculate proportion of intraspecific variation from C and N variance
SiteVar$propintraspecific_N<-SiteVar$site_intraspe_var_N/(SiteVar$site_interspe_var_N+SiteVar$site_intraspe_var_N)
SiteVar$propintraspecific_C<-SiteVar$site_intraspe_var_C/(SiteVar$site_interspe_var_C+SiteVar$site_intraspe_var_C)
SiteVar$propintraspecific_Total<-(SiteVar$site_intraspe_var_N+SiteVar$site_intraspe_var_C)/(SiteVar$site_interspe_var_N+SiteVar$site_intraspe_var_N+SiteVar$site_interspe_var_C+SiteVar$site_intraspe_var_C)

