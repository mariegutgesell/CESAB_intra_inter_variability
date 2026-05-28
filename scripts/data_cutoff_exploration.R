##Testing different thresholds 



library(tidyverse)
library(ggplot2)
library(purrr)

##Read in min1 dataset -- this has minimum 1 replicate per species, and at least 2-3 species (I think at least 3 species is still a reasonable cutoff)
##Want to look at:
## Number of food webs retained at different cut-offs
## Number of species retained at different cut-offs
## Average number of individuals per species 



##Step 1: get data to 1 year of sampling 
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


##Check for and remove temporal replication within site -- want to ensure only 1 year of data
DataFish <- DataFish %>%
  # mutate(as.Date(collection_date)) %>%
  mutate(is_range = str_detect(collection_date, "^\\d{4}-\\d{4}$")) %>% ##indicate if sample collection date is a range of years
  mutate(year = str_extract(collection_date, "^\\d{4}") %>% as.numeric()) ##create column of year

##calculate the number of sampling years
site_num_sampling_years <- DataFish %>%
  filter(is_range == FALSE) %>% ##remove any sites that only had collection date listed as year ranges 
  select(FWB_id, collection_site_id, year) %>%
  unique() %>%
  group_by(FWB_id, collection_site_id) %>%
  count() %>%
  rename(number_sampling_years = "n")

##number of food webs with more than 1 year of smapling 
num_sites_more1yr <- site_num_sampling_years %>%
  filter(number_sampling_years > 1)

DataFish <- left_join(DataFish, site_num_sampling_years, by = c("FWB_id", "collection_site_id"))


datafish_summary_sp <- DataFish %>%
  group_by(FWB_id, collection_site_id, year, number_sampling_years, scientific_name) %>%
  count() %>%
  rename(num_ind_per_species = "n")


datafish_summary_site <- datafish_summary_sp %>%
  group_by(FWB_id, collection_site_id, year, number_sampling_years) %>%
  summarise(num_species = n_distinct(scientific_name),
            mean_num_ind_per_sample = mean(num_ind_per_species),
            min_num_ind_per_sample = min(num_ind_per_species))


##Annual resolution data - decision tree
annual_site_list <- datafish_summary_site %>%
  # apply decision tree within each site
  group_by(FWB_id, collection_site_id) %>%
  arrange(
    desc(num_species),               # most species first
    desc(mean_num_ind_per_sample),   # then largest mean sample size
    desc(year),                      # then most recent
    .by_group = TRUE
  ) %>%
  
  # keep top-ranked row
  slice(1) %>%
  ungroup() %>%
  filter(!is.na(year)) %>% ##remove sites with year as NA 
  mutate(site_year_code = paste(FWB_id, year, sep = "_"))

##945 sites 


##Select the site and year from datafish
datafish_annual <- DataFish %>%
  mutate(site_year_code = paste(FWB_id, year, sep = "_")) %>%
  filter(site_year_code %in% annual_site_list$site_year_code) %>%
  left_join(datafish_summary_site, by = c("FWB_id", "collection_site_id", "year", "number_sampling_years")) 

##okay lets make a min 2 and min3 species 
datafish_annual_min2sp <- datafish_annual %>%
  filter(num_species >=2)  

num_sites_min2sp <- datafish_annual_min2sp %>%
  ungroup() %>%
  select(FWB_id, collection_site_id, year) %>%
  distinct()
##849 food webs 

datafish_summary_sp_min2sp <- datafish_annual_min2sp %>%
  group_by(FWB_id, collection_site_id, year, site_year_code, number_sampling_years, scientific_name) %>%
  count() %>%
  rename(num_ind_per_species = "n")

head(datafish_summary_sp_min2sp)

datafish_annual_min3sp <- datafish_annual %>%
  filter(num_species >=3) 

num_sites_min3sp <- datafish_annual_min3sp %>%
  ungroup() %>%
  select(FWB_id, collection_site_id, year) %>%
  distinct()
##681 food webs 

datafish_summary_sp_min3sp <- datafish_annual_min3sp%>%
  group_by(FWB_id, collection_site_id, year, site_year_code, number_sampling_years, scientific_name) %>%
  count() %>%
  rename(num_ind_per_species = "n")


##Calculate thresholds:
site_thresholds_2sp <- datafish_summary_sp_min2sp %>%
  group_by(FWB_id, collection_site_id, year, site_year_code) %>%
  summarise(
    total_species = n_distinct(scientific_name),
    species_ge3 = sum(num_ind_per_species >= 3, na.rm = TRUE),
    prop_ge3 = species_ge3 / total_species,
    .groups = "drop"
  )

cutoffs <- seq(0.0, 1, by = 0.05)

species_datasets_2sp <- map(
  cutoffs,
  ~ datafish_annual_min2sp %>%
    semi_join(
      site_thresholds_2sp %>%
        filter(prop_ge3 >= .x),
      by = "site_year_code"
    )
)

names(species_datasets_2sp) <- paste0(cutoffs * 100, "pct")

##Do for 3 species 
site_thresholds_3sp <- datafish_summary_sp_min3sp %>%
  group_by(FWB_id, collection_site_id, year, site_year_code) %>%
  summarise(
    total_species = n_distinct(scientific_name),
    species_ge3 = sum(num_ind_per_species >= 3, na.rm = TRUE),
    prop_ge3 = species_ge3 / total_species,
    .groups = "drop"
  )

cutoffs <- seq(0.0, 1, by = 0.05)

species_datasets_3sp <- map(
  cutoffs,
  ~ datafish_annual_min3sp %>%
    semi_join(
      site_thresholds_3sp %>%
        filter(prop_ge3 >= .x),
      by = "site_year_code"
    )
)

names(species_datasets_3sp) <- paste0(cutoffs * 100, "pct")





##Now, create function that for each dataset, summarizes table: number of food webs, number of species, mean sample size , min sample size 

summarise_cutoff_dataset <- function(df, cutoff_name) {
  
  overall_summary <- df %>%
    ungroup() %>%
    summarise(
      num_food_webs = n_distinct(FWB_id),
      total_num_species = n_distinct(scientific_name)
    )
  
  df %>%
    ungroup() %>%
    group_by(FWB_id, collection_site_id, year, number_sampling_years, scientific_name) %>%
    count(name = "num_ind_per_species") %>%
    group_by(FWB_id, collection_site_id, year, number_sampling_years) %>%
    summarise(
      num_species = n_distinct(scientific_name),
      mean_num_ind_per_sample = mean(num_ind_per_species, na.rm = TRUE),
      min_num_ind_per_sample = min(num_ind_per_species, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      cutoff = cutoff_name,
      cutoff_num = as.numeric(gsub("pct", "", cutoff_name)),
      num_food_webs = overall_summary$num_food_webs,
      total_num_species = overall_summary$total_num_species
    )
}
  
summary_all_cutoffs_2sp <- imap_dfr(
  species_datasets_2sp,
  ~ summarise_cutoff_dataset(.x, .y)
) %>%
  mutate(min_sp_num = 2)


summary_all_cutoffs_3sp <- imap_dfr(
  species_datasets_3sp,
  ~ summarise_cutoff_dataset(.x, .y)
)%>%
  mutate(min_sp_num = 3)

summary_all_cutoffs <- rbind(summary_all_cutoffs_2sp, summary_all_cutoffs_3sp)

###Make some summary plots 
ggplot(summary_all_cutoffs, aes(x = cutoff_num, y = num_food_webs, group = min_sp_num, color = factor(min_sp_num))) +
  geom_point() +
  geom_line() +
  ylab("Number of Food Webs") +
  xlab("% Cutoff") +
  labs(color = "Minimum Number of Species in Food Web")


ggplot(summary_all_cutoffs, aes(x = cutoff_num, y = total_num_species, group = min_sp_num, color = factor(min_sp_num))) +
  geom_point() +
  geom_line() +
  ylab("Total Number of Unique Taxa") +
  xlab("% Cutoff") +
  labs(color = "Minimum Number of Species in Food Web")

ggplot(summary_all_cutoffs, aes(x = factor(cutoff_num), y = mean_num_ind_per_sample, fill = factor(min_sp_num))) +
  geom_boxplot() +
  ylab("Mean Number of Samples per Species within Food Web") +
  xlab("% Cutoff") +
  labs(fill = "Minimum Number of Species in Food Web")

ggplot(summary_all_cutoffs, aes(x = factor(cutoff_num), y = min_num_ind_per_sample, fill = factor(min_sp_num))) +
  geom_boxplot() +
  ylab("Min Number of Samples per Species within Food Web") +
  xlab("% Cutoff") +
  labs(fill = "Minimum Number of Species in Food Web")

ggplot(summary_all_cutoffs, aes(x = factor(cutoff_num), y = num_species, fill = factor(min_sp_num))) +
  geom_boxplot()+
  ylab("Number of Species per Food Web") +
  xlab("% Cutoff") +
  labs(fill = "Minimum Number of Species in Food Web")



##tomorrow morning: code add to script of 02 essentially creating the site and spvar csvs for all cutoffs above like 75%... and for the 3 species min i think - that way can easily send over ? 
##could also then run on the 75% and see if results are different ... 

##Make site list and save -- this will then be easily used to filter the selected sitevar/sp var 

write.csv(summary_all_cutoffs, "data/cuttoff_site_lists.csv", row.names = FALSE)



##Okay so -- current cutoff, 2 species, 75% -- for those 25% w/ less than 3 replicates, what proportion of total number of individuals in that food web are they ? 
