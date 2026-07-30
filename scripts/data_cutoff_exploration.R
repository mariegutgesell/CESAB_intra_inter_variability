##Code to generate cutoff lists based on thresholds of species and sampling number 


library(tidyverse)
library(ggplot2)
library(purrr)
library(readxl)
##Want to look at:
## Number of food webs retained at different cut-offs
## Number of species retained at different cut-offs
## Average number of individuals per species 

##Step 1: load/clean final individual dataset and filter to get data to 1 year of sampling per site ------------
df_all <- read_excel("data/FINAL_ALLindiv_February2026.xlsx") 

##clean data
##select only fish, where have C/N data and associated scientific name and NA  for species name 
#DataFish<-subset(df_all,!is.na(d15N) & !is.na(d13C) & organism_type=="fish" & fish_species != "NA")
##note: when we remove fish_species is NA we are removing any fish that are identified to only to genus level - this can create false diversity and intravar. because we are removing individuals from different speceis from a community .. 
#test <- subset(df_all,!is.na(d15N) & !is.na(d13C) & organism_type=="fish" & fish_species == "NA")

##Going to try an alternative way: 
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

##confirm that no NA species names are retained
test <- DataFish %>%
  filter(fish_scientific_name == "NA")
test <- DataFish %>%
  filter(is.na(fish_scientific_name))
##NAs come from cases where fish were not idenfied (e.g., mix of rudd and mosquitofish)

##add species-site identifier column
DataFish$sp_site<-paste(DataFish$fish_scientific_name,DataFish$collection_site_id,sep="_")
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

##number of food webs with more than 1 year of sampling 
num_sites_more1yr <- site_num_sampling_years %>%
  filter(number_sampling_years > 1)

##join number of sampling years to main fish df 
DataFish <- left_join(DataFish, site_num_sampling_years, by = c("FWB_id", "collection_site_id"))


##calculate the number of individuals sampled per species within each site-year 
datafish_summary_sp <- DataFish %>%
  group_by(FWB_id, collection_site_id, year, number_sampling_years, fish_scientific_name) %>%
  count() %>%
  rename(num_ind_per_species = "n")

##calculate richness and mean/min of # of individuals sampled per species within each site-year
datafish_summary_site <- datafish_summary_sp %>%
  group_by(FWB_id, collection_site_id, year, number_sampling_years) %>%
  summarise(num_species = n_distinct(fish_scientific_name),
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
##936 when NA sp name are removed (fewer sites because some sites only have genus level ID)

##Select the site and year from datafish
datafish_annual <- DataFish %>%
  mutate(site_year_code = paste(FWB_id, year, sep = "_")) %>%
  filter(site_year_code %in% annual_site_list$site_year_code) %>%
  left_join(datafish_summary_site, by = c("FWB_id", "collection_site_id", "year", "number_sampling_years")) 


##Step 2: generate site lists at different minimum species richness and mean sample size cutoffs ----------

## make a min 2 species dataset
datafish_annual_min2sp <- datafish_annual %>%
  filter(num_species >=2)  

num_sites_min2sp <- datafish_annual_min2sp %>%
  ungroup() %>%
  select(FWB_id, collection_site_id, year) %>%
  distinct()
##849 food webs 
##840 when NA sp names removed

##calculate number of individuals per species within each site
datafish_summary_sp_min2sp <- datafish_annual_min2sp %>%
  group_by(FWB_id, collection_site_id, year, site_year_code, number_sampling_years, fish_scientific_name) %>%
  count() %>%
  rename(num_ind_per_species = "n")

head(datafish_summary_sp_min2sp)

## make a min 3 species dataset
datafish_annual_min3sp <- datafish_annual %>%
  filter(num_species >=3) 

num_sites_min3sp <- datafish_annual_min3sp %>%
  ungroup() %>%
  select(FWB_id, collection_site_id, year) %>%
  distinct()
##681 food webs 
##661 when NA sp names removed

##calculate number of individuals per species within each site
datafish_summary_sp_min3sp <- datafish_annual_min3sp%>%
  group_by(FWB_id, collection_site_id, year, site_year_code, number_sampling_years, fish_scientific_name) %>%
  count() %>%
  rename(num_ind_per_species = "n")


##Calculate thresholds:(i.e., where x% of species within a site have at least 3 individuals sampled)
site_thresholds_2sp <- datafish_summary_sp_min2sp %>%
  group_by(FWB_id, collection_site_id, year, site_year_code) %>%
  summarise(
    total_species = n_distinct(fish_scientific_name),
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
    total_species = n_distinct(fish_scientific_name),
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



##Create function that for each dataset, summarizes table: number of food webs, number of species, mean sample size , min sample size 

summarise_cutoff_dataset <- function(df, cutoff_name) {
  
  overall_summary <- df %>%
    ungroup() %>%
    summarise(
      num_food_webs = n_distinct(FWB_id),
      total_num_species = n_distinct(fish_scientific_name)
    )
  
  df %>%
    ungroup() %>%
    group_by(FWB_id, collection_site_id, year, number_sampling_years, fish_scientific_name) %>%
    count(name = "num_ind_per_species") %>%
    group_by(FWB_id, collection_site_id, year, number_sampling_years) %>%
    summarise(
      num_species = n_distinct(fish_scientific_name),
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


##Make site list and save -- this will then be easily used to filter the selected sitevar/sp var 
write.csv(summary_all_cutoffs, "data/cuttoff_site_lists.csv", row.names = FALSE)

##how does this compare to old cutoffs? ---- 
cutoff_old <- read.csv("data/cuttoff_site_lists_old.csv")

##Okay so -- current cutoff, 2 species, 75% -- for those 25% w/ less than 3 replicates, what proportion of total number of individuals in that food web are they ? 

old_sites <- cutoff_old %>%
  mutate(site_year_code = paste(FWB_id, year, sep = "_")) %>%
  filter(min_sp_num == 2) %>% ##select minimum species number
  filter(cutoff == "75pct") %>% ##select cutoff 
  select(site_year_code, num_species, mean_num_ind_per_sample, min_num_ind_per_sample) %>%
  mutate(site_type = "old")
  #453 sites 

new_sites <- summary_all_cutoffs %>%
  mutate(site_year_code = paste(FWB_id, year, sep = "_")) %>%
  filter(min_sp_num == 2) %>% ##select minimum species number
  filter(cutoff == "75pct") %>%##select cutoff 
  select(site_year_code, num_species, mean_num_ind_per_sample, min_num_ind_per_sample) %>%
  mutate(site_type = "new")

  ##455 sites 


##check overlap 
site_overlap <- full_join(old_sites, new_sites, by = c("site_year_code", "num_species", "mean_num_ind_per_sample", "min_num_ind_per_sample", "site_type"))


test <- site_overlap %>%
  group_by(site_year_code) %>%
  count()
##okay, so yes i think this makes sense, that 

##Looking at mixed id levels:
##note: need to read in and filter to selected sites first, just keeping code here so don't lose it 
##Look at within these sites, how many genus/family/order are there? 
fam_ord_sites <- DataFish_final %>%
  filter(fish_name_level %in% c("family", "order")) %>%
  select(FWB_id) %>%
  distinct()
##31 sites that have family and order level ID 

genus_sites <- DataFish_final %>%
  filter(fish_name_level %in% c("genus")) %>%
  select(FWB_id) %>%
  distinct()
##69 sites have genus level ID 

mixed_taxon_sites <- DataFish_final %>%
  filter(fish_name_level %in% c("genus", "family", "order")) %>%
  select(FWB_id) %>%
  distinct()
##Create df of sites that contain mix of taxonomic id levels 
mixed_taxon_df <- DataFish_final %>%
  filter(FWB_id %in% mixed_taxon_sites$FWB_id)

write.csv(mixed_taxon_df, "data/mixed_taxon_sites_df.csv")
richness_by_level <- DataFish_final %>%
  group_by(FWB_id, fish_name_level) %>%
  summarise(
    richness = n_distinct(fish_scientific_name),
    .groups = "drop"
  )

total_richness <- DataFish_final %>%
  group_by(FWB_id) %>%
  summarise(
    total_richness = n_distinct(fish_scientific_name),
    .groups = "drop"
  )

richness_prop <- richness_by_level %>%
  left_join(total_richness, by = "FWB_id") %>%
  mutate(
    prop_richness = richness / total_richness,
    percent_richness = 100 * prop_richness
  )

richness_prop

richness_prop %>%
  filter(fish_name_level %in% c("genus", "family", "order")) %>%
  group_by(fish_name_level) %>%
  summarise(
    n_sites = n(),
    mean_percent = mean(percent_richness),
    median_percent = median(percent_richness),
    min_percent = min(percent_richness),
    max_percent = max(percent_richness)
  )


