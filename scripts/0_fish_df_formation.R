##Generate site list and set up df for analysis 

library(tidyverse)
library(ggplot2)
library(purrr)
library(readxl)


##Step 1: load/clean final individual dataset and filter to get data to 1 year of sampling per site ------------
df_all <- read_excel("data/FINAL_ALLindiv_February2026.xlsx") 

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

##1) read in Julians sorting and clean up accordingly
julien_site_sorting <- read_excel("data/Sites_Julien_27.07.2026.xlsx") %>%
  select(FWB_id, Use, Reason, Action, Action_description) %>%
  distinct(FWB_id, .keep_all = TRUE)
head(julien_site_sorting)

##number of sites with mixed taxonomy levels that were included
julien_site_sorting %>%
  filter(Use == "YES") %>%
  select(FWB_id) %>%
  distinct() %>%
  count()

actions <- julien_site_sorting %>% 
  select(Action_description) %>%
  distinct()

#sites_remove <- julien_site_sorting %>%
#  filter(Use == "NO")

DataFish_2 <- DataFish %>%
  left_join(julien_site_sorting, by = "FWB_id") %>%
  mutate(Use = coalesce(Use, "YES")) %>% ##if Use is NA (i.e., sites not in sorting table as all species-level) then keep site
  filter(Use != "NO") %>% ##remove sites marked for removal
  filter (case_when(
    Action_description == "Remove individuals with NA" ~ !is.na(fish_species) & !is.na(fish_genus) &!is.na(fish_family),
    Action_description %in% c(
      "Luciobarbus steindachneri; Remove individual with NA",
      "Luciobarbus microcephalus; Remove individual with NA",
      "Luciobarbus sclateri; Remove individual with NA"
    ) ~
      !is.na(fish_species) |
      fish_genus == "Luciobarbus",
    
    # All other rows are retained
    TRUE ~ TRUE
  )) %>%
  
  # !(coalesce(Action_description == "Remove individuals with NA", FALSE) & (is.na(fish_species) | is.na(fish_genus) |is.na(fish_family)))) %>% #remove an individual if any taxonomic field is NA 
  mutate(
    fish_scientific_name = case_when(
      
      # Keep genus-level identification and add "sp."
      Action_description == "Keep genus and add sp." & fish_name_level == "genus" ~ str_squish(paste(fish_genus, "sp.")),
      
      # Keep family-level identification and add "sp."
      Action_description == "Keep family and add sp." & fish_name_level == "family" ~ str_squish(paste(fish_family, "sp.")),
      
      ##Add species level identification for specific cases 
      Action_description == "Luciobarbus steindachneri" & fish_genus == "Luciobarbus" ~ "Luciobarbus steindachneri",
      
      Action_description == "Luciobarbus steindachneri; Remove individual with NA" & fish_genus == "Luciobarbus" ~ "Luciobarbus steindachneri",
      Action_description == "Luciobarbus microcephalus; Remove individual with NA"& fish_genus == "Luciobarbus" ~ "Luciobarbus microcephalus",
      Action_description == "Luciobarbus sclateri; Remove individual with NA"& fish_genus == "Luciobarbus" ~ "Luciobarbus sclateri",
      # Otherwise retain the existing name
      TRUE ~ fish_scientific_name ##for all species level(i.e., fish_name_level == species) it should already be there and all good, nothing to do 
    ))



##add species-site identifier column
DataFish_2$sp_site<-paste(DataFish_2$fish_scientific_name,DataFish_2$collection_site_id,sep="_")
##Assign 1 - number of fish
DataFish_2$num<-1



##Check for and remove temporal replication within site -- want to ensure only 1 year of data
DataFish_2 <- DataFish_2 %>%
  # mutate(as.Date(collection_date)) %>%
  mutate(is_range = str_detect(collection_date, "^\\d{4}-\\d{4}$")) %>% ##indicate if sample collection date is a range of years
  mutate(year = str_extract(collection_date, "^\\d{4}") %>% as.numeric()) ##create column of year

##calculate the number of sampling years
site_num_sampling_years <- DataFish_2 %>%
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
DataFish_2 <- left_join(DataFish_2, site_num_sampling_years, by = c("FWB_id", "collection_site_id"))


##calculate the number of individuals sampled per species within each site-year 
datafish_summary_sp <- DataFish_2 %>%
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
datafish_annual <- DataFish_2 %>%
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
write.csv(summary_all_cutoffs, "data/cuttoff_site_lists.csv")

##Step 3: Create fish df for selected cutoffs/minimum species, and for sites with multiple levels of taxonomic resolution adjust based on sorting by Julian ------------
##Read in cutoff site list, and then select the cutoff/criteria 
#cutoff_site_list_all <- read.csv("data/cuttoff_site_lists.csv")  


cutoff_site_list <- summary_all_cutoffs %>%
  mutate(site_year_code = paste(FWB_id, year, sep = "_")) %>%
  filter(min_sp_num == 2) %>% ##select minimum species number
  filter(cutoff == "75pct") ##select cutoff 
##436 sites 

##select site/years from chosen cutoff 
DataFish_final <- datafish_annual %>%
  filter(site_year_code %in% cutoff_site_list$site_year_code)





##Correcting to final sites:



##Calculate the proporton of total individuals in a food web that the species with less than 3 replicates make up (i.e., calculating rare species)
##Okay so -- current cutoff, 2 species, 75% -- for those 25% w/ less than 3 replicates, what proportion of total number of individuals in that food web are they ? 

cutoff_site_list_25 <- cutoff_site_list %>%
  filter(min_num_ind_per_sample < 3)

datafish_25 <- DataFish_final %>%
  mutate(site_year_code = paste(FWB_id, year, sep = "_")) %>%
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





##checking into Luciobarbus steindachneri
ls1 <- DataFish %>%
  filter(fish_scientific_name == "Luciobarbus sclateri") %>% 
  select(FWB_id, fish_scientific_name) %>%
  distinct()


ls2 <- DataFish_2 %>%
  filter(fish_scientific_name == "Luciobarbus sclateri") %>% 
  select(FWB_id, fish_scientific_name) %>%
  distinct()

ls3 <- DataFish_final %>%
  filter(fish_scientific_name == "Luciobarbus sclateri") %>% 
  select(FWB_id, fish_scientific_name) %>%
  distinct()

#Luciobarbus sclateri

##test - went from 437 to 436 which site is missing?
#sites_437 <- read.csv("data/SiteVar_2sp_75cutoff.csv")  %>%
#  select(FWB_id) %>%
#  distinct()
#sites_436 <- cutoff_site_list

#test <- anti_join(sites_437, sites_436, by = "FWB_id")
##FWB_0757 -- why did this site disappear? because only 70% of its species (less than 75% have more than 3 individuals)