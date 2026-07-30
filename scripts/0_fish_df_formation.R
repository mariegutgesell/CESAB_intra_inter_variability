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

##add species-site identifier column
DataFish$sp_site<-paste(DataFish$fish_scientific_name,DataFish$collection_site_id,sep="_")
##Assign 1 - number of fish
DataFish$num<-1



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
write.csv(summary_all_cutoffs, "data/cuttoff_site_lists.csv")

##Step 3: Create fish df for selected cutoffs/minimum species, and for sites with multiple levels of taxonomic resolution adjust based on sorting by Julian ------------
##Read in cutoff site list, and then select the cutoff/criteria 
#cutoff_site_list_all <- read.csv("data/cuttoff_site_lists.csv")  


cutoff_site_list <- summary_all_cutoffs %>%
  mutate(site_year_code = paste(FWB_id, year, sep = "_")) %>%
  filter(min_sp_num == 2) %>% ##select minimum species number
  filter(cutoff == "75pct") ##select cutoff 
##453 sites 

##select site/years from chosen cutoff 
DataFish_2 <- datafish_annual %>%
  filter(site_year_code %in% cutoff_site_list$site_year_code)

##Correcting to final sites:
##1) read in Julians sorting and clean up accordingly
julien_site_sorting <- read_excel("data/Sites_Julien_27.07.2026.xlsx") %>%
  select(FWB_id, Use, Reason, Action, Action_description) %>%
  distinct(FWB_id, .keep_all = TRUE)
head(julien_site_sorting)


actions <- julien_site_sorting %>% 
  select(Action_description) %>%
  distinct()

#sites_remove <- julien_site_sorting %>%
#  filter(Use == "NO")

DataFish_final <- DataFish_2 %>%
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
  









