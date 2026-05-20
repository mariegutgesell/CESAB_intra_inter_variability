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
library(betareg)
library(broom)

##Load intraspecific variability at site level data - 3 levels of sampling replication 
SiteVar_min1 <- read.csv("data/analysis_csvs/SiteVar_min1.csv")
SiteVar_min3 <- read.csv("data/analysis_csvs/SiteVar_min3.csv")
SiteVar_min5 <- read.csv("data/analysis_csvs/SiteVar_min5.csv")

##checking duplicate site names/lat longs -- don't worry about for now, come back to once using final data
##sometimes have same lat/long with 2 names, and sometimes have same site name with different lat/longs 
latlong1 <- SiteVar_min1 %>%
  select(collection_decimal_latitude, collection_decimal_longitude) %>%
  unique() 

##find which site ids have duplicate lat/longs
dup_latlong <- SiteVar_min1 %>%
  select(FWB_id, collection_site_id, collection_decimal_latitude, collection_decimal_longitude) %>%
  group_by(collection_decimal_latitude, collection_decimal_longitude) %>%
  filter(n() > 1) %>%
  arrange(collection_decimal_latitude, collection_decimal_longitude)
##39 cases with duplicate lat/longs -- some are sampled wet and dry season (may need to only select 1?), some have duplicate site names (see below), others not clear, do we need to think about this if they are truly independent sites? Or just not worry about it? 
##the issue is, is that the environmental data is joined by the latlong - actually this shouldnt matter, because both have the same environmental data - just need to double check the isotope data itself isn't getting duplicated ... checking back, i dont think it is so I think this is okay
##note: these issues disappear when we use minimum 3 or 5 samples

sites1 <- SiteVar_min5 %>%
  select(collection_site_id, waterbody_type) %>%
  group_by(collection_site_id, waterbody_type) %>%
  count()
##2 sites with the same name, one says lentic the other lotic - but have the same lat/long and environmental data .. 
##One-Mile Waterhole
##Twelve Mile Lagoon 
##again, these are not issues when using the min 3/5 samples


##Environmental Variables
##productivity = TP
##type = ecosystem type (lake vs. stream)
##human ftp = hft
##size = size_z_scored
##flow variability = hydro_dis_z_scored
##climate = temp , climate zone 
##Data transformations 
SiteVar_min1 <- SiteVar_min1 %>%
  mutate(site_nbspe_log = log(site_nbspe),
         site_nbspe_log_z_scored = as.numeric(scale(site_nbspe_log,center = TRUE, scale = TRUE)),
         TP_z_scored = as.numeric(scale(TP, center = TRUE, scale = TRUE)),
         hft_z_scored = as.numeric(scale(hft, center = TRUE, scale = TRUE))
         )%>%
  mutate(data_group = "min_1")

SiteVar_min3 <- SiteVar_min3 %>%
  mutate(site_nbspe_log = log(site_nbspe),
         site_nbspe_log_z_scored = as.numeric(scale(site_nbspe_log,center = TRUE, scale = TRUE)),
         TP_z_scored = as.numeric(scale(TP, center = TRUE, scale = TRUE)),
         hft_z_scored = as.numeric(scale(hft, center = TRUE, scale = TRUE))
         )%>%
  mutate(data_group = "min_3")
SiteVar_min5 <- SiteVar_min5 %>% 
  mutate(site_nbspe_log = log(site_nbspe),
         site_nbspe_log_z_scored = as.numeric(scale(site_nbspe_log,center = TRUE, scale = TRUE)),
         TP_z_scored = as.numeric(scale(TP, center = TRUE, scale = TRUE)),
         hft_z_scored = as.numeric(scale(hft, center = TRUE, scale = TRUE))
         ) %>%
  mutate(data_group = "min_5")

##bind all together 
SiteVar_all <- rbind(SiteVar_min1, SiteVar_min3, SiteVar_min5)

##potential random effects: species richness, C:N, tissue type 
##Make correlation plot of numerical explanatory variables
exp_var_min1 <- SiteVar_min1%>%
  filter(Size != "-Inf") %>%
  ungroup() %>%
  select(TP_z_scored, size_z_scored, hydro_dis_z_scored, hft_z_scored, site_nbspe_log_z_scored)# %>%
#  unique() #


exp_var_min3 <- SiteVar_min3%>%
  filter(Size != "-Inf") %>%
  ungroup() %>%
  select(TP_z_scored, size_z_scored, hydro_dis_z_scored, hft_z_scored, site_nbspe_log_z_scored)

exp_var_min5 <- SiteVar_min5%>%
  filter(Size != "-Inf") %>%
  ungroup() %>%
  select(TP_z_scored, size_z_scored, hydro_dis_z_scored, hft_z_scored, site_nbspe_log_z_scored)

duplicated_exp_variables <- SiteVar_min5 %>%
  select(collection_site_id,collection_decimal_longitude, collection_decimal_latitude,
         waterbody_type, Climate_zone_e2, TP,
         size_z_scored, hydro_dis_z_scored, hft, site_nbspe_log) %>%
  filter(
    duplicated(select(., waterbody_type, Climate_zone_e2, TP,
                      size_z_scored, hydro_dis_z_scored, hft, site_nbspe_log)) |
      duplicated(select(., waterbody_type, Climate_zone_e2, TP,
                        size_z_scored, hydro_dis_z_scored, hft, site_nbspe_log), fromLast = TRUE)
  )

##so there are a fair amount of duplicated variables, likely because sites are close togehter  - need to check this carefully when have final data
##also, some points are in the ocean, again check if this issue remains with final data 
sites_coord <- st_as_sf(duplicated_exp_variables, coords = c("collection_decimal_longitude", "collection_decimal_latitude"), crs = 4326)
site_map <- mapview(sites_coord, map.types = "Esri.NatGeoWorldMap", legend = FALSE, cex = 2)
site_map

corr_mat_min1 <- cor(exp_var_min1, use = "pairwise.complete.obs")
corr_mat_min3 <- cor(exp_var_min3, use = "pairwise.complete.obs")
corr_mat_min5 <- cor(exp_var_min5, use = "pairwise.complete.obs")

# visualize
corrplot(corr_mat_min1, method = "color", type = "upper",  addCoef.col = "black",
         tl.col = "black", tl.srt = 45)

corrplot(corr_mat_min3, method = "color", type = "upper",  addCoef.col = "black",
         tl.col = "black", tl.srt = 45)
corrplot(corr_mat_min5, method = "color", type = "upper",  addCoef.col = "black",
         tl.col = "black", tl.srt = 45)


##Look at distributions of proportion intraspecific variation
ggplot(SiteVar_all, aes(x = propintraspecific_C)) +
  geom_histogram() +
  facet_wrap(~data_group)

ggplot(SiteVar_all, aes(x = propintraspecific_N)) +
  geom_histogram() +
  facet_wrap(~data_group)

ggplot(SiteVar_all, aes(x = TP_z_scored)) +
  geom_histogram() +
  facet_wrap(~data_group)


#SiteVar_all$Climate_zone_e2 <- ordered(SiteVar_all$Climate_zone_e2, 
#                                       levels = c(
#                                         "Cold and wet/mesic",
#                                         "Cool and moist",
#                                         "Cool temperate and dry/xeric",
#                                         "Warm temperate",
#                                         "Hot and moist",
#                                         "Hot and dry"
#                                       ))
ggplot(SiteVar_all, aes(x = data_group, y = TP_z_scored, fill = Climate_zone_e2)) +
  geom_boxplot() 

##Write function that conducts beta regressions or ANOVAs over all 3 datasets 
beta_reg_univariate <- function(dat,
                                res = "propintraspecific_C",
                                exp = "site_nbspe_log") {
  
  # keep only needed cols and remove NAs
  df <- dat %>%
    select(data_group, all_of(res), all_of(exp)) %>%
    filter(
      !is.na(.data[[res]]),
      !is.na(.data[[exp]])
    )
  
  # fit one beta regression per data group
  models <- df %>%
    group_by(data_group) %>%
    group_split() %>%
    map(function(d) {
      formula_i <- as.formula(paste(res, "~", exp))
      mod <- betareg(formula_i, data = d)
      
      tibble(
        data_group = unique(d$data_group),
        model = list(mod),
        n = nrow(d),
        pseudo_r2 = summary(mod)$pseudo.r.squared,
        estimate = coef(summary(mod))$mean[exp, "Estimate"],
        p_value = coef(summary(mod))$mean[exp, "Pr(>|z|)"]
      )
    }) %>%
    bind_rows()
  
  # prediction data for each group
  pred_df <- df %>%
    group_by(data_group) %>%
    summarise(
      x_min = min(.data[[exp]], na.rm = TRUE),
      x_max = max(.data[[exp]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      newdat = list(tibble(!!exp := seq(x_min, x_max, length.out = 100)))
    ) %>%
    left_join(models, by = "data_group") %>%
    mutate(
      pred = list(predict(model, newdata = newdat, type = "response")),
      newdat = list(bind_cols(newdat, tibble(pred = pred)))
    ) %>%
    select(data_group, newdat) %>%
    tidyr::unnest(newdat)
  
  # plot
  p <- ggplot(df, aes(x = .data[[exp]], y = .data[[res]], color = data_group)) +
    geom_point(alpha = 0.5) +
    geom_line(
      data = pred_df,
      aes(x = .data[[exp]], y = pred, color = data_group),
      linewidth = 1.2
    ) +
    theme_classic() +
    labs(
      x = exp,
      y = res,
      color = "Dataset"
    )
  
  list(
    models = models,
    predictions = pred_df,
    plot = p
  )
}

aov_univariate <- function(dat,
                           res = "propintraspecific_C",
                           exp = "Climate_zone_e2",
                           order_levels = NULL) {
  
  df <- dat %>%
    select(data_group, all_of(res), all_of(exp)) %>%
    filter(!is.na(.data[[res]]), !is.na(.data[[exp]]))
  
  # optionally order factor levels
  if (!is.null(order_levels)) {
    df[[exp]] <- factor(df[[exp]], levels = order_levels, ordered = TRUE)
  }
  
  # run ANOVA per dataset
  models <- df %>%
    group_by(data_group) %>%
    group_split() %>%
    map(function(d) {
      formula_i <- as.formula(paste(res, "~", exp))
      mod <- aov(formula_i, data = d)
      
      tibble(
        data_group = unique(d$data_group),
        model = list(mod),
        anova = list(broom::tidy(mod)),
        tukey = list(TukeyHSD(mod)[[exp]])
      )
    }) %>%
    bind_rows()
  
  # combine Tukey results nicely
  tukey_df <- models %>%
    select(data_group, tukey) %>%
    unnest(tukey) %>%
    rownames_to_column("comparison")
  
  # plot (all datasets together)
  p <- ggplot(df, aes(x = .data[[exp]], y = .data[[res]], fill = data_group)) +
    geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.8)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = exp,
      y = res,
      fill = "Dataset"
    )
  p2 <- ggplot(df, aes(x = data_group, y = .data[[res]], fill = .data[[exp]])) +
    geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.8)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      x = "Dataset",
      y = res,
      fill = exp  
    )
  
  list(
    models = models,
    tukey = tukey_df,
    plot = p,
    plot2 = p2
  )
}

##Univariate Beta Regression Results - Proportion Intraspecific Variability Contribution & Environmental Variables (and ) -------------------------
##Species Richness
sp_richness_C <- beta_reg_univariate(SiteVar_all, res = "propintraspecific_C", exp = "site_nbspe_log")
sp_richness_C$models
sp_richness_C$plot

sp_richness_N <- beta_reg_univariate(SiteVar_all, res = "propintraspecific_N", exp = "site_nbspe_log")
sp_richness_N$models
sp_richness_N$plot


##Total Phosphorous
tp_C <- beta_reg_univariate(SiteVar_all, res = "propintraspecific_C", exp = "TP")
tp_C$models
tp_C$plot

tp_N <- beta_reg_univariate(SiteVar_all, res = "propintraspecific_N", exp = "TP")
tp_N$models
tp_N$plot

##Human Footprint
hft_C <- beta_reg_univariate(SiteVar_all, res = "propintraspecific_C", exp = "hft")
hft_C$models
hft_C$plot


hft_N <- beta_reg_univariate(SiteVar_all, res = "propintraspecific_N", exp = "hft")
hft_N$models
hft_N$plot

##Flow Variability
fv_C <- beta_reg_univariate(SiteVar_all, res = "propintraspecific_C", exp = "hydro_dis_z_scored")
fv_C$models
fv_C$plot


fv_N <- beta_reg_univariate(SiteVar_all, res = "propintraspecific_N", exp = "hydro_dis_z_scored")
fv_N$models
fv_N$plot

##Ecosystem Size
es_C <- beta_reg_univariate(SiteVar_all, res = "propintraspecific_C", exp = "size_z_scored")
es_C$models
es_C$plot

es_N <- beta_reg_univariate(SiteVar_all, res = "propintraspecific_N", exp = "size_z_scored")
es_N$models
es_N$plot

##Climate category
climate_levels <- c(
  "Cold and wet/mesic",
  "Cool and moist",
  "Cool temperate and dry/xeric",
  "Warm temperate",
  "Hot and moist",
  "Hot and dry"
)

out_C_climate <- aov_univariate(SiteVar_all, res = "propintraspecific_C", exp = "Climate_zone_e2",order_levels = climate_levels)
out_C_climate$plot
out_C_climate$plot2
out_C_climate$tukey

out_N_climate <- aov_univariate(SiteVar_all, res = "propintraspecific_N", exp = "Climate_zone_e2",order_levels = climate_levels)
out_N_climate$plot
out_N_climate$plot2
out_N_climate$tukey

##waterbody type
out_C_type <- aov_univariate(SiteVar_all, res = "propintraspecific_C", exp = "waterbody_type")
out_C_type$plot
out_C_type$plot2
out_C_type$tukey

out_N_type <- aov_univariate(SiteVar_all, res = "propintraspecific_N", exp = "waterbody_type")
out_N_type$plot
out_N_type$plot2
out_N_type$tukey

##Note: ANOVA assumes normal residuals, may need to transform proportions (logit) or use beta regressions w/ categorical predictor

##Does the effect of species richness differ my climate or habitat type?
sp_rich_climate <- aov_univariate(SiteVar_all, res = "site_nbspe_log", exp = "Climate_zone_e2", order_levels = climate_levels)
sp_rich_climate$plot
sp_rich_climate$plot2
sp_rich_climate$tukey

sp_rich_type <- aov_univariate(SiteVar_all, res = "site_nbspe_log", exp = "waterbody_type")
sp_rich_type$plot
sp_rich_type$plot2
sp_rich_type$tukey


###Multiple Regression Models ------------------------
##Look at series of potential regression models - hypothesis driven
##does the environmental variable affect the response, after accounting for species richness and climate? is there an overall relationship between environmental driver and contribution fo intrasp. var to total var? 
m1 <- betareg(y ~ exp + richness + climate)
#does the effect of richness on the response differ among climate zones, while still including an overall efffect of exp? does species richness matter differently in tropical vs. temperate vs. cold systems?
m2 <- betareg(y ~ exp + richness * climate)

##does the effect of the environmental variable depend on species richness? does environmental variable affect intraspecific variability differently in species rich vs. species poor communities? 
m3 <- betareg(y ~ exp * richness + climate)

m4 <- betareg()

##these next two are not really needed...not really hypothesis driven unless we specifically want to exclude climate? 
#does environmental variable explain response after accountinf for richness? ignoring cliamte differences, do TP and richness explain variation? 
m4 <- betareg(y ~ exp + richness)

#does the relationship between environmental variable and the response depend on richness, without accounting for cliamte? is there an itneraction between TP and richness, if we ignore climatic context? 
m5 <- betareg(y ~ exp + richness)
##do all explanatory and species richness etc variables need to be scaled? i think if including them in one model..yes... 

##Checking Model Assumptions:
##1) ensure response variable is 0 < y < 1 (can't have exactly 0 or 1)

##2) ensure no colinearity between predictors 
##function to check colinearity:
check_beta_collinearity <- function(dat,
                                    res,
                                    exp,
                                    richness = "site_nbspe_log_z_scored",
                                    climate = "Climate_zone_e2") {
  
  cor_by_group <- dat %>%
    group_by(data_group) %>%
    summarise(
      cor_exp_richness = cor(
        .data[[exp]],
        .data[[richness]],
        use = "complete.obs"
      ),
      .groups = "drop"
    )
  
  vif_by_group <- dat %>%
    select(data_group, all_of(c(res, exp, richness, climate))) %>%
    filter(if_all(everything(), ~ !is.na(.x))) %>%
    group_split(data_group) %>%
    purrr::map_df(function(d) {
      
      form <- as.formula(
        paste(res, "~", exp, "*", climate, "+", richness)
      )
      
      mod <- lm(form, data = d)
      
      vif_out <- as.data.frame(car::vif(mod))
      vif_out$predictor <- rownames(vif_out)
      
      vif_out %>%
        tibble::as_tibble() %>%
        rename(
          degrees_freedom = Df,
          adjusted_vif = `GVIF^(1/(2*Df))`
        ) %>%
        mutate(
          data_group = unique(d$data_group),
          vif_flag = case_when(
            adjusted_vif < 2 ~ "low",
            adjusted_vif < 5 ~ "moderate",
            adjusted_vif < 10 ~ "high",
            TRUE ~ "problematic"
          )
        ) %>%
        select(data_group, predictor, degrees_freedom, adjusted_vif, vif_flag)
    })
  
  list(
    correlation = cor_by_group,
    vif = vif_by_group
  )
}

##3) need to have enough sites for each climate zone .. 
site_climate_summary <- SiteVar_all %>%
  group_by(data_group, Climate_zone_e2) %>%
  count()
#for min  = 5, don't have enough sites in each climate category ... so i don't think climate interactions are interpretable 

##Function to run model candidate types,
run_beta_models <- function(dat,
                            res,
                            exp,
                            richness = "site_nbspe_log_z_scored",
                            climate = "Climate_zone_e2") {
  
  df <- dat %>%
    select(data_group, all_of(c(res, exp, richness, climate))) %>%
    filter(if_all(everything(), ~ !is.na(.x)))
  
  safe_betareg <- purrr::possibly(
    function(formula, data) betareg::betareg(as.formula(formula), data = data),
    otherwise = NULL
  )
  
  fit_one_group <- function(d) {
    
    formulas <- list(
      m1_main = paste(res, "~", exp, "+", richness, "+", climate),
      m2_richness_x_climate = paste(res, "~", exp, "+", richness, "*", climate),
      m3_exp_x_richness = paste(res, "~", exp, "*", richness, "+", climate) ,
      m4_exp_x_climate = paste(res, "~", exp, "*", climate, "+", richness) ,
      m5_exp_x_richness_noclimate = paste(res, "~", exp, "*", richness),
      m6_exp_richness_noclimate = paste(res, "~", exp, "+", richness)
    )
    
    models <- purrr::map(formulas, ~ safe_betareg(.x, d))
    
    tibble(
      data_group = unique(d$data_group),
      response = res,
      explanatory = exp,
      data = list(d),
      formulas = list(formulas),
      models = list(models)
    )
  }
  
  df %>%
    group_by(data_group) %>%
    group_split() %>%
    map_dfr(fit_one_group)
}



select_beta_models_aic <- function(model_results, delta_cutoff = 2) {
  
  model_results %>%
    mutate(
      aic_table = map(models, function(mods) {
        
        tibble(
          model_name = names(mods),
          model = mods
        ) %>%
          filter(!map_lgl(model, is.null)) %>%
          mutate(AIC = map_dbl(model, AIC)) %>%
          arrange(AIC) %>%
          mutate(
            delta_AIC = AIC - min(AIC),
            weight = exp(-0.5 * delta_AIC) / sum(exp(-0.5 * delta_AIC)),
            supported = delta_AIC < delta_cutoff
          )
      }),
      
      top_model_name = map_chr(aic_table, ~ {
        if (nrow(.x) == 0) NA_character_ else .x$model_name[1]
      }),
      
      top_model = map(aic_table, ~ {
        if (nrow(.x) == 0) NULL else .x$model[[1]]
      }),
      
      supported_model_names = map(aic_table, ~ {
        if (nrow(.x) == 0) character(0) else .x$model_name[.x$supported]
      }),
      
      n_supported_models = map_int(supported_model_names, length),
      
      model_selection_result = case_when(
        n_supported_models == 0 ~ "no models fitted",
        n_supported_models == 1 ~ "clear top model",
        n_supported_models > 1 ~ "model uncertainty"
      )
    )
}

plot_beta_model_predictions <- function(selected,
                                        exp,
                                        model_name,
                                        richness = "site_nbspe_log_z_scored",
                                        climate = "Climate_zone_e2",
                                        n_points = 100) {
  
  make_predictions <- function(d, mods, res) {
    
    mod <- mods[[model_name]]
    
    if (is.null(mod)) {
      return(tibble())
    }
    
    newdat <- tibble(
      !!exp := seq(
        min(d[[exp]], na.rm = TRUE),
        max(d[[exp]], na.rm = TRUE),
        length.out = n_points
      ),
      !!richness := median(d[[richness]], na.rm = TRUE),
      !!climate := names(sort(table(d[[climate]]), decreasing = TRUE))[1]
    )
    
    newdat$pred <- predict(mod, newdata = newdat, type = "response")
    newdat$data_group <- unique(d$data_group)
    newdat$response <- res
    newdat$explanatory <- exp
    newdat$model_name <- model_name
    
    newdat
  }
  
  pred_df <- pmap_dfr(
    list(
      d = selected$data,
      mods = selected$models,
      res = selected$response
    ),
    make_predictions
  )
  
  plot_df <- selected %>%
    select(data) %>%
    unnest(data)
  
  res <- unique(selected$response)
  
  p <- ggplot(plot_df, aes(x = .data[[exp]], y = .data[[res]], color = data_group)) +
    geom_point(alpha = 0.45) +
    geom_line(
      data = pred_df,
      aes(x = .data[[exp]], y = pred, color = data_group),
      linewidth = 1.2
    ) +
    theme_classic() +
    labs(
      x = exp,
      y = res,
      color = "Dataset",
      title = paste("Model predictions:", model_name, "|", res, "~", exp)
    )
  
  list(
    predictions = pred_df,
    plot_data = plot_df,
    plot = p
  )
}

##function to be able to show/plot interactions - how relationships change between climate zones 
plot_beta_climate_interaction <- function(selected,
                                          exp,
                                          model_name,
                                          richness = "site_nbspe_log_z_scored",
                                          climate = "Climate_zone_e2",
                                          n_points = 100) {
  
  make_predictions <- function(d, mods, res) {
    
    mod <- mods[[model_name]]
    
    if (is.null(mod)) return(tibble())
    
    newdat <- tidyr::expand_grid(
      !!exp := seq(
        min(d[[exp]], na.rm = TRUE),
        max(d[[exp]], na.rm = TRUE),
        length.out = n_points
      ),
      !!climate := unique(d[[climate]])
    ) %>%
      mutate(
        !!richness := median(d[[richness]], na.rm = TRUE),
        data_group = unique(d$data_group),
        response = res,
        explanatory = exp,
        model_name = model_name
      )
    
    newdat$pred <- predict(mod, newdata = newdat, type = "response")
    
    newdat
  }
  
  pred_df <- purrr::pmap_dfr(
    list(
      d = selected$data,
      mods = selected$models,
      res = selected$response
    ),
    make_predictions
  )
  
  plot_df <- selected %>%
    select(data) %>%
    unnest(data)
  
  res <- unique(selected$response)
  
  p <- ggplot(plot_df, aes(x = .data[[exp]], y = .data[[res]])) +
    geom_point(aes(color = .data[[climate]]), alpha = 0.45) +
    geom_line(
      data = pred_df,
      aes(x = .data[[exp]], y = pred, color = .data[[climate]]),
      linewidth = 1.1
    ) +
    facet_wrap(~ data_group) +
    theme_classic() +
    labs(
      x = exp,
      y = res,
      color = "Climate zone",
      title = paste("Climate interaction predictions:", model_name)
    )
  
  list(
    predictions = pred_df,
    plot_data = plot_df,
    plot = p
  )
}

##for now, lets just focus on min1 and min3, not enough samples to look at effects of climate in min 5 
SiteVar_all <- SiteVar_all %>%
  filter(data_group != "min_5")

##TP vs. dC -----------
#check colinearity 
diag <- check_beta_collinearity(
  SiteVar_all,
  res = "propintraspecific_C",
  exp = "TP_z_scored"
)

diag$correlation
diag$vif

mods <- run_beta_models(
  SiteVar_all,
  res = "propintraspecific_C",
  exp = "TP_z_scored"
)

# inspect all models first
mods$models[[1]]
summary(mods$models[[1]]$m1_main)
summary(mods$models[[1]]$m2_richness_x_climate)
summary(mods$models[[1]]$m3_exp_x_richness)
summary(mods$models[[1]]$m4_exp_x_climate)
summary(mods$models[[1]]$m5_exp_x_richness_noclimate)
summary(mods$models[[1]]$m6_exp_richness_noclimate)

##min3
mods$models[[1]]
summary(mods$models[[2]]$m1_main)
summary(mods$models[[2]]$m2_richness_x_climate)
summary(mods$models[[2]]$m3_exp_x_richness)
summary(mods$models[[2]]$m4_exp_x_climate)
summary(mods$models[[2]]$m5_exp_x_richness_noclimate)
summary(mods$models[[2]]$m6_exp_richness_noclimate)

##min5
mods$models[[3]]
summary(mods$models[[3]]$m1_main)
#summary(mods$models[[3]]$m2_richness_x_climate)
summary(mods$models[[3]]$m3_exp_x_richness)
#summary(mods$models[[1]]$m4_exp_x_climate)

# then select by AIC
selected <- select_beta_models_aic(mods)

selected %>%
  select(
    data_group,
    response,
    explanatory,
    top_model_name,
    supported_model_names,
    n_supported_models,
    model_selection_result
  )

selected %>%
  select(data_group, response, explanatory, aic_table) %>%
  unnest(aic_table)

p1 <- plot_beta_model_predictions(selected, exp = "TP_z_scored", model_name = "m1_main")
p2 <- plot_beta_model_predictions(selected, exp = "TP_z_scored", model_name = "m2_richness_x_climate")
p3 <- plot_beta_model_predictions(selected, exp = "TP_z_scored", model_name = "m3_exp_x_richness")
p4 <- plot_beta_model_predictions(selected, exp = "TP_z_scored", model_name = "m4_exp_x_climate")
p5 <- plot_beta_model_predictions(selected, exp = "TP_z_scored", model_name = "m5_exp_x_richness_noclimate")
p6 <- plot_beta_model_predictions(selected, exp = "TP_z_scored", model_name = "m6_exp_richness_noclimate")

p1$plot
p2$plot
p3$plot
#p4$plot
p5$plot
p6$plot

p_clim_exp <- plot_beta_climate_interaction(
  selected,
  exp = "TP_z_scored",
  model_name = "m4_exp_x_climate"
)

p_clim_exp$plot

p_rich_clim <- plot_beta_climate_interaction(
  selected,
  exp = "site_nbspe_log_z_scored",
  model_name = "m2_richness_x_climate",
  richness = "TP_z_scored" # held at median here
)
p_rich_clim$plot



##TP vs. dN -----------
#check colinearity 
#check colinearity 
diag <- check_beta_collinearity(
  SiteVar_all,
  res = "propintraspecific_N",
  exp = "TP_z_scored"
)

diag$correlation

diag$vif


mods <- run_beta_models(
  SiteVar_all,
  res = "propintraspecific_N",
  exp = "TP_z_scored"
)

# inspect all models first
mods$models[[1]]
summary(mods$models[[1]]$m1_main)
summary(mods$models[[1]]$m2_richness_x_climate)
summary(mods$models[[1]]$m3_exp_x_richness)
summary(mods$models[[1]]$m4_exp_x_climate)
summary(mods$models[[1]]$m5_exp_x_richness_noclimate)
summary(mods$models[[1]]$m6_exp_richness_noclimate)

##min3
mods$models[[1]]
summary(mods$models[[2]]$m1_main)
summary(mods$models[[2]]$m2_richness_x_climate)
summary(mods$models[[2]]$m3_exp_x_richness)
summary(mods$models[[2]]$m4_exp_x_climate)
summary(mods$models[[2]]$m5_exp_x_richness_noclimate)
summary(mods$models[[2]]$m6_exp_richness_noclimate)

##min5
mods$models[[3]]
summary(mods$models[[3]]$m1_main)
#summary(mods$models[[3]]$m2_richness_x_climate)
summary(mods$models[[3]]$m3_exp_x_richness)
#summary(mods$models[[1]]$m4_exp_x_climate)

# then select by AIC
selected <- select_beta_models_aic(mods)

selected %>%
  select(
    data_group,
    response,
    explanatory,
    top_model_name,
    supported_model_names,
    n_supported_models,
    model_selection_result
  )

selected %>%
  select(data_group, response, explanatory, aic_table) %>%
  unnest(aic_table)

p1 <- plot_beta_model_predictions(selected, exp = "TP_z_scored", model_name = "m1_main")
p2 <- plot_beta_model_predictions(selected, exp = "TP_z_scored", model_name = "m2_richness_x_climate")
p3 <- plot_beta_model_predictions(selected, exp = "TP_z_scored", model_name = "m3_exp_x_richness")
p4 <- plot_beta_model_predictions(selected, exp = "TP_z_scored", model_name = "m4_exp_x_climate")
p5 <- plot_beta_model_predictions(selected, exp = "TP_z_scored", model_name = "m5_exp_x_richness_noclimate")
p6 <- plot_beta_model_predictions(selected, exp = "TP_z_scored", model_name = "m6_exp_richness_noclimate")


p1$plot
p2$plot
p3$plot
#p4$plot
p5$plot
p6$plot

p_clim_exp <- plot_beta_climate_interaction(
  selected,
  exp = "TP_z_scored",
  model_name = "m4_exp_x_climate"
)

p_clim_exp$plot

p_rich_clim <- plot_beta_climate_interaction(
  selected,
  exp = "site_nbspe_log_z_scored",
  model_name = "m2_richness_x_climate",
  richness = "TP_z_scored" # held at median here
)
p_rich_clim$plot





##hft vs. dC -----------
#check colinearity 
diag <- check_beta_collinearity(
  SiteVar_all,
  res = "propintraspecific_C",
  exp = "hft_z_scored"
)

diag$correlation
diag$vif


mods <- run_beta_models(
  SiteVar_all,
  res = "propintraspecific_C",
  exp = "hft_z_scored"
)

# inspect all models first
mods$models[[1]]
summary(mods$models[[1]]$m1_main)
summary(mods$models[[1]]$m2_richness_x_climate)
summary(mods$models[[1]]$m3_exp_x_richness)
summary(mods$models[[1]]$m4_exp_x_climate)
summary(mods$models[[1]]$m5_exp_x_richness_noclimate)
summary(mods$models[[1]]$m6_exp_richness_noclimate)


##min3
mods$models[[1]]
summary(mods$models[[2]]$m1_main)
summary(mods$models[[2]]$m2_richness_x_climate)
summary(mods$models[[2]]$m3_exp_x_richness)
summary(mods$models[[2]]$m4_exp_x_climate)
summary(mods$models[[2]]$m5_exp_x_richness_noclimate)
summary(mods$models[[2]]$m6_exp_richness_noclimate)

##min3
mods$models[[3]]
summary(mods$models[[3]]$m1_main)
summary(mods$models[[3]]$m2_richness_x_climate)
summary(mods$models[[3]]$m3_exp_x_richness)
summary(mods$models[[3]]$m4_exp_x_climate)

# then select by AIC
selected <- select_beta_models_aic(mods)

selected %>%
  select(
    data_group,
    response,
    explanatory,
    top_model_name,
    supported_model_names,
    n_supported_models,
    model_selection_result
  )

selected %>%
  select(data_group, response, explanatory, aic_table) %>%
  unnest(aic_table)

p1 <- plot_beta_model_predictions(selected, exp = "hft_z_scored", model_name = "m1_main")
p2 <- plot_beta_model_predictions(selected, exp = "hft_z_scored", model_name = "m2_richness_x_climate")
p3 <- plot_beta_model_predictions(selected, exp = "hft_z_scored", model_name = "m3_exp_x_richness")
p4 <- plot_beta_model_predictions(selected, exp = "hft_z_scored", model_name = "m4_exp_x_climate")
p5 <- plot_beta_model_predictions(selected, exp = "hft_z_scored", model_name = "m5_exp_x_richness_noclimate")
p6 <- plot_beta_model_predictions(selected, exp = "hft_z_scored", model_name = "m6_exp_richness_noclimate")


p1$plot
p2$plot
p3$plot
p4$plot
p5$plot
p6$plot

p_clim_exp <- plot_beta_climate_interaction(
  selected,
  exp = "hft_z_scored",
  model_name = "m4_exp_x_climate"
)

p_clim_exp$plot

p_rich_clim <- plot_beta_climate_interaction(
  selected,
  exp = "site_nbspe_log_z_scored",
  model_name = "m2_richness_x_climate",
  richness = "hft_z_scored" # held at median here
)
p_rich_clim$plot


##hft vs. dN -----------
#check colinearity 
#check colinearity 
diag <- check_beta_collinearity(
  SiteVar_all,
  res = "propintraspecific_N",
  exp = "hft_z_scored"
)

diag$correlation
diag$vif


mods <- run_beta_models(
  SiteVar_all,
  res = "propintraspecific_N",
  exp = "hft_z_scored"
)

# inspect all models first
mods$models[[1]]
summary(mods$models[[1]]$m1_main)
summary(mods$models[[1]]$m2_richness_x_climate)
summary(mods$models[[1]]$m3_exp_x_richness)
summary(mods$models[[1]]$m4_exp_x_climate)
summary(mods$models[[1]]$m5_exp_x_richness_noclimate)
summary(mods$models[[1]]$m6_exp_richness_noclimate)



##min3
mods$models[[2]]
summary(mods$models[[2]]$m1_main)
summary(mods$models[[2]]$m2_richness_x_climate)
summary(mods$models[[2]]$m3_exp_x_richness)
summary(mods$models[[2]]$m4_exp_x_climate)
summary(mods$models[[2]]$m5_exp_x_richness_noclimate)
summary(mods$models[[2]]$m6_exp_richness_noclimate)

##min3
mods$models[[3]]
summary(mods$models[[3]]$m1_main)
summary(mods$models[[3]]$m2_richness_x_climate)
summary(mods$models[[3]]$m3_exp_x_richness)
summary(mods$models[[3]]$m4_exp_x_climate)

# then select by AIC
selected <- select_beta_models_aic(mods)

selected %>%
  select(
    data_group,
    response,
    explanatory,
    top_model_name,
    supported_model_names,
    n_supported_models,
    model_selection_result
  )

selected %>%
  select(data_group, response, explanatory, aic_table) %>%
  unnest(aic_table)

p1 <- plot_beta_model_predictions(selected, exp = "hft_z_scored", model_name = "m1_main")
p2 <- plot_beta_model_predictions(selected, exp = "hft_z_scored", model_name = "m2_richness_x_climate")
p3 <- plot_beta_model_predictions(selected, exp = "hft_z_scored", model_name = "m3_exp_x_richness")
p4 <- plot_beta_model_predictions(selected, exp = "hft_z_scored", model_name = "m4_exp_x_climate")
p5 <- plot_beta_model_predictions(selected, exp = "hft_z_scored", model_name = "m5_exp_x_richness_noclimate")
p6 <- plot_beta_model_predictions(selected, exp = "hft_z_scored", model_name = "m6_exp_richness_noclimate")

p1$plot
p2$plot
p3$plot
p4$plot
p5$plot
p6$plot

p_clim_exp <- plot_beta_climate_interaction(
  selected,
  exp = "hft_z_scored",
  model_name = "m4_exp_x_climate"
)

p_clim_exp$plot

p_rich_clim <- plot_beta_climate_interaction(
  selected,
  exp = "site_nbspe_log_z_scored",
  model_name = "m2_richness_x_climate",
  richness = "hft_z_scored" # held at median here
)
p_rich_clim$plot



##size vs. dC -----------
#check colinearity 
diag <- check_beta_collinearity(
  SiteVar_all,
  res = "propintraspecific_C",
  exp = "size_z_scored"
)

diag$correlation
diag$vif

mods <- run_beta_models(
  SiteVar_all,
  res = "propintraspecific_C",
  exp = "size_z_scored"
)

# inspect all models first
mods$models[[1]]
summary(mods$models[[1]]$m1_main)
summary(mods$models[[1]]$m2_richness_x_climate)
summary(mods$models[[1]]$m3_exp_x_richness)
summary(mods$models[[1]]$m4_exp_x_climate)
summary(mods$models[[1]]$m5_exp_x_richness_noclimate)
summary(mods$models[[1]]$m6_exp_richness_noclimate)


##min3
mods$models[[2]]
summary(mods$models[[2]]$m1_main)
summary(mods$models[[2]]$m2_richness_x_climate)
summary(mods$models[[2]]$m3_exp_x_richness)
summary(mods$models[[2]]$m4_exp_x_climate)
summary(mods$models[[2]]$m5_exp_x_richness_noclimate)
summary(mods$models[[2]]$m6_exp_richness_noclimate)

##min3
mods$models[[3]]
summary(mods$models[[3]]$m1_main)
summary(mods$models[[3]]$m2_richness_x_climate)
summary(mods$models[[3]]$m3_exp_x_richness)
summary(mods$models[[3]]$m4_exp_x_climate)
summary(mods$models[[3]]$m5_exp_x_richness_noclimate)
summary(mods$models[[3]]$m6_exp_richness_noclimate)

# then select by AIC
selected <- select_beta_models_aic(mods)

selected %>%
  select(
    data_group,
    response,
    explanatory,
    top_model_name,
    supported_model_names,
    n_supported_models,
    model_selection_result
  )

selected %>%
  select(data_group, response, explanatory, aic_table) %>%
  unnest(aic_table)

p1 <- plot_beta_model_predictions(selected, exp = "size_z_scored", model_name = "m1_main")
p2 <- plot_beta_model_predictions(selected, exp = "size_z_scored", model_name = "m2_richness_x_climate")
p3 <- plot_beta_model_predictions(selected, exp = "size_z_scored", model_name = "m3_exp_x_richness")
p4 <- plot_beta_model_predictions(selected, exp = "size_z_scored", model_name = "m4_exp_x_climate")
p5 <- plot_beta_model_predictions(selected, exp = "size_z_scored", model_name = "m5_exp_x_richness_noclimate")
p6 <- plot_beta_model_predictions(selected, exp = "size_z_scored", model_name = "m6_exp_richness_noclimate")


p1$plot
p2$plot
p3$plot
p4$plot
p5$plot
p6$plot

p_clim_exp <- plot_beta_climate_interaction(
  selected,
  exp = "size_z_scored",
  model_name = "m4_exp_x_climate"
)


p_clim_exp$plot

p_rich_clim <- plot_beta_climate_interaction(
  selected,
  exp = "site_nbspe_log_z_scored",
  model_name = "m2_richness_x_climate",
  richness = "size_z_scored" # held at median here
)
p_rich_clim$plot

##size vs. dN -----------
#check colinearity 
diag <- check_beta_collinearity(
  SiteVar_all,
  res = "propintraspecific_N",
  exp = "size_z_scored"
)

diag$correlation
diag$vif

mods <- run_beta_models(
  SiteVar_all,
  res = "propintraspecific_N",
  exp = "size_z_scored"
)

# inspect all models first
mods$models[[1]]
summary(mods$models[[1]]$m1_main)
summary(mods$models[[1]]$m2_richness_x_climate)
summary(mods$models[[1]]$m3_exp_x_richness)
summary(mods$models[[1]]$m4_exp_x_climate)
summary(mods$models[[1]]$m5_exp_x_richness_noclimate)
summary(mods$models[[1]]$m6_exp_richness_noclimate)


##min3
mods$models[[2]]
summary(mods$models[[2]]$m1_main)
summary(mods$models[[2]]$m2_richness_x_climate)
summary(mods$models[[2]]$m3_exp_x_richness)
summary(mods$models[[2]]$m4_exp_x_climate)
summary(mods$models[[2]]$m5_exp_x_richness_noclimate)
summary(mods$models[[2]]$m6_exp_richness_noclimate)

##min3
mods$models[[3]]
summary(mods$models[[3]]$m1_main)
summary(mods$models[[3]]$m2_richness_x_climate)
summary(mods$models[[3]]$m3_exp_x_richness)
summary(mods$models[[3]]$m4_exp_x_climate)
summary(mods$models[[3]]$m5_exp_x_richness_noclimate)
summary(mods$models[[3]]$m6_exp_richness_noclimate)

# then select by AIC
selected <- select_beta_models_aic(mods)

selected %>%
  select(
    data_group,
    response,
    explanatory,
    top_model_name,
    supported_model_names,
    n_supported_models,
    model_selection_result
  )

selected %>%
  select(data_group, response, explanatory, aic_table) %>%
  unnest(aic_table)

p1 <- plot_beta_model_predictions(selected, exp = "size_z_scored", model_name = "m1_main")
p2 <- plot_beta_model_predictions(selected, exp = "size_z_scored", model_name = "m2_richness_x_climate")
p3 <- plot_beta_model_predictions(selected, exp = "size_z_scored", model_name = "m3_exp_x_richness")
p4 <- plot_beta_model_predictions(selected, exp = "size_z_scored", model_name = "m4_exp_x_climate")
p5 <- plot_beta_model_predictions(selected, exp = "size_z_scored", model_name = "m5_exp_x_richness_noclimate")
p6 <- plot_beta_model_predictions(selected, exp = "size_z_scored", model_name = "m6_exp_richness_noclimate")

p1$plot
p2$plot
p3$plot
p4$plot
p5$plot
p6$plot

p_clim_exp <- plot_beta_climate_interaction(
  selected,
  exp = "size_z_scored",
  model_name = "m4_exp_x_climate"
)

p_clim_exp$plot

p_rich_clim <- plot_beta_climate_interaction(
  selected,
  exp = "site_nbspe_log_z_scored",
  model_name = "m2_richness_x_climate",
  richness = "size_z_scored" # held at median here
)
p_rich_clim$plot



##discharge vs. dC -----------
#check colinearity 
diag <- check_beta_collinearity(
  SiteVar_all,
  res = "propintraspecific_C",
  exp = "hydro_dis_z_scored"
)

diag$correlation
diag$vif

mods <- run_beta_models(
  SiteVar_all,
  res = "propintraspecific_C",
  exp = "hydro_dis_z_scored"
)

# inspect all models first
mods$models[[1]]
summary(mods$models[[1]]$m1_main)
summary(mods$models[[1]]$m2_richness_x_climate)
summary(mods$models[[1]]$m3_exp_x_richness)
summary(mods$models[[1]]$m4_exp_x_climate)
summary(mods$models[[1]]$m5_exp_x_richness_noclimate)
summary(mods$models[[1]]$m6_exp_richness_noclimate)


##min3
mods$models[[1]]
summary(mods$models[[2]]$m1_main)
summary(mods$models[[2]]$m2_richness_x_climate)
summary(mods$models[[2]]$m3_exp_x_richness)
summary(mods$models[[2]]$m4_exp_x_climate)
summary(mods$models[[2]]$m5_exp_x_richness_noclimate)
summary(mods$models[[2]]$m6_exp_richness_noclimate)

##min3
mods$models[[3]]
summary(mods$models[[3]]$m1_main)
summary(mods$models[[3]]$m2_richness_x_climate)
summary(mods$models[[3]]$m3_exp_x_richness)
summary(mods$models[[3]]$m4_exp_x_climate)
summary(mods$models[[3]]$m5_exp_x_richness_noclimate)
summary(mods$models[[3]]$m6_exp_richness_noclimate)

# then select by AIC
selected <- select_beta_models_aic(mods)

selected %>%
  select(
    data_group,
    response,
    explanatory,
    top_model_name,
    supported_model_names,
    n_supported_models,
    model_selection_result
  )

selected %>%
  select(data_group, response, explanatory, aic_table) %>%
  unnest(aic_table)

p1 <- plot_beta_model_predictions(selected, exp = "hydro_dis_z_scored", model_name = "m1_main")
p2 <- plot_beta_model_predictions(selected, exp = "hydro_dis_z_scored", model_name = "m2_richness_x_climate")
p3 <- plot_beta_model_predictions(selected, exp = "hydro_dis_z_scored", model_name = "m3_exp_x_richness")
p4 <- plot_beta_model_predictions(selected, exp = "hydro_dis_z_scored", model_name = "m4_exp_x_climate")
p5 <- plot_beta_model_predictions(selected, exp = "hydro_dis_z_scored", model_name = "m5_exp_x_richness_noclimate")
p6 <- plot_beta_model_predictions(selected, exp = "hydro_dis_z_scored", model_name = "m6_exp_richness_noclimate")


p1$plot
p2$plot
p3$plot
#p4$plot
p5$plot
p6$plot

p_clim_exp <- plot_beta_climate_interaction(
  selected,
  exp = "hydro_dis_z_scored",
  model_name = "m4_exp_x_climate"
)

p_clim_exp$plot

p_rich_clim <- plot_beta_climate_interaction(
  selected,
  exp = "site_nbspe_log_z_scored",
  model_name = "m2_richness_x_climate",
  richness = "hydro_dis_z_scored" # held at median here
)
p_rich_clim$plot

##discharge vs. dN -----------
#check colinearity 
diag <- check_beta_collinearity(
  SiteVar_all,
  res = "propintraspecific_N",
  exp = "hydro_dis_z_scored"
)

diag$correlation
diag$vif


mods <- run_beta_models(
  SiteVar_all,
  res = "propintraspecific_N",
  exp = "hydro_dis_z_scored"
)

# inspect all models first
mods$models[[1]]
summary(mods$models[[1]]$m1_main)
summary(mods$models[[1]]$m2_richness_x_climate)
summary(mods$models[[1]]$m3_exp_x_richness)
summary(mods$models[[1]]$m4_exp_x_climate)
summary(mods$models[[1]]$m5_exp_x_richness_noclimate)
summary(mods$models[[1]]$m6_exp_richness_noclimate)


##min3
mods$models[[2]]
summary(mods$models[[2]]$m1_main)
summary(mods$models[[2]]$m2_richness_x_climate)
summary(mods$models[[2]]$m3_exp_x_richness)
summary(mods$models[[2]]$m4_exp_x_climate)
summary(mods$models[[2]]$m5_exp_x_richness_noclimate)
summary(mods$models[[2]]$m6_exp_richness_noclimate)

##min3
mods$models[[3]]
summary(mods$models[[3]]$m1_main)
summary(mods$models[[3]]$m2_richness_x_climate)
summary(mods$models[[3]]$m3_exp_x_richness)
summary(mods$models[[3]]$m4_exp_x_climate)
summary(mods$models[[3]]$m5_exp_x_richness_noclimate)
summary(mods$models[[3]]$m6_exp_richness_noclimate)

# then select by AIC
selected <- select_beta_models_aic(mods)

selected %>%
  select(
    data_group,
    response,
    explanatory,
    top_model_name,
    supported_model_names,
    n_supported_models,
    model_selection_result
  )

selected %>%
  select(data_group, response, explanatory, aic_table) %>%
  unnest(aic_table)

p1 <- plot_beta_model_predictions(selected, exp = "hydro_dis_z_scored", model_name = "m1_main")
p2 <- plot_beta_model_predictions(selected, exp = "hydro_dis_z_scored", model_name = "m2_richness_x_climate")
p3 <- plot_beta_model_predictions(selected, exp = "hydro_dis_z_scored", model_name = "m3_exp_x_richness")
p4 <- plot_beta_model_predictions(selected, exp = "hydro_dis_z_scored", model_name = "m4_exp_x_climate")
p5 <- plot_beta_model_predictions(selected, exp = "hydro_dis_z_scored", model_name = "m5_exp_x_richness_noclimate")
p6 <- plot_beta_model_predictions(selected, exp = "hydro_dis_z_scored", model_name = "m6_exp_richness_noclimate")


p1$plot
p2$plot
p3$plot
p4$plot
p5$plot
p6$plot

p_clim_exp <- plot_beta_climate_interaction(
  selected,
  exp = "hydro_dis_z_scored",
  model_name = "m4_exp_x_climate"
)

p_clim_exp$plot

p_rich_clim <- plot_beta_climate_interaction(
  selected,
  exp = "site_nbspe_log_z_scored",
  model_name = "m2_richness_x_climate",
  richness = "hydro_dis_z_scored" # held at median here
)
p_rich_clim$plot




##Function to run all models:-----------
run_all_beta_model_sets <- function(dat,
                                    responses,
                                    drivers,
                                    richness = "site_nbspe_log_z_scored",
                                    climate = "Climate_zone_e2",
                                    delta_cutoff = 2) {
  
  crossing(
    response = responses,
    driver = drivers
  ) %>%
    mutate(
      mods = map2(
        response,
        driver,
        ~ run_beta_models(
          dat = dat,
          res = .x,
          exp = .y,
          richness = richness,
          climate = climate
        )
      ),
      
      selected = map(
        mods,
        ~ select_beta_models_aic(.x, delta_cutoff = delta_cutoff)
      )
    )
}

all_model_outputs <- run_all_beta_model_sets(
  dat = SiteVar_all,
  responses = c("propintraspecific_C", "propintraspecific_N"),
  drivers = c(
    "TP_z_scored",
    "hft_z_scored",
    "size_z_scored",
    "hydro_dis_z_scored"
  )
)

all_aic_summary <- all_model_outputs %>%
  select(selected) %>%
  unnest(selected) %>%
  rename(driver = explanatory) %>%
  select(
    response,
    driver,
    data_group,
    top_model_name,
    supported_model_names,
    n_supported_models,
    model_selection_result,
    aic_table
  ) %>%
  unnest(aic_table)


summarise_beta_model_coefficients <- function(all_model_outputs) {
  
  all_model_outputs %>%
    select(selected) %>%
    unnest(selected) %>%
    rename(driver = explanatory) %>%
    select(response, driver, data_group, models, aic_table) %>%
    mutate(
      coef_table = map2(models, aic_table, function(mods, aic_tab) {
        
        purrr::imap_dfr(mods, function(mod, mod_name) {
          
          if (is.null(mod)) return(tibble())
          
          coef_df <- broom::tidy(mod) %>%
            mutate(model_name = mod_name)
          
          fit_df <- aic_tab %>%
            filter(.data$model_name == mod_name) %>%
            select(model_name, AIC, delta_AIC, weight, supported)
          
          coef_df %>%
            left_join(fit_df, by = "model_name")
        })
      })
    ) %>%
    select(response, driver, data_group, coef_table) %>%
    unnest(coef_table)
}

coef_summary <- summarise_beta_model_coefficients(all_model_outputs)


write.csv(coef_summary, "outputs/beta_regression_summary_table.csv", row.names = FALSE)
##OLD  -----------------------



out <- run_beta_candidate_models(
  SiteVar_all %>% filter(data_group != "min_5"),
  res = "propintraspecific_C",
  exp = "TP"
)

out$aic
out$results %>%
  select(data_group, response, explanatory, best_model_name)

summary(out$results$best_model[[1]])
summary(out$results$best_model[[2]])

out$plot




out$best_model_name

test <- SiteVar_all %>%
  group_by(data_group, Climate_zone_e2) %>%
  count()


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



