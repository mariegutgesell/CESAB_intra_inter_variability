##Relating intrasp. variability contribution to environmental drivers -- for single cutoff 

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
library(mgcv)
library(glmmTMB)
library(DHARMa)
library(marginaleffects)

##read in contribution of site level intraspecific contribution -- from cutoff chosen
SiteVar <- read.csv("data/SiteVar_2sp_75cutoff.csv") ##produced by 2_Intra_vs_Inter_21_01_26.R

##checking duplicate site names/lat longs -- don't worry about for now, come back to once using final data
##sometimes have same lat/long with 2 names, and sometimes have same site name with different lat/longs 
latlong1 <- SiteVar %>%
  select(collection_decimal_latitude, collection_decimal_longitude) %>%
  unique() 

##find which site ids have duplicate lat/longs
dup_latlong <- SiteVar %>%
  select(FWB_id, collection_site_id, collection_decimal_latitude, collection_decimal_longitude) %>%
  group_by(collection_decimal_latitude, collection_decimal_longitude) %>%
  filter(n() > 1) %>%
  arrange(collection_decimal_latitude, collection_decimal_longitude)
##8 cases with duplicate lat/longs -- some are sampled wet and dry season (may need to only select 1?), some have duplicate site names (see below), others not clear, do we need to think about this if they are truly independent sites? Or just not worry about it? 
##the issue is, is that the environmental data is joined by the latlong - actually this shouldnt matter, because both have the same environmental data - just need to double check the isotope data itself isn't getting duplicated ... checking back, i dont think it is so I think this is okay
##note: these issues disappear when we use minimum 3 or 5 samples

sites1 <- SiteVar%>%
  select(collection_site_id, waterbody_type) %>%
  group_by(collection_site_id, waterbody_type) %>%
  count()


##Environmental Variables
##productivity = TP
##type = ecosystem type (lake vs. stream)
##human ftp = hft
##size = size_z_scored -- were these standardized within ecosysetm type?
##flow variability = hydro_dis_z_scored - were these standardized within ecosysetm type?
##climate = temp , climate zone 

##need. to check if TP, sp richness and hft are signficantly different between ecosystem types so if should be standardized within ecosystem type 
##Data transformations 
SiteVar <- SiteVar %>%
  mutate(site_nbspe_log = log(site_nbspe),
         site_nbspe_log_z_scored = as.numeric(scale(site_nbspe_log,center = TRUE, scale = TRUE)),
         TP_z_scored = as.numeric(scale(TP, center = TRUE, scale = TRUE)),
         hft_z_scored = as.numeric(scale(hft, center = TRUE, scale = TRUE))
  ) %>%
  mutate(Climate_zone_2cat = case_when(
    startsWith(Climate_zone_e2, "Co") ~"Cold/Cool",
    startsWith(Climate_zone_e2, "Ho") ~"Warm/Hot",
    startsWith(Climate_zone_e2, "Wa") ~"Warm/Hot",
  )) %>%
  mutate(data_group = "sp2_75pct") 



##potential random effects: species richness, C:N, tissue type 
##Make correlation plot of numerical explanatory variables
exp_var <- SiteVar%>%
  filter(Size != "-Inf") %>%
  ungroup() %>%
  select(TP_z_scored, size_z_scored, hydro_dis_z_scored, hft_z_scored, site_nbspe_log_z_scored)# %>%
#  unique() #



duplicated_exp_variables <- SiteVar %>%
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

corr_mat <- cor(exp_var, use = "pairwise.complete.obs")

# visualize
corrplot(corr_mat, method = "color", type = "upper",  addCoef.col = "black",
         tl.col = "black", tl.srt = 45)



##Look at distributions of proportion intraspecific variation
ggplot(SiteVar, aes(x = propintraspecific_C)) +
  geom_histogram() 

ggplot(SiteVar, aes(x = site_intraspe_var_C)) +
  geom_histogram() 

ggplot(SiteVar, aes(x = propintraspecific_N)) +
  geom_histogram() 

ggplot(SiteVar, aes(x = site_intraspe_var_N)) +
  geom_histogram() 

ggplot(SiteVar, aes(x = TP_z_scored)) +
  geom_histogram()

ggplot(SiteVar, aes(x = propintraspecific_Total)) +
  geom_histogram() 

#SiteVar$Climate_zone_e2propintraspecific_Total#SiteVar$Climate_zone_e2 <- ordered(SiteVar$Climate_zone_e2, 
#                                       levels = c(
#                                         "Cold and wet/mesic",
#                                         "Cool and moist",
#                                         "Cool temperate and dry/xeric",
#                                         "Warm temperate",
#                                         "Hot and moist",
#                                         "Hot and dry"
#                                       ))
ggplot(SiteVar, aes(x = data_group, y = TP_z_scored, fill = Climate_zone_2cat)) +
  geom_boxplot() 


ggplot(SiteVar, aes(x = data_group, y = propintraspecific_C, fill = Climate_zone_2cat)) +
  geom_boxplot() 

ggplot(SiteVar, aes(x = data_group, y = propintraspecific_N, fill = Climate_zone_2cat)) +
  geom_boxplot() 

ggplot(SiteVar, aes(x = site_nbspe_log_z_scored, y = propintraspecific_N, color = Climate_zone_2cat)) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(SiteVar, aes(x = site_nbspe_log_z_scored, y = propintraspecific_C, color = Climate_zone_2cat)) +
  geom_point() +
  geom_smooth(method = "lm")



###Generalized linear model 


##CARBON 
glmm2_C <- glmmTMB(propintraspecific_C ~ log(site_nbspe) * Climate_zone_2cat * waterbody_type +(1 |FWB_id), data = SiteVar, family = beta_family(link = "logit"))
summary(glmm2_C)
Anova(glmm2_C)

glmm3_C <- glmmTMB(propintraspecific_C ~  log(site_nbspe)*Climate_zone_2cat + log(site_nbspe)*waterbody_type  + (1 |FWB_id), data = SiteVar, family = beta_family(link = "logit"))
summary(glmm3_C)
Anova(glmm3_C)

glmm4_C <- glmmTMB(propintraspecific_C ~  log(site_nbspe) + Climate_zone_2cat + log(site_nbspe)*waterbody_type  + (1 |FWB_id), data = SiteVar, family = beta_family(link = "logit"))
summary(glmm4_C)
Anova(glmm4_C)

glmm5_C <- glmmTMB(propintraspecific_C ~  site_nbspe_log , data = SiteVar, family = beta_family(link = "logit"))
summary(glmm5_C)
Anova(glmm5_C)

##final model -- that retains all three fixed effects
glmm6_C <- glmmTMB(propintraspecific_C ~  site_nbspe_log * Climate_zone_2cat + waterbody_type, data = SiteVar, family = beta_family(link = "logit"))
summary(glmm6_C)
Anova(glmm6_C)

min(SiteVar$site_nbspe_log)
max(SiteVar$site_nbspe_log)

nd_C <- datagrid(site_nbspe_log = seq(0.5, 4, 0.1), Climate_zone_2cat = "Warm/Hot", waterbody_type = "lotic", model = glmm6_C)
pred_C <- predict(glmm6_C, newdata = nd_C, type = "link", se.fit = TRUE)
#pred_C <- predict(glmm5_C, newdata = nd_C, type = "response", se.fit = TRUE)
nd_C <- nd_C %>%
  mutate(
    prediction = plogis(pred_C$fit),
    #se = pred_C$se.fit,
    lower = plogis(pred_C$fit - 1.96 * pred_C$se.fit),
    upper = plogis(pred_C$fit + 1.96 * pred_C$se.fit)
  )

c_site <- ggplot() +
  geom_point(data = SiteVar, aes(x = site_nbspe_log, y = propintraspecific_C, color= Climate_zone_2cat) ) +
  scale_color_manual(name = "Climate zone", values = c("#99CCCC","#993333" )) + 
  geom_ribbon(data = nd_C, aes(x = site_nbspe_log,ymin = lower, ymax = upper), alpha = 0.4) +
  geom_line(data = nd_C, aes(x = site_nbspe_log, y = prediction), color = "black")+
 # geom_abline(slope = -0.44131, intercept =  0.83235 , color = "black", linewidth = 2) +
  theme_classic() +
  xlab("Species Richness (log)") +
  ylab("CIV in Carbon") + 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))
c_site

c_violin <- ggplot(SiteVar, aes(y = propintraspecific_C, x = "")) +
  geom_violin() +
  geom_jitter(width = 0.1, alpha = 0.2, size = 1.5) +
  stat_summary(fun = median, geom = "point", size = 5, shape = 18) +
  theme_classic() +
  labs(y = "CIV in Carbon", x = "") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+ 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))
c_violin
##NITROGEN 

glmm2_N <- glmmTMB(propintraspecific_N ~ log(site_nbspe) * Climate_zone_2cat * waterbody_type +(1 |FWB_id), data = SiteVar, family = beta_family(link = "logit"))
summary(glmm2_N)
Anova(glmm2_N)

glmm3_N <- glmmTMB(propintraspecific_N ~ log(site_nbspe) + Climate_zone_2cat + waterbody_type + log(site_nbspe)*Climate_zone_2cat + log(site_nbspe)*waterbody_type + Climate_zone_2cat*waterbody_type + (1 |FWB_id), data = SiteVar, family = beta_family(link = "logit"))
summary(glmm3_N)
Anova(glmm3_N)


glmm4_N <- glmmTMB(propintraspecific_N ~ log(site_nbspe) + Climate_zone_2cat + waterbody_type + log(site_nbspe)*Climate_zone_2cat + Climate_zone_2cat*waterbody_type + (1 |FWB_id), data = SiteVar, family = beta_family(link = "logit"))
summary(glmm4_N)
Anova(glmm4_N)

glmm5_N <- glmmTMB(propintraspecific_N ~  log(site_nbspe)*waterbody_type + log(site_nbspe)*Climate_zone_2cat  + (1 |FWB_id), data = SiteVar, family = beta_family(link = "logit"))
summary(glmm5_N)
Anova(glmm5_N)

glmm6_N <- glmmTMB(propintraspecific_N ~   log(site_nbspe)*Climate_zone_2cat +waterbody_type + (1 |FWB_id), data = SiteVar, family = beta_family(link = "logit"))
summary(glmm6_N)
Anova(glmm6_N)

glmm7_N <- glmmTMB(propintraspecific_N ~   log(site_nbspe)*Climate_zone_2cat, data = SiteVar, family = beta_family(link = "logit"))
summary(glmm7_N)
Anova(glmm7_N)


glmm8_N <- glmmTMB(propintraspecific_N ~   site_nbspe_log*Climate_zone_2cat, data = SiteVar, family = beta_family(link = "logit"))
summary(glmm8_N)
Anova(glmm8_N)
simulationOutput <- simulateResiduals(fittedModel = glmm8_N, plot = FALSE)
plot(simulationOutput)
#DHARMa::residuals(simulationOutput)

##final model -- retain significant interaction and fixed effect of ecosystem type 
glmm9_N <- glmmTMB(propintraspecific_N ~   site_nbspe_log*Climate_zone_2cat + waterbody_type, data = SiteVar, family = beta_family(link = "logit"))
summary(glmm9_N)
Anova(glmm9_N)


nd_N <- datagrid(site_nbspe_log = seq(0.5, 4, 0.1), Climate_zone_2cat = c("Cold/Cool", "Warm/Hot"),waterbody_type = "lotic", model = glmm9_N)
pred_N <- predict(glmm9_N, newdata = nd_N, type = "link", se.fit = TRUE)
#pred_N <- predict(glmm8_N, newdata = nd_N, type = "response", se.fit = TRUE)

nd_N <- nd_N %>%
  mutate(
    prediction = plogis(pred_N$fit),
    lower = plogis(pred_N$fit - 1.96 * pred_N$se.fit),
    upper = plogis(pred_N$fit + 1.96 * pred_N$se.fit)
  )


n_site <- ggplot() +
  geom_point(data = SiteVar, aes(x = site_nbspe_log, y = propintraspecific_N, color= Climate_zone_2cat) ) +
 # scale_color_manual(name = "Climate zone", values = c("#99CCCC","#993333" )) + 
  geom_ribbon(data = nd_N, aes(x = site_nbspe_log,ymin = lower, ymax = upper, fill = Climate_zone_2cat), alpha = 0.4) +
   geom_line(data = nd_N, aes(x = site_nbspe_log, y = prediction, group = Climate_zone_2cat, color = Climate_zone_2cat))+
  scale_color_manual(name ="Climate zone", values = c("#99CCCC","#993333" )) + 
  scale_fill_manual(name ="Climate zone", values = c("#99CCCC","#993333" )) + 
  # geom_abline(slope = -0.44131, intercept =  0.83235 , color = "black", linewidth = 2) +
  theme_classic() +
  xlab("Species Richness (log)") +
  ylab("CIV in Nitrogen")+ 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))
n_site

n_violin <- ggplot(SiteVar, aes(y = propintraspecific_N, x = "")) +
  geom_violin() +
  geom_jitter(width = 0.1, alpha = 0.2, size = 1.5) +
  stat_summary(fun = median, geom = "point", size = 5, shape = 18) +
  theme_classic() +
  labs(y = "CIV in Nitrogen", x = "") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+ 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))
n_violin


c_plot <- ggarrange(c_violin, c_site, nrow = 1, ncol = 2 ,  widths = c(1, 3))
c_plot
n_plot <- ggarrange(n_violin, n_site, nrow = 1, ncol = 2 , widths = c(1, 3))
n_plot

ggarrange(c_plot, n_plot, labels = c("a", "b"),font.label = list(size = 18, face = "bold"), nrow = 2, ncol =1)

###2D analysis -- both isotopes combined
##final model -- retain significant interaction and fixed effect of ecosystem type 
glmm_t1 <- glmmTMB(propintraspecific_Total ~   site_nbspe_log*Climate_zone_2cat * waterbody_type, data = SiteVar, family = beta_family(link = "logit"))
summary(glmm_t1)
Anova(glmm_t1)




glmm_t <- glmmTMB(propintraspecific_Total ~   site_nbspe_log*Climate_zone_2cat + waterbody_type, data = SiteVar, family = beta_family(link = "logit"))
summary(glmm_t)
Anova(glmm_t)

max(SiteVar$site_nbspe_log)
min(SiteVar$site_nbspe_log)

nd_T <- datagrid(site_nbspe_log = seq(0.5, 4, 0.1), Climate_zone_2cat = c("Cold/Cool", "Warm/Hot"),waterbody_type = "lotic", model = glmm_t)
pred_T <- predict(glmm_t, newdata = nd_T, type = "link", se.fit = TRUE)
#pred_N <- predict(glmm8_N, newdata = nd_N, type = "response", se.fit = TRUE)

nd_T <- nd_T %>%
  mutate(
    prediction = plogis(pred_T$fit),
    lower = plogis(pred_T$fit - 1.96 * pred_T$se.fit),
    upper = plogis(pred_T$fit + 1.96 * pred_T$se.fit)
  )


total_site <- ggplot() +
  geom_point(data = SiteVar, aes(x = site_nbspe_log, y = propintraspecific_Total, color= Climate_zone_2cat) ) +
  # scale_color_manual(name = "Climate zone", values = c("#99CCCC","#993333" )) + 
  geom_ribbon(data = nd_T, aes(x = site_nbspe_log,ymin = lower, ymax = upper, fill = Climate_zone_2cat), alpha = 0.4) +
  geom_line(data = nd_T, aes(x = site_nbspe_log, y = prediction, group = Climate_zone_2cat, color = Climate_zone_2cat))+
  scale_color_manual(name ="Climate zone", values = c("#99CCCC","#993333" )) + 
  scale_fill_manual(name ="Climate zone", values = c("#99CCCC","#993333" )) + 
  # geom_abline(slope = -0.44131, intercept =  0.83235 , color = "black", linewidth = 2) +
  theme_classic() +
  xlab("Species Richness (log)") +
  ylab("Total CIV")+ 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))
total_site

total_violin <- ggplot(SiteVar, aes(y = propintraspecific_Total, x = "")) +
  geom_violin() +
  geom_jitter(width = 0.1, alpha = 0.2, size = 1.5) +
  stat_summary(fun = median, geom = "point", size = 5, shape = 18) +
  theme_classic() +
  labs(y = "Total CIV", x = "") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )+ 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))
total_violin




##difference between lotic and lentic?
ggplot(SiteVar, aes(x = waterbody_type, y = propintraspecific_Total)) +
  geom_boxplot()

nd_T_2 <- datagrid(site_nbspe_log = seq(0.5, 4, 0.1), Climate_zone_2cat = c("Cold/Cool"),waterbody_type = c("lotic", "lentic"), model = glmm_t)
pred_T_2 <- predict(glmm_t, newdata = nd_T_2, type = "link", se.fit = TRUE)
#pred_N <- predict(glmm8_N, newdata = nd_N, type = "response", se.fit = TRUE)

nd_T_2 <- nd_T_2 %>%
  mutate(
    prediction = plogis(pred_T_2$fit),
    lower = plogis(pred_T_2$fit - 1.96 * pred_T_2$se.fit),
    upper = plogis(pred_T_2$fit + 1.96 * pred_T_2$se.fit)
  )

total_site_wb <- ggplot() +
  geom_point(data = SiteVar, aes(x = site_nbspe_log, y = propintraspecific_Total, color= waterbody_type) ) +
  # scale_color_manual(name = "Climate zone", values = c("#99CCCC","#993333" )) + 
  geom_ribbon(data = nd_T_2, aes(x = site_nbspe_log,ymin = lower, ymax = upper, fill = waterbody_type), alpha = 0.4) +
  geom_line(data = nd_T_2, aes(x = site_nbspe_log, y = prediction, group = waterbody_type, color = waterbody_type))+
  scale_color_manual(name ="waterbody", values = c("darkgrey","black" )) + 
  scale_fill_manual(name ="waterbody", values = c("darkgrey","black" )) + 
  # geom_abline(slope = -0.44131, intercept =  0.83235 , color = "black", linewidth = 2) +
  theme_classic() +
  xlab("Species Richness (log)") +
  ylab("Total CIV")+ 
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16))
total_site_wb


sp_total <- ggarrange(total_site, total_site_wb, labels = c("b", "c"),font.label = list(size = 18, face = "bold"), nrow = 2, ncol = 1 )

ggarrange(total_violin, sp_total, labels = c("a", ""),font.label = list(size = 18, face = "bold"), nrow = 1, ncol = 2, widths = c(1, 3)  )

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
sp_richness_C <- beta_reg_univariate(SiteVar, res = "propintraspecific_C", exp = "site_nbspe_log")
sp_richness_C$models
sp_richness_C$plot

sp_richness_N <- beta_reg_univariate(SiteVar, res = "propintraspecific_N", exp = "site_nbspe_log")
sp_richness_N$models
sp_richness_N$plot

sp_richness_N <- beta_reg_univariate(SiteVar, res = "propintraspecific_Total", exp = "site_nbspe_log")
sp_richness_N$models
sp_richness_N$plot


##Total Phosphorous
tp_C <- beta_reg_univariate(SiteVar, res = "propintraspecific_C", exp = "TP")
tp_C$models
tp_C$plot

tp_N <- beta_reg_univariate(SiteVar, res = "propintraspecific_N", exp = "TP")
tp_N$models
tp_N$plot

##Human Footprint
hft_C <- beta_reg_univariate(SiteVar, res = "propintraspecific_C", exp = "hft")
hft_C$models
hft_C$plot


hft_N <- beta_reg_univariate(SiteVar, res = "propintraspecific_N", exp = "hft")
hft_N$models
hft_N$plot

##Flow Variability
fv_C <- beta_reg_univariate(SiteVar, res = "propintraspecific_C", exp = "hydro_dis_z_scored")
fv_C$models
fv_C$plot


fv_N <- beta_reg_univariate(SiteVar, res = "propintraspecific_N", exp = "hydro_dis_z_scored")
fv_N$models
fv_N$plot

##Ecosystem Size
es_C <- beta_reg_univariate(SiteVar, res = "propintraspecific_C", exp = "size_z_scored")
es_C$models
es_C$plot

es_N <- beta_reg_univariate(SiteVar, res = "propintraspecific_N", exp = "size_z_scored")
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

out_C_climate <- aov_univariate(SiteVar, res = "propintraspecific_C", exp = "Climate_zone_e2",order_levels = climate_levels)
out_C_climate$plot
out_C_climate$plot2
out_C_climate$tukey

out_N_climate <- aov_univariate(SiteVar, res = "propintraspecific_N", exp = "Climate_zone_e2",order_levels = climate_levels)
out_N_climate$plot
out_N_climate$plot2
out_N_climate$tukey

##waterbody type
out_C_type <- aov_univariate(SiteVar, res = "propintraspecific_C", exp = "waterbody_type")
out_C_type$plot
out_C_type$plot2
out_C_type$tukey

out_N_type <- aov_univariate(SiteVar, res = "propintraspecific_N", exp = "waterbody_type")
out_N_type$plot
out_N_type$plot2
out_N_type$tukey

##Note: ANOVA assumes normal residuals, may need to transform proportions (logit) or use beta regressions w/ categorical predictor

##Does the effect of species richness differ my climate or habitat type?
sp_rich_climate <- aov_univariate(SiteVar, res = "site_nbspe_log", exp = "Climate_zone_e2", order_levels = climate_levels)
sp_rich_climate$plot
sp_rich_climate$plot2
sp_rich_climate$tukey

sp_rich_type <- aov_univariate(SiteVar, res = "site_nbspe_log", exp = "waterbody_type")
sp_rich_type$plot
sp_rich_type$plot2
sp_rich_type$tukey


###Multiple Regression Models ------------------------
##Look at series of potential regression models - hypothesis driven
##does the environmental variable affect the response, after accounting for species richness and climate? is there an overall relationship between environmental driver and contribution fo intrasp. var to total var? 
#m1 <- betareg(y ~ exp + richness + climate)
#does the effect of richness on the response differ among climate zones, while still including an overall efffect of exp? does species richness matter differently in tropical vs. temperate vs. cold systems?
#m2 <- betareg(y ~ exp + richness * climate)

##does the effect of the environmental variable depend on species richness? does environmental variable affect intraspecific variability differently in species rich vs. species poor communities? 
#m3 <- betareg(y ~ exp * richness + climate)

#m4 <- betareg()

##these next two are not really needed...not really hypothesis driven unless we specifically want to exclude climate? 
#does environmental variable explain response after accountinf for richness? ignoring cliamte differences, do TP and richness explain variation? 
#m4 <- betareg(y ~ exp + richness)

#does the relationship between environmental variable and the response depend on richness, without accounting for cliamte? is there an itneraction between TP and richness, if we ignore climatic context? 
#m5 <- betareg(y ~ exp + richness)
##do all explanatory and species richness etc variables need to be scaled? i think if including them in one model..yes... 

##Checking Model Assumptions:
##1) ensure response variable is 0 < y < 1 (can't have exactly 0 or 1)

##2) ensure no colinearity between predictors 
##function to check colinearity:
check_beta_collinearity <- function(dat,
                                    res,
                                    exp,
                                    richness = "site_nbspe_log_z_scored",
                                    climate = "Climate_zone_2cat") {
  
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
site_climate_summary <- SiteVar %>%
  group_by(data_group, Climate_zone_2cat) %>%
  count()
#for min  = 5, don't have enough sites in each climate category ... so i don't think climate interactions are interpretable 

##Function to run model candidate types,
run_beta_models <- function(dat,
                            res,
                            exp,
                            richness = "site_nbspe_log_z_scored",
                            climate = "Climate_zone_2cat") {
  
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
                                        climate = "Climate_zone_2cat",
                                        n_points = 100) {
  
  make_predictions <- function(d, mods, res) {
    
    mod <- mods[[model_name]]
    if (is.null(mod)) return(tibble())
    
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
  
  extract_model_stats <- function(d, mods, res) {
    
    mod <- mods[[model_name]]
    if (is.null(mod)) return(tibble())
    
    coef_tab <- summary(mod)$coefficients$mean
    
    vars_to_show <- rownames(coef_tab) %>%
      stringr::str_subset(paste(c(exp, richness, climate), collapse = "|"))
    
    stats_text <- coef_tab[vars_to_show, , drop = FALSE] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("term") %>%
      mutate(
        label = paste0(
          term,
          ": est = ", round(Estimate, 3),
          ", p = ", signif(`Pr(>|z|)`, 3)
        )
      ) %>%
      pull(label) %>%
      paste(collapse = "\n")
    
    pseudo_r2 <- tryCatch(
      summary(mod)$pseudo.r.squared,
      error = function(e) NA_real_
    )
    
    tibble(
      data_group = unique(d$data_group),
      response = res,
      explanatory = exp,
      model_name = model_name,
      label = paste0(
        stats_text,
        "\nPseudo R² = ", round(pseudo_r2, 3)
      )
    )
  }
  
  pred_df <- purrr::pmap_dfr(
    list(
      d = selected$data,
      mods = selected$models,
      res = selected$response
    ),
    make_predictions
  )
  
  stats_df <- purrr::pmap_dfr(
    list(
      d = selected$data,
      mods = selected$models,
      res = selected$response
    ),
    extract_model_stats
  )
  
  plot_df <- selected %>%
    select(data) %>%
    unnest(data)
  
  res <- unique(selected$response)
  
  stats_df <- stats_df %>%
    mutate(
      x = min(plot_df[[exp]], na.rm = TRUE),
      y = max(plot_df[[res]], na.rm = TRUE) -
        (row_number() - 1) * 0.10 * diff(range(plot_df[[res]], na.rm = TRUE))
    )
  
  p <- ggplot(plot_df, aes(x = .data[[exp]], y = .data[[res]], color = data_group)) +
    geom_point(alpha = 0.45) +
    geom_line(
      data = pred_df,
      aes(x = .data[[exp]], y = pred, color = data_group),
      linewidth = 1.2
    ) +
    geom_label(
      data = stats_df,
      aes(x = x, y = y, label = label, color = data_group),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3,
      show.legend = FALSE
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
    model_stats = stats_df,
    plot = p
  )
}

##function to be able to show/plot interactions - how relationships change between climate zones 
plot_beta_climate_interaction <- function(selected,
                                          exp,
                                          model_name,
                                          richness = "site_nbspe_log_z_scored",
                                          climate = "Climate_zone_2cat",
                                          n_points = 100) {
  
  make_predictions <- function(d, mods, res) {
    
    mod <- mods[[model_name]]
    if (is.null(mod)) return(tibble())
    
    newdat <- tidyr::expand_grid(
      !!exp := seq(min(d[[exp]], na.rm = TRUE),
                   max(d[[exp]], na.rm = TRUE),
                   length.out = n_points),
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
  
  extract_model_stats <- function(d, mods, res) {
    
    mod <- mods[[model_name]]
    if (is.null(mod)) return(tibble())
    
    coef_tab <- summary(mod)$coefficients$mean
    
    vars_to_show <- rownames(coef_tab) %>%
      stringr::str_subset(paste(c(exp, richness, climate), collapse = "|"))
    
    stats_text <- coef_tab[vars_to_show, , drop = FALSE] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("term") %>%
      mutate(
        label = paste0(
          term,
          ": est = ", round(Estimate, 3),
          ", p = ", signif(`Pr(>|z|)`, 3)
        )
      ) %>%
      pull(label) %>%
      paste(collapse = "\n")
    
    pseudo_r2 <- tryCatch(
      summary(mod)$pseudo.r.squared,
      error = function(e) NA_real_
    )
    
    tibble(
      data_group = unique(d$data_group),
      label = paste0(
        stats_text,
        "\nPseudo R² = ", round(pseudo_r2, 3)
      ),
      x = min(d[[exp]], na.rm = TRUE),
      y = max(d[[res]], na.rm = TRUE)
    )
  }
  
  pred_df <- purrr::pmap_dfr(
    list(
      d = selected$data,
      mods = selected$models,
      res = selected$response
    ),
    make_predictions
  )
  
  stats_df <- purrr::pmap_dfr(
    list(
      d = selected$data,
      mods = selected$models,
      res = selected$response
    ),
    extract_model_stats
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
    geom_label(
      data = stats_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3
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
    model_stats = stats_df,
    plot = p
  )
}

##for now, lets just focus on min1 and min3, not enough samples to look at effects of climate in min 5 


##TP vs. dC -----------
#check colinearity 
diag <- check_beta_collinearity(
  SiteVar,
  res = "propintraspecific_C",
  exp = "TP_z_scored"
)

diag$correlation
diag$vif

mods <- run_beta_models(
  SiteVar,
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
  SiteVar,
  res = "propintraspecific_N",
  exp = "TP_z_scored"
)

diag$correlation

diag$vif


mods <- run_beta_models(
  SiteVar,
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
  SiteVar,
  res = "propintraspecific_C",
  exp = "hft_z_scored"
)

diag$correlation
diag$vif


mods <- run_beta_models(
  SiteVar,
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
  SiteVar,
  res = "propintraspecific_N",
  exp = "hft_z_scored"
)

diag$correlation
diag$vif


mods <- run_beta_models(
  SiteVar,
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
  SiteVar,
  res = "propintraspecific_C",
  exp = "size_z_scored"
)

diag$correlation
diag$vif

mods <- run_beta_models(
  SiteVar,
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
  SiteVar,
  res = "propintraspecific_N",
  exp = "size_z_scored"
)

diag$correlation
diag$vif

mods <- run_beta_models(
  SiteVar,
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
  SiteVar,
  res = "propintraspecific_C",
  exp = "hydro_dis_z_scored"
)

diag$correlation
diag$vif

mods <- run_beta_models(
  SiteVar,
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
  SiteVar,
  res = "propintraspecific_N",
  exp = "hydro_dis_z_scored"
)

diag$correlation
diag$vif


mods <- run_beta_models(
  SiteVar,
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
                                    climate = "Climate_zone_2cat",
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
  dat = SiteVar,
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



coef_summary_mod5 <- coef_summary %>%
  filter(model_name == "m5_exp_x_richness_noclimate")
##function to make fig 
plot_one_beta_prediction <- function(all_model_outputs,
                                     response_choice,
                                     driver_choice,
                                     model_name,
                                     richness = "site_nbspe_log_z_scored",
                                     climate = "Climate_zone_2cat",
                                     n_points = 100) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  selected_row <- all_model_outputs %>%
    select(selected) %>%
    unnest(selected) %>%
    mutate(driver = explanatory) %>%
    filter(
      response == response_choice,
      driver == driver_choice
    )
  
  if (nrow(selected_row) != 1) {
    stop(
      paste0(
        "Expected exactly 1 row, but found ", nrow(selected_row),
        ". Check response_choice and driver_choice."
      )
    )
  }
  
  d <- selected_row$data[[1]]
  models <- selected_row$models[[1]]
  mod <- models[[model_name]]
  
  if (is.null(mod)) {
    stop("Model name not found in this model list.")
  }
  
  x_seq <- seq(
    min(d[[driver_choice]], na.rm = TRUE),
    max(d[[driver_choice]], na.rm = TRUE),
    length.out = n_points
  )
  
  newdat <- tibble(x_value = x_seq)
  newdat[[driver_choice]] <- x_seq
  
  if (richness %in% all.vars(formula(mod))) {
    newdat[[richness]] <- median(d[[richness]], na.rm = TRUE)
  }
  
  if (climate %in% all.vars(formula(mod))) {
    newdat[[climate]] <- names(sort(table(d[[climate]]), decreasing = TRUE))[1]
  }
  
  newdat$pred <- predict(mod, newdata = newdat, type = "response")
  
  plot_df <- tibble(
    x_value = d[[driver_choice]],
    y_value = d[[response_choice]]
  )
  
  ggplot(plot_df, aes(x = x_value, y = y_value)) +
    geom_point(alpha = 0.35, size = 1.2) +
    geom_line(
      data = newdat,
      aes(x = x_value, y = pred),
      linewidth = 1.1, 
      color = "black",
      linetype = "dashed",
    ) +
    theme_classic() +
    labs(
      x = driver_choice,
      y = "% CIV"
    )
}
p1 <- plot_one_beta_prediction(
  all_model_outputs = all_model_outputs,
  response_choice = "propintraspecific_C",
  driver_choice = "TP_z_scored",
  model_name = "m5_exp_x_richness_noclimate"
)

p1

p2 <- plot_one_beta_prediction(
  all_model_outputs = all_model_outputs,
  response_choice = "propintraspecific_N",
  driver_choice = "TP_z_scored",
  model_name = "m5_exp_x_richness_noclimate"
)

p2

p3 <- plot_one_beta_prediction(
  all_model_outputs = all_model_outputs,
  response_choice = "propintraspecific_C",
  driver_choice = "hft_z_scored",
  model_name = "m5_exp_x_richness_noclimate"
)

p3

p4 <- plot_one_beta_prediction(
  all_model_outputs = all_model_outputs,
  response_choice = "propintraspecific_N",
  driver_choice = "hft_z_scored",
  model_name = "m5_exp_x_richness_noclimate"
)

p4

p5 <- plot_one_beta_prediction(
  all_model_outputs = all_model_outputs,
  response_choice = "propintraspecific_C",
  driver_choice = "size_z_scored",
  model_name = "m5_exp_x_richness_noclimate"
)

p5

p6 <- plot_one_beta_prediction(
  all_model_outputs = all_model_outputs,
  response_choice = "propintraspecific_N",
  driver_choice = "size_z_scored",
  model_name = "m5_exp_x_richness_noclimate"
)

p6

p7 <- plot_one_beta_prediction(
  all_model_outputs = all_model_outputs,
  response_choice = "propintraspecific_C",
  driver_choice = "hydro_dis_z_scored",
  model_name = "m5_exp_x_richness_noclimate"
)

p7

p8 <- plot_one_beta_prediction(
  all_model_outputs = all_model_outputs,
  response_choice = "propintraspecific_N",
  driver_choice = "hydro_dis_z_scored",
  model_name = "m5_exp_x_richness_noclimate"
)

p8

plot_env <- ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, legend = "none", nrow = 4, ncol = 2, labels = c("a)", "b)", "c)", "d)", "e)", "f)", "g)", "h)"), font.label = list(colour = "black", size = 12))
plot_env




all_model_outputs %>%
  select(selected) %>%
  unnest(selected) %>%
  slice(1) %>%
  pull(models) %>%
  .[[1]] %>%
  names()
