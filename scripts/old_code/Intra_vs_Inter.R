#Load main packages------------
library(ggplot2)
library(readxl)
library(tidyr)
library(dplyr)
library(ggpubr)

#Load data and "normalize" d13C and d15N per site for same "baseline"
##per site: d13C_i=(d13C_i - d13C_min)/(d13C_max-d13C_min)
##per site: d15N_i=(d15N_i - d15N_min)/(d15N_max-d15N_min)

ALLindiv_June2025 <- read_excel("ALLindiv_June2025.xlsx")

DataFish<-subset(ALLindiv_June2025,!is.na(d15N) & !is.na(d13C) & organism_type=="fish" & !is.na(scientific_name))

DataFish$sp_site<-paste(DataFish$Scientific.nameFishBase,DataFish$collection_site_id,sep="_")
DataFish$num<-1

DataFish<-DataFish %>% 
  group_by(collection_site_id) %>% 
  mutate(d15N_norm = (d15N - min(d15N, na.rm = TRUE))/(max(d15N, na.rm = TRUE)-min(d15N, na.rm = TRUE)),
          d13C_norm = (d13C - min(d13C, na.rm = TRUE))/(max(d13C, na.rm = TRUE)-min(d13C, na.rm = TRUE)))


#I. Intraspecific variation per species------------
#some species have var = NA when only one individual sampled
#we ignore such species in average intra var but still include it in inter var

SpVar<-DataFish %>% 
  group_by(sp_site,collection_site_id,Scientific.nameFishBase,Fish.family,Ecosystem_Type,lentic_ecosystem_size_km2,lotic_ecosystem_width_m) %>% 
  summarise(sp_site_mean_N = mean(d15N_norm, na.rm = TRUE),
            sp_site_var_N = var(d15N_norm, na.rm = TRUE),
            sp_site_mean_C = mean(d13C_norm, na.rm = TRUE),
            sp_site_var_C = var(d13C_norm, na.rm = TRUE),
            sp_site_mean_length = mean(collected_sample_total_length, na.rm = TRUE),
            sp_site_var_length = var(collected_sample_total_length, na.rm = TRUE),
            collection_decimal_longitude = mean(collection_decimal_longitude, na.rm = TRUE),
            collection_decimal_latitude = mean(collection_decimal_latitude, na.rm = TRUE),
            sp_site_num = sum(num))

SpVar$VarTot<-SpVar$sp_site_var_N+SpVar$sp_site_var_C


#II. Intraspecific vs. interspecific variation per site------------
SpVar$num<-1
SiteVar <-SpVar %>% 
  group_by(collection_site_id,Ecosystem_Type,lentic_ecosystem_size_km2,lotic_ecosystem_width_m) %>% 
  summarise(site_interspe_var_N = var(sp_site_mean_N, na.rm = TRUE),
            site_intraspe_var_N = mean(sp_site_var_N,na.rm=TRUE),
            site_interspe_var_C = var(sp_site_mean_C, na.rm = TRUE),
            site_intraspe_var_C = mean(sp_site_var_C,na.rm=TRUE),
            site_nbspe = sum(num),
            collection_decimal_longitude = mean(collection_decimal_longitude, na.rm = TRUE),
            collection_decimal_latitude = mean(collection_decimal_latitude, na.rm = TRUE),
            site_mean_sample_id = mean(sp_site_num,na.rm=TRUE),
            site_min_sample_id = min(sp_site_num,na.rm=TRUE))

SiteVar$propintraspecific_N<-SiteVar$site_intraspe_var_N/(SiteVar$site_interspe_var_N+SiteVar$site_intraspe_var_N)
SiteVar$propintraspecific_C<-SiteVar$site_intraspe_var_C/(SiteVar$site_interspe_var_C+SiteVar$site_intraspe_var_C)
SiteVar$propintraspecific_Total<-(SiteVar$site_intraspe_var_N+SiteVar$site_intraspe_var_C)/(SiteVar$site_interspe_var_N+SiteVar$site_intraspe_var_N+SiteVar$site_interspe_var_C+SiteVar$site_intraspe_var_C)

SiteVar<-subset(SiteVar,site_nbspe>=2 & site_mean_sample_id>1)

write.table(SiteVar,file="Intraspecific_contribution_perSite.txt",row.names=FALSE,sep="\t")
write.table(SpVar,file="SpeciesIntraspecific_variance_perSite.txt",row.names=FALSE,sep="\t")

#III. Some plots and preliminary stats

p1<-ggplot(SiteVar,aes(x=Ecosystem_Type,y=propintraspecific_Total,col=Ecosystem_Type))+geom_boxplot()+theme_bw() +ylab("Intraspecific contribution Total")
p1

p1bis<-ggplot(SiteVar,aes(x=abs(collection_decimal_latitude),y=propintraspecific_Total,col=Ecosystem_Type))+geom_point(alpha=0.5) +
  theme_bw()+ 
  geom_smooth(method='lm',formula= y~x)+xlab("abs latitude")+ylab("Intraspecific contribution Total")
p1bis

p5<-ggplot(SiteVar,aes(x=site_nbspe,y=propintraspecific_Total,col=Ecosystem_Type))+geom_point(alpha=0.5) +
  scale_x_continuous(trans='log10')+theme_bw()+ 
  geom_smooth(method='lm',formula= y~x)+xlab("number of species")+ylab("Intraspecific contribution Total")
p5

p7<-ggplot(SiteVar,aes(x=site_mean_sample_id,y=propintraspecific_Total,col=Ecosystem_Type))+geom_point(alpha=0.5) +
  scale_x_continuous(trans='log10')+theme_bw()+ 
  geom_smooth(method='lm',formula= y~x)+xlab("mean number of individuals sampled per species")+ylab("Intraspecific contribution Total")
p7

library(DHARMa)
library(glmmTMB)
library(car)
library(performance)

m1<-glmmTMB(propintraspecific_Total~log(site_nbspe)+log(site_mean_sample_id)+abs(collection_decimal_latitude)*Ecosystem_Type,
            data=SiteVar)
summary(m1)
simulationOutput <- simulateResiduals(fittedModel = m1, n = 250)
plot(simulationOutput)
Anova(m1)

SpVar<-subset(SpVar,VarTot>0 & sp_site_var_length>0)
SpVar$log_VarTot<-log10(SpVar$VarTot)
SpVar$CV_length<-sqrt(SpVar$sp_site_var_length)/SpVar$sp_site_mean_length

m1<-glmmTMB(log_VarTot~log(CV_length)+log(sp_site_mean_length)
            +log(sp_site_num)+abs(collection_decimal_latitude)+(1|Fish.family)+(1|Scientific.nameFishBase)+(1|collection_site_id),data=SpVar)
summary(m1)
simulationOutput <- simulateResiduals(fittedModel = m1, n = 250)
plot(simulationOutput)
Anova(m1)
r2(m1)
