
library(lme4)
library(cati)
library(ade4)
library(adegraphics)


# Select traits and variables for MULTIVARIATE analyses
#traits.fish <- total.sp[,c(4,20,27,29,22,31,35,36,32,33,44,40,42,39)]
traits.fish <- total.sp[,c(20,27,29,22,31,35,36,32,33,40,42,39)]
traits.fish <- traits.fish[complete.cases(traits.fish),]
species <- traits.fish$species
habitat <- traits.fish$group
site <- traits.fish$Site1


#####  GLMM for Plot 1 (inter vs intra by TRAITS) (Albert 2010)

traits <- c('BD','PFL','Throt','ED','EP','GS','GSh','MP','BS')
variables <- c('Interspecific var.', 'Intraspecific var./ site', 'Intraspecific var. unexplained')
m = matrix(, nrow = 9, ncol = 3)


#mod1 <- lmer(log(Length_SL)~species  +(1| Site1), data= total.sp)
#VarF <- var(as.vector(fixef(mod1) %*% t(getME(mod1,"X"))))
#m[1,1] <- VarF/(VarF + VarCorr(mod1)$Site1[1] + attr(VarCorr(mod1), "sc")^2)
#m[1,2] <- VarCorr(mod1)$Site1[1]/(VarF + VarCorr(mod1)$Site1[1] + attr(VarCorr(mod1), #"sc")^2)
#m[1,3] <- attr(VarCorr(mod1), "sc")^2/(VarF + VarCorr(mod1)$Site1[1] + attr(VarCorr(mod1), "sc")^2)

mod2 <- lmer(BD~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod2) %*% t(getME(mod2,"X"))))
m[1,1] <- VarF/(VarF + VarCorr(mod2)$Site1[1] + attr(VarCorr(mod2), "sc")^2)
m[1,2] <- VarCorr(mod2)$Site1[1]/(VarF + VarCorr(mod2)$Site1[1] + attr(VarCorr(mod2), "sc")^2)
m[1,3] <- attr(VarCorr(mod2), "sc")^2/(VarF + VarCorr(mod2)$Site1[1] + attr(VarCorr(mod2), "sc")^2)

mod3 <- lmer(PFL~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod3) %*% t(getME(mod3,"X"))))
m[2,1] <- VarF/(VarF + VarCorr(mod3)$Site1[1] + attr(VarCorr(mod3), "sc")^2)
m[2,2] <- VarCorr(mod3)$Site1[1]/(VarF + VarCorr(mod3)$Site1[1] + attr(VarCorr(mod3), "sc")^2)
m[2,3] <- attr(VarCorr(mod3), "sc")^2/(VarF + VarCorr(mod3)$Site1[1] + attr(VarCorr(mod3), "sc")^2)

mod4 <- lmer(throt~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod4) %*% t(getME(mod4,"X"))))
m[3,1] <- VarF/(VarF + VarCorr(mod4)$Site1[1] + attr(VarCorr(mod4), "sc")^2)
m[3,2] <- VarCorr(mod4)$Site1[1]/(VarF + VarCorr(mod4)$Site1[1] + attr(VarCorr(mod4), "sc")^2)
m[3,3] <- attr(VarCorr(mod4), "sc")^2/(VarF + VarCorr(mod4)$Site1[1] + attr(VarCorr(mod4), "sc")^2)

mod5 <- lmer(ED~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod5) %*% t(getME(mod5,"X"))))
m[4,1] <- VarF/(VarF + VarCorr(mod5)$Site1[1] + attr(VarCorr(mod5), "sc")^2)
m[4,2] <- VarCorr(mod5)$Site1[1]/(VarF + VarCorr(mod5)$Site1[1] + attr(VarCorr(mod5), "sc")^2)
m[4,3] <- attr(VarCorr(mod5), "sc")^2/(VarF + VarCorr(mod5)$Site1[1] + attr(VarCorr(mod5), "sc")^2)


mod6 <- lmer(EP~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod6) %*% t(getME(mod6,"X"))))
m[5,1] <- VarF/(VarF + VarCorr(mod6)$Site1[1] + attr(VarCorr(mod6), "sc")^2)
m[5,2] <- VarCorr(mod6)$Site1[1]/(VarF + VarCorr(mod6)$Site1[1] + attr(VarCorr(mod6), "sc")^2)
m[5,3] <- attr(VarCorr(mod6), "sc")^2/(VarF + VarCorr(mod6)$Site1[1] + attr(VarCorr(mod6), "sc")^2)

mod7 <- lmer(log(GS)~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod7) %*% t(getME(mod7,"X"))))
m[6,1] <- VarF/(VarF + VarCorr(mod7)$Site1[1] + attr(VarCorr(mod7), "sc")^2)
m[6,2] <- VarCorr(mod7)$Site1[1]/(VarF + VarCorr(mod7)$Site1[1] + attr(VarCorr(mod7), "sc")^2)
m[6,3] <- attr(VarCorr(mod7), "sc")^2/(VarF + VarCorr(mod7)$Site1[1] + attr(VarCorr(mod7), "sc")^2)

mod8 <- lmer(log(GSh)~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod8) %*% t(getME(mod8,"X"))))
m[7,1] <- VarF/(VarF + VarCorr(mod8)$Site1[1] + attr(VarCorr(mod8), "sc")^2)
m[7,2] <- VarCorr(mod8)$Site1[1]/(VarF + VarCorr(mod8)$Site1[1] + attr(VarCorr(mod8), "sc")^2)
m[7,3] <- attr(VarCorr(mod8), "sc")^2/(VarF + VarCorr(mod8)$Site1[1] + attr(VarCorr(mod8), "sc")^2)

mod9 <- lmer(MP~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod9) %*% t(getME(mod9,"X"))))
m[8,1] <- VarF/(VarF + VarCorr(mod9)$Site1[1] + attr(VarCorr(mod9), "sc")^2)
m[8,2] <- VarCorr(mod9)$Site1[1]/(VarF + VarCorr(mod9)$Site1[1] + attr(VarCorr(mod9), "sc")^2)
m[8,3] <- attr(VarCorr(mod9), "sc")^2/(VarF + VarCorr(mod9)$Site1[1] + attr(VarCorr(mod9), "sc")^2)

mod10 <- lmer(BS~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod10) %*% t(getME(mod10,"X"))))
m[9,1] <- VarF/(VarF + VarCorr(mod10)$Site1[1] + attr(VarCorr(mod10), "sc")^2)
m[9,2] <- VarCorr(mod10)$Site1[1]/(VarF + VarCorr(mod10)$Site1[1] + attr(VarCorr(mod10), "sc")^2)
m[9,3] <- attr(VarCorr(mod10), "sc")^2/(VarF + VarCorr(mod10)$Site1[1] + attr(VarCorr(mod10), "sc")^2)

#mod11 <- lmer(Rel_Gut_Length~species+(1| Site1), data= total.sp)
#VarF <- var(as.vector(fixef(mod10) %*% t(getME(mod10,"X"))))
#m[11,1] <- VarF/(VarF + VarCorr(mod11)$Site1[1] + attr(VarCorr(mod11), "sc")^2)
#m[11,2] <- VarCorr(mod11)$Site1[1]/(VarF + VarCorr(mod11)$Site1[1] + attr(VarCorr(#mod11), "sc")^2)
#m[11,3] <- attr(VarCorr(mod11), "sc")^2/(VarF + VarCorr(mod11)$Site1[1] + attr(VarCorr(mod11), "sc")^2)

colnames(m) <- variables
rownames(m) <- traits

# Plot Single traits
par(xpd=TRUE)
mycol=c('grey45', 'grey75', 'grey95')
barplot(t(m), las=1, ylab="% of variance", xlim=c(0,15),cex.lab=1.2, col=mycol, legend=TRUE, args.legend=list(x='bottomright', legend = rev(colnames(m)),bty='n',  fill=rev(mycol),y.intersp = 1.5,inset=c(-0.0001,0)))
mtext('Traits', side=1, line=3, adj=0.4, cex=1.2)
lines(c(-0.5, 11.3),c(0.5,0.5), lwd=1.4, lty=2)

# Body size traits plots
plot1 <- ggplot(total.sp, aes(x=Length_SL, y=Rel_Gut_Length, col=species))+
#geom_point(size=2,shape=21, col=1,)+
geom_smooth(method='lm',fullrange=TRUE,alpha=0.4) +
guides(fill=FALSE) +
theme_classic(base_size = 18)


#+
#scale_x_continuous("Relative gape")+
#scale_y_continuous("Condition factor (Kn)")+
#scale_colour_manual(name='Season',breaks=c("LD16","ED","DR","LD17"), labels=c("Late dry #16","Early dry", "Dry", "Late dry 17"),values=c('#A020F0','#CD4F39','#008B00','##4F94CD'))+
#scale_fill_manual(values=c('#A020F0','#CD4F39','#008B00','#4F94CD'))+
#scale_linetype_manual(values = 1:4)+
#theme(legend.position= 'top')+
#guides(color=guide_legend(override.aes=list(fill=NA)))

# Sex - traits plots

total.sp$Sex[is.na(total.sp$Sex)] <- 'IND'
total.spSex <- total.sp[total.sp$Sex=='M' | total.sp$Sex=='F',]

ggplot(total.spSex, aes(x=group, y=Rel_Gut_Length, col=species))+
#geom_point(size=2,shape=21, col=1,)+
geom_boxplot(fullrange=TRUE,alpha=0.4)+
theme_bw(base_size = 18)

## Boxplot
ggplot(total.sp, aes(x=group, y=EP, col=species))+
#geom_point(size=2,shape=21, col=1,)+
geom_boxplot(fullrange=TRUE,alpha=0.4)+
theme_bw(base_size = 18)





mod1 <- lmer(BD~1+(1| group), data= total.sp)

mod2 <- lmer(BD~species+(1| group), data= total.sp)



VarF <- var(as.vector(fixef(mod2) %*% t(getME(mod2,"X"))))


VarF/(VarF + VarCorr(mod2)$group[1] + attr(VarCorr(mod2), "sc")^2)

VarCorr(mod2)$Site_Code[1]/(VarF + VarCorr(mod2)$Site_Code[1] + attr(VarCorr(mod2), "sc")^2)

attr(VarCorr(mod2), "sc")^2/(VarF + VarCorr(mod2)$Site_Code[1] + attr(VarCorr(mod2), "sc")^2)




# Multi traits

# Tutorial about the colorful plot here: https://cran.r-project.org/web/packages/adegraphics/vignettes/adegraphics.html#customizing-a-graph

(traitpca <- dudi.pca(scale(traits.fish[,1:9]), scannf = FALSE, nf=3))
g1 <- s.corcircle(traitpca$co, plot = FALSE, grid=FALSE)
g2 <- s.label(traitpca$li, plot = FALSE)
ADEgS(list(g1, g2))

s.class(traitpca$li, species, chullSize=0,starSize=1,ellipseSize=0, ppoints.cex=1,xlim = c(-6, 7), ylim = c(-5, 5),col = TRUE)
s.class(traitpca$li, site, xlim = c(-4, 8), ylim = c(-4, 4),col = TRUE)
s.class(traitpca$li, habitat, xlim = c(-4, 8), ylim = c(-4, 4),col = TRUE)

fviz_pca_var(traitpca,axes=c(1,2),col.var="contrib", title="Traits")+
     scale_color_gradient2(low="white", mid="blue", 
           high="red", midpoint=4) + theme_minimal()

# Between-species Analysis

betsp <- bca(traitpca, species, scannf = FALSE)
dim(betsp$tab)

betsp$tab[1:3, 1:5]
mean(traitpca$tab$Temp[site == "COAR"])

betsp

plot(betsp, row.pellipses.col = adegpar()$ppalette$quali(1))

(rtbetsp <- randtest(betsp))
plot(rtbetsp)


### Between sites across species

m2 = matrix(, nrow = 15, ncol = 3)
colnames(m2) <- variables
rownames(m2) <- levels(species)
m2[,1] <- rep(betsp$ratio, 15)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$species=='AMNPER',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, traits.fish[traits.fish$species=='AMNPER',]$Site1, scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$species=='AMNPER',]$Site1, scannf = FALSE)
m2[1,2] <- betgroup$ratio * (1-betsp$ratio)
m2[1,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$species=='CRASTE',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$species=='CRASTE',]$Site1), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$species=='CRASTE',]$Site1, scannf = FALSE)
m2[2,2] <- betgroup$ratio * (1-betsp$ratio)
m2[2,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$species=='GLOAPR',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$species=='GLOAPR',]$Site1), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$species=='GLOAPR',]$Site1, scannf = FALSE)
m2[3,2] <- betgroup$ratio * (1-betsp$ratio)
m2[3,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$species=='HEPFUL',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$species=='HEPFUL',]$Site1), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$species=='HEPFUL',]$Site1, scannf = FALSE)
m2[4,2] <- betgroup$ratio * (1-betsp$ratio)
m2[4,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$species=='LATCAL',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$species=='LATCAL',]$Site1), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$species=='LATCAL',]$Site1, scannf = FALSE)
m2[5,2] <- betgroup$ratio * (1-betsp$ratio)
m2[5,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$species=='LEIUNI',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$species=='LEIUNI',]$Site1), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$species=='LEIUNI',]$Site1, scannf = FALSE)
m2[6,2] <- betgroup$ratio * (1-betsp$ratio)
m2[6,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$species=='MELAUS',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$species=='MELAUS',]$Site1), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$species=='MELAUS',]$Site1, scannf = FALSE)
m2[7,2] <- betgroup$ratio * (1-betsp$ratio)
m2[7,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$species=='MELSPL',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$species=='MELSPL',]$Site1), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$species=='MELSPL',]$Site1, scannf = FALSE)
m2[8,2] <- betgroup$ratio * (1-betsp$ratio)
m2[8,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$species=='MOGMOG',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$species=='MOGMOG',]$Site1), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$species=='MOGMOG',]$Site1, scannf = FALSE)
m2[9,2] <- betgroup$ratio * (1-betsp$ratio)
m2[9,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$species=='NEMERE',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$species=='NEMERE',]$Site1), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$species=='NEMERE',]$Site1, scannf = FALSE)
m2[10,2] <- betgroup$ratio * (1-betsp$ratio)
m2[10,3] <- witgroup$ratio * (1-betsp$ratio)


neoate.df <-as.data.frame(scale(traits.fish[traits.fish$species=='NEOATE',1:9]))
is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))
neoate.df$throt[is.nan(neoate.df$throt)] <- 0
(traitpca <- dudi.pca(neoate.df, scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$species=='NEOATE',]$Site1), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$species=='NEOATE',]$Site1, scannf = FALSE)
m2[11,2] <- betgroup$ratio * (1-betsp$ratio)
m2[11,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$species=='NEOGRA',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$species=='NEOGRA',]$Site1), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$species=='NEOGRA',]$Site1, scannf = FALSE)
m2[12,2] <- betgroup$ratio * (1-betsp$ratio)
m2[12,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$species=='OXYLIN',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$species=='OXYLIN',]$Site1), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$species=='OXYLIN',]$Site1, scannf = FALSE)
m2[13,2] <- betgroup$ratio * (1-betsp$ratio)
m2[13,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$species=='PLAORD',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$species=='PLAORD',]$Site1), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$species=='PLAORD',]$Site1, scannf = FALSE)
m2[14,2] <- betgroup$ratio * (1-betsp$ratio)
m2[14,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$species=='TOXCHA',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$species=='TOXCHA',]$Site1), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$species=='TOXCHA',]$Site1, scannf = FALSE)
m2[15,2] <- betgroup$ratio * (1-betsp$ratio)
m2[15,3] <- witgroup$ratio * (1-betsp$ratio)

# Plot Single species
par(xpd=TRUE)
mycol=c('grey75', 'grey95')

rotate_x <- function(data, labels_vec, rot_angle) {
     plt <- barplot(t(data), col=mycol, xaxt="n", las=2, ylab="% of variance")
     text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1.1,1.1), xpd =TRUE, cex=0.8) 
 }
 rotate_x(m2[,2:3], row.names(m2), 45)

m2 <- transform(m2, LH = c('Periodic', 'Opportunistic', 'Opportunistic', 'Periodic', 'Periodic', 'Opportunistic', 'Opportunistic', 'Opportunistic', 'Opportunistic', 'Periodic', 'Periodic', 'Equilibrium', 'Periodic', 'Periodic', 'Opportunistic'))

m2a <- cbind(m2, aus_dat4[c(3,4,6,7,9,10,13,14,16,17,18,19, 21,11,25) , 25:27])

m2$LH <- factor(m2$LH, levels = c("Opportunistic", "Periodic", "Equilibrium"))
boundaries <- boxplot(m2[,2] ~ m2$LH , col="skyblue2", ylim=c(0.03, 0.12), ylab="% Site variability")
# Now you can type boundaries$stats to get the boundaries of the boxes

# Add sample size on top
nbLH <- nlevels(m2$LH)
text( x=c(1:nbLH), y=boundaries$stats[nrow(boundaries$stats),] +0.003, paste("n = ",table(m2$LH),sep=""))

summary(lm(sqrt(m2[,2]) ~ m2$LH))

##### BY WEIGHT
plot(m2a[,2]~m2a$wgtOpp)
plot(m2a[,2]~m2a$wgtPer)
plot(m2a[,2]~m2a$wgtEqu)

mOp <- lm(m2a[,2]~m2a$wgtOpp)
mPe <- lm(m2a[,2]~m2a$wgtPer)
mEq <- lm(m2a[,2]~m2a$wgtEqu)

summary(mOp)
summary(mPe)
summary(mEq)

### Between species across sites


variables <- c('Interspecific var.', 'Intraspecific var./ site', 'Intraspecific var. unexplained')

m2 = matrix(, nrow = 14, ncol = 3)
colnames(m2) <- variables[c(2,1,3)]
rownames(m2) <- c('SHCR', 'BSBC', 'SHFR', 'GCFR', 'DRAR', 'COAR', 'BDFR', 'SRAR', 'MHMR', 'GPHR', 'GJKR','OODR', 'AHMR', 'DMAR')

levels(site)

betsp <- bca(traitpca, site, scannf = FALSE)
m2[,1] <- rep(betsp$ratio, 14)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$Site1=='SHCR',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$Site1=='SHCR',]$species), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$Site1=='SHCR',]$species, scannf = FALSE)
m2[1,2] <- betgroup$ratio * (1-betsp$ratio)
m2[1,3] <- witgroup$ratio * (1-betsp$ratio)
m2[1,2] <- 1 - (m2[1,1]+m2[1,3])


(traitpca <- dudi.pca(scale(traits.fish[traits.fish$Site1=='BSBC',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$Site1=='BSBC',]$species), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$Site1=='BSBC',]$species, scannf = FALSE)
m2[2,2] <- betgroup$ratio * (1-betsp$ratio)
m2[2,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$Site1=='SHFR',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$Site1=='SHFR',]$species), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$Site1=='SHFR',]$species, scannf = FALSE)
m2[3,2] <- betgroup$ratio * (1-betsp$ratio)
m2[3,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$Site1=='GCFR',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$Site1=='GCFR',]$species), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$Site1=='GCFR',]$species, scannf = FALSE)
m2[4,2] <- betgroup$ratio * (1-betsp$ratio)
m2[4,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$Site1=='DRAR',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$Site1=='DRAR',]$species), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$Site1=='DRAR',]$species, scannf = FALSE)
m2[5,2] <- betgroup$ratio * (1-betsp$ratio)
m2[5,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$Site1=='COAR',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$Site1=='COAR',]$species), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$Site1=='COAR',]$species, scannf = FALSE)
m2[6,2] <- betgroup$ratio * (1-betsp$ratio)
m2[6,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$Site1=='BDFR',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$Site1=='BDFR',]$species), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$Site1=='BDFR',]$species, scannf = FALSE)
m2[7,2] <- betgroup$ratio * (1-betsp$ratio)
m2[7,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$Site1=='SRAR',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$Site1=='SRAR',]$species), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$Site1=='SRAR',]$species, scannf = FALSE)
m2[8,2] <- betgroup$ratio * (1-betsp$ratio)
m2[8,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$Site1=='MHMR',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$Site1=='MHMR',]$species), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$Site1=='MHMR',]$species, scannf = FALSE)
m2[9,2] <- betgroup$ratio * (1-betsp$ratio)
m2[9,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$Site1=='GPHR',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$Site1=='GPHR',]$species), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$Site1=='GPHR',]$species, scannf = FALSE)
m2[10,2] <- betgroup$ratio * (1-betsp$ratio)
m2[10,3] <- witgroup$ratio * (1-betsp$ratio)


(traitpca <- dudi.pca(scale(traits.fish[traits.fish$Site1=='GJKR',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$Site1=='GJKR',]$species), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$Site1=='GJKR',]$species, scannf = FALSE)
m2[11,2] <- betgroup$ratio * (1-betsp$ratio)
m2[11,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$Site1=='OODR',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$Site1=='OODR',]$species), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$Site1=='OODR',]$species, scannf = FALSE)
m2[12,2] <- betgroup$ratio * (1-betsp$ratio)
m2[12,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$Site1=='AHMR',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$Site1=='AHMR',]$species), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$Site1=='AHMR',]$species, scannf = FALSE)
m2[13,2] <- betgroup$ratio * (1-betsp$ratio)
m2[13,3] <- witgroup$ratio * (1-betsp$ratio)

(traitpca <- dudi.pca(scale(traits.fish[traits.fish$Site1=='DMAR',1:9]), scannf = FALSE, nf=3))
betgroup <- bca(traitpca, droplevels(traits.fish[traits.fish$Site1=='DMAR',]$species), scannf = FALSE)
witgroup <- wca(traitpca, traits.fish[traits.fish$Site1=='DMAR',]$species, scannf = FALSE)
m2[14,2] <- betgroup$ratio * (1-betsp$ratio)
m2[14,3] <- witgroup$ratio * (1-betsp$ratio)

# Plot Single species
par(xpd=TRUE)
mycol=c('grey45', 'grey95')

rotate_x <- function(data, labels_vec, rot_angle) {
     plt <- barplot(t(data), col=mycol, xaxt="n", las=2, ylab="% of variance", ylim=c(0,1))
     text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1.1,1.1), xpd =TRUE, cex=0.8) 
 }
 rotate_x(m2[,2:3], row.names(m2), 45)

 comb <- unique(traits.fish[,c("Site1", "group")])

m2 <- transform(m2, Group = comb$group[match(rownames(m2), comb$Site1)])

boxplot (m2[,3]~m2$Group, col='grey', notch=T)

boundaries <- boxplot(m2[,3] ~ m2$Group , col="#69b3a2",  ylab="% Intraspecific variability", las=1)
# Now you can type boundaries$stats to get the boundaries of the boxes

m2$perc.intra <- m2[,3] / (m2[,3]+m2[,2])
m2 <- m2[order(rownames(m2)),]
m2$cat.area <- c(5811.56,371.56,146.93,84.44,1713.50,25.59,1044.68,9250.67,165.55,2838.60,36327.60,456.72,1685.69,712.05)
m2$strahler <- c(6,4,3,4,5,2,5,6,3,6,7,4,5,5)

plot(perc.intra*100~log(cat.area), data=m2)
plot(perc.intra*100~strahler, data=m2)

summary(lm(perc.intra*100~log(cat.area), data=m2))
summary(lm(perc.intra*100~strahler, data=m2))
summary(lm(perc.intra*100~Group, data=m2))

summary(lm(perc.intra*100~Group+log(cat.area)+strahler, data=m2))

summary(lm(log(cat.area)~Group, data=m2))
summary(lm(strahler~Group, data=m2))


# Add sample size on top
nbGroup <- nlevels(m2$Group)
text( x=c(1:nbGroup), y=boundaries$stats[nrow(boundaries$stats),] +0.005, paste("n = ",table(m2$Group),sep=""))

summary(lm(sqrt(m2[,3]) ~ m2$Group))
mod.river<- aov(sqrt(m2[,3]) ~ m2$Group)

TukeyHSD(mod.river)

# Between-group Analysis

betgroup <- bca(traitpca, habitat, scannf = FALSE)
dim(betgroup$tab)

betgroup$tab[1:3, 1:5]
mean(traitpca$tab$Temp[sites == "S1"])

betgroup

plot(betgroup, row.pellipses.col = adegpar()$ppalette$quali(6))

(rtbetgroup <- randtest(betgroup))
plot(rtbetgroup)

# Within-Class Analysis

witsp <- wca(traitpca, species, scannf = FALSE)
g1 <- s.corcircle(witsp$co, plot = FALSE)
g2 <- s.class(witsp$li, species, col = TRUE, plot = FALSE)
ADEgS(list(g1, g2))

plot(witsp, row.pellipses.col = adegpar()$ppalette$quali(6))


 sum(traitpca$eig)

sum(betsp$eig) + sum(witsp$eig)

sum(betsp$eig) / sum(traitpca$eig)

sum(witsp$eig) / sum(traitpca$eig)

s.class(witsp$li, group, col = TRUE)

####  Each Species PCA

(sppca <- dudi.pca(traits.fish[traits.fish$species=="TOXCHA",2:10], scannf = FALSE, nf=3))
g1 <- s.corcircle(sppca$co, plot = FALSE)
g1
g2 <- s.label(sppca$li, plot = FALSE)
ADEgS(list(g1, g2))



#### Average PCA

sp.only <- with(traits.fish, aggregate(traits.fish[,2:10], list(species), mean))
rownames(sp.only) <- sp.only[,1]


sp_group <- with(traits.fish, aggregate(traits.fish[,2:10], list(species, group), mean))
colnames(sp_group)[1]<- 'species'
colnames(sp_group)[2]<- 'habitat'
habitat2 <- sp_group$habitat

group.only <- with(traits.fish, aggregate(traits.fish[,1:11], list(group), mean))
rownames(group.only) <- group.only[,1]

sp_site <- with(traits.fish, aggregate(traits.fish[,1:11], list(species, Site1), mean))
colnames(sp_site)[1]<- 'species'
colnames(sp_site)[2]<- 'site'
species2 <- sp_site$species 


(sppca <- dudi.pca(scale(sp.only[,2:10]), scannf = FALSE, nf=3))
g1 <- s.corcircle(sppca$co, plot = FALSE)
g2 <- s.label(sppca$li, plot = FALSE)
ADEgS(list(g1, g2))



(sphabpca <- dudi.pca(scale(sp_site[,3:13]), scannf = FALSE, nf=3))
g1 <- s.corcircle(sphabpca$co, plot = FALSE)
g2 <- s.label(sphabpca$li, plot = FALSE)
ADEgS(list(g1, g2))

betsite <- bca(sphabpca, species2, scannf = FALSE)
dim(betsite$tab)

betsite$tab[1:3, 1:5]
mean(traitpca$tab$Temp[sites == "S1"])

betsite

plot(betsite, row.pellipses.col = adegpar()$ppalette$quali(6))

(rtbetgroup <- randtest(betgroup))
plot(rtbetgroup)


s.class(sphabpca$li, species2, xlim = c(-6, 7), ylim = c(-5, 5),col = TRUE)

(traitpca <- dudi.pca(scale(sp_site[,4:12]), scannf = FALSE, nf=3))
species <- sp_site$species
s.class(traitpca$li, species, xlim = c(-6, 7), ylim = c(-5, 5),col = TRUE)

  #### PLOTS CATI PACKAGE

library(cati)
library("mice")
library("hypervolume")
library("e1071")


  par(mfrow = c(4,4), cex = 0.5)
plotDistri(traits.fish[,1:10], species, habitat,
           ylim.cex = 3, plot.ask = F, multipanel = T, leg = F)

  par(mfrow = c(4,4), cex = 0.5)
plotDistri(traits.fish[,1:10], habitat, species, 
           ylim.cex = 3, plot.ask = F, multipanel = T, leg = F)

par(mfrow = c(1,3))
plotDistri(as.matrix(traits.fish[,"BD"]), habitat, species,
          ylim.cex = 8, plot.ask = F, multipanel = F, leg = T, cex.leg = 0.5)

# T-stats

res.fish <- Tstats(traits.fish[,1:10], ind.plot = habitat, sp = species,
         nperm = 99, print = FALSE)


res.fish

barplot(res.fish, ylim = c(0,3.5))

plot(res.fish)


 levelplot(t(ses(res.fish$Tstats$T_IP.IC, res.fish$Tstats$T_IP.IC_nm)$ses),border = "black")





#####  GLMM (Albert 2010)  ###  for sites

traits <- c('TLgth','BD','PFL','Throt','ED','EP','GS','GSh','MP','BS')
variables <- c('SHCR','BSBC','SHFR','GPHR','SRAR','MHMR','GCFR','BDFR','DRAR','AHMR','COAR','DMAR','GJKR','OODR')
m = matrix(, nrow = 10, ncol = 3)


mod1 <- lmer(log(Length_SL)~species  +(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod1) %*% t(getME(mod1,"X"))))
m[1,1] <- VarF/(VarF + VarCorr(mod1)$Site1[1] + attr(VarCorr(mod1), "sc")^2)
m[1,2] <- VarCorr(mod1)$Site1[1]/(VarF + VarCorr(mod1)$Site1[1] + attr(VarCorr(mod1), "sc")^2)
m[1,3] <- attr(VarCorr(mod1), "sc")^2/(VarF + VarCorr(mod1)$Site1[1] + attr(VarCorr(mod1), "sc")^2)

mod2 <- lmer(BD~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod2) %*% t(getME(mod2,"X"))))
m[2,1] <- VarF/(VarF + VarCorr(mod2)$Site1[1] + attr(VarCorr(mod2), "sc")^2)
m[2,2] <- VarCorr(mod2)$Site1[1]/(VarF + VarCorr(mod2)$Site1[1] + attr(VarCorr(mod2), "sc")^2)
m[2,3] <- attr(VarCorr(mod2), "sc")^2/(VarF + VarCorr(mod2)$Site1[1] + attr(VarCorr(mod2), "sc")^2)

mod3 <- lmer(PFL~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod3) %*% t(getME(mod3,"X"))))
m[3,1] <- VarF/(VarF + VarCorr(mod3)$Site1[1] + attr(VarCorr(mod3), "sc")^2)
m[3,2] <- VarCorr(mod3)$Site1[1]/(VarF + VarCorr(mod3)$Site1[1] + attr(VarCorr(mod3), "sc")^2)
m[3,3] <- attr(VarCorr(mod3), "sc")^2/(VarF + VarCorr(mod3)$Site1[1] + attr(VarCorr(mod3), "sc")^2)

mod4 <- lmer(throt~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod4) %*% t(getME(mod4,"X"))))
m[4,1] <- VarF/(VarF + VarCorr(mod4)$Site1[1] + attr(VarCorr(mod4), "sc")^2)
m[4,2] <- VarCorr(mod4)$Site1[1]/(VarF + VarCorr(mod4)$Site1[1] + attr(VarCorr(mod4), "sc")^2)
m[4,3] <- attr(VarCorr(mod4), "sc")^2/(VarF + VarCorr(mod4)$Site1[1] + attr(VarCorr(mod4), "sc")^2)

mod5 <- lmer(ED~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod5) %*% t(getME(mod5,"X"))))
m[5,1] <- VarF/(VarF + VarCorr(mod5)$Site1[1] + attr(VarCorr(mod5), "sc")^2)
m[5,2] <- VarCorr(mod5)$Site1[1]/(VarF + VarCorr(mod5)$Site1[1] + attr(VarCorr(mod5), "sc")^2)
m[5,3] <- attr(VarCorr(mod5), "sc")^2/(VarF + VarCorr(mod5)$Site1[1] + attr(VarCorr(mod5), "sc")^2)


mod6 <- lmer(EP~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod6) %*% t(getME(mod6,"X"))))
m[6,1] <- VarF/(VarF + VarCorr(mod6)$Site1[1] + attr(VarCorr(mod6), "sc")^2)
m[6,2] <- VarCorr(mod6)$Site1[1]/(VarF + VarCorr(mod6)$Site1[1] + attr(VarCorr(mod6), "sc")^2)
m[6,3] <- attr(VarCorr(mod6), "sc")^2/(VarF + VarCorr(mod6)$Site1[1] + attr(VarCorr(mod6), "sc")^2)

mod7 <- lmer(log(GS)~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod7) %*% t(getME(mod7,"X"))))
m[7,1] <- VarF/(VarF + VarCorr(mod7)$Site1[1] + attr(VarCorr(mod7), "sc")^2)
m[7,2] <- VarCorr(mod7)$Site1[1]/(VarF + VarCorr(mod7)$Site1[1] + attr(VarCorr(mod7), "sc")^2)
m[7,3] <- attr(VarCorr(mod7), "sc")^2/(VarF + VarCorr(mod7)$Site1[1] + attr(VarCorr(mod7), "sc")^2)

mod8 <- lmer(log(GSh)~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod8) %*% t(getME(mod8,"X"))))
m[8,1] <- VarF/(VarF + VarCorr(mod8)$Site1[1] + attr(VarCorr(mod8), "sc")^2)
m[8,2] <- VarCorr(mod8)$Site1[1]/(VarF + VarCorr(mod8)$Site1[1] + attr(VarCorr(mod8), "sc")^2)
m[8,3] <- attr(VarCorr(mod8), "sc")^2/(VarF + VarCorr(mod8)$Site1[1] + attr(VarCorr(mod8), "sc")^2)

mod9 <- lmer(MP~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod9) %*% t(getME(mod9,"X"))))
m[9,1] <- VarF/(VarF + VarCorr(mod9)$Site1[1] + attr(VarCorr(mod9), "sc")^2)
m[9,2] <- VarCorr(mod9)$Site1[1]/(VarF + VarCorr(mod9)$Site1[1] + attr(VarCorr(mod9), "sc")^2)
m[9,3] <- attr(VarCorr(mod9), "sc")^2/(VarF + VarCorr(mod9)$Site1[1] + attr(VarCorr(mod9), "sc")^2)

mod10 <- lmer(BS~species+(1| Site1), data= total.sp)
VarF <- var(as.vector(fixef(mod10) %*% t(getME(mod10,"X"))))
m[10,1] <- VarF/(VarF + VarCorr(mod10)$Site1[1] + attr(VarCorr(mod10), "sc")^2)
m[10,2] <- VarCorr(mod10)$Site1[1]/(VarF + VarCorr(mod10)$Site1[1] + attr(VarCorr(mod10), "sc")^2)
m[10,3] <- attr(VarCorr(mod10), "sc")^2/(VarF + VarCorr(mod10)$Site1[1] + attr(VarCorr(mod10), "sc")^2)

colnames(m) <- variables
rownames(m) <- traits

##### Mahalanobis distance

coord <- traitpca$li
maha.mat <- cbind(traits.fish$species, coord)
names(maha.mat)[1] <- "species"




par(mfrow=c(3,3))


head(maha.mat)
amn.per <- maha.mat[maha.mat$species=='AMNPER', 2:4]
cra.ste <- maha.mat[maha.mat$species=='CRASTE', 2:4]
glo.apr <- maha.mat[maha.mat$species=='GLOAPR', 2:4]
hep.ful <- maha.mat[maha.mat$species=='HEPFUL', 2:4]
lat.cal <- maha.mat[maha.mat$species=='LATCAL', 2:4]
lei.uni <- maha.mat[maha.mat$species=='LEIUNI', 2:4]
mel.aus <- maha.mat[maha.mat$species=='MELAUS', 2:4]
mel.spl <- maha.mat[maha.mat$species=='MELSPL', 2:4]
mog.mog <- maha.mat[maha.mat$species=='MOGMOG', 2:4]
nem.ere <- maha.mat[maha.mat$species=='NEMERE', 2:4]
neo.ate <- maha.mat[maha.mat$species=='NEOATE', 2:4]
neo.gra <- maha.mat[maha.mat$species=='NEOGRA', 2:4]
oxy.lin <- maha.mat[maha.mat$species=='OXYLIN', 2:4]
pla.ord <- maha.mat[maha.mat$species=='PLAORD', 2:4]
tox.cha <- maha.mat[maha.mat$species=='TOXCHA', 2:4]

library(StatMatch)

#amnper
amnamn <- mahalanobis.dist(data.x=amn.per, data.y=amn.per)
amncra <- mahalanobis.dist(data.x=amn.per, data.y=cra.ste)
amnglo <- mahalanobis.dist(data.x=amn.per, data.y=glo.apr)
amnhep <- mahalanobis.dist(data.x=amn.per, data.y=hep.ful)
amnlat <- mahalanobis.dist(data.x=amn.per, data.y=lat.cal)
amnlei <- mahalanobis.dist(data.x=amn.per, data.y=lei.uni)
amnaus <- mahalanobis.dist(data.x=amn.per, data.y=mel.aus)
amnspl <- mahalanobis.dist(data.x=amn.per, data.y=mel.spl)
amnmog <- mahalanobis.dist(data.x=amn.per, data.y=mog.mog)
amnnem <- mahalanobis.dist(data.x=amn.per, data.y=nem.ere)
amnate <- mahalanobis.dist(data.x=amn.per, data.y=neo.ate)
amngra <- mahalanobis.dist(data.x=amn.per, data.y=neo.gra)
amnlin <- mahalanobis.dist(data.x=amn.per, data.y=oxy.lin)
amnpla <- mahalanobis.dist(data.x=amn.per, data.y=pla.ord)
amntox <- mahalanobis.dist(data.x=amn.per, data.y=tox.cha)

plot (density(amnamn), col=1,  las=1, xlim=c(0,6), ylim=c(0,0.8), xlab='', main='', ylab="Kernel density estimate", lwd=1.5, lty=4)
lines(density(amncra), col=2)
lines(density(amnglo), col=3)
lines(density(amnhep), col=4)
lines(density(amnlat), col=5)
lines(density(amnlei), col=6)
lines(density(amnaus), col=7)
lines(density(amnspl), col=8)
lines(density(amnmog), col=9)
lines(density(amnnem), col=10)
lines(density(amnate), col=11)
lines(density(amngra), col=12)
lines(density(amnlin), col=13)
lines(density(amnpla), col=14)
lines(density(amntox), col=15)
mtext(expression(paste("A) ",italic("Amniataba percoides"))),line=0.2, adj=0, padj=-0.1 )
mtext("Mahalanobis distance between individuals",1 , line=2)

#craste
craamn <- mahalanobis.dist(data.x=cra.ste, data.y=amn.per)
cracra <- mahalanobis.dist(data.x=cra.ste, data.y=cra.ste)
craglo <- mahalanobis.dist(data.x=cra.ste, data.y=glo.apr)
crahep <- mahalanobis.dist(data.x=cra.ste, data.y=hep.ful)
cralat <- mahalanobis.dist(data.x=cra.ste, data.y=lat.cal)
cralei <- mahalanobis.dist(data.x=cra.ste, data.y=lei.uni)
craaus <- mahalanobis.dist(data.x=cra.ste, data.y=mel.aus)
craspl <- mahalanobis.dist(data.x=cra.ste, data.y=mel.spl)
cramog <- mahalanobis.dist(data.x=cra.ste, data.y=mog.mog)
cranem <- mahalanobis.dist(data.x=cra.ste, data.y=nem.ere)
craate <- mahalanobis.dist(data.x=cra.ste, data.y=neo.ate)
cragra <- mahalanobis.dist(data.x=cra.ste, data.y=neo.gra)
cralin <- mahalanobis.dist(data.x=cra.ste, data.y=oxy.lin)
crapla <- mahalanobis.dist(data.x=cra.ste, data.y=pla.ord)
cratox <- mahalanobis.dist(data.x=cra.ste, data.y=tox.cha)

plot (density(craamn), col=1,  las=1, xlim=c(0,6), ylim=c(0,0.8), xlab='', main='', ylab="Kernel density estimate")
lines(density(cracra), col=2, lwd=1.5, lty=4)
lines(density(craglo), col=3)
lines(density(crahep), col=4)
lines(density(cralat), col=5)
lines(density(cralei), col=6)
lines(density(craaus), col=7)
lines(density(craspl), col=8)
lines(density(cramog), col=9)
lines(density(cranem), col=10)
lines(density(craate), col=11)
lines(density(cragra), col=12)
lines(density(cralin), col=13)
lines(density(crapla), col=14)
lines(density(cratox), col=15)
mtext(expression(paste("B) ",italic("Craterocephalus stercusmuscarum"))),line=0.2, adj=0, padj=-0.1 )
mtext("Mahalanobis distance between individuals",1 , line=2)

#gloapr
gloamn <- mahalanobis.dist(data.x=glo.apr, data.y=amn.per)
glocra <- mahalanobis.dist(data.x=glo.apr, data.y=cra.ste)
gloglo <- mahalanobis.dist(data.x=glo.apr, data.y=glo.apr)
glohep <- mahalanobis.dist(data.x=glo.apr, data.y=hep.ful)
glolat <- mahalanobis.dist(data.x=glo.apr, data.y=lat.cal)
glolei <- mahalanobis.dist(data.x=glo.apr, data.y=lei.uni)
gloaus <- mahalanobis.dist(data.x=glo.apr, data.y=mel.aus)
glospl <- mahalanobis.dist(data.x=glo.apr, data.y=mel.spl)
glomog <- mahalanobis.dist(data.x=glo.apr, data.y=mog.mog)
glonem <- mahalanobis.dist(data.x=glo.apr, data.y=nem.ere)
gloate <- mahalanobis.dist(data.x=glo.apr, data.y=neo.ate)
glogra <- mahalanobis.dist(data.x=glo.apr, data.y=neo.gra)
glolin <- mahalanobis.dist(data.x=glo.apr, data.y=oxy.lin)
glopla <- mahalanobis.dist(data.x=glo.apr, data.y=pla.ord)
glotox <- mahalanobis.dist(data.x=glo.apr, data.y=tox.cha)

plot (density(gloamn), col=1,  las=1, xlim=c(0,6), ylim=c(0,0.8), xlab='', main='', ylab="Kernel density estimate")
lines(density(glocra), col=2)
lines(density(gloglo), col=3, lwd=1.5, lty=4)
lines(density(glohep), col=4)
lines(density(glolat), col=5)
lines(density(glolei), col=6)
lines(density(gloaus), col=7)
lines(density(glospl), col=8)
lines(density(glomog), col=9)
lines(density(glonem), col=10)
lines(density(gloate), col=11)
lines(density(glogra), col=12)
lines(density(glolin), col=13)
lines(density(glopla), col=14)
lines(density(glotox), col=15)
mtext(expression(paste("C) ",italic("Glossamia aprion"))),line=0.2, adj=0, padj=-0.1 )
mtext("Mahalanobis distance between individuals",1 , line=2)

#hepful
hepamn <- mahalanobis.dist(data.x=hep.ful, data.y=amn.per)
hepcra <- mahalanobis.dist(data.x=hep.ful, data.y=cra.ste)
hepglo <- mahalanobis.dist(data.x=hep.ful, data.y=glo.apr)
hephep <- mahalanobis.dist(data.x=hep.ful, data.y=hep.ful)
heplat <- mahalanobis.dist(data.x=hep.ful, data.y=lat.cal)
heplei <- mahalanobis.dist(data.x=hep.ful, data.y=lei.uni)
hepaus <- mahalanobis.dist(data.x=hep.ful, data.y=mel.aus)
hepspl <- mahalanobis.dist(data.x=hep.ful, data.y=mel.spl)
hepmog <- mahalanobis.dist(data.x=hep.ful, data.y=mog.mog)
hepnem <- mahalanobis.dist(data.x=hep.ful, data.y=nem.ere)
hepate <- mahalanobis.dist(data.x=hep.ful, data.y=neo.ate)
hepgra <- mahalanobis.dist(data.x=hep.ful, data.y=neo.gra)
heplin <- mahalanobis.dist(data.x=hep.ful, data.y=oxy.lin)
heppla <- mahalanobis.dist(data.x=hep.ful, data.y=pla.ord)
heptox <- mahalanobis.dist(data.x=hep.ful, data.y=tox.cha)

plot (density(hepamn), col=1,  las=1, xlim=c(0,6), ylim=c(0,0.8), xlab='', main='', ylab="Kernel density estimate")
lines(density(hepcra), col=2)
lines(density(hepglo), col=3)
lines(density(hephep), col=4, lwd=1.5, lty=4)
lines(density(heplat), col=5)
lines(density(heplei), col=6)
lines(density(hepaus), col=7)
lines(density(hepspl), col=8)
lines(density(hepmog), col=9)
lines(density(hepnem), col=10)
lines(density(hepate), col=11)
lines(density(hepgra), col=12)
lines(density(heplin), col=13)
lines(density(heppla), col=14)
lines(density(heptox), col=15)
mtext(expression(paste("D) ",italic("Hephaestus fuliginosus"))),line=0.2, adj=0, padj=-0.1 )
mtext("Mahalanobis distance between individuals",1 , line=2)

#latcal
latamn <- mahalanobis.dist(data.x=lat.cal, data.y=amn.per)
latcra <- mahalanobis.dist(data.x=lat.cal, data.y=cra.ste)
latglo <- mahalanobis.dist(data.x=lat.cal, data.y=glo.apr)
lathep <- mahalanobis.dist(data.x=lat.cal, data.y=hep.ful)
latlat <- mahalanobis.dist(data.x=lat.cal, data.y=lat.cal)
latlei <- mahalanobis.dist(data.x=lat.cal, data.y=lei.uni)
lataus <- mahalanobis.dist(data.x=lat.cal, data.y=mel.aus)
latspl <- mahalanobis.dist(data.x=lat.cal, data.y=mel.spl)
latmog <- mahalanobis.dist(data.x=lat.cal, data.y=mog.mog)
latnem <- mahalanobis.dist(data.x=lat.cal, data.y=nem.ere)
latate <- mahalanobis.dist(data.x=lat.cal, data.y=neo.ate)
latgra <- mahalanobis.dist(data.x=lat.cal, data.y=neo.gra)
latlin <- mahalanobis.dist(data.x=lat.cal, data.y=oxy.lin)
latpla <- mahalanobis.dist(data.x=lat.cal, data.y=pla.ord)
lattox <- mahalanobis.dist(data.x=lat.cal, data.y=tox.cha)

plot (density(latamn), col=1,  las=1, xlim=c(0,6), ylim=c(0,0.8), xlab='', main='', ylab="Kernel density estimate")
lines(density(latcra), col=2)
lines(density(latglo), col=3)
lines(density(lathep), col=4)
lines(density(latlat), col=5, lwd=1.5, lty=4)
lines(density(latlei), col=6)
lines(density(lataus), col=7)
lines(density(latspl), col=8)
lines(density(latmog), col=9)
lines(density(latnem), col=10)
lines(density(latate), col=11)
lines(density(latgra), col=12)
lines(density(latlin), col=13)
lines(density(latpla), col=14)
lines(density(lattox), col=15)
mtext(expression(paste("E) ",italic("Lates calcarifer"))),line=0.2, adj=0, padj=-0.1 )
mtext("Mahalanobis distance between individuals",1 , line=2)

#leiuni
leiamn <- mahalanobis.dist(data.x=lei.uni, data.y=amn.per)
leicra <- mahalanobis.dist(data.x=lei.uni, data.y=cra.ste)
leiglo <- mahalanobis.dist(data.x=lei.uni, data.y=glo.apr)
leihep <- mahalanobis.dist(data.x=lei.uni, data.y=hep.ful)
leilat <- mahalanobis.dist(data.x=lei.uni, data.y=lat.cal)
leilei <- mahalanobis.dist(data.x=lei.uni, data.y=lei.uni)
leiaus <- mahalanobis.dist(data.x=lei.uni, data.y=mel.aus)
leispl <- mahalanobis.dist(data.x=lei.uni, data.y=mel.spl)
leimog <- mahalanobis.dist(data.x=lei.uni, data.y=mog.mog)
leinem <- mahalanobis.dist(data.x=lei.uni, data.y=nem.ere)
leiate <- mahalanobis.dist(data.x=lei.uni, data.y=neo.ate)
leigra <- mahalanobis.dist(data.x=lei.uni, data.y=neo.gra)
leilin <- mahalanobis.dist(data.x=lei.uni, data.y=oxy.lin)
leipla <- mahalanobis.dist(data.x=lei.uni, data.y=pla.ord)
leitox <- mahalanobis.dist(data.x=lei.uni, data.y=tox.cha)

plot (density(leiamn), col=1,  las=1, xlim=c(0,6), ylim=c(0,0.8), xlab='', main='', ylab="Kernel density estimate")
lines(density(leicra), col=2)
lines(density(leiglo), col=3)
lines(density(leihep), col=4)
lines(density(leilat), col=5)
lines(density(leilei), col=6, lwd=1.5, lty=4)
lines(density(leiaus), col=7)
lines(density(leispl), col=8)
lines(density(leimog), col=9)
lines(density(leinem), col=10)
lines(density(leiate), col=11)
lines(density(leigra), col=12)
lines(density(leilin), col=13)
lines(density(leipla), col=14)
lines(density(leitox), col=15)
mtext(expression(paste("E) ",italic("Leiopotherapon unicolor"))),line=0.2, adj=0, padj=-0.1 )
mtext("Mahalanobis distance between individuals",1 , line=2)

#melaus
ausamn <- mahalanobis.dist(data.x=mel.aus, data.y=amn.per)
auscra <- mahalanobis.dist(data.x=mel.aus, data.y=cra.ste)
ausglo <- mahalanobis.dist(data.x=mel.aus, data.y=glo.apr)
aushep <- mahalanobis.dist(data.x=mel.aus, data.y=hep.ful)
auslat <- mahalanobis.dist(data.x=mel.aus, data.y=lat.cal)
auslei <- mahalanobis.dist(data.x=mel.aus, data.y=lei.uni)
ausaus <- mahalanobis.dist(data.x=mel.aus, data.y=mel.aus)
ausspl <- mahalanobis.dist(data.x=mel.aus, data.y=mel.spl)
ausmog <- mahalanobis.dist(data.x=mel.aus, data.y=mog.mog)
ausnem <- mahalanobis.dist(data.x=mel.aus, data.y=nem.ere)
ausate <- mahalanobis.dist(data.x=mel.aus, data.y=neo.ate)
ausgra <- mahalanobis.dist(data.x=mel.aus, data.y=neo.gra)
auslin <- mahalanobis.dist(data.x=mel.aus, data.y=oxy.lin)
auspla <- mahalanobis.dist(data.x=mel.aus, data.y=pla.ord)
austox <- mahalanobis.dist(data.x=mel.aus, data.y=tox.cha)

plot (density(ausamn), col=1,  las=1, xlim=c(0,6), ylim=c(0,0.8), xlab='', main='', ylab="Kernel density estimate")
lines(density(auscra), col=2)
lines(density(ausglo), col=3)
lines(density(aushep), col=4)
lines(density(auslat), col=5)
lines(density(auslei), col=6)
lines(density(ausaus), col=7, lwd=1.5, lty=4)
lines(density(ausspl), col=8)
lines(density(ausmog), col=9)
lines(density(ausnem), col=10)
lines(density(ausate), col=11)
lines(density(ausgra), col=12)
lines(density(auslin), col=13)
lines(density(auspla), col=14)
lines(density(austox), col=15)
mtext(expression(paste("E) ",italic("Melanotaenia australis"))),line=0.2, adj=0, padj=-0.1 )
mtext("Mahalanobis distance between individuals",1 , line=2)

#melspl
splamn <- mahalanobis.dist(data.x=mel.spl, data.y=amn.per)
splcra <- mahalanobis.dist(data.x=mel.spl, data.y=cra.ste)
splglo <- mahalanobis.dist(data.x=mel.spl, data.y=glo.apr)
splhep <- mahalanobis.dist(data.x=mel.spl, data.y=hep.ful)
spllat <- mahalanobis.dist(data.x=mel.spl, data.y=lat.cal)
spllei <- mahalanobis.dist(data.x=mel.spl, data.y=lei.uni)
splaus <- mahalanobis.dist(data.x=mel.spl, data.y=mel.aus)
splspl <- mahalanobis.dist(data.x=mel.spl, data.y=mel.spl)
splmog <- mahalanobis.dist(data.x=mel.spl, data.y=mog.mog)
splnem <- mahalanobis.dist(data.x=mel.spl, data.y=nem.ere)
splate <- mahalanobis.dist(data.x=mel.spl, data.y=neo.ate)
splgra <- mahalanobis.dist(data.x=mel.spl, data.y=neo.gra)
spllin <- mahalanobis.dist(data.x=mel.spl, data.y=oxy.lin)
splpla <- mahalanobis.dist(data.x=mel.spl, data.y=pla.ord)
spltox <- mahalanobis.dist(data.x=mel.spl, data.y=tox.cha)

plot (density(splamn), col=1,  las=1, xlim=c(0,6), ylim=c(0,0.8), xlab='', main='', ylab="Kernel density estimate")
lines(density(splcra), col=2)
lines(density(splglo), col=3)
lines(density(splhep), col=4)
lines(density(spllat), col=5)
lines(density(spllei), col=6)
lines(density(splaus), col=7)
lines(density(splspl), col=8, lwd=1.5, lty=4)
lines(density(splmog), col=9)
lines(density(splnem), col=10)
lines(density(splate), col=11)
lines(density(splgra), col=12)
lines(density(spllin), col=13)
lines(density(splpla), col=14)
lines(density(spltox), col=15)
mtext(expression(paste("E) ",italic("Melanotaenia splendida"))),line=0.2, adj=0, padj=-0.1 )
mtext("Mahalanobis distance between individuals",1 , line=2)

#mogmog
mogamn <- mahalanobis.dist(data.x=mog.mog, data.y=amn.per)
mogcra <- mahalanobis.dist(data.x=mog.mog, data.y=cra.ste)
mogglo <- mahalanobis.dist(data.x=mog.mog, data.y=glo.apr)
moghep <- mahalanobis.dist(data.x=mog.mog, data.y=hep.ful)
moglat <- mahalanobis.dist(data.x=mog.mog, data.y=lat.cal)
moglei <- mahalanobis.dist(data.x=mog.mog, data.y=lei.uni)
mogaus <- mahalanobis.dist(data.x=mog.mog, data.y=mel.aus)
mogspl <- mahalanobis.dist(data.x=mog.mog, data.y=mel.spl)
mogmog <- mahalanobis.dist(data.x=mog.mog, data.y=mog.mog)
mognem <- mahalanobis.dist(data.x=mog.mog, data.y=nem.ere)
mogate <- mahalanobis.dist(data.x=mog.mog, data.y=neo.ate)
moggra <- mahalanobis.dist(data.x=mog.mog, data.y=neo.gra)
moglin <- mahalanobis.dist(data.x=mog.mog, data.y=oxy.lin)
mogpla <- mahalanobis.dist(data.x=mog.mog, data.y=pla.ord)
mogtox <- mahalanobis.dist(data.x=mog.mog, data.y=tox.cha)

plot (density(mogamn), col=1,  las=1, xlim=c(0,6), ylim=c(0,0.8), xlab='', main='', ylab="Kernel density estimate")
lines(density(mogcra), col=2)
lines(density(mogglo), col=3)
lines(density(moghep), col=4)
lines(density(moglat), col=5)
lines(density(moglei), col=6)
lines(density(mogaus), col=7)
lines(density(mogspl), col=8)
lines(density(mogmog), col=9, lwd=1.5, lty=4)
lines(density(mognem), col=10)
lines(density(mogate), col=11)
lines(density(moggra), col=12)
lines(density(moglin), col=13)
lines(density(mogpla), col=14)
lines(density(mogtox), col=15)
mtext(expression(paste("E) ",italic("Mogurnda Mogurnda"))),line=0.2, adj=0, padj=-0.1 )
mtext("Mahalanobis distance between individuals",1 , line=2)

#nemere
nemamn <- mahalanobis.dist(data.x=nem.ere, data.y=amn.per)
nemcra <- mahalanobis.dist(data.x=nem.ere, data.y=cra.ste)
nemglo <- mahalanobis.dist(data.x=nem.ere, data.y=glo.apr)
nemhep <- mahalanobis.dist(data.x=nem.ere, data.y=hep.ful)
nemlat <- mahalanobis.dist(data.x=nem.ere, data.y=lat.cal)
nemlei <- mahalanobis.dist(data.x=nem.ere, data.y=lei.uni)
nemaus <- mahalanobis.dist(data.x=nem.ere, data.y=mel.aus)
nemspl <- mahalanobis.dist(data.x=nem.ere, data.y=mel.spl)
nemmog <- mahalanobis.dist(data.x=nem.ere, data.y=mog.mog)
nemnem <- mahalanobis.dist(data.x=nem.ere, data.y=nem.ere)
nemate <- mahalanobis.dist(data.x=nem.ere, data.y=neo.ate)
nemgra <- mahalanobis.dist(data.x=nem.ere, data.y=neo.gra)
nemlin <- mahalanobis.dist(data.x=nem.ere, data.y=oxy.lin)
nempla <- mahalanobis.dist(data.x=nem.ere, data.y=pla.ord)
nemtox <- mahalanobis.dist(data.x=nem.ere, data.y=tox.cha)
plot (density(nemamn), col=1,  las=1, xlim=c(0,6), ylim=c(0,0.8), xlab='', main='', ylab="Kernel density estimate")
lines(density(nemcra), col=2)
lines(density(nemglo), col=3)
lines(density(nemhep), col=4)
lines(density(nemlat), col=5)
lines(density(nemlei), col=6)
lines(density(nemaus), col=7)
lines(density(nemspl), col=8)
lines(density(nemmog), col=9)
lines(density(nemnem), col=10, lwd=1.5, lty=4)
lines(density(nemate), col=11)
lines(density(nemgra), col=12)
lines(density(nemlin), col=13)
lines(density(nempla), col=14)
lines(density(nemtox), col=15)
mtext(expression(paste("E) ",italic("Nematolosa erebi"))),line=0.2, adj=0, padj=-0.1 )
mtext("Mahalanobis distance between individuals",1 , line=2)

#neogra
graamn <- mahalanobis.dist(data.x=neo.gra, data.y=amn.per)
gracra <- mahalanobis.dist(data.x=neo.gra, data.y=cra.ste)
graglo <- mahalanobis.dist(data.x=neo.gra, data.y=glo.apr)
grahep <- mahalanobis.dist(data.x=neo.gra, data.y=hep.ful)
gralat <- mahalanobis.dist(data.x=neo.gra, data.y=lat.cal)
gralei <- mahalanobis.dist(data.x=neo.gra, data.y=lei.uni)
graaus <- mahalanobis.dist(data.x=neo.gra, data.y=mel.aus)
graspl <- mahalanobis.dist(data.x=neo.gra, data.y=mel.spl)
gramog <- mahalanobis.dist(data.x=neo.gra, data.y=mog.mog)
granem <- mahalanobis.dist(data.x=neo.gra, data.y=nem.ere)
graate <- mahalanobis.dist(data.x=neo.gra, data.y=neo.ate)
gragra <- mahalanobis.dist(data.x=neo.gra, data.y=neo.gra)
gralin <- mahalanobis.dist(data.x=neo.gra, data.y=oxy.lin)
grapla <- mahalanobis.dist(data.x=neo.gra, data.y=pla.ord)
gratox <- mahalanobis.dist(data.x=neo.gra, data.y=tox.cha)
plot (density(graamn), col=1,  las=1, xlim=c(0,6), ylim=c(0,0.8), xlab='', main='', ylab="Kernel density estimate")
lines(density(gracra), col=2)
lines(density(graglo), col=3)
lines(density(grahep), col=4)
lines(density(gralat), col=5)
lines(density(gralei), col=6)
lines(density(graaus), col=7)
lines(density(graspl), col=8)
lines(density(gramog), col=9)
lines(density(granem), col=10)
lines(density(graate), col=11)
lines(density(gragra), col=12, lwd=1.5, lty=4)
lines(density(gralin), col=13)
lines(density(grapla), col=14)
lines(density(gratox), col=15)
mtext(expression(paste("E) ",italic("Neoarius graeffei"))),line=0.2, adj=0, padj=-0.1 )
mtext("Mahalanobis distance between individuals",1 , line=2)

#oxylin
oxyamn <- mahalanobis.dist(data.x=oxy.lin, data.y=amn.per)
oxycra <- mahalanobis.dist(data.x=oxy.lin, data.y=cra.ste)
oxyglo <- mahalanobis.dist(data.x=oxy.lin, data.y=glo.apr)
oxyhep <- mahalanobis.dist(data.x=oxy.lin, data.y=hep.ful)
oxylat <- mahalanobis.dist(data.x=oxy.lin, data.y=lat.cal)
oxylei <- mahalanobis.dist(data.x=oxy.lin, data.y=lei.uni)
oxyaus <- mahalanobis.dist(data.x=oxy.lin, data.y=mel.aus)
oxyspl <- mahalanobis.dist(data.x=oxy.lin, data.y=mel.spl)
oxymog <- mahalanobis.dist(data.x=oxy.lin, data.y=mog.mog)
oxynem <- mahalanobis.dist(data.x=oxy.lin, data.y=nem.ere)
oxyate <- mahalanobis.dist(data.x=oxy.lin, data.y=neo.ate)
oxygra <- mahalanobis.dist(data.x=oxy.lin, data.y=neo.gra)
oxylin <- mahalanobis.dist(data.x=oxy.lin, data.y=oxy.lin)
oxypla <- mahalanobis.dist(data.x=oxy.lin, data.y=pla.ord)
oxytox <- mahalanobis.dist(data.x=oxy.lin, data.y=tox.cha)
plot (density(oxyamn), col=1,  las=1, xlim=c(0,6), ylim=c(0,0.8), xlab='', main='', ylab="Kernel density estimate")
lines(density(oxycra), col=2)
lines(density(oxyglo), col=3)
lines(density(oxyhep), col=4)
lines(density(oxylat), col=5)
lines(density(oxylei), col=6)
lines(density(oxyaus), col=7)
lines(density(oxyspl), col=8)
lines(density(oxymog), col=9)
lines(density(oxynem), col=10)
lines(density(oxyate), col=11)
lines(density(oxygra), col=12)
lines(density(oxylin), col=13, lwd=1.5, lty=4)
lines(density(oxypla), col=14)
lines(density(oxytox), col=15)
mtext(expression(paste("E) ",italic("Oxyeleotris lineolata"))),line=0.2, adj=0, padj=-0.1 )
mtext("Mahalanobis distance between individuals",1 , line=2)

#plaord
plaamn <- mahalanobis.dist(data.x=pla.ord, data.y=amn.per)
placra <- mahalanobis.dist(data.x=pla.ord, data.y=cra.ste)
plaglo <- mahalanobis.dist(data.x=pla.ord, data.y=glo.apr)
plahep <- mahalanobis.dist(data.x=pla.ord, data.y=hep.ful)
plalat <- mahalanobis.dist(data.x=pla.ord, data.y=lat.cal)
plalei <- mahalanobis.dist(data.x=pla.ord, data.y=lei.uni)
plaaus <- mahalanobis.dist(data.x=pla.ord, data.y=mel.aus)
plaspl <- mahalanobis.dist(data.x=pla.ord, data.y=mel.spl)
plamog <- mahalanobis.dist(data.x=pla.ord, data.y=mog.mog)
planem <- mahalanobis.dist(data.x=pla.ord, data.y=nem.ere)
plaate <- mahalanobis.dist(data.x=pla.ord, data.y=neo.ate)
plagra <- mahalanobis.dist(data.x=pla.ord, data.y=neo.gra)
plalin <- mahalanobis.dist(data.x=pla.ord, data.y=oxy.lin)
plapla <- mahalanobis.dist(data.x=pla.ord, data.y=pla.ord)
platox <- mahalanobis.dist(data.x=pla.ord, data.y=tox.cha)
plot (density(plaamn), col=1,  las=1, xlim=c(0,6), ylim=c(0,0.8), xlab='', main='', ylab="Kernel density estimate")
lines(density(placra), col=2)
lines(density(plaglo), col=3)
lines(density(plahep), col=4)
lines(density(plalat), col=5)
lines(density(plalei), col=6)
lines(density(plaaus), col=7)
lines(density(plaspl), col=8)
lines(density(plamog), col=9)
lines(density(planem), col=10)
lines(density(plaate), col=11)
lines(density(plagra), col=12)
lines(density(plalin), col=13)
lines(density(plapla), col=14, lwd=1.5, lty=4)
lines(density(platox), col=15)
mtext(expression(paste("E) ",italic("Planiliza ordensis"))),line=0.2, adj=0, padj=-0.1 )
mtext("Mahalanobis distance between individuals",1 , line=2)

#toxcha
toxamn <- mahalanobis.dist(data.x=tox.cha, data.y=amn.per)
toxcra <- mahalanobis.dist(data.x=tox.cha, data.y=cra.ste)
toxglo <- mahalanobis.dist(data.x=tox.cha, data.y=glo.apr)
toxhep <- mahalanobis.dist(data.x=tox.cha, data.y=hep.ful)
toxlat <- mahalanobis.dist(data.x=tox.cha, data.y=lat.cal)
toxlei <- mahalanobis.dist(data.x=tox.cha, data.y=lei.uni)
toxaus <- mahalanobis.dist(data.x=tox.cha, data.y=mel.aus)
toxspl <- mahalanobis.dist(data.x=tox.cha, data.y=mel.spl)
toxmog <- mahalanobis.dist(data.x=tox.cha, data.y=mog.mog)
toxnem <- mahalanobis.dist(data.x=tox.cha, data.y=nem.ere)
toxate <- mahalanobis.dist(data.x=tox.cha, data.y=neo.ate)
toxgra <- mahalanobis.dist(data.x=tox.cha, data.y=neo.gra)
toxlin <- mahalanobis.dist(data.x=tox.cha, data.y=oxy.lin)
toxpla <- mahalanobis.dist(data.x=tox.cha, data.y=pla.ord)
toxtox <- mahalanobis.dist(data.x=tox.cha, data.y=tox.cha)
plot (density(toxamn), col=1,  las=1, xlim=c(0,6), ylim=c(0,0.8), xlab='', main='', ylab="Kernel density estimate")
lines(density(toxcra), col=2)
lines(density(toxglo), col=3)
lines(density(toxhep), col=4)
lines(density(toxlat), col=5)
lines(density(toxlei), col=6)
lines(density(toxaus), col=7)
lines(density(toxspl), col=8)
lines(density(toxmog), col=9)
lines(density(toxnem), col=10)
lines(density(toxate), col=11)
lines(density(toxgra), col=12)
lines(density(toxlin), col=13)
lines(density(toxpla), col=14)
lines(density(toxtox), col=15, lwd=1.5, lty=4)
mtext(expression(paste("E) ",italic("Toxodes chatareous"))),line=0.2, adj=0, padj=-0.1 )
mtext("Mahalanobis distance between individuals",1 , line=2)

### 10-fold cross validation
library(caret)
flds <- createFolds(1:4302, 10)

fold.1 <- flds$Fold01
fold.2 <- flds$Fold02
fold.3 <- flds$Fold03
fold.4 <- flds$Fold04
fold.5 <- flds$Fold05
fold.6 <- flds$Fold06
fold.7 <- flds$Fold07
fold.8 <- flds$Fold08
fold.9 <- flds$Fold09
fold.10 <- flds$Fold10

####  Linear discriminant analysis

library(MASS)

x<-1:4302
store <- list ()

for (i in 1:10) {

samp.i <- x[-fold.10]
rest.i <- fold.10

model.final <- lda(species[samp.i]~., traits.fish[samp.i,1:9])
new.data.i<-traits.fish[rest.i,]

pred.class.i <- predict(model.final, list(BD=new.data.i$BD, PFL=new.data.i$PFL, throt=new.data.i$throt, ED=new.data.i$ED, EP=new.data.i$EP, GS=new.data.i$GS, GSh=new.data.i$GSh, MP=new.data.i$MP, BS=new.data.i$BS), type="class")

amnper.pred.i <- table(pred.class.i$class[new.data.i$species=='AMNPER'])
craste.pred.i <- table(pred.class.i$class[new.data.i$species=='CRASTE'])
gloapr.pred.i <- table(pred.class.i$class[new.data.i$species=='GLOAPR'])
hepful.pred.i <- table(pred.class.i$class[new.data.i$species=='HEPFUL'])
latcal.pred.i <- table(pred.class.i$class[new.data.i$species=='LATCAL'])
leiuni.pred.i <- table(pred.class.i$class[new.data.i$species=='LEIUNI'])
melaus.pred.i <- table(pred.class.i$class[new.data.i$species=='MELAUS'])
melspl.pred.i <- table(pred.class.i$class[new.data.i$species=='MELSPL'])
mogmog.pred.i <- table(pred.class.i$class[new.data.i$species=='MOGMOG'])
nemere.pred.i <- table(pred.class.i$class[new.data.i$species=='NEMERE'])
neoate.pred.i <- table(pred.class.i$class[new.data.i$species=='NEOATE'])
neogra.pred.i <- table(pred.class.i$class[new.data.i$species=='NEOGRA'])
oxylin.pred.i <- table(pred.class.i$class[new.data.i$species=='OXYLIN'])
plaord.pred.i <- table(pred.class.i$class[new.data.i$species=='PLAORD'])
toxcha.pred.i <- table(pred.class.i$class[new.data.i$species=='TOXCHA'])

lda.i <- rbind(amnper.pred.i,craste.pred.i,gloapr.pred.i,hepful.pred.i,latcal.pred.i,leiuni.pred.i,melaus.pred.i,melspl.pred.i,mogmog.pred.i,nemere.pred.i,neoate.pred.i,neogra.pred.i,oxylin.pred.i,plaord.pred.i,toxcha.pred.i)
rownames(lda.i) <- colnames(lda.i)
adds <- rowSums(lda.i)
add.mat <- cbind(adds,adds, adds, adds, adds, adds,adds,adds,adds,adds,adds,adds,adds,adds,adds)
for.plot.i <- lda.i/add.mat
store[[10]] <- for.plot.i
}



####  Linear discriminant analysis
library(MASS)

x<-1:4302

store <- list ()

for (i in 1:9999) {
samp.i <- sample(4302,2151)
rest.i <- x[-samp.i] 
#rest.i <- sample(113,28)

model.final <- lda(species[samp.i]~., traits.fish[samp.i,1:9])
new.data.i<-traits.fish[rest.i,]

pred.class.i <- predict(model.final, list(BD=new.data.i$BD, PFL=new.data.i$PFL, throt=new.data.i$throt, ED=new.data.i$ED, EP=new.data.i$EP, GS=new.data.i$GS, GSh=new.data.i$GSh, MP=new.data.i$MP, BS=new.data.i$BS), type="class")

amnper.pred.i <- table(pred.class.i$class[new.data.i$species=='AMNPER'])
craste.pred.i <- table(pred.class.i$class[new.data.i$species=='CRASTE'])
gloapr.pred.i <- table(pred.class.i$class[new.data.i$species=='GLOAPR'])
hepful.pred.i <- table(pred.class.i$class[new.data.i$species=='HEPFUL'])
latcal.pred.i <- table(pred.class.i$class[new.data.i$species=='LATCAL'])
leiuni.pred.i <- table(pred.class.i$class[new.data.i$species=='LEIUNI'])
melaus.pred.i <- table(pred.class.i$class[new.data.i$species=='MELAUS'])
melspl.pred.i <- table(pred.class.i$class[new.data.i$species=='MELSPL'])
mogmog.pred.i <- table(pred.class.i$class[new.data.i$species=='MOGMOG'])
nemere.pred.i <- table(pred.class.i$class[new.data.i$species=='NEMERE'])
neoate.pred.i <- table(pred.class.i$class[new.data.i$species=='NEOATE'])
neogra.pred.i <- table(pred.class.i$class[new.data.i$species=='NEOGRA'])
oxylin.pred.i <- table(pred.class.i$class[new.data.i$species=='OXYLIN'])
plaord.pred.i <- table(pred.class.i$class[new.data.i$species=='PLAORD'])
toxcha.pred.i <- table(pred.class.i$class[new.data.i$species=='TOXCHA'])

lda.i <- rbind(amnper.pred.i,craste.pred.i,gloapr.pred.i,hepful.pred.i,latcal.pred.i,leiuni.pred.i,melaus.pred.i,melspl.pred.i,mogmog.pred.i,nemere.pred.i,neoate.pred.i,neogra.pred.i,oxylin.pred.i,plaord.pred.i,toxcha.pred.i)
rownames(lda.i) <- colnames(lda.i)
adds <- rowSums(lda.i)
add.mat <- cbind(adds,adds, adds, adds, adds, adds,adds,adds,adds,adds,adds,adds,adds,adds,adds)
for.plot.i <- lda.i/add.mat
store[[i]] <- for.plot.i
}

for.plot<- Reduce("+", store) / length(store)


library("lattice")
jet.colors <-colorRampPalette(c("#FFFFFF","#BFBFBF","#7F7F7F","#404040","#000000"))
jet.colors1 <-colorRampPalette(c("#FFFFFF","#BFBFBF","#7F7F7F","#404040","#000000"))
levelplot( t(for.plot),col.regions=jet.colors, scale=list(x=list(rot=45)), xlab="Observed", ylab="Predicted")


predictions <- table(pred.class.i$class,new.data.i$species)
sensitivity(predictions, "GLOAPR")
specificity(predictions, "GLOAPR")

library(dplyr)
traits.fish2 <- traits.fish%>% group_by(species) %>% slice(21:n())
traits.fish3 <-traits.fish %>% group_by(species) %>% filter(row_number()==1:20)


library(MASS)
model1 <- lda(species~., traits.fish[,1:9])
pred <- predict (model1, newdata=traits.fish[,1:9])
amnper.pred <- table(pred$class[1:619])
craste.pred <- table(pred$class[620:732])
gloapr.pred <- table(pred$class[733:1120])
hepful.pred <- table(pred$class[1121:1251])
latcal.pred <- table(pred$class[1252:1450])
leiuni.pred <- table(pred$class[1451:1998])
melaus.pred <- table(pred$class[1999:2306])
melspl.pred <- table(pred$class[2307:2511])
mogmog.pred <- table(pred$class[2512:2768])
nemere.pred <- table(pred$class[3291:3849])
neoate.pred <- table(pred$class[3850:4136])
neogra.pred <- table(pred$class[4137:4254])
oxylin.pred <- table(pred$class[2769:2960])
plaord.pred <- table(pred$class[2961:3117])
toxcha.pred <- table(pred$class[3118:3290])
###
library(MASS)
model1 <- lda(traits.fish2$species~., traits.fish2[,2:10])
pred <- predict (model1, newdata=traits.fish3[,2:10])
amnper.pred <- table(pred$class[1:20])
craste.pred <- table(pred$class[21:40])
gloapr.pred <- table(pred$class[41:60])
hepful.pred <- table(pred$class[61:80])
latcal.pred <- table(pred$class[81:100])
leiuni.pred <- table(pred$class[101:120])
melaus.pred <- table(pred$class[121:140])
melspl.pred <- table(pred$class[141:160])
mogmog.pred <- table(pred$class[161:180])
oxylin.pred <- table(pred$class[181:200])
plaord.pred <- table(pred$class[201:220])
toxcha.pred <- table(pred$class[221:240])
nemere.pred <- table(pred$class[241:260])
neoate.pred <- table(pred$class[261:280])
neogra.pred <- table(pred$class[281:300])

lda <- rbind(amnper.pred,craste.pred,gloapr.pred,hepful.pred,latcal.pred,leiuni.pred,melaus.pred,melspl.pred,mogmog.pred,nemere.pred,neoate.pred,neogra.pred,oxylin.pred,plaord.pred,toxcha.pred)
rownames(lda) <- colnames(lda)
adds <- rowSums(lda)
add.mat <- cbind(adds,adds, adds, adds, adds, adds,adds,adds,adds,adds,adds,adds,adds,adds,adds)
for.plot <- lda/add.mat


jet.colors <-colorRampPalette(c("#FFFFFF","#BFBFBF","#7F7F7F","#404040","#000000"))
jet.colors1 <-colorRampPalette(c("#FFFFFF","#BFBFBF","#7F7F7F","#404040","#000000"))
jet.colors2 <- colorRampPalette(n=256, start='green', end='red', alpha=1)

library(RColorBrewer)

bluecols <- brewer.pal(9, 'YlOrRd')
newcol <- colorRampPalette(bluecols, bias=1.7)
ncols <- 200
bluecols2 <- newcol(ncols)#apply the function to get 100 colours


library("lattice")
levelplot( t(for.plot),col.regions=bluecols2, scale=list(x=list(rot=45)), xlab="Observed", ylab="Predicted")

library(caret)
confusionMatrix(species,pred$class )
