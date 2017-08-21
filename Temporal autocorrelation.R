###Temporal autocorrelation with modeling.

library(reshape2)
library(lme4)
library(nlme)
library(car)

# Loading the data ----
plants <- read.csv("data/BE.plants08.16.csv", header =T)
lui <- read.csv("data/LUI06_15.csv", header =T)

lui.only <- lui[-c(2:23)]## remove 2006 and 2007

plant.only <- plants[-c(1:5)]

# Select 20 most common species ----

#Region 1 (AEG)----
region.aeg <- plants[grep("AEG",plants$EP_PlotID),]
top15.aeg <- rev(sort(apply(region.aeg[,-c(1:5)], 2, sum,na.rm=T)))[1:15]
top15.aeg.short <- c("Alo_pra", "Poa_tri", "Tri_fla", "Tri_rep", "Fes_rub", "Dac_glo", "Bro_ere", "Ran_acr",
                     "Tri_pra", "Tar_sp", "Lol_per", "Gal_mol", "Hel_pub", "Poa_pra", "Her_sph") 
                     
#"Arr_ela","Ant_syl", "Fes_pra", "Pla_lan", "Ant_odo")

plants2.aeg <- region.aeg[,match(names(top15.aeg), names(region.aeg))]
plant.only <- region.aeg[-c(1:5)]
Rest <- apply(plant.only[, -match(names(top15.aeg), names(plant.only))], 1, sum,na.rm=T)
plants3 <- cbind(plants2.aeg, Rest)
names(plants3)[1:15] <- top15.aeg.short
#Transform the abundance data. Divide by the sum to make it a proportion of max 1 and transform with logit
plants3 <- logit(plants3/apply(plants3,1,sum))
#Summarize all data
region.aeg.final <- cbind(region.aeg[c(1:4)], plants3)
#Do the same with the LUI
lui.aeg <- lui.only[grep("AEG",lui.only$Plot),]
lui.aeg <- lui.aeg[-c(1)]


#Region 2 (HEG)----
region.heg <- plants[grep("HEG",plants$EP_PlotID),]
top15.heg <- rev(sort(apply(region.heg[,-c(1:5)], 2, sum,na.rm=T)))[1:20]
top15.heg.short <- c("Poa_pra", "Tar_sp", "Dac_glo", "Alo_pra", "Fes_pra", "Lol_per", "Tri_rep", "Arr_ela", 
                     "Fes_rub", "Poa_tri", "Tri_pra", "Ely_rep", "Bra_pin", "Tri_fla", "Bro_hor")

#"Bro_ere", "Achi_mil", "Phl_pra", "Pla_lan", "Cre_bie")

plants2.heg <- region.heg[,match(names(top15.heg), names(region.heg))]
plant.only <- region.heg[-c(1:5)]
Rest <- apply(plant.only[, -match(names(top15.heg), names(plant.only))], 1, sum,na.rm=T)
plants3 <- cbind(plants2.heg, Rest)
names(plants3)[1:15] <- top15.heg.short
#Transform the abundance data to logit
plants3 <- plants3/apply(plants3,1,sum)
#Summarize all data
region.heg.final <- cbind(region.heg[c(1:4)], plants3)
#Do the same with the LUI
lui.heg <- lui.only[grep("HEG",lui.only$Plot),]
lui.heg <- lui.heg[-c(1)]


#Region 3 (SEG)----
region.seg <- plants[grep("SEG",plants$EP_PlotID),]
top15.seg <- rev(sort(apply(region.seg[,-c(1:5)], 2, sum,na.rm=T)))[1:15]
top15.seg.short <- c("Poa_tri", "Poa_pra", "Tri_rep", "Dac_glo", "Lol_per", "Tar_sp", "Ely_rep", "Arr_ela", 
                     "Alo_pra", "Ran_rep", "Fes_rub", "Hol_lan", "Bro_hor", "Car_hir", "Pha_aru") 
                     
                     #"Achi_mil", "Fes_pra", "Pla_lan", "Cir_ole", "Urt_dio")
plants2.seg <- region.seg[,match(names(top15.seg), names(region.seg))]
plant.only <- region.seg[-c(1:5)]
Rest <- apply(plant.only[, -match(names(top15.seg), names(plant.only))], 1, sum,na.rm=T)
plants3 <- cbind(plants2.seg, Rest)
names(plants3)[1:15] <- top15.seg.short
#Transform the abundance data to logit
plants3 <- logit(plants3/apply(plants3,1,sum))
#Summarize all data
region.seg.final <- cbind(region.seg[c(1:4)], plants3)
#Do the same with the LUI
lui.seg <- lui.only[grep("SEG",lui.only$Plot),]
lui.seg <- lui.seg[-c(1)]


#hold the final objects to start modeling 
rm(list= ls()[!(ls() %in% c('region.aeg.final', 'lui.aeg', 'top15.aeg.short', 'region.heg.final', 'lui.heg',
                            'top15.heg.short', 'region.seg.final', 'lui.seg', 'top15.seg.short'))])

#Just get the data
region.aeg.final2 <- region.aeg.final[-c(1:4)]
region.heg.final2 <- region.heg.final[-c(1:4)]
region.seg.final2 <- region.seg.final[-c(1:4)]

#Analyses for Region AEG with lme4----

pyear <- split(region.aeg.final2, region.aeg.final$Year) # A list to separate between years.
pchange <- list()

for(i in 1:(length(pyear)-1)){
  xx <- pyear[[i+1]]
  names(xx) <- paste(names(xx), "_t1",sep="")
  xx2 <- cbind("LUI" = lui.aeg[,i], xx)
  pchange[[i]] <- cbind(xx2, pyear[[i]])
}

pchange.all <- do.call("rbind", pchange)
Year_t1 <- 2009:2016
pchange.all2 <- data.frame("Plot" = rep(region.aeg.final$Useful_EP_PlotID, 8), "Year_t1" = rep(Year_t1, each =50), "Year_t" = rep(2008:2015, each = 50), "Year_number" = rep(1:8, each = 50), pchange.all) # 8years, 50locations per year within site
yy <- names(pchange.all2)[grep("t1", names(pchange.all2))]
yy2 <- gsub("_t1", "", yy)

### fit the model for each of the 20 species + the rest
mlist <- list()

for(i in 1:16){
  mlist[[i]] <- lmer(as.formula(paste(yy[i], "~ ", paste(top15.aeg.short, collapse="+"),"+Rest",
"+(0 + Year_number|Plot) + (1|Plot) + (1|Year_t)")), data= pchange.all2)
}

coef.list <- lapply(mlist, function(x)summary(x)$coef)

alpha.mat.aeg <- matrix(nrow=length(yy)-1, ncol=length(yy)-1)
gamma.mat.aeg <- matrix(nrow=length(yy)-1, ncol=length(yy)-1)

for(i in 1:length(coef.list)){
  
  cc <- coef.list[[i]][,1] ## extract coefficients
  
  cc2 <- cc[-c(1:2)]
  cc3 <- cc2[-grep("LUI", names(cc2))]
  cc4 <- cc2[grep("LUI", names(cc2))]
  
  alpha.mat.aeg[i,] <- cc3
  gamma.mat.aeg[i,] <- cc4
}

row.names(alpha.mat.aeg) <- yy2[-1]
colnames(alpha.mat.aeg) <- yy2[-1]
row.names(gamma.mat.aeg) <- yy2[-1]
colnames(gamma.mat.aeg) <- yy2[-1]


#Analyses for Region HEG with lme4----

pyear <- split(region.heg.final2, region.heg.final$Year) # A list to separate between years.
pchange <- list()

for(i in 1:(length(pyear)-1)){
  xx <- pyear[[i+1]]
  names(xx) <- paste(names(xx), "_t1",sep="")
  xx2 <- cbind("LUI" = lui.heg[,i], xx)
  pchange[[i]] <- cbind(xx2, pyear[[i]])
}

pchange.all <- do.call("rbind", pchange)
Year_t1 <- 2009:2016
pchange.all2 <- data.frame("Plot" = rep(region.heg.final$Useful_EP_PlotID, 8), "Year_t1" = rep(Year_t1, each =50), "Year_t" = rep(2008:2015, each = 50), "Year_number" = rep(1:8, each = 50), pchange.all) # 8years, 50locations per year within site
yy <- names(pchange.all2)[grep("t1", names(pchange.all2))]
yy2 <- gsub("_t1", "", yy)

### fit the model for each of the 20 species + the rest
mlist <- list()

for(i in 1:16){
  mlist[[i]] <- lmer(as.formula(paste(yy[i], "~ LUI + LUI*(", paste(top15.heg.short, collapse="+"),"+Rest)",
                                      "+(0 + Year_number|Plot) + (1|Plot) + (1|Year_t)")), data= pchange.all2)
}

coef.list <- lapply(mlist, function(x)summary(x)$coef)

alpha.mat.heg <- matrix(nrow=length(yy)-1, ncol=length(yy)-1)
gamma.mat.heg <- matrix(nrow=length(yy)-1, ncol=length(yy)-1)

for(i in 1:length(coef.list)){
  
  cc <- coef.list[[i]][,1] ## extract coefficients
  
  cc2 <- cc[-c(1:2)]
  cc3 <- cc2[-grep("LUI", names(cc2))]
  cc4 <- cc2[grep("LUI", names(cc2))]
  
  alpha.mat.heg[i,] <- cc3
  gamma.mat.heg[i,] <- cc4
}

row.names(alpha.mat.heg) <- yy2[-1]
colnames(alpha.mat.heg) <- yy2[-1]
row.names(gamma.mat.heg) <- yy2[-1]
colnames(gamma.mat.heg) <- yy2[-1]

#Analyses for Region SEG with lme4----

pyear <- split(region.seg.final2, region.seg.final$Year) # A list to separate between years.
pchange <- list()

for(i in 1:(length(pyear)-1)){
  xx <- pyear[[i+1]]
  names(xx) <- paste(names(xx), "_t1",sep="")
  xx2 <- cbind("LUI" = lui.seg[,i], xx)
  pchange[[i]] <- cbind(xx2, pyear[[i]])
}

pchange.all <- do.call("rbind", pchange)
Year_t1 <- 2009:2016
pchange.all2 <- data.frame("Plot" = rep(region.seg.final$Useful_EP_PlotID, 8), "Year_t1" = rep(Year_t1, each =50), "Year_t" = rep(2008:2015, each = 50), "Year_number" = rep(1:8, each = 50), pchange.all) # 8years, 50locations per year within site
yy <- names(pchange.all2)[grep("t1", names(pchange.all2))]
yy2 <- gsub("_t1", "", yy)

### fit the model for each of the 20 species + the rest
mlist <- list()

for(i in 1:16){
  mlist[[i]] <- lmer(as.formula(paste(yy[i], "~ LUI + LUI*(", paste(top15.seg.short, collapse="+"),"+Rest)",
                                      "+(0 + Year_number|Plot) + (1|Plot) + (1|Year_t)")), data= pchange.all2)
}

coef.list <- lapply(mlist2, function(x)summary(x)$coef)

alpha.mat.seg <- matrix(nrow=length(yy)-1, ncol=length(yy)-1)
gamma.mat.seg <- matrix(nrow=length(yy)-1, ncol=length(yy)-1)

for(i in 1:length(coef.list)){
  
  cc <- coef.list[[i]][,1] ## extract coefficients
  
  cc2 <- cc[-c(1:2)]
  cc3 <- cc2[-grep("LUI", names(cc2))]
  cc4 <- cc2[grep("LUI", names(cc2))]
  
  alpha.mat.seg[i,] <- cc3
  gamma.mat.seg[i,] <- cc4
}

row.names(alpha.mat.seg) <- yy2[-1]
colnames(alpha.mat.seg) <- yy2[-1]
row.names(gamma.mat.seg) <- yy2[-1]
colnames(gamma.mat.seg) <- yy2[-1]


#### Example with lme from package lmer 

pchange.all2$Plot <-paste(pchange.all2$Plot, pchange.all2$Year_number  , sep="_")

mlist2 <- list()

for(i in 1:16){
mlist2[[i]] <- lme(as.formula(paste(yy[i], "~ LUI + LUI*(", paste(top15.seg.short, collapse="+"),"+Rest)")),random= ~1|Plot, data=pchange.all2, 
              correlation = corAR1(), na.action=na.omit, control=lmeControl(returnObject=TRUE)) ##Add the control to help convergence
}
coef.list <- lapply(mlist2, function(x)summary(x)$coef)



#Hold the final data for structural equation analyses. 





#### Example with lme from package lmer 
mdl.ac <- lme(Alo_pra ~ Poa_tri+Tri_fla+Tri_rep+Fes_rub+Dac_glo+Bro_ere+Ran_acr+Tri_pra+Tar_sp+Lol_per+Gal_mol+Hel_pub+Poa_pra+Her_sph+Arr_ela+Ant_syl+Fes_pra+Pla_lan+Ant_odo+Rest, random = ~1 | EP_PlotID, data=region.aeg.final, 
              correlation = corAR1(form = ~ Year), na.action=na.omit)


summary(mdl.ac)
plot(fitted(mdl.ac),residuals(mdl.ac))
abline(h=0,lty=3)
qqnorm(mdl.ac)
acf(residuals(mdl.ac,type="p"))



