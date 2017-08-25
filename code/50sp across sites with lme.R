###MODELING INTERACTION COEFFICIENTS AND LUI EFFECTS WITH LME4.

library(reshape2)
library(nlme)
library(car)
#library(boot) only for the inv.logit

plants <- read.csv("data/BE.plants08.16.csv", header =T)

lui <- read.csv("data/LUI06_15.csv", header =T)

lui.only <- lui[, grep("LUI", names(lui))]
lui.only2 <- lui.only[, -c(1:2)] ## remove 2006 and 2007

plant.only <- plants[-c(1:5)]

### top50
top50 <- rev(sort(apply(plants[,-c(1:5)], 2, mean,na.rm=T)))[1:51] 
top50.short <- c("Poa_tri", "Poa_pra", "Alo_pra", "Dac_glo", "Tri_rep", "Tar_off", "Lol_per", "Arr_ela", 
                 "Fes_rub", "Fes_pra", "Tri_fla", "Ely_rep", "Tri_pra", "Bro_ere", "Ran_rep", "Bro_hor", 
                 "Ran_acr", "Pla_lan","Ach_mil", "Gal_mol", "Her_sph", "Ant_syl", "Hol_lan", "Hel_pub",
                 "Ant_odo", "Bra_pin", "Car_hir", "Ver_cha", "Rum_ace", "Fes_ovi", "Phl_pra", "Pha_aru",
                 "Des_ces", "Agr_sto", "Cyn_cri", "Cir_ole", "Cer_hol", "Pla_med", "Cre_bie", "Urt_dio",
                 "Thy_pul", "Lol_mul", "Cir_arv", "Lot_cor", "Ran_bul", "Tri_dub", "Med_lup", "Leo_his",
                 "Car_car", "Vic_sep", "Pru_sp")

plants2 <- plants[,match(names(top50), names(plants))]

Rest <- apply(plant.only[, -match(names(top50), names(plant.only))], 1, sum,na.rm=T)

plants3 <- cbind(plants2, Rest)
plants3 <- plants3/apply(plants3,1,sum) # to rscale to no more than 100%
names(plants3)[1:51] <- top50.short

pyear <- split(plants3, plants$Year)
pchange <- list()
pchange.logit <- list()
pchange.log <- list()
pchange.relative <- list()

# This is to check which is the best transformation.
#1. differences
for(i in 1:(length(pyear)-1)){
  xx <- pyear[[i+1]]-pyear[[i]] #absolute change
  names(xx) <- paste(names(xx), "_delta",sep="")
  xx2 <- cbind("LUI" = lui.only2[,i], xx)
  pchange[[i]] <- cbind(xx2, pyear[[i]])
}

pchange.all <- do.call("rbind", pchange)
Year_change <- paste(2008:2015, 2009:2016, sep="to")
pchange.all2 <- data.frame("Plot" = rep(unique(plants$Plot), 8), "Site" = strtrim(unique(plants$Plot), 1), "Year_change" = rep(Year_change, each =150), "Yeart" = rep(1:8, each = 150), pchange.all) 
yy <- names(pchange.all2)[grep("delta", names(pchange.all2))]

plot(pchange.all2$Poa_tri, pchange.all2$Poa_tri_delta, xlab="Proportion cover Poa_tri_year_t",
     ylab = "Change in proportion cover Poa_tri_yeart+1 - yeart", main="No transformation")

#2. logit transformation
for(i in 1:(length(pyear)-1)){
  xx.logit <- logit(pyear[[i+1]]) - logit(pyear[[i]]) #logit transformation 
  names(xx.logit) <- paste(names(xx.logit), "_delta",sep="")
  xx2.logit <- cbind("LUI" = lui.only2[,i], xx.logit)
  pchange.logit[[i]] <- cbind(xx2.logit, pyear[[i]])
}

pchange.all.logit <- do.call("rbind", pchange.logit)
Year_change <- paste(2008:2015, 2009:2016, sep="to")
pchange.all2.logit <- data.frame("Plot" = rep(unique(plants$Plot), 8), "Site" = strtrim(unique(plants$Plot), 1), "Year_change" = rep(Year_change, each =150), "Yeart" = rep(1:8, each = 150), pchange.all.logit) 
yy <- names(pchange.all2.logit)[grep("delta", names(pchange.all2.logit))]

plot(pchange.all2.logit$Poa_tri, pchange.all2.logit$Poa_tri_delta, xlab="Proportion cover Poa_tri_logit(year_t)",
     ylab = "Change in proportion cover Poa_tri_logit(yeart+1) - logit(yeart)", main="Logit transformation")

#3. log response ratio transformation
for(i in 1:(length(pyear)-1)){
  xx.log <- log(pyear[[i+1]]/pyear[[i]]) # log response ratio
  names(xx.log) <- paste(names(xx.log), "_delta",sep="")
  xx2.log <- cbind("LUI" = lui.only2[,i], xx.log)
  pchange.log[[i]] <- cbind(xx2.log, pyear[[i]])
}

pchange.all.log <- do.call("rbind", pchange.log)
Year_change <- paste(2008:2015, 2009:2016, sep="to")
pchange.all2.log <- data.frame("Plot" = rep(unique(plants$Plot), 8), "Site" = strtrim(unique(plants$Plot), 1), "Year_change" = rep(Year_change, each =150), "Yeart" = rep(1:8, each = 150), pchange.all.log) 
yy <- names(pchange.all2.log)[grep("delta", names(pchange.all2.log))]

plot(pchange.all2.log$Poa_tri, pchange.all2.log$Poa_tri_delta, xlab="Proportion cover Poa_tri_year_t",
     ylab = "Change in proportion cover Poa_tri_log(yeart+1 / yeart)", main="Log response ratio transformation")

#4. relative change in cover transformation
for(i in 1:(length(pyear)-1)){
  xx.relative <- (pyear[[i+1]] - pyear[[i]]) / (max(pyear[[i+1]],pyear[[i]])) # relative change cover
  names(xx.relative) <- paste(names(xx.relative), "_delta",sep="")
  xx2.relative <- cbind("LUI" = lui.only2[,i], xx.relative)
  pchange.relative[[i]] <- cbind(xx2.relative, pyear[[i]])
}

pchange.all.relative <- do.call("rbind", pchange.relative)
Year_change <- paste(2008:2015, 2009:2016, sep="to")
pchange.all2.relative <- data.frame("Plot" = rep(unique(plants$Plot), 8), "Site" = strtrim(unique(plants$Plot), 1), "Year_change" = rep(Year_change, each =150), "Yeart" = rep(1:8, each = 150), pchange.all.relative) 
yy <- names(pchange.all2.relative)[grep("delta", names(pchange.all2.relative))]

plot(pchange.all2.relative$Poa_tri, pchange.all2.relative$Poa_tri_delta, xlab="Proportion cover Poa_tri_year_t",
     ylab = "Change in proportion cover Poa_tri_(yeart+1 - yeart)/max(yeart+1, yeart)", main="Relative change in cover transformation")

#perform the modelling with lme and temporal autocorrelation

str(lCtr <- lmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-6, niterEM = 250, msMaxEval = 200))

##running all the models
mlist <- list()
for(i in 1:51){
  mlist[[i]] <- lme(as.formula(paste(yy[i], "~ Site + LUI*(", paste(top50.short, collapse="+"),")")), data= pchange.all2,
                    random=~1|Plot, control=lCtr, correlation=corAR1(form=~Yeart), method='REML', na.action=na.omit)
}

mlist.logit <- list()
for(i in 1:51){
  mlist.logit[[i]] <- lme(as.formula(paste(yy[i], "~ Site + LUI*(", paste(top50.short, collapse="+"),")")), data= pchange.all2.logit,
                    random=~1|Plot, control=lCtr, correlation=corAR1(form=~Yeart), method='REML', na.action=na.omit)
}

mlist.log <- list()
for(i in 1:51){
  mlist.log[[i]] <- lme(as.formula(paste(yy[i], "~ Site + LUI*(", paste(top50.short, collapse="+"),")")), data= pchange.all2.log,
                          random=~1|Plot, control=lCtr, correlation=corAR1(form=~Yeart), method='REML', na.action=na.omit)
}

mlist.relative <- list()
for(i in 1:51){
  mlist.relative[[i]] <- lme(as.formula(paste(yy[i], "~ Site + LUI*(", paste(top50.short, collapse="+"),")")), data= pchange.all2.relative,
                        random=~1|Plot, control=lCtr, correlation=corAR1(form=~Yeart), method='REML', na.action=na.omit)
}

##Check by AIC which transformation is best
for (i in 1:51){
  print(AIC(mlist[[i]],mlist.logit[[i]], mlist.relative[[i]]))
}

##According to AIC raw differences in proportion of cover are the best model supported across all species.
##Lets save these coefficients

h<- lme(Poa_tri_delta ~ Site + Poa_tri,data= pchange.all2,random=~1|Plot, control=lCtr, correlation=corAR1(form=~Yeart), method='REML', na.action=na.omit)

predictors <- expand.grid(Site= c("A","H","S"),
                          Sp= "Poa_tri")
predict(h, newdata=predictors, levels=0)


coef.list <- lapply(mlist, function(x)summary(x)$coef$fixed)

#the miscelaneous species "REST" is not saved
inter.mat <- matrix(nrow=length(yy)-1, ncol=length(yy)-1)
lui.mat <- matrix(nrow=length(yy)-1, ncol=length(yy)-1)
intrinsic.site.lui <- matrix(nrow=length(yy)-1, ncol=4)

for(i in 1:length(coef.list)){
  
  cc <- coef.list[[i]] ## extract coefficients
  
  cc2 <- cc[-c(1:4)]
  cc3 <- cc2[c(1:51)]
  cc4 <- cc2[52:102]
  cc5 <- cc[1:4]
  
  inter.mat[i,] <- cc3
  lui.mat[i,] <- cc4
  intrinsic.site.lui[i,]<- cc5
}

yy2 <- gsub("_delta", "", yy)
yy2 <- yy2[-c(52)] # to remove name "rest"
row.names(inter.mat) <- yy2
colnames(inter.mat) <- yy2
row.names(lui.mat) <- yy2
colnames(lui.mat) <- yy2
row.names(intrinsic.site.lui) <- yy2
colnames(intrinsic.site.lui) <- c("Intrinsic", "SiteH", "SiteS", "LUI")

#check whether diagonal values of all species are negtive, if not these species needs to be removed in further analyses
#of structural stability as can not handle facilitation effects for intraspecific effects. 
diag(inter.mat)>0
#in three species this are the cases

#Save all these matrices and then go to the structural stability approach in another R file

write.csv(inter.mat, "results/interaction_matrix_lme_average_50.csv")
write.csv(lui.mat, "results/lui_matrix_lme_average_50.csv")
write.csv(intrinsic.site.lui, "results/intrinsic_site_lui_average_lme_50.csv")

##To plot 
## http://myweb.uiowa.edu/pbreheny/publications/visreg.pdf
visreg(mlist[[1]], "Cer_hol", by="LUI", overlay=T, partial=T, type="contrast")

#We also calculate errors ###This needs to be asked to ERIC. Basically is find a way to propagate error

coef.list.error <- lapply(mlist, function(x)summary(x)$tTable[,2]) #column two correspond to std error.

#the miscelaneous species "REST" is not saved
inter.mat.error <- matrix(nrow=length(yy)-1, ncol=length(yy)-1)
lui.mat.error <- matrix(nrow=length(yy)-1, ncol=length(yy)-1)
intrinsic.site.lui.error <- matrix(nrow=length(yy)-1, ncol=4)

for(i in 1:length(coef.list.error)){
  
  cc <- coef.list.error[[i]] ## extract coefficients
  
  cc2 <- cc[-c(1:4)]
  cc3 <- cc2[c(1:51)]
  cc4 <- cc2[52:102]
  cc5 <- cc[1:4]
  
  inter.mat.error[i,] <- cc3
  lui.mat.error[i,] <- cc4
  intrinsic.site.lui.error[i,]<- cc5
}

row.names(inter.mat.error) <- yy2
colnames(inter.mat.error) <- yy2
row.names(lui.mat.error) <- yy2
colnames(lui.mat.error) <- yy2
row.names(intrinsic.site.lui.error) <- yy2
colnames(intrinsic.site.lui.error) <- c("Intrinsic", "SiteH", "SiteS", "LUI")

write.csv(inter.mat.error, "results/interaction_matrix_lme_std_error_50.csv")
write.csv(lui.mat.error, "results/lui_matrix_lme_std_error_50.csv")
write.csv(intrinsic.site.lui.error, "results/intrinsic_site_lui_std_error_lme_50.csv")
#From here we can calculate latter 95CI intervals multiplying by ±1.96. 

## calculate response to LUI of each species (Put as a Table)----
source("code/lui.niche.R")
lui.mean <- read.table("data/LUI.08.15.mean.txt", header =T)

library(SDMTools)
plants4 <- aggregate(plants2, list(plants$Useful_EP_PlotID), mean,na.rm=T)
names(plants4)[1] <- "Plot"
plants5 <- melt(plants4)
plants6 <- merge(plants5, lui.mean[,c(1,5)], by = "Plot", all.x =T)
names(plants6)[2:3] <- c("Species", "Abundance")

ln <- lui.niche(Data = plants6, plot.col = "Plot", gradient.col = "LUI", abundance.col = "Abundance", runs = 1000)

fates <- ln[[3]][["fate"]][,2]
ses <- ln[[3]][["ses"]][,1]

lui.resp <- data.frame("Species" = unique(plants6$Species), fates, ses)
lui.resp <- lui.resp[match(names(top50), lui.resp$Species),]


