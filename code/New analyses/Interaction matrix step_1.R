###MODELING INTERACTION COEFFICIENTS AND LUI EFFECTS WITH LMER 
## IN ORDER TO ACCOUNT FOR TEMPORAL AUTOCORRELATION.

library(reshape2)
library(nlme)
library(car)

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

# This is to prepare the database in order to conduct the analyses where t+1 will be compared to t
for(i in 1:(length(pyear)-1)){
  xx <- pyear[[i+1]]
  names(xx) <- paste(names(xx), "_delta",sep="")
  xx2 <- cbind("LUI" = lui.only2[,i], xx)
  pchange[[i]] <- cbind(xx2, pyear[[i]])
}

pchange.all <- do.call("rbind", pchange)
Year_change <- paste(2008:2015, 2009:2016, sep="to")
pchange.all2 <- data.frame("Plot" = rep(unique(plants$Plot), 8), "Site" = strtrim(unique(plants$Plot), 1), "Year_change" = rep(Year_change, each =150), "Yeart" = rep(1:8, each = 150), pchange.all) 
yy <- names(pchange.all2)[grep("delta", names(pchange.all2))]

#Let's plot a couple of examples to see how the data looks like 
par(mfrow=c(2,2))
plot(pchange.all2$Poa_tri, pchange.all2$Poa_pra_delta, xlab="P. cover Poa_tri_year_t",
     ylab = "P. cover Poa_pra_yeart+1")
plot(pchange.all2$Poa_pra, pchange.all2$Poa_pra_delta, xlab="P. cover Poa_pra_year_t",
     ylab = "P. cover Poa_pra_yeart+1")
plot(pchange.all2$Alo_pra, pchange.all2$Alo_pra_delta, xlab="P. cover Alo_pra_year_t",
     ylab = "P. cover Alo_pra_yeart+1")
plot(pchange.all2$Dac_glo, pchange.all2$Dac_glo_delta, xlab="P. cover Dac_glo_year_t",
     ylab = "P. cover Dac_glo_yeart+1")

#perform the modelling with lme and temporal autocorrelation

lCtr <- lmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-6, niterEM = 250, msMaxEval = 200)

##running all the models for the 51 sps. Site does not improve the model, so I have removed it. 
## I haven't also included weight factors because does not either improve the model. 
mlist <- list()
for(i in 1:51){
  mlist[[i]] <- lme(as.formula(paste(yy[i], " ~ LUI*(", paste(top50.short, collapse="+"),")")), data= pchange.all2, 
                random=~1|Plot/Site, control=lCtr, correlation=corAR1(form=~Yeart), method='REML', na.action=na.omit)
}

coef.list <- lapply(mlist, function(x)summary(x)$coef$fixed)

#the miscelaneous species "REST" is not saved
inter.mat <- matrix(nrow=length(yy)-1, ncol=length(yy)-1) #matrix of species interactions.
lui.mat <- matrix(nrow=length(yy)-1, ncol=length(yy)-1) #matrix of how LUI modify pairwise interactions.
intrinsic.site.lui <- matrix(nrow=length(yy)-1, ncol=2) #matrix of intrinsic ability to growth and how LUI modifies it. 

for(i in 1:length(coef.list)){
  
  cc <- coef.list[[i]] ## extract coefficients
  
  cc2 <- cc[-c(1:2)]
  cc3 <- cc2[c(1:51)]
  cc4 <- cc2[52:102]
  cc5 <- cc[1:2]
  
  inter.mat[i,] <- cc3
  lui.mat[i,] <- cc4
  intrinsic.site.lui[i,]<- cc5
}

yy2 <- gsub("_delta", "", yy)
yy2 <- yy2[-c(52)] # to remove name "rest"
diag(inter.mat)<- diag(inter.mat) -1 #because assuming that the 45 degrees slope means no change.
row.names(inter.mat) <- yy2
colnames(inter.mat) <- yy2
row.names(lui.mat) <- yy2
colnames(lui.mat) <- yy2
row.names(intrinsic.site.lui) <- yy2
colnames(intrinsic.site.lui) <- c("Intrinsic", "LUI")

#check whether diagonal values of all species are negative, if not these species needs to be removed in further analyses
#of structural stability as the method does not handle facilitation effects for intraspecific effects. 
diag(inter.mat)>0
#in three species this are the case but might change to more or less depending of the LUI's effect.

#Save all these matrices and then go to the structural stability approach in another R file

write.csv(inter.mat, "new results/interaction_matrix_lme_average_50.csv")
write.csv(lui.mat, "new results/lui_matrix_lme_average_50.csv")
write.csv(intrinsic.site.lui, "new results/intrinsic_site_lui_average_lme_50.csv")

#We also calculate errors----

coef.list.error <- lapply(mlist, function(x)summary(x)$tTable[,2]) #column two correspond to std error.

#the miscelaneous species "REST" is not saved
inter.mat.error <- matrix(nrow=length(yy)-1, ncol=length(yy)-1)
lui.mat.error <- matrix(nrow=length(yy)-1, ncol=length(yy)-1)
intrinsic.site.lui.error <- matrix(nrow=length(yy)-1, ncol=2)

for(i in 1:length(coef.list.error)){
  
  cc <- coef.list.error[[i]] ## extract coefficients
  
  cc2 <- cc[-c(1:2)]
  cc3 <- cc2[c(1:51)]
  cc4 <- cc2[52:102]
  cc5 <- cc[1:2]
  
  inter.mat.error[i,] <- cc3
  lui.mat.error[i,] <- cc4
  intrinsic.site.lui.error[i,]<- cc5
}

row.names(inter.mat.error) <- yy2
colnames(inter.mat.error) <- yy2
row.names(lui.mat.error) <- yy2
colnames(lui.mat.error) <- yy2
row.names(intrinsic.site.lui.error) <- yy2
colnames(intrinsic.site.lui.error) <- c("Intrinsic", "LUI")

write.csv(inter.mat.error, "new results/interaction_matrix_lme_std_error_50.csv")
write.csv(lui.mat.error, "new results/lui_matrix_lme_std_error_50.csv")
write.csv(intrinsic.site.lui.error, "new results/intrinsic_site_lui_std_error_lme_50.csv")
#From here we can calculate latter 95CI intervals multiplying by ±1.96. 

## calculate response to LUI of each species (Put as a Table)----
##This serve to understand whether species are winners, neutrals or loosers to LUI
source("code/New analyses/lui.niche.R")
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


