###MODELING INTERACTION COEFFICIENTS AND LUI EFFECTS WITH LME4.

library(reshape2)
library(lme4)
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
plants3 <- logit(plants3/apply(plants3,1,sum))
names(plants3)[1:51] <- top50.short

#Pha_aru is removed because only occur at one of the three sites, we can not therefore evaluate the effect of site
#Same occur with Ver_cha
plants3 <- plants3[,-c(28,32)]

pyear <- split(plants3, plants$Year)
pchange <- list()


for(i in 1:(length(pyear)-1)){
  xx <- pyear[[i+1]]
  names(xx) <- paste(names(xx), "_t1",sep="")
  xx2 <- cbind("LUI" = lui.only2[,i], xx)
  pchange[[i]] <- cbind(xx2, pyear[[i]])
}

pchange.all <- do.call("rbind", pchange)
Year_t1 <- 2009:2016
pchange.all2 <- data.frame("Plot" = rep(unique(plants$Plot), 8),"Site" = rep(1:3, each=50), "Year_t1" = rep(Year_t1, each =150), "Year_t" = rep(2008:2015, each = 150), "Year_number" = rep(1:8, each = 150), pchange.all) # 8years, 150locations per year across site
yy <- names(pchange.all2)[grep("t1", names(pchange.all2))]
yy2 <- gsub("_t1", "", yy)

top50.short <- c("Poa_tri", "Poa_pra", "Alo_pra", "Dac_glo", "Tri_rep", "Tar_off", "Lol_per", "Arr_ela", 
                 "Fes_rub", "Fes_pra", "Tri_fla", "Ely_rep", "Tri_pra", "Bro_ere", "Ran_rep", "Bro_hor", 
                 "Ran_acr", "Pla_lan","Ach_mil", "Gal_mol", "Her_sph", "Ant_syl", "Hol_lan", "Hel_pub",
                 "Ant_odo", "Bra_pin", "Car_hir", #"Ver_cha", 
                 "Rum_ace", "Fes_ovi", "Phl_pra", #"Pha_aru",
                 "Des_ces", "Agr_sto", "Cyn_cri", "Cir_ole", "Cer_hol", "Pla_med", "Cre_bie", "Urt_dio",
                 "Thy_pul", "Lol_mul", "Cir_arv", "Lot_cor", "Ran_bul", "Tri_dub", "Med_lup", "Leo_his",
                 "Car_car", "Vic_sep", "Pru_sp")

## calculate response to LUI of each species
source("code\\lui.niche.R")
lui.mean <- read.table("data\\LUI.08.15.mean.txt", header =T)


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


### fit the models

mlist <- list()
for(i in 1:50){
  mlist[[i]] <- lmer(as.formula(paste(yy[i], "~ (Site+LUI)*(", paste(top50.short, collapse="+"),"+Rest)","+(0+Year_number|Plot)+(1|Plot)+(1|Year_t)")), 
                     data= pchange.all2, REML=FALSE)
}

#There is warning when the model run, to check whether this is an issue do the following and see whether values are very low e.g. 0.001.
#If the values are as low as happen in our case then, we are ok to follow. https://stats.stackexchange.com/questions/97929/lmer-model-fails-to-converge
for(i in 1:50){
relgrad <- with(mlist[[i]]@optinfo$derivs,solve(Hessian,gradient))
print(max(abs(relgrad)))
}
      
coef.list <- lapply(mlist, function(x)summary(x)$coef)


inter.mat <- matrix(nrow=length(yy)-1, ncol=length(yy)-1)
lui.mat <- matrix(nrow=length(yy)-1, ncol=length(yy)-1)
site.mat <- matrix(nrow=length(yy)-1, ncol=length(yy)-1) #It has one dimension less because the species Pha_aru does not occur across sites
intrinsic.site.lui <- matrix(nrow=length(yy)-1, ncol=3)

for(i in 1:length(coef.list)){
  
  cc <- coef.list[[i]][,1] ## extract coefficients
  
  cc2 <- cc[-c(1:3)]
  cc3 <- cc2[c(1:50)]
  cc4 <- cc2[grep("LUI", names(cc2))]
  cc5 <- cc2[grep("Site", names(cc2))]
  cc6 <- cc[1:3]
  
  inter.mat[i,] <- cc3
  lui.mat[i,] <- cc4
  site.mat[i,] <- cc5
  intrinsic.site.lui[i,]<- cc6
}

yy2 <-yy2[-1] # to remove name "year"

row.names(inter.mat) <- yy2
colnames(inter.mat) <- yy2
row.names(lui.mat) <- yy2
colnames(lui.mat) <- yy2
row.names(site.mat) <- yy2
colnames(site.mat) <- yy2
row.names(intrinsic.site.lui) <- yy2
colnames(intrinsic.site.lui) <- c("Intrinsic", "Site", "LUI")

#Save all these matrices and then go to the structural stability approach in another R file

write.csv(inter.mat, "results/interaction_matrix_50.csv")
write.csv(lui.mat, "results/lui_matrix_50.csv")
write.csv(site.mat, "results/site_matrix_50.csv")
write.csv(intrinsic.site.lui, "results/intrinsic_site_lui_50.csv")


