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
plants3 <- logit(plants3/apply(plants3,1,sum)) 
names(plants3)[1:51] <- top50.short

pyear <- split(plants3, plants$Year)
pchange <- list()


for(i in 1:(length(pyear)-1)){
  xx <- pyear[[i+1]]-pyear[[i]]
  names(xx) <- paste(names(xx), "_delta",sep="")
  xx2 <- cbind("LUI" = lui.only2[,i], xx)
  pchange[[i]] <- cbind(xx2, pyear[[i]])
}

pchange.all <- do.call("rbind", pchange)
Year_change <- paste(2008:2015, 2009:2016, sep="to")
pchange.all2 <- data.frame("Plot" = rep(unique(plants$Plot), 8), "Site" = strtrim(unique(plants$Plot), 1), "Year_change" = rep(Year_change, each =150), "Yeart" = rep(1:8, each = 150), pchange.all) 
yy <- names(pchange.all2)[grep("delta", names(pchange.all2))]
yy2 <- gsub("_delta", "", yy)

str(lCtr <- lmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-6, niterEM = 250, msMaxEval = 200))

##running al the models
mlist <- list()
for(i in 1:51){
  mlist[[i]] <- lme(as.formula(paste(yy[i], "~ Site + LUI*(", paste(top50.short, collapse="+"),"+Rest)")), data= pchange.all2,
                    random=~1|Plot, control=lCtr, correlation=corAR1(form=~Yeart), method='REML', na.action=na.omit)
}

coef.list <- lapply(mlist, function(x)summary(x)$coef$fixed)

#the miscelaneous species "REST" is not saved
inter.mat <- matrix(nrow=length(yy)-1, ncol=length(yy)-1)
lui.mat <- matrix(nrow=length(yy)-1, ncol=length(yy)-1)
intrinsic.site.lui <- matrix(nrow=length(yy)-1, ncol=4)

for(i in 1:length(coef.list)){
  
  cc <- coef.list[[i]] ## extract coefficients
  
  cc2 <- cc[-c(1:4)]
  cc3 <- cc2[c(1:51)]
  cc4 <- cc2[53:103]
  cc5 <- cc[1:4]
  
  inter.mat[i,] <- cc3
  lui.mat[i,] <- cc4
  intrinsic.site.lui[i,]<- cc5
}

yy2 <- yy2[-c(52)] # to remove name "rest"
row.names(inter.mat) <- yy2
colnames(inter.mat) <- yy2
row.names(lui.mat) <- yy2
colnames(lui.mat) <- yy2
row.names(intrinsic.site.lui) <- yy2
colnames(intrinsic.site.lui) <- c("Intrinsic", "SiteH", "SiteS", "LUI")

#Save all these matrices and then go to the structural stability approach in another R file

write.csv(inter.mat, "results/interaction_matrix_lme_50.csv")
write.csv(lui.mat, "results/lui_matrix_lme_50.csv")
write.csv(intrinsic.site.lui, "results/intrinsic_site_lui_lme_50.csv")

## calculate response to LUI of each species
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




###OLD STUFF IN CASE IT IS USEFUL, if it is not at the end then delete. 
pchange.all <- do.call("rbind", pchange)
Year_t1 <- 2009:2016
pchange.all2 <- data.frame("Plot" = rep(unique(plants$Plot), 8),"Site" = strtrim(unique(plants$Plot), 1), "Year_t1" = rep(Year_t1, each =150), "Year_t" = rep(2008:2015, each = 150), "Year_number" = rep(1:8, each = 150), pchange.all) # 8years, 150locations per year across site
yy <- names(pchange.all2)[grep("t1", names(pchange.all2))]
yy2 <- gsub("_t1", "", yy)
pchange.all2<- pchange.all2[complete.cases(pchange.all2), ]
##Either divide or substract to check which is better to get a stationary dynamic. 

#Poa_tri----
Poa_tri_diff <- (pchange.all2$Poa_tri_t1 - pchange.all2$Poa_tri)

Poa_tri<-cbind(pchange.all2[,1:6], Poa_tri_diff, pchange.all2[,59:110])
Poa_tri<- Poa_tri[complete.cases(Poa_tri[,7]), ]
plot.ts(Poa_tri$Year_t1, Poa_tri$Poa_tri_diff) #check that at year 1 is stationary then in CorARA d=1
pacf(Poa_tri$Poa_tri_diff) #higer values are at lag 4 and 5. then q is around these values 


x2 <- lme(Poa_tri_diff ~ Site + LUI*(Poa_tri+Poa_pra+Alo_pra+Dac_glo+Tri_rep+Tar_off+Lol_per+Arr_ela+Fes_rub+Fes_pra+Tri_fla+
    Ely_rep+Tri_pra+Bro_ere+Ran_rep+Bro_hor+Ran_acr+Pla_lan+Ach_mil+Gal_mol+Her_sph+Ant_syl+Hol_lan+Hel_pub+Ant_odo+
    Bra_pin+Car_hir+Ver_cha+Rum_ace+Fes_ovi+Phl_pra+Pha_aru+Des_ces+Agr_sto+Cyn_cri+Cir_ole+Cer_hol+Pla_med+Cre_bie+
      Urt_dio+Thy_pul+Lol_mul+Cir_arv+Lot_cor+Ran_bul+Tri_dub+Med_lup+Leo_his+Car_car+Vic_sep+Pru_sp+Rest), 
    random=~1|Plot, data= Poa_tri, control=lCtr, correlation=corARMA(form=~Year_number, p=1, q=1), method='REML', na.action=na.omit)

plot(Poa_tri$Poa_tri,Poa_tri$Poa_tri_diff)

coef_Poa_tri<-summary(x2)$tTable

#Poa_pra----
Poa_pra_diff <- (pchange.all2$Poa_pra_t1 - pchange.all2$Poa_pra)

Poa_pra<-cbind(pchange.all2[,1:6], Poa_pra_diff, pchange.all2[,59:110])
Poa_pra<- Poa_pra[complete.cases(Poa_pra[,7]), ]
plot.ts(Poa_pra$Year_t1, Poa_pra$Poa_pra_diff) #check that at year 1 is stationary then in CorARA d=1
pacf(Poa_pra$Poa_pra_diff) #higer values are at lag 4 and 5. then q is around these values 

str(lCtr <- lmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-6, niterEM = 250, msMaxEval = 200))

x <- lme(Poa_pra_diff ~ Site + LUI*(Poa_tri+Poa_pra+Alo_pra+Dac_glo+Tri_rep+Tar_off+Lol_per+Arr_ela+Fes_rub+Fes_pra+Tri_fla+
                                       Ely_rep+Tri_pra+Bro_ere+Ran_rep+Bro_hor+Ran_acr+Pla_lan+Ach_mil+Gal_mol+Her_sph+Ant_syl+Hol_lan+Hel_pub+Ant_odo+
                                       Bra_pin+Car_hir+Ver_cha+Rum_ace+Fes_ovi+Phl_pra+Pha_aru+Des_ces+Agr_sto+Cyn_cri+Cir_ole+Cer_hol+Pla_med+Cre_bie+
                                       Urt_dio+Thy_pul+Lol_mul+Cir_arv+Lot_cor+Ran_bul+Tri_dub+Med_lup+Leo_his+Car_car+Vic_sep+Pru_sp+Rest), 
          random=~1|Plot, data= Poa_pra, control=lCtr, correlation=corARMA(form=~Year_number, p=1, q=0), method='REML', na.action=na.omit)

plot(Poa_pra$Poa_pra,Poa_pra$Poa_pra_diff)

coef_Poa_pra<-summary(x2)$tTable


#Alo_pra----
Alo_pra_diff <- (pchange.all2$Alo_pra_t1 - pchange.all2$Alo_pra)

Alo_pra<-cbind(pchange.all2[,1:6], Alo_pra_diff, pchange.all2[,59:110])
Alo_pra<- Alo_pra[complete.cases(Alo_pra[,7]), ]
plot.ts(Alo_pra$Year_t1, Alo_pra$Alo_pra_diff) #check that at year 1 is stationary then in CorARA d=1
pacf(Alo_pra$Alo_pra_diff) #higer values are at lag 4 and 5. then q is around these values 

str(lCtr <- lmeControl(maxIter = 500, msMaxIter = 500, tolerance = 1e-6, niterEM = 250, msMaxEval = 200))

x <- lme(Alo_pra_diff ~ Site + LUI*(Poa_tri+Poa_pra+Alo_pra+Dac_glo+Tri_rep+Tar_off+Lol_per+Arr_ela+Fes_rub+Fes_pra+Tri_fla+
                                       Ely_rep+Tri_pra+Bro_ere+Ran_rep+Bro_hor+Ran_acr+Pla_lan+Ach_mil+Gal_mol+Her_sph+Ant_syl+Hol_lan+Hel_pub+Ant_odo+
                                       Bra_pin+Car_hir+Ver_cha+Rum_ace+Fes_ovi+Phl_pra+Pha_aru+Des_ces+Agr_sto+Cyn_cri+Cir_ole+Cer_hol+Pla_med+Cre_bie+
                                       Urt_dio+Thy_pul+Lol_mul+Cir_arv+Lot_cor+Ran_bul+Tri_dub+Med_lup+Leo_his+Car_car+Vic_sep+Pru_sp+Rest), 
          random=~1|Plot, data= Alo_pra, control=lCtr, correlation=corARMA(form=~Year_number, p=1, q=0), method='REML', na.action=na.omit)

plot(Alo_pra$Alo_pra,Alo_pra$Alo_pra_diff)

coef_Alo_pra<-summary(x2)$tTable






xx <- profilelike.lme(formula = Poa_pra_t1 ~ Site + LUI*(Poa_pra+Poa_tri), random = ~ 1 | Plot, 
                      correlation=corAR1(form=~Year_number), data=pchange.all2, subject = Plot,
                      profile.theta="Poa_pra", method="ML", lo.theta=1, hi.theta=5, length=500, round=2)


x1 <- lme(Poa_tri_t1 ~ Site + LUI*Poa_tri, random=~1|Plot, data= pchange.all2, correlation=corAR1(form=~Year_number), method='ML', na.action=na.omit)
x2 <- lme(Poa_tri_t1 ~ Site + LUI*Poa_tri, random=~1|Plot, data= pchange.all2, correlation=corARMA(form=~Year_number, p=1, q=1), method='REML', na.action=na.omit)
x3 <- lme(Poa_tri_t1 ~ Site + LUI*Poa_tri, random=~1|Plot, data= pchange.all2, correlation= corCAR1(0.2, form = ~Year_number), method='REML', na.action=na.omit)

