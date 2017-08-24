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
pchange.all2 <- data.frame("Plot" = rep(unique(plants$Plot), 8),"Site" = strtrim(unique(plants$Plot), 1), "Year_t1" = rep(Year_t1, each =150), "Year_t" = rep(2008:2015, each = 150), "Year_number" = rep(1:8, each = 150), pchange.all) # 8years, 150locations per year across site
yy <- names(pchange.all2)[grep("t1", names(pchange.all2))]
yy2 <- gsub("_t1", "", yy)

lmeCtlList <-lmeControl(maxIter=2,  msMaxIter=3, tolerance=1e-35, niter=4,  
                        msTol=1e-5, nlmStepMax=5,msVerbose=TRUE  ,returnObject=TRUE ) 

x <- lme(Poa_pra_t1 ~ Site + LUI*(Poa_tri+Poa_pra+Alo_pra+Dac_glo+Tri_rep+Tar_off+Lol_per+Arr_ela+Fes_rub+Fes_pra+Tri_fla+
      Ely_rep+Tri_pra+Bro_ere+Ran_rep+Bro_hor+Ran_acr+Pla_lan+Ach_mil+Gal_mol+Her_sph+Ant_syl+Hol_lan+Hel_pub+Ant_odo+
      Bra_pin+Car_hir+Ver_cha+Rum_ace+Fes_ovi+Phl_pra+Pha_aru+Des_ces+Agr_sto+Cyn_cri+Cir_ole+Cer_hol+Pla_med+Cre_bie+
      Urt_dio+Thy_pul+Lol_mul+Cir_arv+Lot_cor+Ran_bul+Tri_dub+Med_lup+Leo_his+Car_car+Vic_sep+Pru_sp+Rest),
      random=~1|Plot, data= pchange.all2, correlation=corAR1(form=~Year_number), method='REML', na.action=na.omit, control=lmeCtlList)

x1 <- lme(Poa_tri_t1 ~ Site + LUI*Poa_tri, random=~1|Plot, data= pchange.all2, correlation=corAR1(form=~Year_number), method='ML', na.action=na.omit)
x2 <- lme(Poa_tri_t1 ~ Site + LUI*Poa_tri, random=~1|Plot, data= pchange.all2, correlation=corARMA(form=~Year_number, p=1, q=1), method='REML', na.action=na.omit)
x3 <- lme(Poa_tri_t1 ~ Site + LUI*Poa_tri, random=~1|Plot, data= pchange.all2, correlation= corCAR1(0.2, form = ~Year_number), method='REML', na.action=na.omit)

