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

# Select 20 most common speceis ----

#Region 1 (AEG)
region.aeg <- plants[grep("AEG",plants$EP_PlotID),]
top20.aeg <- rev(sort(apply(region.aeg[,-c(1:5)], 2, sum,na.rm=T)))[1:20]
top20.aeg.short <- c("Alo_pra", "Poa_tri", "Tri_fla", "Tri_rep", "Fes_rub", "Dac_glo", "Bro_ere", "Ran_acr",
                     "Tri_pra", "Tar_sp", "Lol_per", "Gal_mol", "Hel_pub", "Poa_pra", "Her_sph", "Arr_ela",
                     "Ant_syl", "Fes_pra", "Pla_lan", "Ant_odo")
plants2.aeg <- region.aeg[,match(names(top20.aeg), names(region.aeg))]
plant.only <- region.aeg[-c(1:5)]
Rest <- apply(plant.only[, -match(names(top20.aeg), names(plant.only))], 1, sum,na.rm=T)
plants3 <- cbind(plants2.aeg, Rest)
names(plants3)[1:20] <- top20.aeg.short
#Transform the abundance data. Divide by 100 to make it a proportion and transform with logit
plants3 <- logit(plants3/100) # we need to divided by the sum to reach a max of 100 across plots
#Summarize all data
region.aeg.final <- cbind(region.aeg[c(1:4)], plants3)
#Do the same with the LUI
lui.aeg <- lui.only[grep("AEG",lui.only$Plot),]


#Region 2 (HEG)
region.heg <- plants[grep("HEG",plants$EP_PlotID),]
top20.heg <- rev(sort(apply(region.heg[,-c(1:5)], 2, sum,na.rm=T)))[1:20]
top20.heg.short <- c("Poa_pra", "Tar_sp", "Dac_glo", "Alo_pra", "Fes_pra", "Lol_per", "Tri_rep", "Arr_ela", 
                     "Fes_rub", "Poa_tri", "Tri_pra", "Ely_rep", "Bra_pin", "Tri_fla", "Bro_hor", "Bro_ere",
                     "Achi_mil", "Phl_pra", "Pla_lan", "Cre_bie")
plants2.heg <- region.heg[,match(names(top20.heg), names(region.heg))]
plant.only <- region.heg[-c(1:5)]
Rest <- apply(plant.only[, -match(names(top20.heg), names(plant.only))], 1, sum,na.rm=T)
plants3 <- cbind(plants2.heg, Rest)
names(plants3)[1:20] <- top20.heg.short
#Transform the abundance data to logit
plants3 <- logit(plants3/100)
#Summarize all data
region.heg.final <- cbind(region.heg[c(1:4)], plants3)
#Do the same with the LUI
lui.heg <- lui.only[grep("HEG",lui.only$Plot),]



#Region 3 (SEG)
region.seg <- plants[grep("SEG",plants$EP_PlotID),]
top20.seg <- rev(sort(apply(region.seg[,-c(1:5)], 2, sum,na.rm=T)))[1:20]
top20.seg.short <- c("Poa_tri", "Poa_pra", "Tri_rep", "Dac_glo", "Lol_per", "Tar_sp", "Ely_rep", "Arr_ela", 
                     "Alo_pra", "Ran_rep", "Fes_rub", "Hol_lan", "Bro_hor", "Car_hir", "Pha_aru", "Achi_mil",
                     "Fes_pra", "Pla_lan", "Cir_ole", "Urt_dio")
plants2.seg <- region.seg[,match(names(top20.seg), names(region.seg))]
plant.only <- region.seg[-c(1:5)]
Rest <- apply(plant.only[, -match(names(top20.seg), names(plant.only))], 1, sum,na.rm=T)
plants3 <- cbind(plants2.seg, Rest)
names(plants3)[1:20] <- top20.seg.short
#Transform the abundance data to logit
plants3 <- logit(plants3/100)
#Summarize all data
region.seg.final <- cbind(region.seg[c(1:4)], plants3)
#Do the same with the LUI
lui.seg <- lui.only[grep("SEG",lui.only$Plot),]


#hold the final objects to start modeling 
rm(list= ls()[!(ls() %in% c('region.aeg.final', 'lui.aeg', 'top20.aeg.short', 'region.heg.final', 'lui.heg',
                            'top20.heg.short', 'region.seg.final', 'lui.seg', 'top20.seg.short'))])


## Create the matrices of year t+1 as a function of year t. 

pyear <- split(region.aeg.final, region.aeg.final$Year) # A list to separate between years.

pchange <- list()

for(i in 1:(length(pyear)-1)){
  
  xx <- pyear[[i+1]]
  
  names(xx) <- paste(names(xx), "t_1",sep="")
  
  xx2 <- cbind("LUI" = lui.aeg[,i], xx)
  
  pchange[[i]] <- cbind(xx2, pyear[[i]])
}




#### Example with lme from package lmer 
mdl.ac <- lme(Alo_pra ~ Poa_tri+Tri_fla+Tri_rep+Fes_rub+Dac_glo+Bro_ere+Ran_acr+Tri_pra+Tar_sp+Lol_per+Gal_mol+Hel_pub+Poa_pra+Her_sph+Arr_ela+Ant_syl+Fes_pra+Pla_lan+Ant_odo+Rest, random = ~1 | EP_PlotID, data=region.aeg.final, 
              correlation = corAR1(form = ~ Year), na.action=na.omit)


summary(mdl.ac)
plot(fitted(mdl.ac),residuals(mdl.ac))
abline(h=0,lty=3)
qqnorm(mdl.ac)
acf(residuals(mdl.ac,type="p"))



