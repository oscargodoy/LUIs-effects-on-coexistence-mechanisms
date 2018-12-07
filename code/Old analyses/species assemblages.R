library(reshape)

## Assessing which species co-occur at each site for the 150 sites.
plants <- read.csv("data/BE.plants08.16.csv", header =T)

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
plants3 <- plants3/apply(plants3,1,sum)
names(plants3)[1:51] <- top50.short

#Pha_aru and Ver_cha are removed because only occur at one of the three sites, we can not therefore evaluate the effect of site
plants3 <- plants3[,-c(28,32,52)]
##Convert NA to Zeros to avoid future problems when summarizing
plants3[is.na(plants3)] <- 0
##Add the code to sum by year.
plants3 <- cbind(plants[,2], plants3)
names(plants3)[1] <- "Year"

#Sum over years to know which species are together 
pyear <- split(plants3, plants3$Year)
plants4 <- Reduce('+', pyear)
plants5 <- plants4[,-1]

#Put the code again and change values to 1 and 0, present or not.
plants5[plants5>0] <- 1

#Some metics to know the number of coexisting species within plots
mean(rowSums(plants5))
max(rowSums(plants5))
min(rowSums(plants5))

plants5 <- cbind(plants[1:150,c(1,4)], plants5)
write.csv(plants5, "results/presence_ausence_plots.csv")

