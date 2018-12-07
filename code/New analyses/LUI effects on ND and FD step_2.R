# This code is taken from:

#R-code of "A structural approach for understanding multispecies coexistence" by:
#Serguei Saavedra, Rudolf P. Rohr, Jordi Bascompte, Oscar Godoy, Nathan J.B. Kraft,
#and Jonathan M. Levine. Published in: Ecological Monographs

# Important! Serguei Saavedra and Rudolf Rohr confirmed that the framework only works with intraspecific positive values meaning all species experience some sort of intraspecific competition.
# So after sorting species alpha matrices needs to be multiply by -1

rm(list=ls())
source('code/toolbox_coexistence.R')
source('code/toolbox_figure.R')

#1.Obtaining combinations to calculate niche and fitness diff.----
#loading average values
intrinsic <- read.table("new results/intrinsic_site_lui_average_lme_50.csv", header=T, sep=",", row.names=1)
alpha <- as.matrix(read.table("new results/interaction_matrix_lme_average_50.csv", header=T, sep=",", row.names=1))
lui_modify_alpha <- as.matrix(read.table("new results/lui_matrix_lme_average_50.csv", header=T, sep=",", row.names=1))
top50.short <- c("Poa_tri", "Poa_pra", "Alo_pra", "Dac_glo", "Tri_rep", "Tar_off", "Lol_per", "Arr_ela", 
                 "Fes_rub", "Fes_pra", "Tri_fla", "Ely_rep", "Tri_pra", "Bro_ere", "Ran_rep", "Bro_hor", 
                 "Ran_acr", "Pla_lan","Ach_mil", "Gal_mol", "Her_sph", "Ant_syl", "Hol_lan", "Hel_pub",
                 "Ant_odo", "Bra_pin", "Car_hir", "Ver_cha", "Rum_ace", "Fes_ovi", "Phl_pra", "Pha_aru",
                 "Des_ces", "Agr_sto", "Cyn_cri", "Cir_ole", "Cer_hol", "Pla_med", "Cre_bie", "Urt_dio",
                 "Thy_pul", "Lol_mul", "Cir_arv", "Lot_cor", "Ran_bul", "Tri_dub", "Med_lup", "Leo_his",
                 "Car_car", "Vic_sep", "Pru_sp")

#Simulate x values of LUI at equal intervals. 
lui <- seq(from = 0.001, to = 4.5, by = 0.500) # we obtain 36 points to model. Is this range of intervals ok?

lui_intrinsic <- list()
lui_alpha <- list()
xx<- matrix(nrow=51, ncol=2, NA)
colnames(xx) <- c("intrinsic", "lui_value")
row.names(xx)<- paste(top50.short, sep=",")

#to obtain all value of intrinsic abundance at different LUI values
for(i in 1:(length(lui))){
  for(j in 1:length(row.names(intrinsic))){
    xx[j,1] <- intrinsic[j,1] + intrinsic[j,2]*lui[i]
    xx[j,2] <- lui[i]
    }
  lui_intrinsic[[i]] <- xx
}

#to obtain species with positive intrinsic at each site.
lui_intrinsic_positive <-  list()
zz <- list()
for(i in 1:(length(lui_intrinsic))){
  zz[[i]] <- row.names(subset(lui_intrinsic[[i]], lui_intrinsic[[i]][,1]>0))
  lui_intrinsic_positive[[i]] <- lui_intrinsic[[i]][which(rownames(lui_intrinsic[[i]]) %in% zz[[i]]),]
}

#to obtain interaction matrices at different LUI values
#subset the interaction matrices at different LUIs to those species showing positive growth rates.
hh <- list()

for(i in 1:(length(lui))){
    lui_alpha[[i]] <- alpha + lui_modify_alpha*lui[i]
    lui_alpha[[i]] <- lui_alpha[[i]][which(rownames(lui_alpha[[i]]) %in% row.names(lui_intrinsic_positive[[i]])),
                                     which(colnames(lui_alpha[[i]]) %in% row.names(lui_intrinsic_positive[[i]]))]
    hh[[i]] <- row.names(subset(lui_alpha[[i]], diag(lui_alpha[[i]])<0)) #we hold those with negative intra
    lui_alpha[[i]] <- lui_alpha[[i]][which(rownames(lui_alpha[[i]]) %in% hh[[i]]),
                                     which(colnames(lui_alpha[[i]]) %in% hh[[i]])]
}

#get back to the intrinsic and reduce the species list to those with negative intras
#also multiply by 100%
for(i in 1:(length(lui_intrinsic_positive))){
  lui_intrinsic_positive[[i]] <- lui_intrinsic_positive[[i]][which(rownames(lui_intrinsic_positive[[i]]) %in% rownames(lui_alpha[[i]])),]
}

#Obtain a list of combination of all species all combination of 3, 4, 5 and so on species until all species with positive abundance.
#This is super high demanding so limited to 3 so far
combos <- list()
combos_lui <- list()
for(i in 1:length(lui_alpha)){
  for(j in 2:3){
    combos[[j]] <- t(combn(row.names(lui_alpha[[i]]),j))
  }
  combos_lui[[i]] <- combos
}
  
#saveRDS(combos, 'new results/combos_lui.rds') In case the file is too big to compute every time.
#combos_lui <- readRDS("combos_lui.rds")

#2. Compute niche and fitness diff----
results <- list()
results_coex <- list()
col_results <- c("omega", "theta", "overlap", "differential", "feasibility", "lui")


for(i in 1:length(combos_lui)){
  
  for(j in 2:3){ #from 2 to 3 because there are combos up to three species
    results_combo <- matrix (nrow = nrow(combos_lui[[i]][[j]]), ncol= length(col_results))
    row.names(results_combo) <- apply(combos_lui[[i]][[j]], 1, paste, collapse=".") 
    colnames(results_combo) <- col_results
    
    for (h in 1:(length(combos_lui[[i]][[j]])/j)) {
      ll <- combos_lui[[i]][[j]][h,]
      mm <- lui_alpha[[i]][which(rownames(lui_alpha[[i]]) %in% ll), which(colnames(lui_alpha[[i]]) %in% ll)]
      mm <- -1*mm #this is because intras has to be positive. 
      ii <- subset(lui_intrinsic_positive[[i]], rownames(lui_intrinsic_positive[[i]]) %in% ll)
      om <- 10^Omega(mm) # this is niche differences
      th <- theta(mm, ii[,1]) # this is fitness differences
      co <- compute_overlap(mm,1000)
      oo <- co$overlap # this is community overlap
      di <- co$Omega - co$Omega_all # this is community differential
      fe <- test_feasibility(mm,ii[,1]) #this is whether all speceis can coexist. 
      results_combo[h,1] <- om
      results_combo[h,2] <- th
      results_combo[h,3] <- oo
      results_combo[h,4] <- di
      results_combo[h,5] <- fe
      results_combo[h,6] <- unique(lui_intrinsic_positive[[i]][,2])
    }
    results[[j]] <- results_combo
  }
  results_coex[[i]] <- results
}

#3. Plot results----

library(ggplot2)
library(gridExtra)

#Omega with LUI two species
omega <- list()
lui_omega <-list()

for(i in 1:length(results_coex)){
  omega[[i]] <-results_coex[[i]][[2]][,1]
  lui_omega[[i]] <-results_coex[[i]][[2]][,6]
}
omega <- unlist(omega)
lui_omega <- unlist(lui_omega)
omega_2 <- as.data.frame(cbind(omega, lui_omega))
colnames <- c("omega", "lui")
om_2 <- ggplot(omega_2,aes(x=lui_omega, y=omega)) + geom_point(alpha = 0.3) + labs(x = "LUI", y= "Omega", title = "Structural Niche Diff. 2 sps", tag = "A") + geom_smooth(method = "lm")

#Omega with LUI three species
omega <- list()
lui_omega <-list()

for(i in 1:length(results_coex)){
  omega[[i]] <-results_coex[[i]][[3]][,1]
  lui_omega[[i]] <-results_coex[[i]][[3]][,6]
}
omega <- unlist(omega)
lui_omega <- unlist(lui_omega)
omega_3 <- as.data.frame(cbind(omega, lui_omega))
colnames <- c("omega", "lui")
om_3 <- ggplot(omega_3,aes(x=lui_omega, y=omega)) + geom_point(alpha = 0.3) + labs(x = "LUI", y= "Omega", title = "Structural Niche Diff. 3 sps", tag = "B") + geom_smooth(method = "lm")

#Theta with LUI two species
theta <- list()
lui_theta <-list()

for(i in 1:length(results_coex)){
  theta[[i]] <-results_coex[[i]][[2]][,2]
  lui_theta[[i]] <-results_coex[[i]][[2]][,6]
}
theta <- log(unlist(theta))
lui_theta <- unlist(lui_theta)
theta_2 <- as.data.frame(cbind(theta, lui_theta))
colnames <- c("theta", "lui")
th_2 <- ggplot(theta_2,aes(x=lui_theta, y=theta)) + geom_point(alpha = 0.3) + labs(x = "LUI", y= "Log. transf. Theta", title = "Structural Fitness differences. 2 sps", tag = "C") + geom_smooth(method = "lm")

#Theta with LUI three species
theta <- list()
lui_theta <-list()

for(i in 1:length(results_coex)){
  theta[[i]] <-results_coex[[i]][[3]][,2]
  lui_theta[[i]] <-results_coex[[i]][[3]][,6]
}
theta <- log(unlist(theta))
lui_theta <- unlist(lui_theta)
theta_3 <- as.data.frame(cbind(theta, lui_theta))
colnames <- c("theta", "lui")
th_3 <- ggplot(theta_3,aes(x=lui_theta, y=theta)) + geom_point(alpha = 0.3) + labs(x = "LUI", y= "Log. transf. Theta", title = "Structural Fitness differences. 3 sps", tag = "D") + geom_smooth(method = "lm")

#put all together. 
grid.arrange(om_2, om_3, th_2, th_3, nrow=2, ncol=2)

#Plot how community overlap and community differential change with LUI.
#These are metrics of multispecies assemblages, so only for 3 species. 

#Overlap with LUI three species
overlap <- list()
lui_overlap <-list()

for(i in 1:length(results_coex)){
  overlap[[i]] <-results_coex[[i]][[3]][,3]
  lui_overlap[[i]] <-results_coex[[i]][[3]][,6]
}
overlap <- unlist(overlap)
lui_overlap <- unlist(lui_overlap)
overlap_3 <- as.data.frame(cbind(overlap, lui_overlap))
colnames <- c("overlap", "lui")
ov_3 <- ggplot(overlap_3,aes(x=lui_overlap, y=overlap)) + geom_point(alpha = 0.3) + labs(x = "LUI", y= "Overlap", title = "Community-pair Overlap", tag = "A") + geom_smooth(method = "lm")

#Differential with LUI three species
differential <- list()
lui_differential <-list()

for(i in 1:length(results_coex)){
  differential[[i]] <-results_coex[[i]][[3]][,4]
  lui_differential[[i]] <-results_coex[[i]][[3]][,6]
}
differential <- unlist(differential)
lui_differential <- unlist(lui_differential)
differential_3 <- as.data.frame(cbind(differential, lui_differential))
colnames <- c("differential", "lui")
di_3 <- ggplot(differential_3,aes(x=lui_differential, y=differential)) + geom_point(alpha = 0.3) + labs(x = "LUI", y= "Differential", title = "Community-pair Differential", tag = "B") + geom_smooth(method = "lm")

#put all together. 
grid.arrange(ov_3, di_3, ncol=2)

#Finally plot the number of species triplets and pairs that coexist across LUI

#Feasibility with LUI two species
feasibility <- list()
lui_feasibility <-list()

for(i in 1:length(results_coex)){
  feasibility[[i]] <-results_coex[[i]][[2]][,5]
  lui_feasibility[[i]] <-results_coex[[i]][[2]][,6]
}
feasibility <- unlist(feasibility)
lui_feasibility <- unlist(lui_feasibility)
feasibility_2 <- as.data.frame(cbind(feasibility, lui_feasibility))
colnames <- c("feasibility", "lui")
fe_2 <- ggplot(feasibility_2,aes(x=lui_feasibility, y=feasibility)) + geom_point(alpha = 0.3) + labs(x = "LUI", y= "Feasibility", title = " Feasible 2 sps", tag = "A") + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

#Feasibility with LUI three species
feasibility <- list()
lui_feasibility <-list()

for(i in 1:length(results_coex)){
  feasibility[[i]] <-results_coex[[i]][[3]][,5]
  lui_feasibility[[i]] <-results_coex[[i]][[3]][,6]
}
feasibility <- unlist(feasibility)
lui_feasibility <- unlist(lui_feasibility)
feasibility_3 <- as.data.frame(cbind(feasibility, lui_feasibility))
colnames <- c("feasibility", "lui")
fe_3 <- ggplot(feasibility_3,aes(x=lui_feasibility, y=feasibility)) + geom_point(alpha = 0.3) + labs(x = "LUI", y= "Feasibility", title = "Feasible 3 sps", tag = "B") + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

#put all together. 
grid.arrange(fe_2, fe_3, ncol=2)


