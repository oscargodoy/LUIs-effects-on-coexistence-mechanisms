
library(reshape2)
library(lme4)

plants <- read.csv("data/BE.plants08.16.csv", header =T)
lui <- read.csv("data/LUI06_15.csv", header =T)

lui.only <- lui[, grep("LUI", names(lui))]
lui.only2 <- lui.only[, -c(1:2)] ## remove 2006 and 2007

plant.only <- plants[-c(1:5)]

### top20
top20 <- rev(sort(apply(plants[,-c(1:5)], 2, sum,na.rm=T)))[1:20]
top20.short <- c("Poa_tri", "Poa_pra", "Alo_pra", "Dac_glo", "Tri_rep", "Tar_off", "Lol_per", "Arr_ela", "Fes_rub", "Fes_pra", "Tri_fla", "Ely_rep", "Tri_pra", "Bro_ere", "Ran_rep", "Bro_hor", "Ran_acr", "Pla_lan","Ach_mil", "Gal_mol")
 
plants2 <- plants[,match(names(top20), names(plants))]

Rest <- apply(plant.only[, -match(names(top20), names(plant.only))], 1, sum,na.rm=T)
 
plants3 <- cbind(plants2, Rest)
names(plants3)[1:20] <- top20.short

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

pchange.all2 <- data.frame("Plot" = rep(unique(plants$Plot), 8), "Year_change" = rep(Year_change, each =150), "Yeart" = rep(1:8, each = 150), pchange.all) 

yy <- names(pchange.all2)[grep("delta", names(pchange.all2))]
yy2 <- gsub("_delta", "", yy)

### fit the model for each of the 20 species + the rest

mlist <- list()

for(i in 1:21){

mlist[[i]] <- lmer(as.formula(paste(yy[i], "~ LUI + LUI*(", paste(top20.short, collapse="+"),"+Rest)","+(0+Yeart|Plot)+(1|Plot)+(1|Year_change)")), data= pchange.all2)

}

coef.list <- lapply(mlist, function(x)summary(x)$coef)

alpha.mat <- matrix(nrow=length(yy), ncol=length(yy))
gamma.mat <- matrix(nrow=length(yy), ncol=length(yy))

for(i in 1:length(coef.list)){

cc <- coef.list[[i]][,1] ## extract coefficients

cc2 <- cc[-c(1:2)]
cc3 <- cc2[-grep("LUI", names(cc2))]
cc4 <- cc2[grep("LUI", names(cc2))]

alpha.mat[i,] <- cc3
gamma.mat[i,] <- cc4
}

row.names(alpha.mat) <- yy2
colnames(alpha.mat) <- yy2
row.names(gamma.mat) <- yy2
colnames(gamma.mat) <- yy2

niche_diff<- matrix(NA, nrow= length(yy), ncol= length(yy))
rownames(niche_diff)<- yy2
colnames(niche_diff)<- yy2
fitness_diff<- matrix(NA, nrow= length(yy), ncol= length(yy))
rownames(fitness_diff)<- yy2
colnames(fitness_diff)<- yy2

for(i in 1:length(yy)){
for(j in 1:length(yy)){

niche_diff[i,j] <- 1-(sqrt((alpha.mat[i,j]*alpha.mat[j,i])/(alpha.mat[i,i]*alpha.mat[j,j])))
fitness_diff[i,j] <- sqrt((alpha.mat[j,j]*alpha.mat[j,i])/(alpha.mat[i,i]*alpha.mat[i,j]))
}
}

#### select the largest of the fitness diffs for each pair
fit.select <- function(fitdiffs, species){

fd1 <- fitdiffs
fd1[lower.tri(fd1)] <- NA
fd2 <- fitdiffs
fd2[upper.tri(fd2)] <- NA
fd3 <- t(fd2)

fd4 <- cbind(fd1[upper.tri(fd1)], fd3[upper.tri(fd3)])
ff <- unlist(apply(fd4, 1, function(x) x[which(x==max(x, na.rm=T))]))

ff2 <- rep(NA, nrow(fd3))
ff2[apply(fd4, 1, function(x)all(complete.cases(x)))] <- ff

ss <- paste(species[combn(21,2)[1,]], species[combn(21,2)[2,]], sep=".")

ff3 <- data.frame(ss, ff2)

return(ff3)
}

fd.selected <- fit.select(fitness_diff, yy2)

hist(niche_diff[niche_diff>0],xlab="Niche differences", col ="indianred",main="Histogram of pairwise niche differences")
hist(fd.selected[,2],xlab="Fitness differences", col ="indianred",main="Histogram of pairwise fitness differences")

plot(as.vector(niche_diff[lower.tri(niche_diff)]), ff3[,2], pch=16, col="indianred", xlab = "Pairwise niche differences", ylab="Pairwise fitness differences")
n <- seq(0,1,len=100)
f <- 1/(1-n)
lines(n, f)

### effect of LUI on niche and fitness diffs

LUI <- seq(0.5, 3.5, len=10)
niche_diff_lui <- rep(list(matrix(NA, nrow= length(yy), ncol= length(yy))),length(LUI))
fitness_diff_lui <- rep(list(matrix(NA, nrow= length(yy), ncol= length(yy))),length(LUI))

for(h in 1:length(LUI)){
for(i in 1:length(yy)){
for(j in 1:length(yy)){

niche_diff_lui[[h]][i,j] <- 1-(sqrt(((alpha.mat[i,j]+gamma.mat[i,j]*LUI[h])*(alpha.mat[j,i]+gamma.mat[j,i]*LUI[h]))/
((alpha.mat[i,i]+gamma.mat[i,i]*LUI[h])*(alpha.mat[j,j]+gamma.mat[j,j]*LUI[h]))))

fitness_diff_lui[[h]][i,j] <- sqrt(((alpha.mat[j,j]+gamma.mat[j,j]*LUI[h])*(alpha.mat[j,i]+gamma.mat[j,i]*LUI[h]))/
((alpha.mat[i,i]+gamma.mat[i,i]*LUI[h])*(alpha.mat[i,j]+gamma.mat[i,j]*LUI[h])))
}
}}

ndlui <- sapply(niche_diff_lui, function(x)mean(x[lower.tri(x)], na.rm=T))

fd.lui <- sapply(fd.sel.lui, function(x) (mean(x[,2],na.rm=T)))
fd.sel.lui <- lapply(fitness_diff_lui, fit.select, yy2)

plot(LUI, fd.lui, pch=16)

nd <- niche_diff_lui[[10]]

plot(as.vector(nd[lower.tri(nd)]), fd.sel.lui[[10]][,2], pch=16, col="indianred", xlab = "Pairwise niche differences", ylab="Pairwise fitness differences")
n <- seq(0,1,len=100)
f <- 1/(1-n)
lines(n, f)





