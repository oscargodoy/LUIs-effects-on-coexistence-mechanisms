# This code is taken from:

#R-code of "A structural approach for understanding multispecies coexistence" by:
#Serguei Saavedra, Rudolf P. Rohr, Jordi Bascompte, Oscar Godoy, Nathan J.B. Kraft,
#and Jonathan M. Levine. Published in: Ecological Monographs

# Important! Serguei Saavedra and Rudolf Rohr confirmed that the framework only works with intraspecific positive values meaning all species experience some sort of intraspecific competition.
# So after sorting species alpha matrices needs to be multiply by -1

rm(list=ls())
source('code/toolbox_coexistence.R')
source('code/toolbox_figure.R')

#loading average values
alpha <- as.matrix(read.table("results/interaction_matrix_lme_average_50.csv", header=T, sep=",", row.names=1))
intrinsic <- read.table("results/intrinsic_site_lui_average_lme_50.csv", header=T, sep=",", row.names=1)

#calculate average intrinsic growth rates across three sites
intrinsic$SiteA <- intrinsic$Intrinsic
intrinsic$SiteH <- intrinsic$Intrinsic + intrinsic$SiteH
intrinsic$SiteS <- intrinsic$Intrinsic + intrinsic$SiteS
lui <- intrinsic$LUI
intrinsic <- intrinsic[,-c(1,4)]
intrinsic$Site_average <- apply(intrinsic,1,mean)
intrinsic$lui <-lui

#loading std.error values
alpha.error <- as.matrix(read.table("results/interaction_matrix_lme_std_error_50.csv", header=T, sep=",", row.names=1))
alpha_95_low <- alpha * (-1.96*alpha.error)
diag(alpha_95_low)>0 
# Result: Same as before, three sps in the alpha matrix needs to be removed Urt_dio, Car_car, Ely_rep.
alpha_95_high <- alpha * (1.96*alpha.error)  
diag(alpha_95_high)>0
#This is the opposite these three species have now negative values. 
# as CI 95% can not be evaluated because of positive values yet they are symmetrical we will evaluate only at 95% of one of the sides.
rm(alpha_95_low)

#Coexistence no LUI----

#We need to perform some steps

#Step1: Which species can not coexist without LUI (negative intrinsic growth rates)

sp_no_lui <- row.names(subset(intrinsic, Site_average>0))

#subset in the alpha matrix and in the intrinsic growth rate matrix

alpha<- alpha[which(rownames(alpha) %in% sp_no_lui),which(colnames(alpha) %in% sp_no_lui)]
alpha_95_high<- alpha_95_high[which(rownames(alpha_95_high) %in% sp_no_lui),which(colnames(alpha_95_high) %in% sp_no_lui)]

#Step2: Which species can not be incorporated into the models due to facilitative intra effects
diag(alpha)>0 
# Result: Two sps in the alpha matrix needs to be removed at average values because they have intra facilitation: Urt_dio, Car_car
names<-c("Car_car", "Urt_dio")
alpha<- alpha[-which(rownames(alpha) %in% names),-which(colnames(alpha) %in% names)]
alpha_95_high<- alpha_95_high[-which(rownames(alpha_95_high) %in% names),-which(colnames(alpha_95_high) %in% names)]

# A total set of 31 species to work with NO LUI, the number os species varies as LUI increases.

#Three species
combos3 <- t(combn(rownames(alpha),3))
#combos4 <- t(combn(rownames(alpha),4))

nd_fd3 <- matrix(nrow=dim(combos3)[1], ncol=6)
row.names(nd_fd3) <- apply(combos3,1,paste,collapse=".") 
colnames(nd_fd3) <- c("omega", "theta", "feasibility","coexisting pairs1:2","coexisting pairs1:3","coexisting pairs2:3") 

#multipy by -1
alpha <- alpha*-1
alpha_95_high <- alpha_95_high*-1

#this loop can not calculate outcomes for two combinations (number 654 and 807)
for(i in 1:dim(combos3)[1]){
  alpha2 <- alpha[combos3[i,], combos3[i,]]
  nd_fd3[i,1] <- exp(Omega(alpha2))
  intrinsic2 <- intrinsic[combos3[i,],] ## This needs to be a vector of Site Average, now is a matrix !!!!!!
  nd_fd3[i,2] <- theta(alpha2,intrinsic2)
  nd_fd3[i,3] <- test_feasibility(alpha2,intrinsic2)
  x <- test_feasibility_pairs(alpha2,intrinsic2)
  nd_fd3[i,4:6] <- x$feasibility
  #nd_fd3[i,2] <- r_centroid(alpha2)
}

#Remove NA cases
nd_fd3 <- as.data.frame(nd_fd3[complete.cases(nd_fd3), ])

#Remove -Inf Omega values that occur when the intras has negative and positive values
nd_fd3 <- subset(nd_fd3, omega > -Inf)
hist(nd_fd3$omega)
summary(nd_fd3$omega)

# Replace negative niche differences with niche overlap=0
nd_fd3.2 <- nd_fd3
nd_fd3.2<-within(nd_fd3.2, omega[omega<0] <- 0)
hist(nd_fd3.2$omega)
summary(nd_fd3.2$omega)

par(mfrow=c(1,2))
hist(nd_fd3.2$omega) # these are the niche differences
hist(nd_fd3.2$theta) # these are the fitness differences

#Separate data.frame between the coexisting and the non-coexisting triplets
coex_triplets <- subset(nd_fd3.2, feasibility ==1)
non_coex_triplets <- subset(nd_fd3.2, feasibility ==0)

#Which of these triplets coexist due to intransitivity
intransitivity <- subset(coex_triplets, coex_triplets$`coexisting pairs1:2`==0 & coex_triplets$`coexisting pairs1:3`==0 & coex_triplets$`coexisting pairs2:3`==0)

#Which of these triplets coexist when also some species pair coexist (diffuse competition)
diffuse <- subset(coex_triplets, !(coex_triplets$`coexisting pairs1:2`==0 & coex_triplets$`coexisting pairs1:3`==0 & coex_triplets$`coexisting pairs2:3`==0))
diffuse <- subset(diffuse, !(diffuse$`coexisting pairs1:2`==1 & diffuse$`coexisting pairs1:3`==1 & diffuse$`coexisting pairs2:3`==1))

#Which of these triplets coexist when also the three pairs coexist. 
all_coex_pairs<- subset(coex_triplets, coex_triplets$`coexisting pairs1:2`==1 & coex_triplets$`coexisting pairs1:3`==1 & coex_triplets$`coexisting pairs2:3`==1)








#How the intras in the matrix are modified when include low and high LUI
#Load the data without the three species removed
lui <- as.matrix(read.table("results/lui_matrix_lme_average_50.csv", header=T, sep=",", row.names=1))
lui<- lui[-which(rownames(lui) %in% names),-which(colnames(lui) %in% names)]
lui.error <- as.matrix(read.table("results/lui_matrix_lme_std_error_50.csv", header=T, sep=",", row.names=1))
lui.error<- lui.error[-which(rownames(lui.error) %in% names),-which(colnames(lui.error) %in% names)]

#evaluate the effect of LUI
#low LUI set at 0.5
  #Average and error
alpha_low_lui <- alpha + (lui*0.5) #always is + no -
alpha_low_lui_95_high <- alpha_95_high - (lui.error*0.5*1.96)

#High LUI set at 3.5
#Average and error
alpha_high_lui<- alpha + (lui*3.5)
alpha_low_lui_95_high <- alpha_95_high - (lui.error*3.5*1.96)


#Check that with low and high LUI we are ok of all values in the diagonal being competition (negative values so far), remember to change to positive before performing analyses
diag(alpha_low_lui)<0
diag(alpha_high_lui)<0 


#calculate all possible combinations of three and four species.
combos3 <- t(combn(rownames(alpha),3))
combos4 <- t(combn(rownames(alpha),4))

#No LUI with 3sp----

#All combination of three species


#This is to remove cases with no niche differences
#row_sub = apply(nd_fd3.2, 1, function(row) all(row !=0 ))
#nd_fd3.2<- nd_fd3.2[row_sub,]

#Low LUI with 3sp----
rm(list=ls())
source('code/toolbox_coexistence.R')
source('code/toolbox_figure.R')

#3-species interaction matrix of figure 5 and 6
alpha <- as.matrix(read.table("results/interaction_matrix_50.csv", header=T, sep=",", row.names=1))
intrinsic <- as.matrix(read.table("results/intrinsic.csv", header=T, sep=",", row.names=1))



for(i in 2378:dim(combos4)[1]){
  alpha2 <- alpha[combos4[i,], combos4[i,]]
  nd_fd4[i,1] <- Omega(alpha2)
  intrinsic2 <- intrinsic[combos4[i,],]
  nd_fd4[i,2] <- theta(alpha2,intrinsic2)
  # nd_fd3[i,2] <- r_centroid(alpha2)
  #nd_fd3[i,4] <- test_feasibility(alpha2,intrinsic2)
  #nd_fd3[i,5] <- test_feasibility_pairs(alpha2,intrinsic2)
}

 
## This is only code to be removed if not used. 
#All combination of four species
combos4 <- t(combn(rownames(alpha),4))
nd_fd4 <- matrix(nrow=dim(combos4)[1], ncol=3)
row.names(nd_fd4) <- apply(combos4,1,paste,collapse=".") 
colnames(nd_fd4) <- c("omega", "theta", "feasibility") 

# It serves to check other properties of the coexisting species. 
#"centroid", "feasibility_pairs")
#structural niche difference (Omega)
Omega(alpha)

#centroid of the feasibility domain
r_c <- r_centroid(alpha)

#structural fitness difference (theta), for a given vector of intrinsic growth rates (r = (1,1,1)), 
r <- c(1,1,0.5)
theta(alpha,r)

#test if the system (alpah,r) is feasible
test_feasibility(alpha,r)

#test the feasibility of all the pairs in the system (alpah,r)
test_feasibility_pairs(alpha,r)

#compute the overlap betwen the feasibility domain and the domain of coexistence of all the pairs (with 10000 randomizations)
compute_overlap(alpha,10000)


#draw the 3D cone of feasibility (like figure 4)
cone_3sp(alpha,'cone_3D.pdf',c('Sp1','Sp2','Sp3'))

#draw the projection of the feasibility of a 3-species system 
#on the simplex (like figure 5b and 6c)
projection_3sp(alpha,'projection_3D.pdf',c('Sp1','Sp2','Sp3'))

#draw the projection of the feasibility of a 3-species system 
#on the simplex, and overlap the domaine fo feasibility of all pairs (like figure 6d)
projection_3sp_with_pairwise(alpha,'projection_3D_with_pairs.pdf',c('Sp1','Sp2','Sp3'))

#4-species interaction matrix of figure 5 and the projection of its fesability
#domain on the simplex
alpha <- matrix(c(1, 0.5, 0.05, 0.2, 0.4, 1, 0.5, 0.2, 0.3, 0.6, 1, 0.2, 0.2, 0.2, 0.2, 1),4,4)
projection_4sp(alpha,'projection_4D.pdf',c('Sp1','Sp2','Sp3','Sp4'))
