##### Niche model for empirical data #####
##### Nico Blüthgen, TU Darmstadt, modified by Eric Allan    #####

## function takes a species, plot matrix in long format
## the names of columns containing info on species, plots, abundances of the species and the gradient to be tested can be specified

lui.niche <- function(Data, plot.col="Plot", species.col = "Species", gradient.col = "LUI", abundance.col = "Abundance", runs=1000){

require(SDMTools)  
    
Nplots <- length(unique(Data[,plot.col]))  # How many plots?
Nspp <- length(unique(Data[,species.col])) # How many species?
Gradient <- sort(Data[,gradient.col])      # Gradient values 
Spec <- unique(Data[,species.col])         # Species names

observed <- matrix(ncol = 3, nrow = Nspp)  # define empty matrices

colnames(observed) <- c("Noccur", "AWM", "AWSD")

null.means <- matrix(nrow = Nspp, ncol = 4)
colnames(null.means) <- c("rmean.npos.occ", "rmean.npos.ab", "rmean.nbrdth.occ", "rmean.nbrdth.ab") 

null.sds <- matrix(nrow = Nspp, ncol = 4)
colnames(null.sds) <- c("rsd.npos.occ", "rsd.npos.ab", "rsd.nbrdth.occ", "rsd.nbrdth.ab") 

null.ps <- matrix(nrow = Nspp, ncol = 4)
colnames(null.ps) <- c("pval.npos.occ", "pval.npos.ab", "pval.nbrdth.occ", "pval.nbrdth.ab") 
  
for (i in 1:Nspp) {  # loop through all species
  
  SpData <- subset(Data, Data[,species.col]==Spec[i])
  SpData2 <- na.omit(SpData)                   ## remove NA plots
  observed[i , 1] <- nrow(SpData2)              # on how many plots does species i occur?
  Adist <- SpData2[,abundance.col]              # abundance distribution of species i
  observed[i , 2] <- weighted.mean(SpData2[,gradient.col], SpData2[,abundance.col])  # niche optimum of species i
  observed[i , 3] <- wt.sd(SpData2[,gradient.col], SpData2[,gradient.col])           # niche breadth of species i
  
  
  ### null models for niche model
  nulls <- matrix(nrow=runs, ncol = 4)
  colnames(nulls) <- c("mean.occ", "mean.ab", "sd.occ", "sd.ab")

  
  for (j in 1:runs) { 
    rand.grad <- sample(Gradient, nrow(SpData2))         # randomised gradient values
    nulls[j , 1] <- mean(rand.grad)                      # Occurence-based null model mean LUI (niche position)
    nulls[j , 2] <- weighted.mean(rand.grad, Adist)      # Abundance-based null model abundance weighted mean LUI 
    nulls[j , 3] <- sd(rand.grad)                        # Occurence-based null model sd LUI (niche breadth)
    nulls[j , 4] <- wt.sd(rand.grad, Adist)              # Abundance-based null model abundance weighted sd LUI
  }

  null.means[i , ] <- apply(nulls, 2, mean)
  null.sds[i , ] <- apply(nulls, 2, sd)

  ### p-values for niche null model
  oo <- rep(rep(observed[i , 2:3], each = 2))        ## observed values to compare with random
  oo2 <- matrix(rep(oo, nrow(nulls)),nrow = nrow(nulls), byrow=T)
  Praw <- 1*(nulls >= oo2)                    # how often are null model means greater than or equal to observed value?
  Praw <- apply(Praw, 2, sum) / nrow(nulls)
  nullp0 <- ifelse(Praw > 0.5, 1-Praw, Praw)  # P-value calculated from the non-overlap
  null.ps[i,] <- nullp0 
}

oo.np <- cbind(observed[,2], observed[,2]) ## rep the observed niche pos 2x to compare with random expectations
oo.nb <- cbind(observed[,3],observed[,3])  ## observed niche breadth
oo.all <- cbind(oo.np, oo.nb)
  
null.ses <- (oo.all-null.means)/null.sds

signp <- 1*(null.ps[,1:2]<=0.05) ## convert to significances (1,0) only for test of niche position
signb <- 1*(null.ps[,3:4]<=0.05) ## significances of niche breadth

signp[signp == 1 & oo.np > null.means[,1:2]] <- "winner" ## niche pos sig greater than random
signp[signp == 1 & oo.np < null.means[,1:2]] <- "loser" ## niche pos sig less than random
signp[signp == 0 & oo.nb < null.means[,3:4] & signb <= 0.05] <- "mid_specialist"
signp[signp==0] <- "neutral"


Info <- paste("Dataset contains",Nplots,"plots and",Nspp,"species.",runs,"randomizations performed each null model.")

Citation <- paste("Reference: Chisté MN, Mody K, Gossner MM, Simons NK, Köhler G, Weisser WW, Blüthgen N (2016) Losers, winners and opportunists: how grassland land-use intensity affects orthopteran communities. Ecosphere, DOI 10.1002/ecs2.1545")

Results <- list(observed, null.means, null.sds, null.ses, null.ps, signp)
lapply(Results, function(x){
  rownames(x) <- Spec
  return(x)})

names(Results) <- c("observed", "null_means", "null_sds", "ses", "pvals", "fate")
to.return <- list(Info,Citation,Results)

return(to.return)
}

