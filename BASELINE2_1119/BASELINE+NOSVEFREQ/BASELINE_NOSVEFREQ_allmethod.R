#================================
# Embedded MRP with applications
#================================
#-- last updated: aug 2021 --
#-- uses public use baseline data
#-- just (Z), no X (svefreq2), use MRP to analyze
set.seed(51)
options(mc.cores = 4)
#-- load required packages
require(dplyr)
require(foreign)
require(MASS)
require(polyapost)
require(mvtnorm)
require(rstan)
require(LaplacesDemon)

#-- helper functions
`%notin%` <- Negate(`%in%`)
normalize <- function(vec){
  # does exactly what it's called
  return(vec/sum(vec))
}
resmatfun <- function(x, id=seq(1,ncol(x))){
  estx <- colMeans(x,na.rm=T)
  varx <- apply(x,2,var,na.rm=T)
  CIlower <- apply(x,2,quantile,0.025,na.rm=T)
  CIupper <- apply(x,2,quantile,0.975,na.rm=T)
  resmat <- data.frame(
    Estimate = estx,
    SE = sqrt(varx),
    CIlower = CIlower,
    CIupper = CIupper
  )
  rownames(resmat) <- id
  return(resmat)
}

#-- source cleaned data and group indices
# setwd("~/Desktop/fwbb_mrp/EMRP PAPER BACKUP/realdata-output/FINALCODE_applications/")
source("../BASELINE.R")

# initialize containers for results
sim <- 5000
staniters = 10000
# sim <- 1000
# staniters = 2000
M <- length(unique(samp$sampled_x1_label))
N <- sum(acs$perwt)

#== redefine groupids in terms of Z
x1tbl <- samp %>% group_by(J_cell) %>% summarise(z = mean(sampled_x1_label))

grp1z <-grp2z <- grp3z <- grp4z <- grp5z <- grp6z<- grp7z <- rep(0,M)
grp1z[unique(x1tbl$z[grp1id==1])] <- 1
grp2z[unique(x1tbl$z[grp2id==1])] <- 1
grp3z[unique(x1tbl$z[grp3id==1])] <- 1
grp4z[unique(x1tbl$z[grp4id==1])] <- 1
grp5z[unique(x1tbl$z[grp5id==1])] <- 1
grp6z[unique(x1tbl$z[grp6id==1])] <- 1
grp7z[unique(x1tbl$z[grp7id==1])] <- 1

#== MRP
wmrp_posterior <- rep(0,sim)
wmrp_posterior1 <- rep(0, sim)
wmrp_posterior2 <- rep(0, sim)
wmrp_posterior3 <- rep(0, sim)
wmrp_posterior4 <- rep(0, sim)
wmrp_posterior5 <- rep(0, sim)
wmrp_posterior6 <- rep(0, sim)
wmrp_posterior7 <- rep(0, sim)



#===== (1) Run multilevel regression via Stan=====
# parameters for stan model
Ma <- nlevels(samp$age3) # combinations of x1: age, sex, educat, povgap
Mb <- nlevels(samp$sex2)
Mc <- nlevels(samp$race3)
Md <- nlevels(samp$educat4)
Me <- nlevels(samp$povgap4)
Mab <- nlevels(samp$abcat)
Mbc <- nlevels(samp$bccat)
Mcd <- nlevels(samp$cdcat)
Mdx <- nlevels(samp$dfcat)
L <- nlevels(samp$svefreq2)
n <- nrow(samp)

y <- as.numeric(samp$foodindsev)-1

grptbl <- samp %>%
  group_by(age3,sex2,race3,educat4,povgap4,abcat,bccat,cdcat) %>%
  summarise(n(),zcts = mean(x1_popcts), zlab = mean(sampled_x1_label))
grptbl$sampJ <- 1:nrow(grptbl)
samp <- left_join(samp,grptbl)
J_stan <- nrow(grptbl) # if fitting with sampgroup
agroup <- as.numeric(grptbl$age3)
bgroup <- as.numeric(grptbl$sex2) -1
cgroup <- as.numeric(grptbl$race3) 
dgroup <- as.numeric(grptbl$educat4)
egroup <- as.numeric(grptbl$povgap4)
abgroup <- as.numeric(grptbl$abcat)
bcgroup <- as.numeric(grptbl$bccat)
cdgroup <- as.numeric(grptbl$cdcat)
cell_label <- samp$sampled_x1_label
zpopcts <- grptbl$zcts

stanfit <- stan(file = "outcome_model_BASELINE.stan",
                data = c("Ma", "Mb", "Mc", "Md", "Me", "L", "Mab","Mbc","Mcd",
                         "J_stan", "n", "y", 
                         "agroup", "bgroup","cgroup", "dgroup", "egroup",
                         "abgroup","bcgroup","cdgroup","cell_label"),
                iter=staniters,
                warmup=staniters-sim/4,control=list(adapt_delta=0.99),chains=4)

# extract estimates of cell means from the stan model
stanpars <-
  rstan::extract(object=stanfit, permuted = TRUE)
cellmeans_stan <- stanpars$cellmean
write.csv(stanpars, "mrppars.csv")
wts_wmrp <- normalize(zpopcts)
wts_wmrp1<- normalize(wts_wmrp*grp1z)
wts_wmrp2<- normalize(wts_wmrp*grp2z)
wts_wmrp3<- normalize(wts_wmrp*grp3z)
wts_wmrp4<- normalize(wts_wmrp*grp4z)
wts_wmrp5<- normalize(wts_wmrp*grp5z)
wts_wmrp6<- normalize(wts_wmrp*grp6z)
wts_wmrp7<- normalize(wts_wmrp*grp7z)

#===== (2) PS estimator =====
for(s in 1:sim){
  
  current_cellmean <- cellmeans_stan[s,]
  wmrp_posterior[s] <- crossprod(wts_wmrp, current_cellmean)
  wmrp_posterior1[s] <- crossprod(wts_wmrp1, current_cellmean)
  wmrp_posterior2[s] <- crossprod(wts_wmrp2, current_cellmean)
  wmrp_posterior3[s] <- crossprod(wts_wmrp3, current_cellmean)
  wmrp_posterior4[s] <- crossprod(wts_wmrp4, current_cellmean)
  wmrp_posterior5[s] <- crossprod(wts_wmrp5, current_cellmean)
  wmrp_posterior6[s] <- crossprod(wts_wmrp6, current_cellmean)
  wmrp_posterior7[s] <- crossprod(wts_wmrp7, current_cellmean)
  #************************************************************
  
  print(s)
}

#-- (4) Print results
resmatid <- c("Overall", "Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7")
allgrp_wmrp <- cbind(wmrp_posterior, wmrp_posterior1, wmrp_posterior2, wmrp_posterior3, wmrp_posterior4, wmrp_posterior5, wmrp_posterior6, wmrp_posterior7)

resmat_wmrp <- resmatfun(allgrp_wmrp,resmatid)
write.csv(resmat_wmrp, "BASELINE_mrp.csv")
