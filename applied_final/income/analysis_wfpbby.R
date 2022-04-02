#=====================================#
# WFPBB direct imputation applications
#=====================================#
set.seed(51)
#---- SETUP -----
## load packages
require(dplyr)
require(polyapost)
require(LaplacesDemon)

#-- source cleaned data and group indices
source("BASELINE.R")
N <- sum(acs$perwt)
M <- length(unique(samp$sampled_x1_label))
n <- nrow(samp)
L <- 5000
# L <- 1000

Tfact <- 30
F_draw <- 20
fpbb_temp <- rep(0, ncol=n*Tfact)

#-- helper functions
`%notin%` <- Negate(`%in%`)
normalize <- function(vec){
  res <- vec
  res[!is.na(vec)] <- vec[!is.na(vec)]/sum(vec[!is.na(vec)])
  return(res)
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

# initialize containers for results
ybar <- rep(0, L) # overall
ybar1 <- rep(0, L) # group 1
ybar2 <- rep(0, L) # group 2
ybar3 <- rep(0, L) # group 3
ybar4 <- rep(0, L) # group 4
ybar5 <- rep(0, L) # group 5
ybar6 <- rep(0, L) # group 6
ybar7 <- rep(0, L) # group 7

ybarbb <- rep(0, F_draw) # overall
ybar1bb <- rep(0, F_draw) # group 1
ybar2bb <- rep(0, F_draw) # group 2
ybar3bb <- rep(0, F_draw) # group 3
ybar4bb <- rep(0, F_draw) # group 4
ybar5bb <- rep(0, F_draw) # group 5
ybar6bb <- rep(0, F_draw) # group 6
ybar7bb <- rep(0, F_draw) # group 7



#----- BEGIN SIMULATIONS ----
## beginning of runtime
t0 <- Sys.time()

samp$foodindsev <- as.numeric(samp$foodindsev)-1

#------ DRAW SYNTHETIC POPULATIONS ----
for(l in 1:L){
  
  ## (1) Bootstrap from parent sample
  bootstrap<-BayesianBootstrap(c(1:n),1)
  resample<-rmultinom(1,n,bootstrap)
  
  twt<-N*normalize(samp$wts*resample)
  unique_ind <- c(1:n)[twt!=0]
  twt <- twt[twt!=0]
  
  ## (2) WFPBB with new weights and unique indices
  for(f in 1:F_draw){
    fpbb_temp <- wtpolyap(ysamp = unique_ind, wts = twt, k = (n*Tfact)-length(unique_ind))
    tempdf <- samp[fpbb_temp, c("foodindsev", "grp1id","grp2id","grp3id","grp4id","grp5id","grp6id","grp7id")]
    ybarbb[f] <- mean(tempdf$foodindsev)
    ybar1bb[f] <- mean(tempdf$foodindsev[tempdf$grp1id==1])
    ybar2bb[f] <- mean(tempdf$foodindsev[tempdf$grp2id==1])
    ybar3bb[f] <- mean(tempdf$foodindsev[tempdf$grp3id==1])
    ybar4bb[f] <- mean(tempdf$foodindsev[tempdf$grp4id==1])
    ybar5bb[f] <- mean(tempdf$foodindsev[tempdf$grp5id==1])
    ybar6bb[f] <- mean(tempdf$foodindsev[tempdf$grp6id==1])
    ybar7bb[f] <- mean(tempdf$foodindsev[tempdf$grp7id==1])
  }
  
  ybar[l] <- mean(ybarbb)
  ybar1[l] <- mean(ybar1bb)
  ybar2[l] <- mean(ybar2bb)
  ybar3[l] <- mean(ybar3bb)
  ybar4[l] <- mean(ybar4bb)
  ybar5[l] <- mean(ybar5bb)
  ybar6[l] <- mean(ybar6bb)
  ybar7[l] <- mean(ybar7bb)
}  
#-- (4) Print results
resmatid <- c("Overall","Group 1","Group 2","Group 3","Group 4","Group 5","Group 6","Group 7")

allgrp_wfpbby <- cbind(ybar, ybar1, ybar2, ybar3, ybar4, ybar5, ybar6, ybar7)
resmat_wfpbby <- resmatfun(allgrp_wfpbby,id = resmatid)
write.csv(resmat_wfpbby, "BASELINE_wfpbby.csv")
