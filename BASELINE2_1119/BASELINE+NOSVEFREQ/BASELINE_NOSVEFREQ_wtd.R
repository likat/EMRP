#================================
# Embedded MRP with applications
#================================
#-- last updated: july 21, 2020 --
#-- estimation of survey weights via wfpbb

set.seed(51)

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
#-- source cleaned data and group indices
# setwd("~/Desktop/fwbb_mrp/EMRP PAPER BACKUP/realdata-output/to-run/")

source("../BASELINE.R")
# setwd("~/Desktop/fwbb_mrp/EMRP PAPER BACKUP/realdata-output/BASELINE+NOSVEFREQ/")

# get the joint sample frequency of X1 x X2
svywts_srbi <- samp %>% 
  group_by(age3, sex2, race3, educat4, povgap4) %>% 
  summarise(
    srbi_mkcts = n(),
    ybarj = mean(as.numeric(foodindsev)-1),
    grp1id = mean(povgap4==1), #n=884
    grp2id = mean(povgap4 == 2), #n=273
    grp3id = mean(povgap4 == 3), # n=324
    grp4id = mean(povgap4 == 4), # n=202
    grp5id = mean(sex2 == 1), # n=218
    grp6id = mean(age3 == 1), # n=212
    grp7id = mean(educat4 == 4)) # n=123
J <- nrow(svywts_srbi)
ybarj <- svywts_srbi$ybarj
grp1id <- svywts_srbi$grp1id
grp2id <- svywts_srbi$grp2id
grp3id <- svywts_srbi$grp3id
grp4id <- svywts_srbi$grp4id
grp5id <- svywts_srbi$grp5id
grp6id <- svywts_srbi$grp6id
grp7id <- svywts_srbi$grp7id

svywts_srbi$J_cell <- 1:J
samp$J_cell <- NULL
samp <- samp %>% left_join(svywts_srbi)
samp$foodindsev <- as.numeric(samp$foodindsev)-1

n <- nrow(samp)


# initialize containers for results
sim <- 5000
# sim <- 2000
F_draw = 20
Tfact=30
Njhatmat <- matrix(0,nrow=sim, ncol=J)
srbi_overall <- rep(0,sim)
srbi1 <- rep(0, sim)
srbi2 <- rep(0, sim)
srbi3 <- rep(0, sim)
srbi4 <- rep(0, sim)
srbi5 <- rep(0, sim)
srbi6 <- rep(0, sim)
srbi7 <- rep(0, sim)
srbi_overallbb <- rep(0,F_draw)
srbi1bb <- rep(0, F_draw)
srbi2bb <- rep(0, F_draw)
srbi3bb <- rep(0, F_draw)
srbi4bb <- rep(0, F_draw)
srbi5bb <- rep(0, F_draw)
srbi6bb <- rep(0, F_draw)
srbi7bb <- rep(0, F_draw)


samp <- samp[,c("J_cell", "wts", "foodindsev", "srbi_mkcts", 
                "grp1id", "grp2id", "grp3id", "grp4id","grp5id", "grp6id", "grp7id")]
#-- (2) Estimate Nj via WFPBB
for(s in 1:sim){
  
  # draw from posterior of Nmk
  bootstrap<-BayesianBootstrap(c(1:n),1)
  resample<-rmultinom(1,n,bootstrap)
  bootind <- rep(1:n, resample)
  bootdf <- samp[bootind,]
  boottbl <- bootdf %>%
    group_by(J_cell) %>%
    summarise(bootwts = sum(wts))
  # bootjcellcts <- bootdf %>% group_by(age4,sex2,hispan2,educat4,povgap5,svefreq5) %>% summarise(mk = mean(mk),bootcts=n())
  bootsamp_size <- nrow(boottbl)
  wts_new <- N*normalize(boottbl$bootwts)
  boot_cats <- boottbl$J_cell
  ## find Nj via wfpbb with the bootstrapped categories and new weights
  Nj_hat <- matrix(0, nrow = F_draw, ncol = J)
  bootdf$Njest <- NULL
  srbiwts <- 0
  for(f in 1:F_draw){
    bootdf$Njest <- NULL
    temp <- wtpolyap(ysamp = boot_cats, wts = wts_new, k = (n*Tfact)-bootsamp_size)
    Nj_hat[f, boot_cats] <- table(c(1:J)[temp])
    # boottbl$Njest <- table(temp)
    # bootdf <- bootdf %>% left_join(boottbl, by = "J_cell")
    # srbiwts <- bootdf$Njest/bootdf$srbi_mkcts
    # srbiwts <- (Nj_hat[f,]/srbi_mkcts)[bootdf$J_cell]
    # srbi_overallbb[f] <- sum(normalize(srbiwts)*bootdf$foodindsev,na.rm=T)
    # srbi1bb[f] <- sum(normalize(srbiwts*bootdf$grp1id)*bootdf$foodindsev,na.rm=T)
    # srbi2bb[f] <- sum(normalize(srbiwts*bootdf$grp2id)*bootdf$foodindsev,na.rm=T)
    # srbi3bb[f] <- sum(normalize(srbiwts*bootdf$grp3id)*bootdf$foodindsev,na.rm=T)
    # srbi4bb[f] <- sum(normalize(srbiwts*bootdf$grp4id)*bootdf$foodindsev,na.rm=T)
    # srbi5bb[f] <- sum(normalize(srbiwts*bootdf$grp5id)*bootdf$foodindsev,na.rm=T)
    # srbi6bb[f] <- sum(normalize(srbiwts*bootdf$grp6id)*bootdf$foodindsev,na.rm=T)
    # srbi7bb[f] <- sum(normalize(srbiwts*bootdf$grp7id)*bootdf$foodindsev,na.rm=T)
    srbiwts <- (Nj_hat[f,])
    srbi_overall[s] <- srbi_overall[s]+ sum(normalize(srbiwts)*ybarj,na.rm=T)/F_draw
    srbi1[s] <- srbi1[s]+sum(normalize(srbiwts*grp1id)*ybarj,na.rm=T)/F_draw
    srbi2[s] <- srbi2[s]+sum(normalize(srbiwts*grp2id)*ybarj,na.rm=T)/F_draw
    srbi3[s] <- srbi3[s]+sum(normalize(srbiwts*grp3id)*ybarj,na.rm=T)/F_draw
    srbi4[s] <- srbi4[s]+sum(normalize(srbiwts*grp4id)*ybarj,na.rm=T)/F_draw
    srbi5[s] <- srbi5[s]+sum(normalize(srbiwts*grp5id)*ybarj,na.rm=T)/F_draw
    srbi6[s] <- srbi6[s]+sum(normalize(srbiwts*grp6id)*ybarj,na.rm=T)/F_draw
    srbi7[s] <- srbi7[s]+sum(normalize(srbiwts*grp7id)*ybarj,na.rm=T)/F_draw
    Nj_hat[f, Nj_hat[f,]==0] <- NA
  }
  # boottbl$Njest <- colMeans(Nj_hat)[boot_cats]
  # Njhatmat[s,] <- colMeans(Nj_hat)
  
  # bootdf <- bootdf %>% left_join(boottbl, by = "J_cell")
  # bootdf$srbiwts <- bootdf$Njest / bootdf$srbi_mkcts
  # SRBI estimate
  # srbi_overall[s] <- mean(srbi_overallbb,na.rm=T)
  # srbi1[s] <- mean(srbi1bb,na.rm=T)
  # srbi2[s] <- mean(srbi2bb,na.rm=T)
  # srbi3[s] <- mean(srbi3bb,na.rm=T)
  # srbi4[s] <- mean(srbi4bb,na.rm=T)
  # srbi5[s] <- mean(srbi5bb,na.rm=T)
  # srbi6[s] <- mean(srbi6bb,na.rm=T)
  # srbi7[s] <- mean(srbi7bb,na.rm=T)
  # wtsest[s,] <- colMeans(Nj_hat,na.rm=T)
  print(s)
}

#-- (4) Print results
resmatid <- c("Overall","Group 1","Group 2","Group 3","Group 4","Group 5","Group 6","Group 7")
allgrp_srbi <- cbind(srbi_overall, srbi1, srbi2, srbi3, srbi4,srbi5, srbi6, srbi7)
# allgrp_agency <- cbind(agency_overall, agency1, agency2, agency3, agency4, agency5, agency6,agency7)

resmat_srbi <- resmatfun(allgrp_srbi, id = resmatid)
# resmat_agency <- resmatfun(allgrp_agency, id = resmatid)

write.csv(resmat_srbi,"BASELINE_NOSVEFREQ_wtd.csv")
# write.csv(resmat_agency,"agency_wtd.csv")
