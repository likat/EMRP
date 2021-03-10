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
source("robinhood_cleaned_jan20.R")

# get the joint sample frequency of X1 x X2
svywts_srbi <- samp %>% 
  filter(samptype=="SRBI") %>% 
  group_by(J_cell, samptype) %>% 
  summarise(srbi_mkcts = n())
svywts_agency <- samp %>% 
  filter(samptype=="Agency") %>% 
  group_by(J_cell, samptype) %>% 
  summarise(agency_mkcts = n())
agencycells <- svywts_agency$J_cell
samp <- samp %>% left_join(svywts_srbi) %>% left_join(svywts_agency)


n <- nrow(samp)

# initialize containers for results
sim <- 1000
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

agency_overall <- rep(0,sim)
agency1 <- rep(0, sim)
agency2 <- rep(0, sim)
agency3 <- rep(0, sim)
agency4 <- rep(0, sim)
agency5 <- rep(0, sim)
agency6 <- rep(0, sim)
agency7 <- rep(0, sim)
agency_overallbb <- rep(0,F_draw)
agency1bb <- rep(0, F_draw)
agency2bb <- rep(0, F_draw)
agency3bb <- rep(0, F_draw)
agency4bb <- rep(0, F_draw)
agency5bb <- rep(0, F_draw)
agency6bb <- rep(0, F_draw)
agency7bb <- rep(0, F_draw)

samp <- samp[,c("J_cell", "samptype","wts", "foodindsev", "srbi_mkcts", "agency_mkcts",
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
  wts_new <- boottbl$bootwts
  boot_cats <- boottbl$J_cell
  ## find Nj via wfpbb with the bootstrapped categories and new weights
  Nj_hat <- matrix(0, nrow = F_draw, ncol = J)
  bootdf$Njest <- NULL
  for(f in 1:F_draw){
    bootdf$Njest <- NULL
    temp <- wtpolyap(ysamp = boot_cats, wts = wts_new, k = (n*Tfact)-bootsamp_size)
    # Nj_hat[f, boot_cats] <- table(temp)
    boottbl$Njest <- table(temp)
    bootdf <- bootdf %>% left_join(boottbl, by = "J_cell")
    srbiwts <- bootdf$Njest/bootdf$srbi_mkcts
    agencywts <- bootdf$Njest/bootdf$agency_mkcts
    srbi_overallbb[f] <- sum(normalize(srbiwts)*bootdf$foodindsev,na.rm=T)
    srbi1bb[f] <- sum(normalize(srbiwts*bootdf$grp1id)*bootdf$foodindsev,na.rm=T)
    srbi2bb[f] <- sum(normalize(srbiwts*bootdf$grp2id)*bootdf$foodindsev,na.rm=T)
    srbi3bb[f] <- sum(normalize(srbiwts*bootdf$grp3id)*bootdf$foodindsev,na.rm=T)
    srbi4bb[f] <- sum(normalize(srbiwts*bootdf$grp4id)*bootdf$foodindsev,na.rm=T)
    srbi5bb[f] <- sum(normalize(srbiwts*bootdf$grp5id)*bootdf$foodindsev,na.rm=T)
    srbi6bb[f] <- sum(normalize(srbiwts*bootdf$grp6id)*bootdf$foodindsev,na.rm=T)
    srbi7bb[f] <- sum(normalize(srbiwts*bootdf$grp7id)*bootdf$foodindsev,na.rm=T)
    agency_overallbb[f] <- sum(normalize(agencywts)*bootdf$foodindsev,na.rm=T)
    agency1bb[f] <- sum(normalize(agencywts*bootdf$grp1id)*bootdf$foodindsev,na.rm=T)
    agency2bb[f] <- sum(normalize(agencywts*bootdf$grp2id)*bootdf$foodindsev,na.rm=T)
    agency3bb[f] <- sum(normalize(agencywts*bootdf$grp3id)*bootdf$foodindsev,na.rm=T)
    agency4bb[f] <- sum(normalize(agencywts*bootdf$grp4id)*bootdf$foodindsev,na.rm=T)
    agency5bb[f] <- sum(normalize(agencywts*bootdf$grp5id)*bootdf$foodindsev,na.rm=T)
    agency6bb[f] <- sum(normalize(agencywts*bootdf$grp6id)*bootdf$foodindsev,na.rm=T)
    agency7bb[f] <- sum(normalize(agencywts*bootdf$grp7id)*bootdf$foodindsev,na.rm=T)
    
  }
  # boottbl$Njest <- colMeans(Nj_hat)[boot_cats]
  # Njhatmat[s,] <- colMeans(Nj_hat)
  # 
  # bootdf <- bootdf %>% left_join(boottbl, by = "J_cell")
  # bootdf$srbiwts <- bootdf$Njest / bootdf$srbi_mkcts
  # bootdf$agencywts <- bootdf$Njest / bootdf$agency_mkcts
  ## SRBI estimate
  srbi_overall[s] <- mean(srbi_overallbb,na.rm=T)
  srbi1[s] <- mean(srbi1bb,na.rm=T)
  srbi2[s] <- mean(srbi2bb,na.rm=T)
  srbi3[s] <- mean(srbi3bb,na.rm=T)
  srbi4[s] <- mean(srbi4bb,na.rm=T)
  srbi5[s] <- mean(srbi5bb,na.rm=T)
  srbi6[s] <- mean(srbi6bb,na.rm=T)
  srbi7[s] <- mean(srbi7bb,na.rm=T)
  
  ## Agency estimate
  agency_overall[s] <- mean(agency_overallbb,na.rm=T)
  agency1[s] <- mean(agency1bb,na.rm=T)
  agency2[s] <- mean(agency2bb,na.rm=T)
  agency3[s] <- mean(agency3bb,na.rm=T)
  agency4[s] <- mean(agency4bb,na.rm=T)
  agency5[s] <- mean(agency5bb,na.rm=T)
  agency6[s] <- mean(agency6bb,na.rm=T)
  agency7[s] <- mean(agency7bb,na.rm=T)
  
}

#-- (4) Print results
resmatid <- c("Overall","Group 1","Group 2","Group 3","Group 4","Group 5","Group 6","Group 7")
allgrp_srbi <- cbind(srbi_overall, srbi1, srbi2, srbi3, srbi4,srbi5, srbi6, srbi7)
allgrp_agency <- cbind(agency_overall, agency1, agency2, agency3, agency4, agency5, agency6,agency7)

resmat_srbi <- resmatfun(allgrp_srbi, id = resmatid)
resmat_agency <- resmatfun(allgrp_agency, id = resmatid)

write.csv(resmat_srbi,"srbi_wtd.csv")
write.csv(resmat_agency,"agency_wtd.csv")
