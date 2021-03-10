#================================
# Embedded MRP with applications
#================================
#-- last updated: jan 6, 2021 --
#-- this version uses SRBI data only

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
source("robinhood_cleaned_jan20.R")
samp <- samp[samp$samptype == "SRBI",]

# need to relabel XxZ combinations s.t. they are unique to SRBI (as opposed to SRBI + Agency)
jtbl <- samp %>% group_by(J_cell) %>% summarise(cts = n()) %>% dplyr::select(-cts)
jtbl$Jcellsrbi <- 1:nrow(jtbl)
samp <- left_join(samp, jtbl)
samp <- samp %>% mutate(J_cell = Jcellsrbi) %>% dplyr::select(-Jcellsrbi)
n <- nrow(samp)
Tfact <- 30
grptbl <- samp %>% group_by(age3,sex2,race3,educat4,povgap4,svefreq2, aecat) %>% summarise(
  grp1id = mean(age3==1),#NEW
  grp2id = mean(educat4 %in% c(1,2)), #NEW
  grp3id = mean(povgap4==1),
  grp4id = mean(race3==2),
  grp5id = mean(povgap4 %in% c(1,2) & sex2==1),
  grp6id = mean(svefreq2==1),
  grp7id = mean(sex2==2 & svefreq2==1))
J <- nrow(grptbl)
M <- length(unique(samp$sampled_x1_label))
grp1id <- grptbl$grp1id
grp2id <- grptbl$grp2id
grp3id <- grptbl$grp3id
grp4id <- grptbl$grp4id
grp5id <- grptbl$grp5id
grp6id <- grptbl$grp6id
grp7id <- grptbl$grp7id
samp <- left_join(samp, grptbl)

# initialize containers for results
sim <- 1000
staniters = 10000
# sim <- 10
# staniters = 1000
F_draw = 20
wmrp_Njhatmat <- matrix(0,nrow=sim, ncol=J)
mmrp_Njhatmat <- matrix(0,nrow=sim, ncol=J)
mrp2_Njhatmat <- matrix(0,nrow=sim, ncol=J)

#== WFPBB-MRP
wmrp_posterior <- rep(0,sim)
wmrp_posterior1 <- rep(0, sim)
wmrp_posterior2 <- rep(0, sim)
wmrp_posterior3 <- rep(0, sim)
wmrp_posterior4 <- rep(0, sim)
wmrp_posterior5 <- rep(0, sim)
wmrp_posterior6 <- rep(0, sim)
wmrp_posterior7 <- rep(0, sim)

#== Multinomial-MRP
mmrp_posterior <- mmrp_posterior1 <- mmrp_posterior2 <- mmrp_posterior3 <- mmrp_posterior4 <- mmrp_posterior5 <- mmrp_posterior6 <- mmrp_posterior7 <- rep(0,sim)

#== Two-stage MRP
mrp2_posterior <- mrp2_posterior1 <- mrp2_posterior2 <- mrp2_posterior3 <- mrp2_posterior4 <- mrp2_posterior5 <- mrp2_posterior6 <- mrp2_posterior7 <- rep(0,sim)

#===== (0) Regress X~Z for two-stage MRP====
# parameters for Nj stan model
Ma <- nlevels(samp$age3) # combinations of x1: age, sex, educat, povgap
Mb <- nlevels(samp$sex2)
Mc <- nlevels(samp$race3)
Md <- nlevels(samp$educat4)
Me <- nlevels(samp$povgap4)
L <- nlevels(samp$svefreq2)
Mae <- nlevels(samp$aecat)
n <- nrow(samp)
y <- as.numeric(samp$svefreq2)-1
xj <- samp$J_cell
grptbl <- samp %>%
  group_by(age3,sex2,race3,educat4,povgap4,svefreq2, aecat,x1_popcts) %>%
  summarise(zpopcts = mean(x1_popcts))
J_stan <- nrow(grptbl)
agroup <- as.numeric(grptbl$age3)
bgroup <- as.numeric(grptbl$sex2) -1
cgroup <- as.numeric(grptbl$race3)
dgroup <- as.numeric(grptbl$educat4)
egroup <- as.numeric(grptbl$povgap4)
aegroup <- as.numeric(grptbl$aecat)
xgroup <- as.numeric(grptbl$svefreq2)-1
zpopcts <- grptbl$zpopcts

stanfitNq <- stan(file = "mrp2_Nq_sim.stan",
                data = c("Ma", "Mb", "Mc", "Md", "Me", "L", "Mae", "J_stan", "n", "y", "xj",
                         "agroup", "bgroup", "cgroup", "dgroup", "egroup", "aegroup","xgroup"),
                iter=staniters,warmup=staniters-sim,control=list(adapt_delta=0.99),chains=4)
nq_pars <- 
  extract(stanfitNq,permuted = TRUE, inc_warmup = FALSE,include = TRUE)

# extract cell frequency estimates from MRP
nq_draws <- nq_pars$propmk 

#===== (1) Run multilevel regression via Stan=====
# parameters for stan model
Ma <- nlevels(samp$age3) # combinations of x1: age, sex, educat, povgap
Mb <- nlevels(samp$sex2)
Mc <- nlevels(samp$race3)
Md <- nlevels(samp$educat4)
Me <- nlevels(samp$povgap4)
L <- nlevels(samp$svefreq2)
Mae <- nlevels(samp$aecat)
n <- nrow(samp)

y <- samp$foodindsev

grptbl <- samp %>%
  group_by(age3,sex2,race3,educat4,povgap4,svefreq2, aecat) %>%
  summarise(n())
grptbl$sampJ <- 1:nrow(grptbl)
samp <- left_join(samp,grptbl)
J_stan <- nrow(grptbl) # if fitting with sampgroup
agroup <- as.numeric(grptbl$age3)
bgroup <- as.numeric(grptbl$sex2) -1
cgroup <- as.numeric(grptbl$race3) 
dgroup <- as.numeric(grptbl$educat4)
egroup <- as.numeric(grptbl$povgap4)
xgroup <- as.numeric(grptbl$svefreq2)-1
aegroup <- as.numeric(grptbl$aecat)
# sampgroup <- as.numeric(grptbl$samptype) -1 # for testing samptype sig
cell_label <- samp$sampJ
grp1id <- grp1id[sort(unique(samp$J_cell))]
grp2id <- grp2id[sort(unique(samp$J_cell))]
grp3id <- grp3id[sort(unique(samp$J_cell))]
grp4id <- grp4id[sort(unique(samp$J_cell))]
grp5id <- grp5id[sort(unique(samp$J_cell))]
grp6id <- grp6id[sort(unique(samp$J_cell))]
grp7id <- grp7id[sort(unique(samp$J_cell))]
# xsampgroup <- as.numeric(grptbl$samptype)-1 # for testing samptype sig

#******** use only if testing significance of sample type (Agency or SRBI) ****************
# stanfit <- stan(file = "outcome_model_samp.stan",
#                 data = c("Ma", "Mb", "Mc", "Md", "Me", "L", "Mae", "J_stan", "n", "y", "agroup", "bgroup",
#                          "cgroup", "dgroup", "egroup", "xgroup", "aegroup",
#                          "sampgroup", "cell_label"),
#                 iter=staniters,warmup=staniters-sim,control=list(adapt_delta=0.99),chains=1)
#******************************************************************************************
stanfit <- stan(file = "outcome_model_samp.stan",
                data = c("Ma", "Mb", "Mc", "Md", "Me", "L", "Mae", "J_stan", "n", "y", "agroup", "bgroup",
                         "cgroup", "dgroup", "egroup", "xgroup", "aegroup",
                         "cell_label"),
                iter=staniters,warmup=staniters-sim,control=list(adapt_delta=0.99),chains=4)

# extract estimates of cell means from the stan model
stanpars <-
  extract(stanfit,permuted = TRUE, inc_warmup = FALSE,include = TRUE)
cellmeans_stan <- stanpars$cellmean

#===== (2) Estimate Nj via WFPBB=====
for(s in 1:sim){
  
  # draw from posterior of Nmk
  #- using complete sample
  bootstrap<-BayesianBootstrap(c(1:n),1)
  resample<-rmultinom(1,n,bootstrap)
  bootind <- rep(1:n, resample)
  bootdf <- samp[bootind,] %>% dplyr::select(J_cell, wts, samptype)
  
  #-- WFPBB: Find N_{mk} (joint dist of X1 x svefreq2) from SRBI sample only
    #- although we exclude Agency obs that are not in Z category from SRBI,
    #- we still estimate N_mk if there is Agency obs with mk category unique to Agency
  boottbl <- bootdf %>% 
    group_by(J_cell) %>% 
    summarise(bootwts = sum(wts))
  # bootjcellcts <- bootdf %>% group_by(J_cell) %>% summarise(mk = mean(mk),bootcts=n())
  bootsamp_size <- nrow(boottbl)
  wts_new <- boottbl$bootwts
  boot_cats <- boottbl$J_cell
  ## find Nj via wfpbb with the bootstrapped categories and new weights
  Nmk_hat <- matrix(0, nrow = F_draw, ncol = J)
  
  for(f in 1:F_draw){
    # temp <- wtpolyap(ysamp = boot_cats, wts = wts_new, k = N-bootsamp_size)
    temp <- wtpolyap(ysamp = boot_cats, wts = wts_new, k = (n*Tfact)-bootsamp_size)
    Nmk_hat[f, boot_cats] <- table(temp)
  }
  Nmk_est <- colMeans(Nmk_hat)
  current_cellmean <- cellmeans_stan[s,]
  
  #************************************************************
  #-- WFPBB-MRP
  wmrp_Njhatmat[s,] <- Nmk_est
  wts_wmrp <- normalize(Nmk_est)
  wts_wmrp1<- normalize(wts_wmrp*grp1id)
  wts_wmrp2<- normalize(wts_wmrp*grp2id)
  wts_wmrp3<- normalize(wts_wmrp*grp3id)
  wts_wmrp4<- normalize(wts_wmrp*grp4id)
  wts_wmrp5<- normalize(wts_wmrp*grp5id)
  wts_wmrp6<- normalize(wts_wmrp*grp6id)
  wts_wmrp7<- normalize(wts_wmrp*grp7id)
  
  wmrp_posterior[s] <- crossprod(wts_wmrp, current_cellmean)
  wmrp_posterior1[s] <- crossprod(wts_wmrp1, current_cellmean)
  wmrp_posterior2[s] <- crossprod(wts_wmrp2, current_cellmean)
  wmrp_posterior3[s] <- crossprod(wts_wmrp3, current_cellmean)
  wmrp_posterior4[s] <- crossprod(wts_wmrp4, current_cellmean)
  wmrp_posterior5[s] <- crossprod(wts_wmrp5, current_cellmean)
  wmrp_posterior6[s] <- crossprod(wts_wmrp6, current_cellmean)
  wmrp_posterior7[s] <- crossprod(wts_wmrp7, current_cellmean)
  #************************************************************
  #-- 2-stage MRP
  mrp2_Njhatmat[s,] <- nq_draws[s,] * zpopcts
  
  wts_mrp2 <- normalize(nq_draws[s,] * zpopcts)
  wts_mrp21 <- normalize(wts_mrp2*grp1id)
  wts_mrp22 <- normalize(wts_mrp2*grp2id)
  wts_mrp23 <- normalize(wts_mrp2*grp3id)
  wts_mrp24 <- normalize(wts_mrp2*grp4id) 
  wts_mrp25 <- normalize(wts_mrp2*grp5id)
  wts_mrp26 <- normalize(wts_mrp2*grp6id) 
  wts_mrp27 <- normalize(wts_mrp2*grp7id)
  
  mrp2_posterior[s] <- crossprod(wts_mrp2, current_cellmean)
  mrp2_posterior1[s] <- crossprod(wts_mrp21, current_cellmean)
  mrp2_posterior2[s] <- crossprod(wts_mrp22, current_cellmean)
  mrp2_posterior3[s] <- crossprod(wts_mrp23, current_cellmean)
  mrp2_posterior4[s] <- crossprod(wts_mrp24, current_cellmean)
  mrp2_posterior5[s] <- crossprod(wts_mrp25, current_cellmean)
  mrp2_posterior6[s] <- crossprod(wts_mrp26, current_cellmean)
  mrp2_posterior7[s] <- crossprod(wts_mrp27, current_cellmean)
  #************************************************************
  #-- Multinomial-MRP
  fulltbl <- samp %>% group_by(J_cell,sampled_x1_label) %>% summarise(prob = n()/mean(x1_sampcts), x1popct = mean(x1_popcts))
  nmk_tbl <- fulltbl[,c("J_cell","sampled_x1_label")]
  nmk_tbl$wts <- 0
  for(m in 1:M){
    if(sum(fulltbl$sampled_x1_label==m) ==1){
      nmk_tbl[nmk_tbl$sampled_x1_label==m, "wts"] <- fulltbl$x1popct[fulltbl$sampled_x1_label==m]
    }
    else{
      probs <- fulltbl[fulltbl$sampled_x1_label==m, "prob"] %>% unlist() %>% as.double()
      popcts <- unique(fulltbl$x1popct[fulltbl$sampled_x1_label==m])
      nmk_tbl[nmk_tbl$sampled_x1_label == m, "wts"] <-
        rowMeans(rmultinom(popcts, 1, prob = probs))*popcts
    }
  }
  mmrp_Njhatmat[s,] <- nmk_tbl$wts
  wts_mmrp <- normalize(nmk_tbl$wts)
  wts_mmrp1<- normalize(wts_mmrp*grp1id)
  wts_mmrp2<- normalize(wts_mmrp*grp2id)
  wts_mmrp3<- normalize(wts_mmrp*grp3id)
  wts_mmrp4<- normalize(wts_mmrp*grp4id)
  wts_mmrp5<- normalize(wts_mmrp*grp5id)
  wts_mmrp6<- normalize(wts_mmrp*grp6id)
  wts_mmrp7<- normalize(wts_mmrp*grp7id)
  
  mmrp_posterior[s] <- crossprod(wts_mmrp, current_cellmean)
  mmrp_posterior1[s] <- crossprod(wts_mmrp1, current_cellmean)
  mmrp_posterior2[s] <- crossprod(wts_mmrp2, current_cellmean)
  mmrp_posterior3[s] <- crossprod(wts_mmrp3, current_cellmean)
  mmrp_posterior4[s] <- crossprod(wts_mmrp4, current_cellmean)
  mmrp_posterior5[s] <- crossprod(wts_mmrp5, current_cellmean)
  mmrp_posterior6[s] <- crossprod(wts_mmrp6, current_cellmean)
  mmrp_posterior7[s] <- crossprod(wts_mmrp7, current_cellmean)
  print(s)
}

#-- (4) Print results
resmatid <- c("Overall", "Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7")
allgrp_wmrp <- cbind(wmrp_posterior, wmrp_posterior1, wmrp_posterior2, wmrp_posterior3, wmrp_posterior4, wmrp_posterior5, wmrp_posterior6, wmrp_posterior7)
allgrp_mmrp <- cbind(mmrp_posterior, mmrp_posterior1, mmrp_posterior2, mmrp_posterior3, mmrp_posterior4, mmrp_posterior5, mmrp_posterior6, mmrp_posterior7)
allgrp_mrp2 <- cbind(mrp2_posterior, mrp2_posterior1, mrp2_posterior2, mrp2_posterior3, mrp2_posterior4, mrp2_posterior5, mrp2_posterior6, mrp2_posterior7)

resmat_wmrp <- resmatfun(allgrp_wmrp,resmatid)
resmat_mmrp <- resmatfun(allgrp_mmrp,resmatid)
resmat_mrp2 <- resmatfun(allgrp_mrp2,resmatid)
njmat_wmrp <- resmatfun(wmrp_Njhatmat)
njmat_mmrp <- resmatfun(mmrp_Njhatmat)
njmat_mrp2 <- resmatfun(mrp2_Njhatmat)
write.csv(resmat_wmrp, "srbi_wfpbbmrp.csv")
write.csv(resmat_mmrp, "srbi_mmrp.csv")
write.csv(resmat_mrp2, "srbi_mrp2.csv")
write.csv(njmat_wmrp, "srbi_wfpbbmrp_nj.csv")
write.csv(njmat_mmrp, "srbi_mmrp_nj.csv")
write.csv(njmat_mrp2, "srbi_mrp2_nj.csv")
