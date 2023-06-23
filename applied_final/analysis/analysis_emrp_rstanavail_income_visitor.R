#================================
# Embedded MRP with applications
#================================
#-- uses public use baseline data

### Code for reworking results with different cmat but after we have already fitted the model
type = "income_visitor"
set.seed(51)
# options(mc.cores = 4)
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
source(paste0("BASELINE_", type, ".R"))

# initialize containers for results
sim <- 5000
staniters = 10000
Tfact <- 30
F_draw = 20
M <- length(unique(samp$sampled_x1_label))
N <- sum(acs$perwt)
wmrp_Njhatmat <- matrix(0,nrow=sim, ncol=J)
mmrp_Njhatmat <- matrix(0,nrow=sim, ncol=J)
mrp2_Njhatmat <- matrix(0,nrow=sim, ncol=J)
#== WFPBB-MRP
wmrp_posterior <- wmrp_posterior1 <- wmrp_posterior2 <- wmrp_posterior3 <- wmrp_posterior4 <- rep(0,sim)
wmrp_posterior5 <- wmrp_posterior6 <- wmrp_posterior7 <- wmrp_posterior8 <- wmrp_posterior9 <- wmrp_posterior10 <- rep(0,sim)

#== Multinomial-MRP
mmrp_posterior <- mmrp_posterior1 <- mmrp_posterior2 <- mmrp_posterior3 <- mmrp_posterior4 <- rep(0,sim)
mmrp_posterior5 <- mmrp_posterior6 <- mmrp_posterior7 <- mmrp_posterior8 <- mmrp_posterior9 <- mmrp_posterior10<- rep(0,sim)

#== Two-stage MRP
mrp2_posterior <- mrp2_posterior1 <- mrp2_posterior2 <- mrp2_posterior3 <- mrp2_posterior4 <- rep(0,sim)
mrp2_posterior5 <- mrp2_posterior6 <- mrp2_posterior7 <- mrp2_posterior8 <- mrp2_posterior9 <- mrp2_posterior10<- rep(0,sim)
#===== (0) Regress X~Z for two-stage MRP====
# parameters for Nj stan model
Ma <- nlevels(samp$age3) 
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
y <- as.numeric(samp$svefreq2)-1
xj <- samp$J_cell
grptbl <- samp %>%
  group_by(age3,sex2,race3,educat4,povgap4,svefreq2,abcat,bccat,cdcat,dfcat) %>%
  summarise(zpopcts = mean(x1_popcts))
J_stan <- nrow(grptbl)
agroup <- as.numeric(grptbl$age3)
bgroup <- as.numeric(grptbl$sex2) -1
cgroup <- as.numeric(grptbl$race3)
dgroup <- as.numeric(grptbl$educat4)
egroup <- as.numeric(grptbl$povgap4)
xgroup <- as.numeric(grptbl$svefreq2)-1
abgroup <- as.numeric(grptbl$abcat)
bcgroup <- as.numeric(grptbl$bccat)
cdgroup <- as.numeric(grptbl$cdcat)
dxgroup <- as.numeric(grptbl$dfcat)
zpopcts <- grptbl$zpopcts

nq_draws <- readRDS("nq_stan.rds")
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
  group_by(age3,sex2,race3,educat4,povgap4,svefreq2,abcat,bccat,cdcat,dfcat) %>%
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
abgroup <- as.numeric(grptbl$abcat)
bcgroup <- as.numeric(grptbl$bccat)
cdgroup <- as.numeric(grptbl$cdcat)
dxgroup <- as.numeric(grptbl$dfcat)
cell_label <- samp$sampJ
grp1id <- grp1id[sort(unique(samp$J_cell))]
grp2id <- grp2id[sort(unique(samp$J_cell))]
grp3id <- grp3id[sort(unique(samp$J_cell))]
grp4id <- grp4id[sort(unique(samp$J_cell))]
cellmeans_stan <- readRDS("emrp_stan.rds")
cellmeans_stan <- cellmeans_stan$cellmean

#===== (2) Estimate Nj via WFPBB=====
for(s in 1:sim){
  
  #- using complete sample
  samp$repwts <- NULL
  bootstrap<-BayesianBootstrap(c(1:n),1)
  samp$repwts <- n*normalize(BayesianBootstrap(c(1:n),1))*samp$wts
  
  boottbl <-
    samp %>%
    group_by(J_cell) %>%
    summarise(bootwts = sum(repwts))
  nonzeroind <- boottbl$bootwts != 0
  boot_cats <- boottbl$J_cell[nonzeroind]
  bootsamp_size <- sum(nonzeroind)
  wts_new <- (n*Tfact)*normalize(boottbl$bootwts)
  boot_cats <- boottbl$J_cell
  
  ## find Nj via wfpbb with the bootstrapped categories and new weights
  Nmk_hat <- matrix(0, nrow = F_draw, ncol = J)
  
  for(f in 1:F_draw){
    temp <- wtpolyap(ysamp = boot_cats, wts = wts_new, k = (n*Tfact)-bootsamp_size)
    Nmk_hat[f, boot_cats] <- table(c(1:J)[temp])
  }
  Nmk_est <- colMeans(Nmk_hat,na.rm=T)
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
  wts_wmrp8<- normalize(wts_wmrp*grp8id)
  wts_wmrp9<- normalize(wts_wmrp*grp9id)
  wts_wmrp10<- normalize(wts_wmrp*grp10id)
  
  wmrp_posterior[s] <- crossprod(wts_wmrp, current_cellmean)
  wmrp_posterior1[s] <- crossprod(wts_wmrp1, current_cellmean)
  wmrp_posterior2[s] <- crossprod(wts_wmrp2, current_cellmean)
  wmrp_posterior3[s] <- crossprod(wts_wmrp3, current_cellmean)
  wmrp_posterior4[s] <- crossprod(wts_wmrp4, current_cellmean)
  wmrp_posterior5[s] <- crossprod(wts_wmrp5, current_cellmean)
  wmrp_posterior6[s] <- crossprod(wts_wmrp6, current_cellmean)
  wmrp_posterior7[s] <- crossprod(wts_wmrp7, current_cellmean)
  wmrp_posterior8[s] <- crossprod(wts_wmrp8, current_cellmean)
  wmrp_posterior9[s] <- crossprod(wts_wmrp9, current_cellmean)
  wmrp_posterior10[s] <- crossprod(wts_wmrp10, current_cellmean)
  
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
  wts_mrp28 <- normalize(wts_mrp2*grp8id) 
  wts_mrp29 <- normalize(wts_mrp2*grp9id)
  wts_mrp210 <- normalize(wts_mrp2*grp10id) 
  
  mrp2_posterior[s] <- crossprod(wts_mrp2, current_cellmean)
  mrp2_posterior1[s] <- crossprod(wts_mrp21, current_cellmean)
  mrp2_posterior2[s] <- crossprod(wts_mrp22, current_cellmean)
  mrp2_posterior3[s] <- crossprod(wts_mrp23, current_cellmean)
  mrp2_posterior4[s] <- crossprod(wts_mrp24, current_cellmean)
  mrp2_posterior5[s] <- crossprod(wts_mrp25, current_cellmean)
  mrp2_posterior6[s] <- crossprod(wts_mrp26, current_cellmean)
  mrp2_posterior7[s] <- crossprod(wts_mrp27, current_cellmean)
  mrp2_posterior8[s] <- crossprod(wts_mrp28, current_cellmean)
  mrp2_posterior9[s] <- crossprod(wts_mrp29, current_cellmean)
  mrp2_posterior10[s] <- crossprod(wts_mrp210, current_cellmean)
  
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
        rmultinom(1,popcts,prob = probs)
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
  wts_mmrp8<- normalize(wts_mmrp*grp8id)
  wts_mmrp9<- normalize(wts_mmrp*grp9id)
  wts_mmrp10<- normalize(wts_mmrp*grp10id)
  
  mmrp_posterior[s] <- crossprod(wts_mmrp, current_cellmean)
  mmrp_posterior1[s] <- crossprod(wts_mmrp1, current_cellmean)
  mmrp_posterior2[s] <- crossprod(wts_mmrp2, current_cellmean)
  mmrp_posterior3[s] <- crossprod(wts_mmrp3, current_cellmean)
  mmrp_posterior4[s] <- crossprod(wts_mmrp4, current_cellmean)
  mmrp_posterior5[s] <- crossprod(wts_mmrp5, current_cellmean)
  mmrp_posterior6[s] <- crossprod(wts_mmrp6, current_cellmean)
  mmrp_posterior7[s] <- crossprod(wts_mmrp7, current_cellmean)
  mmrp_posterior8[s] <- crossprod(wts_mmrp8, current_cellmean)
  mmrp_posterior9[s] <- crossprod(wts_mmrp9, current_cellmean)
  mmrp_posterior10[s] <- crossprod(wts_mmrp10, current_cellmean)
  
  print(s)
}

#-- (4) Print results
resmatid <- c("Overall",grplabs)
allgrp_wmrp <- cbind(wmrp_posterior, wmrp_posterior1, wmrp_posterior2, wmrp_posterior3, wmrp_posterior4,
                     wmrp_posterior5, wmrp_posterior6, wmrp_posterior7, wmrp_posterior8, wmrp_posterior9, wmrp_posterior10)
allgrp_mmrp <- cbind(mmrp_posterior, mmrp_posterior1, mmrp_posterior2, mmrp_posterior3, mmrp_posterior4,
                     mmrp_posterior5, mmrp_posterior6, mmrp_posterior7, mmrp_posterior8, mmrp_posterior9, mmrp_posterior10)
allgrp_mrp2 <- cbind(mrp2_posterior, mrp2_posterior1, mrp2_posterior2, mrp2_posterior3, mrp2_posterior4,
                     mrp2_posterior5, mrp2_posterior6, mrp2_posterior7, mrp2_posterior8, mrp2_posterior9, mrp2_posterior10)

resmat_wmrp <- resmatfun(allgrp_wmrp,resmatid)
resmat_mmrp <- resmatfun(allgrp_mmrp,resmatid)
resmat_mrp2 <- resmatfun(allgrp_mrp2,resmatid)
njmat_wmrp <- resmatfun(wmrp_Njhatmat)
njmat_mmrp <- resmatfun(mmrp_Njhatmat)
njmat_mrp2 <- resmatfun(mrp2_Njhatmat)
write.csv(resmat_wmrp, paste0("results/res_", type, "_wfpbbmrp.csv"))
write.csv(resmat_mmrp, paste0("results/res_", type, "_mmrp.csv"))
write.csv(resmat_mrp2, paste0("results/res_", type,"_mrp2.csv"))
write.csv(njmat_wmrp, paste0("results/res_", type, "wfpbbmrp_nj.csv"))
write.csv(njmat_mmrp, paste0("results/res_", type,"mmrp_nj.csv"))
write.csv(njmat_mrp2, paste0("results/res_", type, "mrp2_nj.csv"))
