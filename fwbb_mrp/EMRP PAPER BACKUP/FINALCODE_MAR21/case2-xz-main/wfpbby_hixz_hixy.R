#====================================
# Code for WFPBB simulations
#====================================
# last modified: nov. 27, 2020
#---- SETUP -----
## load packages
require(dplyr)
require(polyapost)
require(LaplacesDemon)

#-- load population generating function:
source("genPop_hixz_hixy.R")
`%notin%` <- Negate(`%in%`)

## define seed, number of simulations ('sim'), 
# number of bootstraps ('L'), 
# number of FPBB draws ('F_draw')
set.seed(50)
sim <- 200
L <- 100
F_draw <- 20
# sim <- 2
# L <- 100
# F_draw <- 20

#-- define parameters from population data
true_mean <- mean(df$Y)
true_mean1 <- mean(df$Y[df$subgroup == 1])
true_mean2 <- mean(df$Y[df$subgroup == 2])
true_mean3 <- mean(df$Y[df$subgroup == 3])
true_mean4 <- mean(df$Y[df$subgroup == 4])
true_mean5 <- mean(df$Y[df$subgroup == 5])
true_mean6 <- mean(df$Y[df$subgroup6 == 1])
true_mean7 <- mean(df$Y[df$subgroup7 == 1])
popmeans <- df %>% group_by(J_cell) %>% summarise(mean(Y)) %>% dplyr::select(-J_cell) %>% as.matrix() %>% as.numeric()
#-- initialize containers
y_bar <- y_bar1 <- y_bar2 <- y_bar3 <- y_bar4 <- y_bar5 <- y_bar6 <- y_bar7 <-  rep(NA, sim)
bbmean <- bbmean1 <- bbmean2 <- bbmean3 <- bbmean4 <- bbmean5 <- bbmean6 <- bbmean7 <-  rep(NA,F_draw)
tempmean <- tempmean1 <- tempmean2 <- tempmean3 <- tempmean4 <- tempmean5 <- tempmean6 <- tempmean7 <- rep(NA,L)
CI <-CI1 <- CI2 <-CI3 <- CI4 <-CI5 <-CI6 <- CI7 <- matrix(nrow = sim, ncol = 2)
CI_covered <- CI_covered1 <-CI_covered2 <-CI_covered3 <-CI_covered4 <-CI_covered5 <-CI_covered6 <-CI_covered7 <-rep(NA,sim)

coverage_rate <- 0
bias <- 0
rMSE <- 0

#-- population counts for sampling
popcts <- as.vector(table(df$J_cell))
p_include <- df$p_include
## CELLMEAN summary
cellmeanbias <- cellmeanrmse <- cellmeanse <- rep(0,J)
CI_cellmean <- rep(NA, J)
CIlength_cellmean <- matrix(NA,nrow=sim, ncol=J)
CI_cellmeancov <- matrix(NA, nrow=sim, ncol=J)

## Nj summary
cellmeans_stan <- matrix(NA,nrow=L,ncol=J)
cellnotpicked <- matrix(NA,nrow=sim, ncol=J)
cellpicked <- rep(0, J)

#----- BEGIN SIMULATIONS ----
## beginning of runtime
t0 <- Sys.time()

for(i in 1:sim){
  
  #-- draw sample 
  df$I <- rbinom(n=N, size = 1, prob = p_include)
  
  # -- create df for sampled data
  sampled <- df[df$I==1,]
  
  #-- sample size
  n <- nrow(sampled)
  
  sampledztbl <- sampled %>% group_by(Z_categorical) %>% summarise(wts = mean(zcts)/n())
  sampled_z <- sampledztbl$Z_categorical
  sampled <- left_join(sampled, sampledztbl)
  sampled_J <- unique(sampled$J_cell)
  # cellmeans_stan <- matrix(NA,nrow=L,ncol=J)
  
  #------ DRAW SYNTHETIC POPULATIONS ----
  for(l in 1:L){
    
    ## (1) Bootstrap from parent sample
    bootstrap<-BayesianBootstrap(c(1:n),1)
    resample<-rmultinom(1,n,bootstrap)
    twt<-sampled$wts*resample
    
    unique_ind <- c(1:n)[twt!=0]
    twt <- twt[twt!=0]
    
    ## (2) WFPBB with new weights and unique indices
    # fpbb_temp <- matrix(0, nrow=F_draw, ncol=N)
    fpbb_temp <- rep(0, N)
    for(f in 1:F_draw){
      fpbb_temp <- wtpolyap(ysamp = unique_ind, wts = twt, k = N-length(unique_ind))
      synthpop <- sampled[fpbb_temp,c("subgroup", "Y", "subgroup6", "subgroup7")]
      bbmean[f] <- mean(synthpop$Y)
      bbmean1[f] <- mean(synthpop$Y[synthpop$subgroup==1])
      bbmean2[f] <- mean(synthpop$Y[synthpop$subgroup==2])
      bbmean3[f] <- mean(synthpop$Y[synthpop$subgroup==3])
      bbmean4[f] <- mean(synthpop$Y[synthpop$subgroup==4])
      bbmean5[f] <- mean(synthpop$Y[synthpop$subgroup==5])
      bbmean6[f] <- mean(synthpop$Y[synthpop$subgroup6==1])
      bbmean7[f] <- mean(synthpop$Y[synthpop$subgroup7==1])
    }
    # bb_cell <- unique(temptbl$J_cell)
    
    ## (3) Collect estimate from the lth synthetic population
    tempmean[l] <- mean(bbmean,na.rm=T)
    tempmean1[l] <- mean(bbmean1,na.rm=T)
    tempmean2[l] <- mean(bbmean2,na.rm=T)
    tempmean3[l] <- mean(bbmean3,na.rm=T)
    tempmean4[l] <- mean(bbmean4,na.rm=T)
    tempmean5[l] <- mean(bbmean5,na.rm=T)
    tempmean6[l] <- mean(bbmean6,na.rm=T)
    tempmean7[l] <- mean(bbmean7,na.rm=T)
    
    print(l)
  }  
  ## collect stats from cellmeans
  # #- cell picked?
  cellpicked[sampled_J] <- cellpicked[sampled_J] + 1
  # cellpicked[sampled_J] <- cellpicked[sampled_J] + 1
  #- bias, SE, rMSE
  # cellmeanbias[sampled_J] <- cellmeanbias[sampled_J]+colMeans(cellmeans_stan[,sampled_J],na.rm=T)-popmeans[sampled_J]
  # cellmeanrmse[sampled_J] <- cellmeanrmse[sampled_J] + cellmeanbias[sampled_J]^2
  
  #- coverage
  # CI_cellmean <- apply(cellmeans_stan[,sampled_J], 2, quantile, c(0.025, 0.975), na.rm=T)
  # CIlength_cellmean[i,sampled_J] <- apply(CI_cellmean,2,diff)
  # CI_cellmeancov[i, sampled_J] <- popmeans[sampled_J] >= CI_cellmean[1,] & popmeans[sampled_J] <= CI_cellmean[2,]
  
  ## (4) Final estimate for parent sample i
  y_bar[i] <- mean(tempmean,na.rm=T)
  y_bar1[i] <- mean(tempmean1,na.rm=T)
  y_bar2[i] <- mean(tempmean2,na.rm=T)
  y_bar3[i] <- mean(tempmean3,na.rm=T)
  y_bar4[i] <- mean(tempmean4,na.rm=T)
  y_bar5[i] <- mean(tempmean5,na.rm=T)
  y_bar6[i] <- mean(tempmean6,na.rm=T)
  y_bar7[i] <- mean(tempmean7,na.rm=T)
  
  #groupvar[i] <- var(tempmean)
  
  ## (5) Confidence Interval
  CI[i,] <- quantile(tempmean, probs = c(0.025, 0.975),na.rm=T)
  CI1[i,] <- quantile(tempmean1, probs = c(0.025, 0.975),na.rm=T)
  CI2[i,] <- quantile(tempmean2, probs = c(0.025, 0.975),na.rm=T)
  CI3[i,] <- quantile(tempmean3, probs = c(0.025, 0.975),na.rm=T)
  CI4[i,] <- quantile(tempmean4, probs = c(0.025, 0.975),na.rm=T)
  CI5[i,] <- quantile(tempmean5, probs = c(0.025, 0.975),na.rm=T)
  CI6[i,] <- quantile(tempmean6, probs = c(0.025, 0.975),na.rm=T)
  CI7[i,] <- quantile(tempmean7, probs = c(0.025, 0.975),na.rm=T)
    sampsub <- unique(sampled$subgroup)
  
  sampsub <- unique(sampled$subgroup)
  nosampsub6 <- 1 %notin% sampled$subgroup6
  nosampsub7 <- 1 %notin% sampled$subgroup7
  if(between(true_mean, CI[i,1], CI[i,2])){CI_covered[i] <- 1}else if(sum(is.na(CI[i,])!=0)){CI_covered[i] <- NA}else{CI_covered[i] <- 0}
  if(1 %notin% sampsub){CI_covered1[i] <- NA}else if(between(true_mean1, CI1[i,1], CI1[i,2])){CI_covered1[i] <- 1} else{CI_covered1[i] <- 0}
  if(2 %notin% sampsub){CI_covered2[i] <- NA}else if(between(true_mean2, CI2[i,1], CI2[i,2])){CI_covered2[i] <- 1} else{CI_covered2[i] <- 0}
  if(3 %notin% sampsub){CI_covered3[i] <- NA}else if(between(true_mean3, CI3[i,1], CI3[i,2])){CI_covered3[i] <- 1} else{CI_covered3[i] <- 0}
  if(4 %notin% sampsub){CI_covered4[i] <- NA}else if(between(true_mean4, CI4[i,1], CI4[i,2])){CI_covered4[i] <- 1} else{CI_covered4[i] <- 0}
  if(5 %notin% sampsub){CI_covered5[i] <- NA}else if(between(true_mean5, CI5[i,1], CI5[i,2])){CI_covered5[i] <- 1} else{CI_covered5[i] <- 0}
  if(nosampsub6){CI_covered6[i] <- NA}else if(between(true_mean6, CI6[i,1], CI6[i,2])){CI_covered6[i] <- 1} else{CI_covered6[i] <- 0}
  if(nosampsub7){CI_covered7[i] <- NA}else if(between(true_mean7, CI7[i,1], CI7[i,2])){CI_covered7[i] <- 1} else{CI_covered7[i] <- 0}
  
  
  ## track progress of code
  print(i)
}

## end of runtime
t1 <- Sys.time()

ybarmat <- cbind(y_bar, y_bar1, y_bar2, y_bar3, y_bar4, y_bar5,y_bar6, y_bar7)
truemeanmat <- c(true_mean,true_mean1,true_mean2,true_mean3,true_mean4,true_mean5,true_mean6,true_mean7)
labs <- c("Overall", "group 1", "group2", "group3", "group4", "group5","group6", "group7")
CImat <- list(CI,CI1,CI2,CI3,CI4,CI5,CI6,CI7)
CIcovmat <- cbind(CI_covered,CI_covered1,CI_covered2,CI_covered3,CI_covered4,CI_covered5,CI_covered6,CI_covered7)

resmatfun <- function(ybarmat, truemean, CI, CIcovered,labs){
  errormat <- t(t(ybarmat)-truemean)
  bias <- colMeans(errormat)
  rMSE <- sqrt(colMeans(errormat^2))
  SE <- t(t(ybarmat)-colMeans(ybarmat))
  SE <- sqrt(colMeans(SE^2))
  CI_length <- vector(length = ncol(ybarmat))
  for(i in 1:ncol(ybarmat)){
    CItemp <- CI[[i]]
    CI_length[i] <- mean(CItemp[,2]-CItemp[,1])
  }
  coverage_rate <- colMeans(CIcovered)
  res <-rbind(rMSE, bias, SE, CI_length, coverage_rate)
  colnames(res) <- labs
  return(res)
}
resultmat <- resmatfun(ybarmat = ybarmat, truemean = truemeanmat, CI = CImat, CIcovered = CIcovmat, labs = labs)


## CELLMEAN Summary
# cellmeanbias <- 
#   cellmeanbias / cellpicked 
# cellmeanrmse <- 
#   sqrt(cellmeanrmse/cellpicked)
# cellmeanse <- 
#   sqrt(cellmeanrmse^2 - cellmeanbias^2)

# print("Cell picking frequency in sample")
# print(cellpicked)
# cellpicked[group1]
# cellpicked[group2]
# cellpicked[group3]
# cellpicked[group4]
# cellpicked[group5]
# 
# print("Cell miss frequency in BB")
# print(cellnotpicked)
# cellnotpicked[group1]
# cellnotpicked[group2]
# cellnotpicked[group3]
# cellnotpicked[group4]
# cellnotpicked[group5]

# print("Bias of cellmeans")
# print(cellmeanbias)
# 
# cellmeanbias[group1]
# cellmeanbias[group2]
# cellmeanbias[group3]
# cellmeanbias[group4]
# cellmeanbias[group5]
# 
# print("rMSE of cellmeans")
# print(cellmeanrmse)
# 
# print("SE of cellmeans")
# print(cellmeanse)
# 
# cellmeanse[group1]
# cellmeanse[group2]
# cellmeanse[group3]
# cellmeanse[group4]
# cellmeanse[group5]
# 
# print("CI length of cellmeans")
# colMeans(CI_cellmean, na.rm=T)
# colMeans(CI_cellmean, na.rm=T)[group1]
# colMeans(CI_cellmean, na.rm=T)[group2]
# colMeans(CI_cellmean, na.rm=T)[group3]
# colMeans(CI_cellmean, na.rm=T)[group4]
# colMeans(CI_cellmean, na.rm=T)[group5]
# 
# print("CI coverage rate of cellmeans")
# colMeans(CI_cellmeancov,na.rm=T)
# colMeans(CI_cellmeancov,na.rm=T)[group1]
# colMeans(CI_cellmeancov,na.rm=T)[group2]
# colMeans(CI_cellmeancov,na.rm=T)[group3]
# colMeans(CI_cellmeancov,na.rm=T)[group4]
# colMeans(CI_cellmeancov,na.rm=T)[group5]
# 
# # PRINT RESULTS
resultmat

# WRITE RESULTS
write.csv(resultmat, "res_wfpbby.csv")
# write.csv(cellmeanbias, "cellmeanbias_wfpbby.csv")
# write.csv(cellmeanse, "cellmeanse_wfpbby.csv")
# write.csv(cellmeanrmse, "cellmeanrmse_wfpbby.csv")
# write.csv(CI_cellmean, "cellmeanCI_wfpbby.csv")
# write.csv(CI_cellmeancov, "cellmeancov_wfpbby.csv")
