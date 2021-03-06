#====================================
# Code for WFPBB-MRP simulations
#====================================
# last modified: dec. 24, 2020
set.seed(50)

#-- load required packages
require(dplyr)
require(MASS)
require(polyapost)
require(rstan)
require(LaplacesDemon)
options(mc.cores =2)

#-- load population generating function:
source("genPop_hixz_hixy.R")
# source("genPop_hixz_loxy.R")
# source("genPop_loxz_hixy.R")
# source("genPop_loxz_loxy.R")

iter = 200 # how many times to run simulation?
L = 100   # number of bootstraps ('L'), 
F_draw = 20  # number of FPBB draws ('F_draw')
sim = 1000 # how many draws from posterior?
staniters = 2000

#-- random integers for naming the stanfit to avoid segfault issues
stanfitrand <- sample(1:iter, iter, replace=F)
modelname <- paste("stanfit",stanfitrand, sep="")

#-- helper fxns
## normalize vectors
normalize <- function(vec){
  return(vec/sum(vec))
}

## boolean operator: are elements in 'vector' not in 'othervector'?
`%notin%` <- Negate(`%in%`)

# population size, number cells
N <- nrow(df)

# what is mean(Y) for the subgroup?
true_mean <- mean(df$Y)
true_mean1 <- mean(df$Y[df$subgroup == 1])
true_mean2 <- mean(df$Y[df$subgroup == 2])
true_mean3 <- mean(df$Y[df$subgroup == 3])
true_mean4 <- mean(df$Y[df$subgroup == 4])
true_mean5 <- mean(df$Y[df$subgroup == 5])
true_mean6 <- mean(df$Y[df$subgroup6 == 1])
true_mean7 <- mean(df$Y[df$subgroup7 == 1])

true_wts <- df %>% group_by(J_cell) %>% summarise(popcts = mean(Nj),wts = mean(Nj)/N, cellmeans = mean(Y))

popmeans <- true_wts$cellmeans
popcts <- true_wts$popcts
true_wts <- normalize(true_wts$wts)

# what are the true fixed effect values?
# true_beta <- c(beta0, alpha_c[2], alpha_x[2], alpha_cx[4])

# Container initialization
## mean estimate
y_bar <- y_bar1 <- y_bar2 <- y_bar3 <- y_bar4 <- y_bar5 <- y_bar6 <- y_bar7 <-  rep(NA, iter)

## general summary
bias <- rMSE <- SE <- 0

## keep track of how many times each cell was selected
cellpicked <- rep(0, J)
cellpickedbb <- rep(0,J)
picked_discount <- rep(0,J) # keep track of cells that miraculously weren't picked in BB step

## posterior of cell mean estimates
mrp_posterior <- mrp_posterior1<- mrp_posterior2 <- mrp_posterior3 <- mrp_posterior4 <-mrp_posterior5 <- mrp_posterior6 <-mrp_posterior7 <- rep(NA, sim)

## CELLMEAN summary
cellmeanbias <- cellmeanrmse <- cellmeanse <- rep(0,J)
CI_cellmean <- rep(NA, J)
CIlength_cellmean <- matrix(NA,nrow=iter, ncol=J)
CI_cellmeancov <- matrix(NA, nrow=iter, ncol=J)

## BETA summary
# betabias <- betarmse <- betase <- rep(0,4)
# CI_beta <- rep(NA, 4)
# CIlength_beta <- matrix(NA,nrow=iter, ncol=4)
# CI_betacov <- matrix(NA, nrow=iter, ncol=4)

## Nj summary
Njbias <- rep(0,J)
Njesttemp <- matrix(NA,nrow=sim, ncol=J)
Njrmse <- rep(0,J)
Njse <- rep(0,J)
CI_Nj <- rep(NA, J)
CIlength_Nj <- matrix(NA,nrow=iter, ncol=J)
CI_Njcov <- matrix(NA, nrow=iter, ncol=J)

## coverage summary
CI <-CI1 <- CI2 <-CI3 <- CI4 <-CI5 <-CI6 <- CI7 <- matrix(nrow = iter, ncol = 2)
CI_covered <- CI_covered1 <-CI_covered2 <-CI_covered3 <-CI_covered4 <-CI_covered5 <-CI_covered6 <-CI_covered7 <-rep(NA,iter)
coverage_rate <-coverage_rate1 <- coverage_rate2 <-coverage_rate3 <- coverage_rate4 <- coverage_rate5 <-coverage_rate6 <- coverage_rate7 <-  0

for(i in 1:iter){
  
  #-- draw sample
  df$I <- rbinom(n=N, size=1, prob = df$p_include)
  # sampled <- 
  #   df[df$I==1,] %>%
  #   dplyr::select(
  #     c(Y, bc, cx, X.0, X.1, J_cell,
  #       Za_categorical, Zb_categorical, Zc_categorical,
  #       Z_categorical, zcts)
  #   )
  sampled <- df[df$I==1, -which(colnames(df)%in% c("Nj","p_include","pxz"))]
  
  # collect weight estimate from each bootstrap iteration
  wts_bootstrapped <- matrix(0,nrow = L, ncol = J)
  
  # re-index Z categories in case we did not sample all of them
  temptbl <- 
    sampled %>% 
    group_by(Za_categorical, Zb_categorical, Zc_categorical) %>% 
    summarise(
      sampled_z = mean(Z_categorical),
      wts = mean(zcts)/n() # empirical inverse probability wt
    )
  sampled_zlength <- 
    nrow(temptbl)
  sampled_z <- 
    temptbl$sampled_z
  temptbl$sampledzlab <- 
    seq(1,sampled_zlength)
  sampled <- 
    left_join(sampled, temptbl)
  rm(list="temptbl")
  ## which J_cells are in our sample?
  temptbl2 <- 
    sampled %>% 
    group_by(J_cell) %>% 
    summarise(
      sampled_J=mean(J_cell)
    )
  sampled_J <- temptbl2$sampled_J
  sampled_Jlength <- length(sampled_J)
  temptbl2$sampledJlab <- seq(1, sampled_Jlength)
  sampled <- left_join(sampled,temptbl2)
  rm(list="temptbl2")
  
  #-- Parameters for stan model
  Ma <- nlevels(sampled$Za_categorical)
  Mb <- nlevels(sampled$Zb_categorical)
  Mc <- nlevels(sampled$Zc_categorical)
  L <- Xcat
  Mbc <- nlevels(sampled$bc)
  Mcx <- nlevels(sampled$cx)
  J_stan <- sampled_Jlength
  n <- nrow(sampled)
  y <- as.vector(sampled$Y)
  cell_label<-sampled$sampledJlab
  
  grptbl <- sampled %>% 
    group_by(J_cell) %>% 
    summarise(zalab = mean(as.numeric(Za_categorical)),
              zblab = mean(as.numeric(Zb_categorical)),
              zclab = mean(as.numeric(Zc_categorical)-1),
              xlab = mean(X.1),
              bclab = mean(as.numeric(bc)),
              cxlab = xlab*zclab
    )
  agroup <- grptbl$zalab
  bgroup <- grptbl$zblab
  cgroup <- grptbl$zclab
  xgroup <- grptbl$xlab 
  bcgroup <- grptbl$bclab
  cxgroup <- grptbl$cxlab
  
  #-- run the model
  assign(modelname[i],stan(file= "mrp_sim.stan",
                           data = c("Ma", "Mb", "Mc","L","Mbc","Mcx","J_stan","n","y","agroup","bgroup","cgroup","xgroup","bcgroup","cxgroup", "cell_label"),
                           iter=staniters, pars=c("cellmean"),warmup = staniters-sim, control=list(adapt_delta=0.99, max_treedepth=13),chains=2))
  
  #-- extract estimates of cell means from the stan model
  # stanpars <-
  #   extract(get(modelname[i]),permuted = TRUE, inc_warmup = FALSE,include = TRUE)
  cellmeans_stan <- extract(get(modelname[i]),permuted = TRUE, inc_warmup = FALSE,include = TRUE)$cellmean
  
  # betas_stan <- stanpars$beta
  # cellmeans_stan <- stanpars$cellmean
  wts_mrp <- rep(0, length(sampled_J))
  
  #- clear containers
  Njesttemp <-  matrix(NA,nrow=sim, ncol=J)
  cellpickedbb <- rep(0,J)
  for(s in 1:sim){
    
    #-- bootstrap resample
    bootstrap<-BayesianBootstrap(c(1:n),1)
    resample<-rmultinom(1,n,bootstrap)
    bootind <- rep(1:n, resample)
    # bootdf <- sampled[bootind,]
    
    boottbl <- sampled[bootind,] %>%
      group_by(J_cell) %>%
      summarise(bootcts = n(), bootwts = sum(wts))
    
    bootsamp_size <- nrow(boottbl)
    wts_new <- boottbl$bootwts
    boot_cats <- boottbl$J_cell
    
    ## find Nj via wfpbb with the bootstrapped categories and new weights
    Nj_hat <- matrix(0, nrow = F_draw, ncol = J)
    
    for(f in 1:F_draw){
      temp <- wtpolyap(ysamp = boot_cats, wts = wts_new, k = N-bootsamp_size)
      Nj_hat[f, boot_cats] <- table(temp)
    }
    
    Nj_est <- colMeans(Nj_hat)
    
    ## collect stats from Nj estimation
    #- bias, SE, rMSE
    cellpickedbb[boot_cats] <- cellpickedbb[boot_cats] +1
    Njesttemp[s,boot_cats] <- Nj_est[boot_cats]
    #- coverage
    
    current_cellmean <- cellmeans_stan[s,]
    wts_mrp <- normalize(Nj_est[sampled_J])  # extract frequencies relevant to estimation of group mean
    wts_mrp1 <- normalize((Nj_est*cmat1)[sampled_J])
    wts_mrp2 <- normalize((Nj_est*cmat2)[sampled_J])
    wts_mrp3 <- normalize((Nj_est*cmat3)[sampled_J])
    wts_mrp4 <- normalize((Nj_est*cmat4)[sampled_J])
    wts_mrp5 <- normalize((Nj_est*cmat5)[sampled_J])
    wts_mrp6 <- normalize((Nj_est*cmat6)[sampled_J])
    wts_mrp7 <- normalize((Nj_est*cmat7)[sampled_J])
    
    mrp_posterior[s] <- crossprod(wts_mrp, current_cellmean)
    mrp_posterior1[s] <- crossprod(wts_mrp1, current_cellmean)
    mrp_posterior2[s] <- crossprod(wts_mrp2, current_cellmean)
    mrp_posterior3[s] <- crossprod(wts_mrp3, current_cellmean)
    mrp_posterior4[s] <- crossprod(wts_mrp4, current_cellmean)
    mrp_posterior5[s] <- crossprod(wts_mrp5, current_cellmean)
    mrp_posterior6[s] <- crossprod(wts_mrp6, current_cellmean)
    mrp_posterior7[s] <- crossprod(wts_mrp7, current_cellmean)
    
    
  }
  y_bar[i] <- mean(mrp_posterior,na.rm=T)
  y_bar1[i] <- mean(mrp_posterior1,na.rm=T)
  y_bar2[i] <- mean(mrp_posterior2,na.rm=T)
  y_bar3[i] <- mean(mrp_posterior3,na.rm=T)
  y_bar4[i] <- mean(mrp_posterior4,na.rm=T)
  y_bar5[i] <- mean(mrp_posterior5,na.rm=T)
  y_bar6[i] <- mean(mrp_posterior6,na.rm=T)
  y_bar7[i] <- mean(mrp_posterior7,na.rm=T)
  
  ## collect stats from cellmeans
  # #- cell picked?
  cellpicked[sampled_J] <- cellpicked[sampled_J] + 1
  #- bias, SE, rMSE
  cellmeanbias[sampled_J] <- cellmeanbias[sampled_J] + colMeans(cellmeans_stan) - popmeans[sampled_J]
  cellmeanrmse[sampled_J] <- cellmeanrmse[sampled_J] + cellmeanbias[sampled_J]^2
  
  #- coverage
  CI_cellmean <- apply(cellmeans_stan, 2, quantile, c(0.025, 0.975))
  CIlength_cellmean[i,sampled_J] <- apply(CI_cellmean,2,diff)
  CI_cellmeancov[i, sampled_J] <- popmeans[sampled_J] >= CI_cellmean[1,] & popmeans[sampled_J] <= CI_cellmean[2,]
  
  ## collect stats from beta
  #- bias, SE, rMSE
  # betabias <- betabias + colMeans(betas_stan) - true_beta
  # betarmse <- betarmse + betabias^2
  
  #- coverage
  # CI_beta <- apply(betas_stan, 2, quantile, c(0.025, 0.975))
  # CIlength_beta[i,] <- apply(CI_beta,2,diff)
  # CI_betacov[i,] <- true_beta >= CI_beta[1,] & true_beta <= CI_beta[2,]
  
  ## collect stats from Nj
  if(0 %in% cellpickedbb[sampled_J]){
    picked_discount[which(cellpickedbb == 0)] <- picked_discount[which(cellpickedbb == 0)]+ 1
    newsampJ <- sampled_J[-which(sampled_J==which(cellpickedbb == 0))]
    Njbias[newsampJ] <- Njbias[newsampJ] + colMeans(Njesttemp[,newsampJ], na.rm=T)-popcts[newsampJ]
    Njrmse[newsampJ] <- Njrmse[newsampJ] + Njbias[newsampJ]^2
    #- coverage
    CI_Nj <- apply(Njesttemp[,newsampJ], 2, quantile, c(0.025, 0.975))
    CIlength_Nj[i,newsampJ] <- apply(CI_Nj,2,diff)
    CI_Nj[i,newsampJ] <- popcts[newsampJ] >= CI_Nj[1,] & popcts[newsampJ] <= CI_Nj[2,]
    # print("bb missed one!")
  }else{
    Njbias[sampled_J] <- Njbias[sampled_J] + colMeans(Njesttemp[,sampled_J], na.rm=T)-popcts[sampled_J]
    Njrmse[sampled_J] <- Njrmse[sampled_J] + Njbias[sampled_J]^2
    #- coverage
    CI_Nj <- apply(Njesttemp[,sampled_J], 2, quantile, c(0.025, 0.975), na.rm=T)
    CIlength_Nj[i,sampled_J] <- apply(CI_Nj,2,diff)
    CI_Njcov[i,sampled_J] <- popcts[sampled_J] >= CI_Nj[1,] & popcts[sampled_J] <= CI_Nj[2,]
  }
  ## collect general stats
  CI_Nj <- apply(Njesttemp[,sampled_J], 2, quantile, c(0.025, 0.975), na.rm=T)
  CIlength_Nj[i,sampled_J] <- apply(CI_Nj,2,diff)
  CI_Njcov[i,sampled_J] <- popcts[sampled_J] >= CI_Nj[1,] & popcts[sampled_J] <= CI_Nj[2,]
  CI[i,] <- quantile(mrp_posterior, c(0.025, 0.975),na.rm=T)
  CI1[i,] <- quantile(mrp_posterior1, c(0.025, 0.975),na.rm=T)
  CI2[i,] <- quantile(mrp_posterior2, c(0.025, 0.975),na.rm=T)
  CI3[i,] <- quantile(mrp_posterior3, c(0.025, 0.975),na.rm=T)
  CI4[i,] <- quantile(mrp_posterior4, c(0.025, 0.975),na.rm=T)
  CI5[i,] <- quantile(mrp_posterior5, c(0.025, 0.975),na.rm=T)
  CI6[i,] <- quantile(mrp_posterior6, c(0.025, 0.975),na.rm=T)
  CI7[i,] <- quantile(mrp_posterior7, c(0.025, 0.975),na.rm=T)
  
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
  #-- remove stanfit object to prevent segfault and conserve memory
  rm(list=c(modelname[i],"cellmeans_stan","Njesttemp",
            "mrp_posterior","mrp_posterior1", "mrp_posterior2","mrp_posterior3","mrp_posterior4","mrp_posterior5", "mrp_posterior6", "mrp_posterior7"))
  

}

## print results
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

# PRINT RESULTS
print("Result Summary")
resultmat  

## CELLMEAN Summary
cellmeanbias <- 
  cellmeanbias / cellpicked 
cellmeanrmse <- 
  sqrt(cellmeanrmse/cellpicked)
cellmeanse <- 
  sqrt(cellmeanrmse^2 - cellmeanbias^2)

# print("========== CELL MEAN SUMMARY ==========")
# 
# print("Cell picking frequency")
# print(cellpicked)
# 
# print("Bias of cellmeans")
# print(cellmeanbias)
# 
# print("rMSE of cellmeans")
# print(cellmeanrmse)
# 
# print("SE of cellmeans")
# print(cellmeanse)
# 
# print("CI length of cellmeans")
# colMeans(CI_cellmean, na.rm=T)
# 
# print("CI coverage rate of cellmeans")
# colMeans(CI_cellmeancov,na.rm=T)
# 
# print("========== Nj SUMMARY ==========")
Njbias <- Njbias / cellpicked
Njrmse <- sqrt(Njrmse/ cellpicked )
Njse <- sqrt(Njrmse^2 - Njbias^2)

# print("Bias of Njs")
# print(Njbias)
# 
# print("rMSE of Njs")
# print(Njrmse)
# 
# print("SE of Njs")
# print(Njse)
# 
# print("CI length of Njs")
# colMeans(CIlength_Nj, na.rm=T)
# 
# print("CI coverage rate of Njs")
# colMeans(CI_Njcov,na.rm=T)

warnings()

write.csv(resultmat, "res_wmrp.csv")
write.csv(cellmeanbias, "cellmeanbias_wmrp.csv")
write.csv(cellmeanse, "cellmeanse_wmrp.csv")
write.csv(cellmeanrmse, "cellmeanrmse_wmrp.csv")
write.csv(Njbias, "Njbias_wmrp.csv")
write.csv(Njse, "Njse_wmrp.csv")
write.csv(Njrmse, "Njrmse_wmrp.csv")
write.csv(colMeans(CIlength_Nj,na.rm=T), "NjCI_wmrp.csv")
write.csv(colMeans(CI_Njcov,na.rm=T),"Njcov_wmrp.csv")