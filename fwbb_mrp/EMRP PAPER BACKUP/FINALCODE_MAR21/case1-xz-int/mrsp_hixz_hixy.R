#====================================
# Code for MrsP simulations
# Adapted from Leemann, Wasserfallen
#====================================
# last modified: nov. 29, 2020

#-- load required packages
require(dplyr)
require(rstan)
options(mc.cores = 2)

#-- load population generating function:
source("genPop_hixz_hixy.R")
# source("genPop_hixz_loxy.R")
# source("genPop_loxz_hixy.R")
# source("genPop_loxz_loxy.R")

#-- parameters for testing the function
set.seed(50)
iter = 200
sim = 1000
staniters = 2000

#-- random integers for naming the stanfit to avoid segfault issues
stanfitrand <- sample(1:iter, iter, replace=F)
modelname <- paste("stanfit",stanfitrand, sep="")

#-- helper function to normalize
normalize <- function(vec){
  return(vec/sum(vec))
}

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

# true_wts <- df %>% group_by(J_cell) %>% summarise(popcts = mean(Nj),wts = mean(Nj)/N, cellmeans = mean(Y))
# popwts <- df %>% group_by(Z_categorical) %>% summarise(popctsz = n())
# popctsz <- popwts$popctsz
# 
# popmeans <- true_wts$cellmeans
# popcts <- true_wts$popcts
# true_wts <- normalize(true_wts$wts)
true_wts <- df %>% group_by(Z_categorical,X.1) %>% summarise(popcts = mean(Nj),zcts = mean(zcts),wts = mean(Nj)/N, cellmeans = mean(Y))
popwts <- df %>% group_by(Z_categorical) %>% summarise(popctsz = n())
popctsz <- popwts$popctsz
# assuming J = 100 for the simulations; if not, need to modify to account for uneven Z x X table when drawing from multinomial
# zskip <- true_wts$popcts ==true_wts$zcts # J-length indicator of whether Z-category has 1 level of X (vs. 2 levels)

popmeans <- true_wts$cellmeans
popcts <- true_wts$popcts
true_wts <- normalize(true_wts$wts)
# initialize containers for results
y_bar <- y_bar1 <- y_bar2 <- y_bar3 <- y_bar4 <- y_bar5 <- y_bar6 <- y_bar7<-rep(NA, iter)
mrp_posterior <- mrp_posterior1<- mrp_posterior2 <- mrp_posterior3 <- mrp_posterior4 <-mrp_posterior5 <- mrp_posterior6 <-mrp_posterior7<-rep(NA, sim)


# Nj Summary
cellpicked <- rep(0, J)
Njbias <- rep(0,J)
Njesttemp <- matrix(NA,nrow=sim, ncol=J)
Njrmse <- rep(0,J)
Njse <- rep(0,J)
CI_Nj <- rep(NA, J)
CIlength_Nj <- matrix(NA,nrow=iter, ncol=J)
CI_Njcov <- matrix(NA, nrow=iter, ncol=J)

CI       <- matrix(nrow = iter, ncol = 2)   # credible interval
CI1       <- matrix(nrow = iter, ncol = 2)   # credible interval
CI2       <- matrix(nrow = iter, ncol = 2)   # credible interval
CI3       <- matrix(nrow = iter, ncol = 2)   # credible interval
CI4       <- matrix(nrow = iter, ncol = 2)   # credible interval
CI5       <- matrix(nrow = iter, ncol = 2)   # credible interval
CI6       <- matrix(nrow = iter, ncol = 2)   # credible interval
CI7       <- matrix(nrow = iter, ncol = 2)   # credible interval

# CI_covered    <- rep(0,iter)    # 1 if covered by CI, 0 else
# CI_covered1    <- rep(0,iter)    # 1 if covered by CI, 0 else
# CI_covered2    <- rep(0,iter)    # 1 if covered by CI, 0 else
# CI_covered3    <- rep(0,iter)    # 1 if covered by CI, 0 else
# CI_covered4    <- rep(0,iter)    # 1 if covered by CI, 0 else
# CI_covered5    <- rep(0,iter)    # 1 if covered by CI, 0 else
# CI_covered6    <- rep(0,iter)    # 1 if covered by CI, 0 else
# CI_covered7    <- rep(0,iter)    # 1 if covered by CI, 0 else
CI_covered    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered1    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered2    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered3    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered4    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered5    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered6    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered7    <- rep(NA,iter)    # 1 if covered by CI, NA else

bias <- 0
rMSE <- 0
SE <- 0


for(i in 1:iter){
  #-- draw sample
  df$I <- rbinom(n=N, size=1, prob = df$p_include)
  sampled <- df[df$I==1,-which(colnames(df)%in% c("Nj","p_include","pxz"))]
  
  sampled_z <- sort(unique(sampled$Z_categorical))
  sampled_J <- sort(unique(sampled$J_cell))
  sampled_Jlength <- length(sampled_J)
  #-- calculate weights
  zrelfreq <- matrix(0,nrow = Zcat, ncol=2)
  zrelfreq[sampled_z,] <- t(table(sampled$X.1, sampled$Z_categorical))*(as.vector(1/table(sampled$Z_categorical)))
  # sampwt_denom <- sampled %>% group_by(Z_categorical) %>% summarise(denomcts=n())
  # sampwt_num <- sampled %>% group_by(Z_categorical,X.1) %>% summarise(numcts=n())
  # fulltbl <- left_join(sampwt_num, sampwt_denom) %>% left_join(popwts) %>% 
  #   mutate(estcts = popctsz / denomcts * numcts,
  #          zrelfreq = numcts/denomcts)
  # sampled_popcts <- fulltbl$popctsz
  
  # sampwt_denom <- sampled %>% group_by(Z_categorical) %>% summarise(denomcts=n())
  # sampwt_num <- sampled %>% group_by(Z_categorical,X.1) %>% summarise(numcts=n())
  # fulltbl <- left_join(sampwt_num, sampwt_denom) %>% left_join(popwts) %>% 
  #   mutate(estcts = popctsz / denomcts * numcts,
  #          zrelfreq = numcts/denomcts)
  # sampled_popcts <- fulltbl$popctsz
  
  # rmultinom_vec <- cbind(fulltbl$Z_categorical, fulltbl$zrelfreq) %>% as.data.frame()
  # colnames(rmultinom_vec) <- c("Z_categorical", "estwts")
  
  # re-index Z categories in case we did not sample all of them
  # 
  # temptbl <- sampled %>% group_by(Za_categorical, Zb_categorical, Zc_categorical) %>% summarise(sampled_z = mean(Z_categorical),wts = mean(zcts)/n())
  # sampled_zlength <- nrow(temptbl)
  # sampled_z <- temptbl$sampled_z
  # temptbl$sampledzlab <- seq(1,sampled_zlength)
  # sampled <- left_join(sampled, temptbl)
  # rm(list="temptbl")
  ## which J_cells are in our sample?
  temptbl2 <- sampled %>% group_by(J_cell) %>% summarise(sampled_J=mean(J_cell))
  cellpicked[sampled_J] <- cellpicked[sampled_J] + 1
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
  
  
  # #-- run the model
  assign(modelname[i],stan(file= "mrp_sim.stan",
                           data = c("Ma", "Mb", "Mc","L","Mbc","Mcx","J_stan","n","y","agroup","bgroup","cgroup","xgroup","bcgroup","cxgroup", "cell_label"),
                           iter=staniters, pars=c("cellmean"),warmup = staniters-sim, control=list(adapt_delta=0.99, max_treedepth=13),chains=2))
  #-- extract estimates of cell means from the stan model
  # stanpars <-
  #   extract(get(modelname[i]),permuted = TRUE, inc_warmup = FALSE,include = TRUE)
  # cellmeans_stan <- stanpars$cellmean
  cellmeans_stan <- extract(get(modelname[i]),permuted = TRUE, inc_warmup = FALSE,include = TRUE)$cellmean
  # wts_mrp <- cbind(fulltbl$Z_categorical, rep(0, sampled_Jlength)) %>% as.data.frame()
  wts_mrp <- matrix(0,nrow = Zcat,ncol=2,byrow=F)
  # colnames(wts_mrp) <- c("Z_categorical", "wts")
  
  Njesttemp <- matrix(NA,nrow=sim, ncol=J)
mrp_posterior <- mrp_posterior1<- mrp_posterior2 <- mrp_posterior3 <- mrp_posterior4 <-mrp_posterior5 <- mrp_posterior6 <-mrp_posterior7<-rep(NA, sim)
  
  for(s in 1:sim){
    # wts_mrp$wts <- 0
    
    for(m in sampled_z){
      
      # if(sum(fulltbl$Z_categorical==m) < 2){
      #   # wts_mrp[wts_mrp$Z_categorical==m, "wts"] <- fulltbl$estcts[fulltbl$Z_categorical==m]
      #   wts_mrp[wts_mrp$Z_categorical==m, "wts"] <- popctsz[m]
      #   # print("if!")
      # }
      wts_mrp[m,] <- rowMeans(rmultinom(popctsz[m], 1, prob = zrelfreq[m,]))*popctsz[m]
      # else{
      #   probs <- fulltbl[fulltbl$Z_categorical == m, "zrelfreq"] %>% unlist() %>% as.double() 
      #   wts_mrp[wts_mrp$Z_categorical == m, "wts"] <- 
      #     rowMeans(rmultinom(popctsz[m], 2-1, prob = probs))*popctsz[m]
      # print("else!")
      #}
    }
    
    
    # current_cellmean <- cellmeans_stan[s,]
    # Njesttemp[s,sampled_J] <- wts_mrp$wts
    # wts <- normalize(wts_mrp$wts)
    # wts1 <- normalize((wts*cmat1[sampled_J]))
    # wts2 <- normalize((wts*cmat2[sampled_J]))
    # wts3 <- normalize((wts*cmat3[sampled_J]))
    # wts4 <- normalize((wts*cmat4[sampled_J]))
    # wts5 <- normalize((wts*cmat5[sampled_J]))
    # 
    # mrp_posterior[s] <- crossprod(wts, current_cellmean)
    # mrp_posterior1[s] <- crossprod(wts1, current_cellmean)
    # mrp_posterior2[s] <- crossprod(wts2, current_cellmean)
    # mrp_posterior3[s] <- crossprod(wts3, current_cellmean)
    # mrp_posterior4[s] <- crossprod(wts4, current_cellmean)
    # mrp_posterior5[s] <- crossprod(wts5, current_cellmean)
    current_cellmean <- cellmeans_stan[s,]
    # Njesttemp[s,sampled_J] <- wts_mrp$wts
    # wts <- normalize(wts_mrp$wts)
    wts <- as.vector(t(wts_mrp))[sampled_J]
    Njesttemp[s,sampled_J] <- wts
    wts <- normalize(wts)
    
    wts1 <- normalize((wts*cmat1[sampled_J]))
    wts2 <- normalize((wts*cmat2[sampled_J]))
    wts3 <- normalize((wts*cmat3[sampled_J]))
    wts4 <- normalize((wts*cmat4[sampled_J]))
    wts5 <- normalize((wts*cmat5[sampled_J]))
    wts6 <- normalize((wts*cmat6[sampled_J]))
    wts7 <- normalize((wts*cmat7[sampled_J]))
    
    
    mrp_posterior[s] <- crossprod(wts, current_cellmean)
    mrp_posterior1[s] <- crossprod(wts1, current_cellmean)
    mrp_posterior2[s] <- crossprod(wts2, current_cellmean)
    mrp_posterior3[s] <- crossprod(wts3, current_cellmean)
    mrp_posterior4[s] <- crossprod(wts4, current_cellmean)
    mrp_posterior5[s] <- crossprod(wts5, current_cellmean)
    mrp_posterior6[s] <- crossprod(wts6, current_cellmean)
    mrp_posterior7[s] <- crossprod(wts7, current_cellmean)
    
  }
  y_bar[i] <- mean(mrp_posterior,na.rm=T)
  y_bar1[i] <- mean(mrp_posterior1,na.rm=T)
  y_bar2[i] <- mean(mrp_posterior2,na.rm=T)
  y_bar3[i] <- mean(mrp_posterior3,na.rm=T)
  y_bar4[i] <- mean(mrp_posterior4,na.rm=T)
  y_bar5[i] <- mean(mrp_posterior5,na.rm=T)
  y_bar6[i] <- mean(mrp_posterior6,na.rm=T)
  y_bar7[i] <- mean(mrp_posterior7,na.rm=T)
  
  # Nj summary
  Njbias[sampled_J] <- Njbias[sampled_J] + colMeans(Njesttemp[,sampled_J], na.rm=T)-popcts[sampled_J]
  Njrmse[sampled_J] <- Njrmse[sampled_J] + Njbias[sampled_J]^2
  #- coverage
  CI_Nj <- apply(Njesttemp[,sampled_J], 2, quantile, c(0.025, 0.975), na.rm=T)
  CIlength_Nj[i,sampled_J] <- apply(CI_Nj,2,diff)
  CI_Njcov[i,sampled_J] <- ((popcts[sampled_J] >= CI_Nj[1,]) & (popcts[sampled_J] <= CI_Nj[2,]))
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
# ybarmat <- cbind(y_bar, y_bar1, y_bar2, y_bar3, y_bar4, y_bar5,y_bar6, y_bar7)
# errormat <- cbind(ybarmat[,1]-true_mean, ybarmat[,2]-true_mean1,ybarmat[,3]-true_mean2,ybarmat[,4]-true_mean3,
#                   ybarmat[,5]-true_mean4, ybarmat[,6]-true_mean5,ybarmat[,7]-true_mean6, ybarmat[,8]-true_mean7)
# bias <- colMeans(errormat)
# rMSE <- sqrt(colMeans(errormat^2))
# SE <- cbind(y_bar - mean(y_bar), y_bar1 - mean(y_bar1),y_bar2 - mean(y_bar2),y_bar3 - mean(y_bar3),y_bar4 - mean(y_bar4),y_bar5 - mean(y_bar5),
#             y_bar6-mean(y_bar6),y_bar7 - mean(y_bar7))
# SE <- sqrt(colMeans(SE^2))
# CI_length <- c(mean(CI[,2]-CI[,1]), mean(CI1[,2]-CI1[,1]),mean(CI2[,2]-CI2[,1]),mean(CI3[,2]-CI3[,1]),
#                mean(CI4[,2]-CI4[,1]),mean(CI5[,2]-CI5[,1]),mean(CI6[,2]-CI6[,1]),mean(CI7[,2]-CI7[,1]))
# coverage_rate <- c(mean(CI_covered),mean(CI_covered1),mean(CI_covered2),mean(CI_covered3),mean(CI_covered4), mean(CI_covered5),
#                    mean(CI_covered6), mean(CI_covered7))
# 
# resultmat <- rbind(rMSE, bias, SE, CI_length, coverage_rate)
# colnames(resultmat) <- c("Overall", "group 1", "group2", "group3", "group4", "group5","group6", "group7")

# CI is a list
#ybarmat, CIcovered is a matrix
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

write.csv(resultmat, "res_mmrp.csv")
write.csv(Njbias, "Njbias_mmrp.csv")
write.csv(Njse, "Njse_mmrp.csv")
write.csv(Njrmse, "Njrmse_mmrp.csv")
write.csv(colMeans(CIlength_Nj, na.rm=T), "NjCI_mmrp.csv")
write.csv(colMeans(CI_Njcov,na.rm=T), "Njcov_mmrp.csv")