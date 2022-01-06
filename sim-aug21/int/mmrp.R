#====================================
# Code for MrsP simulations
# Adapted from Leemann, Wasserfallen
#====================================
# last modified: june 2021

#-- load required packages
require(dplyr)
require(rstan)
# options(mc.cores = 2)

#-- load population generating function:
source("genPop_hixz_hixy.R")

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
true_mean1 <- mean(df$Y[df$J_cell%in% group1])
true_mean2 <- mean(df$Y[df$J_cell%in% group2])
true_mean3 <- mean(df$Y[df$J_cell%in% group3])
true_mean4 <- mean(df$Y[df$J_cell%in% group4])
true_mean5 <- mean(df$Y[df$J_cell%in% group5])

true_wts <- df %>% group_by(Z_categorical,X.1) %>% summarise(popcts = mean(Nj),zcts = mean(zcts),wts = mean(Nj)/N, cellmeans = mean(Y))
popwts <- df %>% group_by(Z_categorical) %>% summarise(popctsz = n())
popctsz <- popwts$popctsz
# assuming J = 100 for the simulations; if not, need to modify to account for uneven Z x X table when drawing from multinomial
zskip <- true_wts$Z_categorical[true_wts$popcts ==true_wts$zcts] # J-length indicator of whether Z-category has 1 level of X (vs. 2 levels)
if(length(zskip)==0){zskip <- 0}
popmeans <- true_wts$cellmeans
popcts <- true_wts$popcts
true_wts <- normalize(true_wts$wts)
# initialize containers for results
y_bar <- y_bar1 <- y_bar2 <- y_bar3 <- y_bar4 <- y_bar5 <-rep(NA, iter)
mrp_posterior <- mrp_posterior1<- mrp_posterior2 <- mrp_posterior3 <- mrp_posterior4 <-mrp_posterior5 <-rep(NA, sim)


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

CI_covered    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered1    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered2    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered3    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered4    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered5    <- rep(NA,iter)    # 1 if covered by CI, NA else

bias <- 0
rMSE <- 0
SE <- 0

for(i in 1:iter){
set.seed(stanfitrand[i])
  #-- draw sample
  df$I <- rbinom(n=N, size=1, prob = df$p_include)
  sampled <- df[df$I==1,-which(colnames(df)%in% c("Nj","p_include","pxz"))]
  
  sampled_z <- sort(unique(sampled$Z_categorical))
  sampled_J <- sort(unique(sampled$J_cell))
  sampled_Jlength <- length(sampled_J)
  #-- calculate weights
  zrelfreq <- matrix(0,nrow = Zcat, ncol=2)
  zrelfreq[sampled_z,] <- t(table(sampled$X.1, sampled$Z_categorical))*(as.vector(1/table(sampled$Z_categorical)))
  ## which J_cells are in our sample?
  temptbl2 <- sampled %>% group_by(J_cell) %>% summarise(sampled_J=mean(J_cell))
  cellpicked[sampled_J] <- cellpicked[sampled_J] + 1
  temptbl2$sampledJlab <- seq(1, sampled_Jlength)
  sampled <- left_join(sampled,temptbl2)
  rm(list="temptbl2")
  
    #-- Parameters for stan model
  Ma <- nlevels(sampled$Za_categorical)
  Mb <- nlevels(sampled$Zb_categorical)
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
    )
  agroup <- grptbl$zalab
  bgroup <- grptbl$zblab
  cgroup <- grptbl$zclab
  xgroup <- grptbl$xlab 

  
  # #-- run the model
  assign(modelname[i],stan(file= "mrp_sim.stan",
                           data = c("Ma", "Mb", "J_stan","n","y","agroup","bgroup","cgroup","xgroup", "cell_label"),
                           iter=staniters, 
                           pars = c("cellmean"),
                           warmup = staniters-sim, control=list(adapt_delta=0.99, max_treedepth=13),chains=1))
  
  #-- extract estimates of cell means from the stan model
  cellmeans_stan <- extract(get(modelname[i]),permuted = TRUE, inc_warmup = FALSE,include = TRUE)$cellmean
  wts_mrp <- matrix(0,nrow = Zcat,ncol=2,byrow=F)

  Njesttemp <- matrix(NA,nrow=sim, ncol=J)
mrp_posterior <- mrp_posterior1<- mrp_posterior2 <- mrp_posterior3 <- mrp_posterior4 <-mrp_posterior5 <-rep(NA, sim)
  
  for(s in 1:sim){
    # wts_mrp$wts <- 0
    
    for(m in sampled_z){
      if(m != zskip){
      wts_mrp[m,] <- rowMeans(rmultinom(popctsz[m], 1, prob = zrelfreq[m,]))*popctsz[m]
      }
      else{
        wts_mrp[m,] <- c(popctsz[m], NA)
      }
      
    }
    
    current_cellmean <- cellmeans_stan[s,]
    wts <- na.omit(as.vector(t(wts_mrp)))[sampled_J]
    Njesttemp[s,sampled_J] <- wts
    wts <- normalize(wts)
    
    wts1 <- normalize((wts*cmat1[sampled_J]))
    wts2 <- normalize((wts*cmat2[sampled_J]))
    wts3 <- normalize((wts*cmat3[sampled_J]))
    wts4 <- normalize((wts*cmat4[sampled_J]))
    wts5 <- normalize((wts*cmat5[sampled_J]))

    
    mrp_posterior[s] <- crossprod(wts, current_cellmean)
    mrp_posterior1[s] <- crossprod(wts1, current_cellmean)
    mrp_posterior2[s] <- crossprod(wts2, current_cellmean)
    mrp_posterior3[s] <- crossprod(wts3, current_cellmean)
    mrp_posterior4[s] <- crossprod(wts4, current_cellmean)
    mrp_posterior5[s] <- crossprod(wts5, current_cellmean)

  }
  y_bar[i] <- mean(mrp_posterior,na.rm=T)
  y_bar1[i] <- mean(mrp_posterior1,na.rm=T)
  y_bar2[i] <- mean(mrp_posterior2,na.rm=T)
  y_bar3[i] <- mean(mrp_posterior3,na.rm=T)
  y_bar4[i] <- mean(mrp_posterior4,na.rm=T)
  y_bar5[i] <- mean(mrp_posterior5,na.rm=T)

  # Nj summary
  Njbias[sampled_J] <- Njbias[sampled_J] + colMeans(Njesttemp[,sampled_J], na.rm=T)-popcts[sampled_J]
  Njrmse[sampled_J] <- Njrmse[sampled_J] + (colMeans(Njesttemp[,sampled_J], na.rm=T)-popcts[sampled_J])^2
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

  # sampsub <- unique(sampled$subgroup)

  if(is.na(between(true_mean, CI[i,1], CI[i,2]))){break}
  if(between(true_mean, CI[i,1], CI[i,2])){CI_covered[i] <- 1}else if(sum(is.na(CI[i,])!=0)){CI_covered[i] <- NA}else{CI_covered[i] <- 0}
  if(between(true_mean1, CI1[i,1], CI1[i,2])){CI_covered1[i] <- 1} else{CI_covered1[i] <- 0}
  if(between(true_mean2, CI2[i,1], CI2[i,2])){CI_covered2[i] <- 1} else{CI_covered2[i] <- 0}
  if(between(true_mean3, CI3[i,1], CI3[i,2])){CI_covered3[i] <- 1} else{CI_covered3[i] <- 0}
  if(between(true_mean4, CI4[i,1], CI4[i,2])){CI_covered4[i] <- 1} else{CI_covered4[i] <- 0}
  if(between(true_mean5, CI5[i,1], CI5[i,2])){CI_covered5[i] <- 1} else{CI_covered5[i] <- 0}
  #-- remove stanfit object to prevent segfault and conserve memory
  rm(list=c(modelname[i],"cellmeans_stan","Njesttemp",
  "mrp_posterior","mrp_posterior1", "mrp_posterior2","mrp_posterior3","mrp_posterior4","mrp_posterior5"))
  
}    

# CI is a list
#ybarmat, CIcovered is a matrix
ybarmat <- cbind(y_bar, y_bar1, y_bar2, y_bar3, y_bar4, y_bar5)
truemeanmat <- c(true_mean,true_mean1,true_mean2,true_mean3,true_mean4,true_mean5)
labs <- c("Overall", "group 1", "group2", "group3", "group4", "group5")
CImat <- list(CI,CI1,CI2,CI3,CI4,CI5)
CIcovmat <- cbind(CI_covered,CI_covered1,CI_covered2,CI_covered3,CI_covered4,CI_covered5)
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

warnings()

write.csv(resultmat, "res_mmrp.csv")
write.csv(Njbias, "Njbias_mmrp.csv")
write.csv(Njse, "Njse_mmrp.csv")
write.csv(Njrmse, "Njrmse_mmrp.csv")
write.csv(colMeans(CIlength_Nj, na.rm=T), "NjCI_mmrp.csv")
write.csv(colMeans(CI_Njcov,na.rm=T), "Njcov_mmrp.csv")
