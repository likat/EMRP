#=========================================
# Code for WFPBB direct impute simulation
#=========================================
# For categorical / discrete Y
# last modified: june 2021

#---- SETUP -----
## load packages
require(dplyr)
require(polyapost)
require(LaplacesDemon)

#-- load population generating function:
source("genPop_hixz_hixy.R")

## helper functions
`%notin%` <- Negate(`%in%`)
#==== FUNCTION FOR ESTIMATING SE  ======
# based on t-approximation
# PARAMETERS:
	# L <real> = number of synthetic populations created for a given sample
	# tempvec[L] <vector> = estimate for each of the synthetic populations
# OUTPUT (res[2] <list>)
	# res[[1]] = SE_est <real> = estimate of SE
	# res[[2]] = CI[2] <vector> = lower and upper endpoints of the CI
estimatedsefun <- function(L,tempvec){
	CI <- rep(NA,2)
	SE_est <- sqrt((1+1/L)*var(tempvec))
	CI[1] <- mean(tempvec) - SE_est*qt(0.975,L-1)
	CI[2] <- mean(tempvec) + SE_est*qt(0.975,L-1)
	res <- list(SE_est,CI)
	return(res)
}

#==== FUNCTION FOR creating objects  =====
# PARAMETERS:
	# basename <string> = what is the base name of the container?
	# subgroups <real> = number of subgroups we are estimating
	# type <string> = matrix or vector?
	# dim <real> or <vector> = <real> if type is vector, <vector>[2] if type is matrix
# OUTPUT 
	# assign subgroups+1 objects to the global environment
containerfun <- function(basename, subgroups=5, type="vector", dims){
	if(type=="vector"){
		assign(basename, rep(NA,dims), envir = .GlobalEnv)
		for(it in 1:subgroups){
			assign(paste(basename,it, sep=""), rep(NA,dims), envir=.GlobalEnv)
			}
	}	
	if(type=="matrix"){
		assign(basename, matrix(nrow = dims[1], ncol=dims[2]), envir = .GlobalEnv)
		for(it in 1:subgroups){
			assign(paste(basename,it, sep=""), matrix(nrow = dims[1], ncol=dims[2]), envir=.GlobalEnv)
			}
	}	
}
## define seed, number of simulations ('sim'), 
# number of bootstraps ('L'), 
# number of FPBB draws ('F_draw')
set.seed(50)
sim <- 200
L <- 1000
F_draw <- 20
grps <- 5 # number of subgroups

stanfitrand <- sample(1:sim, sim, replace=F)
#-- define parameters from population data
true_mean <- mean(df$Y)
true_mean1 <- mean(df$Y[df$J_cell%in% group1])
true_mean2 <- mean(df$Y[df$J_cell%in% group2])
true_mean3 <- mean(df$Y[df$J_cell%in% group3])
true_mean4 <- mean(df$Y[df$J_cell%in% group4])
true_mean5 <- mean(df$Y[df$J_cell%in% group5])

popmeans <- df %>% group_by(J_cell) %>% summarise(mean(Y)) %>% dplyr::select(-J_cell) %>% as.matrix() %>% as.numeric()

#-- initialize containers
containerfun(basename="y_bar", subgroups=grps, type="vector",dims = sim)
containerfun(basename="bbmean", subgroups=grps, type="vector",dims = F_draw)
containerfun(basename="tempmean", subgroups=grps, type="vector",dims = L)
containerfun(basename="CI", subgroups=grps, type="matrix",dims = c(sim,2))
containerfun(basename="CI_est", subgroups=grps, type="matrix",dims = c(sim,2))
containerfun(basename="CI_estcovered", subgroups=grps, type="vector",dims = sim)
containerfun(basename="CI_covered",subgroups=grps, type="vector",dims = sim)
containerfun(basename="SE_est",subgroups=grps, type="vector",dims = sim)
containerfun(basename="bayesvar", subgroups=grps,type="vector",dims=sim)

coverage_rate <- 0
bias <- 0
rMSE <- 0

#-- population counts for sampling
popcts <- as.vector(table(df$J_cell))
p_include <- df$p_include

#----- BEGIN SIMULATIONS ----
t0 <- Sys.time()

for(i in 1:sim){
  
set.seed(stanfitrand[i])  
  #-- draw sample 
  df$I <- rbinom(n=N, size = 1, prob = p_include)
  
  # -- create object for sampled data
  sampled <- df[df$I==1,]
  
  #-- sample size
  n <- nrow(sampled)
  
  #-- create empirical ipw
  sampledztbl <- sampled %>% group_by(Z_categorical) %>% summarise(wts = mean(zcts)/n())
  sampled_z <- sampledztbl$Z_categorical
  sampled <- left_join(sampled, sampledztbl)
  sampled_J <- unique(sampled$J_cell)
  
  #------ DRAW SYNTHETIC POPULATIONS ----
  for(l in 1:L){
    
    ## (1) Bootstrap from parent sample
    bootstrap<-BayesianBootstrap(c(1:n),1)
    resample<-rmultinom(1,n,bootstrap)
    bootind <- rep(1:n, resample)
    boottbl <- sampled[bootind,] %>%
      group_by(Y,J_cell) %>%
      summarise(bootcts = n(), bootwts = sum(wts))

    bootsamp_size <- nrow(boottbl)
    wts_new <- boottbl$bootwts
    boot_cats <- 1:nrow(boottbl)
    
    ## (2) WFPBB with new weights and unique indices
    fpbb_temp <- rep(0, N)
    for(f in 1:F_draw){
      fpbb_temp <- wtpolyap(ysamp = boot_cats, wts = wts_new, k = N-bootsamp_size)
      synthpop <- boottbl[fpbb_temp,]
      bbmean[f] <- mean(synthpop$Y)
      bbmean1[f] <- mean(synthpop$Y[synthpop$J_cell %in% group1])
      bbmean2[f] <- mean(synthpop$Y[synthpop$J_cell %in% group2])
      bbmean3[f] <- mean(synthpop$Y[synthpop$J_cell %in% group3])
      bbmean4[f] <- mean(synthpop$Y[synthpop$J_cell %in% group4])
      bbmean5[f] <- mean(synthpop$Y[synthpop$J_cell %in% group5])
    }

    
    ## (3) Collect estimate from the lth synthetic population
    tempmean[l] <- mean(bbmean,na.rm=T)
    tempmean1[l] <- mean(bbmean1,na.rm=T)
    tempmean2[l] <- mean(bbmean2,na.rm=T)
    tempmean3[l] <- mean(bbmean3,na.rm=T)
    tempmean4[l] <- mean(bbmean4,na.rm=T)
    tempmean5[l] <- mean(bbmean5,na.rm=T)
    
    
    #print(l)
  }  
  ## (4) Final estimate for parent sample i
  y_bar[i] <- mean(tempmean,na.rm=T)
  y_bar1[i] <- mean(tempmean1,na.rm=T)
  y_bar2[i] <- mean(tempmean2,na.rm=T)
  y_bar3[i] <- mean(tempmean3,na.rm=T)
  y_bar4[i] <- mean(tempmean4,na.rm=T)
  y_bar5[i] <- mean(tempmean5,na.rm=T)
   
  ## (5) Confidence Interval
  ## (5a) Empirical CI
  CI[i,] <- quantile(tempmean, probs = c(0.025, 0.975),na.rm=T)
  CI1[i,] <- quantile(tempmean1, probs = c(0.025, 0.975),na.rm=T)
  CI2[i,] <- quantile(tempmean2, probs = c(0.025, 0.975),na.rm=T)
  CI3[i,] <- quantile(tempmean3, probs = c(0.025, 0.975),na.rm=T)
  CI4[i,] <- quantile(tempmean4, probs = c(0.025, 0.975),na.rm=T)
  CI5[i,] <- quantile(tempmean5, probs = c(0.025, 0.975),na.rm=T)  

  if(between(true_mean, CI[i,1], CI[i,2])){CI_covered[i] <- 1}else{CI_covered[i] <- 0}
  if(between(true_mean1, CI1[i,1], CI1[i,2])){CI_covered1[i] <- 1}else{CI_covered1[i] <- 0}
  if(between(true_mean2, CI2[i,1], CI2[i,2])){CI_covered2[i] <- 1}else{CI_covered2[i] <- 0}
  if(between(true_mean3, CI3[i,1], CI3[i,2])){CI_covered3[i] <- 1}else{CI_covered3[i] <- 0}
  if(between(true_mean4, CI4[i,1], CI4[i,2])){CI_covered4[i] <- 1}else{CI_covered4[i] <- 0}
  if(between(true_mean5, CI5[i,1], CI5[i,2])){CI_covered5[i] <- 1}else{CI_covered5[i] <- 0}
  
  ## (5b) Estimated CI
  bayesvar[i] <- var(tempmean)
  bayesvar1[i] <- var(tempmean1)
  bayesvar2[i] <- var(tempmean2)
  bayesvar3[i] <- var(tempmean3)
  bayesvar4[i] <- var(tempmean4)
  bayesvar5[i] <- var(tempmean5)
  CI_est[i,] <- estimatedsefun(L,tempmean)[[2]]
CI_est1[i,] <- estimatedsefun(L,tempmean1)[[2]]
CI_est2[i,] <- estimatedsefun(L,tempmean2)[[2]]
CI_est3[i,] <- estimatedsefun(L,tempmean3)[[2]]
CI_est4[i,] <- estimatedsefun(L,tempmean4)[[2]]
CI_est5[i,] <- estimatedsefun(L,tempmean5)[[2]]
  SE_est[i] <- estimatedsefun(L,tempmean)[[1]]
SE_est1[i] <- estimatedsefun(L,tempmean1)[[1]]
SE_est2[i] <- estimatedsefun(L,tempmean2)[[1]]
SE_est3[i] <- estimatedsefun(L,tempmean3)[[1]]
SE_est4[i] <- estimatedsefun(L,tempmean4)[[1]]
SE_est5[i] <- estimatedsefun(L,tempmean5)[[1]]

  if(between(true_mean, CI_est[i,1], CI_est[i,2])){CI_estcovered[i] <- 1}else{CI_estcovered[i] <- 0}
  if(between(true_mean1, CI_est1[i,1], CI_est1[i,2])){CI_estcovered1[i] <- 1}else{CI_estcovered1[i] <- 0}
  if(between(true_mean2, CI_est2[i,1], CI_est2[i,2])){CI_estcovered2[i] <- 1}else{CI_estcovered2[i] <- 0}
  if(between(true_mean3, CI_est3[i,1], CI_est3[i,2])){CI_estcovered3[i] <- 1}else{CI_estcovered3[i] <- 0}
  if(between(true_mean4, CI_est4[i,1], CI_est4[i,2])){CI_estcovered4[i] <- 1}else{CI_estcovered4[i] <- 0}
  if(between(true_mean5, CI_est5[i,1], CI_est5[i,2])){CI_estcovered5[i] <- 1}else{CI_estcovered5[i] <- 0}
  ## track progress of code
  print(i)
}

## end of runtime
t1 <- Sys.time()

ybarmat <- cbind(y_bar, y_bar1, y_bar2, y_bar3, y_bar4, y_bar5)
truemeanmat <- c(true_mean,true_mean1,true_mean2,true_mean3,true_mean4,true_mean5)
labs <- c("Overall", "group 1", "group2", "group3", "group4", "group5")
CImat <- list(CI,CI1,CI2,CI3,CI4,CI5)
CIcovmat <- cbind(CI_covered,CI_covered1,CI_covered2,CI_covered3,CI_covered4,CI_covered5)
CIestmat <- list(CI_est,CI_est1,CI_est2,CI_est3,CI_est4,CI_est5)
CIestcovmat <- cbind(CI_estcovered,CI_estcovered1,CI_estcovered2,CI_estcovered3,CI_estcovered4,CI_estcovered5)
SEestrow <- c(mean(SE_est),mean(SE_est1),mean(SE_est2),mean(SE_est3),mean(SE_est4),mean(SE_est5))
bayesvarmat <- cbind(bayesvar, bayesvar1, bayesvar2, bayesvar3, bayesvar4, bayesvar5)

resmatfun <- function(ybarmat, truemean, CI, CIcovered,SEestrow,labs,bayesvarmat){
  errormat <- t(t(ybarmat)-truemean)
  bias <- colMeans(errormat)
  rMSE <- sqrt(colMeans(errormat^2))
  bayesse <- sqrt(colMeans(bayesvarmat))
  SE <- sqrt(apply(ybarmat,2, var))
  CI_length <- vector(length = ncol(ybarmat))
  CIestlength <- vector(length =ncol(ybarmat))
  for(i in 1:ncol(ybarmat)){
    CItemp <- CI[[i]]
    CI_length[i] <- mean(CItemp[,2]-CItemp[,1])
    CIesttemp <- CIestmat[[i]]
    CIestlength[i] <- mean(CIesttemp[,2]-CIesttemp[,1])

  }
  coverage_rate <- colMeans(CIcovered)
  estcoverage_rate <- colMeans(CIestcovmat)
  res <-rbind(rMSE, bias, SE, SEestrow,bayesse, CI_length, CIestlength, coverage_rate,estcoverage_rate)
  colnames(res) <- labs
  return(res)
}
resultmat <- resmatfun(ybarmat = ybarmat, truemean = truemeanmat, 
                       CI = CImat, CIcovered = CIcovmat, SEestrow=SEestrow,labs = labs,
                       bayesvarmat=bayesvarmat)

# PRINT RESULTS
resultmat

# # WRITE RESULTS
# write.csv(resultmat, "res_wfpbby.csv")
