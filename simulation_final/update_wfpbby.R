# HELPER: ESTIMATED SE FUNCTION FOR WFPBBY
estimatedsefun <- function(L,tempvec){
  CI <- rep(NA,2)
  SE_est <- sqrt((1+1/L)*var(tempvec))
  CI[1] <- mean(tempvec) - SE_est*qt(0.975,L-1)
  CI[2] <- mean(tempvec) + SE_est*qt(0.975,L-1)
  res <- list(SE_est,CI)
  return(res)
}


# create containers for wfpbby 
containerfun(basename="y_bar_wfpbby", subgroups=grps, type="vector",dims = iter)
containerfun(basename="tempmean", subgroups=grps, type="vector",dims = L)
containerfun(basename="CI_wfpbby", subgroups=grps, type="matrix",dims = c(iter,2))
containerfun(basename="CI_est_wfpbby", subgroups=grps, type="matrix",dims = c(iter,2))
containerfun(basename="CI_estcovered_wfpbby", subgroups=grps, type="vector",dims = iter)
containerfun(basename="CI_covered_wfpbby",subgroups=grps, type="vector",dims = iter)
containerfun(basename="SE_est_wfpbby",subgroups=grps, type="vector",dims = iter)
containerfun(basename="bayesvar_wfpbby", subgroups=grps,type="vector",dims=iter)


# WFPBB DIRECT IMPUTATION OF Y -------
# purpose: update containers with WFPBBY for each iteration in the simulation

update.wfpbby <- function(){
  

#------ DRAW SYNTHETIC POPULATIONS 
for(l in 1:L){
  
  ## (1) Bootstrap from parent sample
  bootstrap<-BayesianBootstrap(c(1:n),1)
  wts_new <- bootstrap*sampled$wts
  ind_nonzero <- wts_new!=0
  wts_new <- wts_new[ind_nonzero]
  wts_new <- N * normalize(wts_new)
  boot_cats <- c(1:n)[ind_nonzero]
  bootsamp_size <- sum(ind_nonzero)
  
  bbmean <- bbmean1 <- bbmean2 <- bbmean3 <- bbmean4 <- bbmean5 <- 0
  ## (2) WFPBB with new weights and unique indices
  fpbb_temp <- rep(0, N)
  for(f in 1:F_draw){
    fpbb_temp <- wtpolyap(ysamp = boot_cats, wts = wts_new, k = N-bootsamp_size)
    synthpop <- sampled[fpbb_temp,]
    bbmean <- bbmean + mean(synthpop$Y)/F_draw
    bbmean1 <- bbmean1 + mean(synthpop$Y[synthpop$J_cell %in% group1])/F_draw
    bbmean2 <- bbmean2 + mean(synthpop$Y[synthpop$J_cell %in% group2])/F_draw
    bbmean3 <- bbmean3 + mean(synthpop$Y[synthpop$J_cell %in% group3])/F_draw
    bbmean4 <- bbmean4 + mean(synthpop$Y[synthpop$J_cell %in% group4])/F_draw
  }
  
  
  ## (3) Collect estimate from the lth synthetic population
  tempmean[l] <- bbmean
  tempmean1[l] <- bbmean1
  tempmean2[l] <- bbmean2
  tempmean3[l] <- bbmean3
  tempmean4[l] <- bbmean4
  #print(l)
}  
## (4) Final estimate for parent sample i
y_bar_wfpbby[i] <<- mean(tempmean,na.rm=T)
y_bar_wfpbby1[i] <<- mean(tempmean1,na.rm=T)
y_bar_wfpbby2[i] <<- mean(tempmean2,na.rm=T)
y_bar_wfpbby3[i] <<- mean(tempmean3,na.rm=T)
y_bar_wfpbby4[i] <<- mean(tempmean4,na.rm=T)

## (5) Confidence Interval
## (5a) Empirical CI
CI_wfpbby[i,] <<- quantile(tempmean, probs = c(0.025, 0.975),na.rm=T)
CI_wfpbby1[i,] <<- quantile(tempmean1, probs = c(0.025, 0.975),na.rm=T)
CI_wfpbby2[i,] <<- quantile(tempmean2, probs = c(0.025, 0.975),na.rm=T)
CI_wfpbby3[i,] <<- quantile(tempmean3, probs = c(0.025, 0.975),na.rm=T)
CI_wfpbby4[i,] <<- quantile(tempmean4, probs = c(0.025, 0.975),na.rm=T)

if(between(true_mean, CI_wfpbby[i,1], CI_wfpbby[i,2])){CI_covered_wfpbby[i] <<- 1}else{CI_covered_wfpbby[i] <<- 0}
if(between(true_mean1, CI_wfpbby1[i,1], CI_wfpbby1[i,2])){CI_covered_wfpbby1[i] <<- 1}else{CI_covered_wfpbby1[i] <<- 0}
if(between(true_mean2, CI_wfpbby2[i,1], CI_wfpbby2[i,2])){CI_covered_wfpbby2[i] <<- 1}else{CI_covered_wfpbby2[i] <<- 0}
if(between(true_mean3, CI_wfpbby3[i,1], CI_wfpbby3[i,2])){CI_covered_wfpbby3[i] <<- 1}else{CI_covered_wfpbby3[i] <<- 0}
if(between(true_mean4, CI_wfpbby4[i,1], CI_wfpbby4[i,2])){CI_covered_wfpbby4[i] <<- 1}else{CI_covered_wfpbby4[i] <<- 0}

## (5b) Estimated CI
bayesvar_wfpbby[i] <<- var(tempmean)
bayesvar_wfpbby1[i] <<- var(tempmean1)
bayesvar_wfpbby2[i] <<- var(tempmean2)
bayesvar_wfpbby3[i] <<- var(tempmean3)
bayesvar_wfpbby4[i] <<- var(tempmean4)
CI_est_wfpbby[i,] <<- estimatedsefun(L,tempmean)[[2]]
CI_est_wfpbby1[i,] <<- estimatedsefun(L,tempmean1)[[2]]
CI_est_wfpbby2[i,] <<- estimatedsefun(L,tempmean2)[[2]]
CI_est_wfpbby3[i,] <<- estimatedsefun(L,tempmean3)[[2]]
CI_est_wfpbby4[i,] <<- estimatedsefun(L,tempmean4)[[2]]
SE_est_wfpbby[i] <<- estimatedsefun(L,tempmean)[[1]]
SE_est_wfpbby1[i] <<- estimatedsefun(L,tempmean1)[[1]]
SE_est_wfpbby2[i] <<- estimatedsefun(L,tempmean2)[[1]]
SE_est_wfpbby3[i] <<- estimatedsefun(L,tempmean3)[[1]]
SE_est_wfpbby4[i] <<- estimatedsefun(L,tempmean4)[[1]]

if(between(true_mean, CI_est_wfpbby[i,1], CI_est_wfpbby[i,2])){CI_estcovered_wfpbby[i] <<- 1}else{CI_estcovered_wfpbby[i] <<- 0}
if(between(true_mean1, CI_est_wfpbby1[i,1], CI_est_wfpbby1[i,2])){CI_estcovered_wfpbby1[i] <<- 1}else{CI_estcovered_wfpbby1[i] <<- 0}
if(between(true_mean2, CI_est_wfpbby2[i,1], CI_est_wfpbby2[i,2])){CI_estcovered_wfpbby2[i] <<- 1}else{CI_estcovered_wfpbby2[i] <<- 0}
if(between(true_mean3, CI_est_wfpbby3[i,1], CI_est_wfpbby3[i,2])){CI_estcovered_wfpbby3[i] <<- 1}else{CI_estcovered_wfpbby3[i] <<- 0}
if(between(true_mean4, CI_est_wfpbby4[i,1], CI_est_wfpbby4[i,2])){CI_estcovered_wfpbby4[i] <<- 1}else{CI_estcovered_wfpbby4[i] <<- 0}
} 
