# WFPBB DIRECT IMPUTATION OF Y -------
# purpose: update containers with WFPBBY for each iteration in the simulation
# input: 
# sample
# output: 
numgrps=5
labelvec <- c(NULL,seq(1,numgrps))

# create containers
containerfun(basename="y_bar_wmrp", subgroups=numgrps, type="vector",dims = iter)
containerfun(basename="wmrp_posterior", subgroups=numgrps, type="vector",dims = sim)
containerfun(basename="CI_wmrp", subgroups=numgrps, type="matrix",dims = c(iter,2))
containerfun(basename="CI_covered_wmrp",subgroups=numgrps, type="vector",dims = iter)
containerfun(basename="bayesvar_wmrp", subgroups=numgrps,type="vector",dims=iter)


wfpbb.mrp <- function(){
  
  containerfun(basename="wmrp_posterior", subgroups=numgrps, type="vector",dims = sim)
  
  for(s in 1:sim){
    
    #-- bootstrap resample
    bootstrap<-BayesianBootstrap(c(1:n),1)
    resample<-rmultinom(1,n,bootstrap)
    bootind <- rep(1:n, resample)
    
    boottbl <- sampled[bootind,] %>%
      group_by(J_cell) %>%
      summarise(bootcts = n(), bootwts = sum(wts))
    
    bootsamp_size <- nrow(boottbl)
    wts_new <- N*normalize(boottbl$bootwts)
    boot_cats <- boottbl$J_cell
    
    ## find Nj via wfpbb with the bootstrapped categories and new weights
    Nj_hat <- matrix(0,nrow = F_draw, ncol = J)
    
    for(f in 1:F_draw){
      temp <- wtpolyap(ysamp = boot_cats, wts = wts_new, k = N-bootsamp_size)
      Nj_hat[f, boot_cats] <- table(c(1:J)[temp]) 
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
    
    wmrp_posterior[s] <- crossprod(wts_mrp, current_cellmean)
    wmrp_posterior1[s] <- crossprod(wts_mrp1, current_cellmean)
    wmrp_posterior2[s] <- crossprod(wts_mrp2, current_cellmean)
    wmrp_posterior3[s] <- crossprod(wts_mrp3, current_cellmean)
    wmrp_posterior4[s] <- crossprod(wts_mrp4, current_cellmean)
    wmrp_posterior5[s] <- crossprod(wts_mrp5, current_cellmean)
    
    
  }
  
  for(k in 1:(numgrps+1)){
    
  }
  y_bar[i] <- mean(wmrp_posterior,na.rm=T)
  y_bar1[i] <- mean(wmrp_posterior1,na.rm=T)
  y_bar2[i] <- mean(wmrp_posterior2,na.rm=T)
  y_bar3[i] <- mean(wmrp_posterior3,na.rm=T)
  y_bar4[i] <- mean(wmrp_posterior4,na.rm=T)
  y_bar5[i] <- mean(wmrp_posterior5,na.rm=T)
  bayesvar[i] <- var(wmrp_posterior,na.rm=T)
  bayesvar1[i] <- var(wmrp_posterior1,na.rm=T)
  bayesvar2[i] <- var(wmrp_posterior2,na.rm=T)
  bayesvar3[i] <- var(wmrp_posterior3,na.rm=T)
  bayesvar4[i] <- var(wmrp_posterior4,na.rm=T)
  bayesvar5[i] <- var(wmrp_posterior5,na.rm=T)
  
  ## collect stats from cellmeans
  # # #- cell picked?
  # cellpicked[sampled_J] <- cellpicked[sampled_J] + 1
  #- bias, SE, rMSE
  cellmeanbias[sampled_J] <- cellmeanbias[sampled_J] + colMeans(cellmeans_stan) - popmeans[sampled_J]
  cellmeanrmse[sampled_J] <- cellmeanrmse[sampled_J] + (colMeans(cellmeans_stan) - popmeans[sampled_J])^2
  
  #- coverage
  CI_cellmean <- apply(cellmeans_stan, 2, quantile, c(0.025, 0.975))
  CIlength_cellmean[i,sampled_J] <- apply(CI_cellmean,2,diff)
  CI_cellmeancov[i, sampled_J] <- ((popmeans[sampled_J] >= CI_cellmean[1,]) & (popmeans[sampled_J] <= CI_cellmean[2,]))
  
  
  ## collect stats from Nj
  if(0 %in% cellpickedbb[sampled_J]){
    picked_discount[which(cellpickedbb == 0)] <- picked_discount[which(cellpickedbb == 0)]+ 1
    newsampJ <- sampled_J[-which(sampled_J==which(cellpickedbb == 0))]
    Njbias[newsampJ] <- Njbias[newsampJ] + colMeans(Njesttemp[,newsampJ], na.rm=T)-popcts[newsampJ]
    Njrmse[newsampJ] <- Njrmse[newsampJ] + (colMeans(Njesttemp[,newsampJ], na.rm=T)-popcts[newsampJ])^2
    #- coverage
    CI_Nj <- apply(Njesttemp[,newsampJ], 2, quantile, c(0.025, 0.975))
    CIlength_Nj[i,newsampJ] <- apply(CI_Nj,2,diff)
    CI_Nj[i,newsampJ] <- ((popcts[newsampJ] >= CI_Nj[1,]) & (popcts[newsampJ] <= CI_Nj[2,]))
    # print("bb missed one!")
  }else{
    Njbias[sampled_J] <- Njbias[sampled_J] + colMeans(Njesttemp[,sampled_J], na.rm=T)-popcts[sampled_J]
    Njrmse[sampled_J] <- Njrmse[sampled_J] + (colMeans(Njesttemp[,sampled_J], na.rm=T)-popcts[sampled_J])^2
    #- coverage
    CI_Nj <- apply(Njesttemp[,sampled_J], 2, quantile, c(0.025, 0.975), na.rm=T)
    CIlength_Nj[i,sampled_J] <- apply(CI_Nj,2,diff)
    CI_Njcov[i,sampled_J] <- ((popcts[sampled_J] >= CI_Nj[1,]) & (popcts[sampled_J] <= CI_Nj[2,]))
  }
  ## collect general stats
  CI_Nj <- apply(Njesttemp[,sampled_J], 2, quantile, c(0.025, 0.975), na.rm=T)
  CIlength_Nj[i,sampled_J] <- apply(CI_Nj,2,diff)
  CI_Njcov[i,sampled_J] <- ((popcts[sampled_J] >= CI_Nj[1,]) & (popcts[sampled_J] <= CI_Nj[2,]))
  CI[i,] <- quantile(wmrp_posterior, c(0.025, 0.975),na.rm=T)
  CI1[i,] <- quantile(wmrp_posterior1, c(0.025, 0.975),na.rm=T)
  CI2[i,] <- quantile(wmrp_posterior2, c(0.025, 0.975),na.rm=T)
  CI3[i,] <- quantile(wmrp_posterior3, c(0.025, 0.975),na.rm=T)
  CI4[i,] <- quantile(wmrp_posterior4, c(0.025, 0.975),na.rm=T)
  CI5[i,] <- quantile(wmrp_posterior5, c(0.025, 0.975),na.rm=T)
  
  if(between(true_mean, CI[i,1], CI[i,2])){CI_covered[i] <- 1}else if(sum(is.na(CI[i,])!=0)){CI_covered[i] <- NA}else{CI_covered[i] <- 0}
  if(between(true_mean1, CI1[i,1], CI1[i,2])){CI_covered1[i] <- 1} else{CI_covered1[i] <- 0}
  if(between(true_mean2, CI2[i,1], CI2[i,2])){CI_covered2[i] <- 1} else{CI_covered2[i] <- 0}
  if(between(true_mean3, CI3[i,1], CI3[i,2])){CI_covered3[i] <- 1} else{CI_covered3[i] <- 0}
  if(between(true_mean4, CI4[i,1], CI4[i,2])){CI_covered4[i] <- 1} else{CI_covered4[i] <- 0}
  if(between(true_mean5, CI5[i,1], CI5[i,2])){CI_covered5[i] <- 1} else{CI_covered5[i] <- 0}
  
}