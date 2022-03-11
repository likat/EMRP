containerfun(basename="y_bar_wmrp", subgroups=grps, type="vector",dims = iter)
containerfun(basename="CI_wmrp", subgroups=grps, type="matrix",dims = c(iter,2))
containerfun(basename="CI_covered_wmrp",subgroups=grps, type="vector",dims = iter)
Njbias_wmrp <- rep(0,J)
Njrmse_wmrp <- rep(0,J)
CIlength_Nj_wmrp <-  rep(0,J)
CI_Njcov_wmrp <-  rep(0,J)

update.wfpbbmrp <- function(){
  # collect weight estimate from each bootstrap iteration
  wts_bootstrapped <- matrix(0,nrow = L, ncol = J)
  wts_mrp <- rep(0, length(sampled_J))
  
  #- clear containers
  Njesttemp <-  matrix(NA,nrow=sim, ncol=J)
    # keep track of cells that weren't picked in BB step
    cellpickedbb <- rep(0,J)
    picked_discount <- rep(0,J)
  mrp_posterior <- mrp_posterior1<- mrp_posterior2 <- mrp_posterior3 <- mrp_posterior4 <-mrp_posterior5 <-rep(NA, sim)
  
  for(s in 1:sim){
    
    #-- bootstrap resample
    bootstrap<-BayesianBootstrap(c(1:n),1)
    resample<-rmultinom(1,n,bootstrap)
    bootind <- rep(1:n, resample)
    
    boottbl <- 
      sampled[bootind,] %>%
      group_by(J_cell) %>%
      summarise(bootcts = n(), bootwts = sum(wts))
    bootsamp_size <- nrow(boottbl)
    wts_new <- N*normalize(boottbl$bootwts)
    boot_cats <- boottbl$J_cell
    
    ## find Nj via wfpbb with the bootstrapped categories and new weights
    Nj_est <- rep(0, J)
    
    for(f in 1:F_draw){
      temp <- wtpolyap(ysamp = boot_cats, wts = wts_new, k = N-bootsamp_size)
      Nj_est[boot_cats] <- Nj_est[boot_cats] + table(c(1:J)[temp])/F_draw
    }
    
    # Nj_est <- colMeans(Nj_hat)
    
    ## collect stats from Nj estimation
    #- bias, SE, rMSE
    cellpickedbb[boot_cats] <- cellpickedbb[boot_cats] +1
    Njesttemp[s,boot_cats] <- Nj_est[boot_cats]
    #- coverage
    
    current_cellmean <- cellmeans_emrp[s,]
    wts_mrp <- normalize(Nj_est[sampled_J])  # extract frequencies relevant to estimation of group mean
    wts_mrp1 <- normalize((Nj_est*cmat1)[sampled_J])
    wts_mrp2 <- normalize((Nj_est*cmat2)[sampled_J])
    wts_mrp3 <- normalize((Nj_est*cmat3)[sampled_J])
    wts_mrp4 <- normalize((Nj_est*cmat4)[sampled_J])
    wts_mrp5 <- normalize((Nj_est*cmat5)[sampled_J])
    
    mrp_posterior[s] <- crossprod(wts_mrp, current_cellmean)
    mrp_posterior1[s] <- crossprod(wts_mrp1, current_cellmean)
    mrp_posterior2[s] <- crossprod(wts_mrp2, current_cellmean)
    mrp_posterior3[s] <- crossprod(wts_mrp3, current_cellmean)
    mrp_posterior4[s] <- crossprod(wts_mrp4, current_cellmean)
    mrp_posterior5[s] <- crossprod(wts_mrp5, current_cellmean)
    
    
  }
  y_bar_wmrp[i] <<- mean(mrp_posterior,na.rm=T)
  y_bar_wmrp1[i] <<- mean(mrp_posterior1,na.rm=T)
  y_bar_wmrp2[i] <<- mean(mrp_posterior2,na.rm=T)
  y_bar_wmrp3[i] <<- mean(mrp_posterior3,na.rm=T)
  y_bar_wmrp4[i] <<- mean(mrp_posterior4,na.rm=T)
  y_bar_wmrp5[i] <<- mean(mrp_posterior5,na.rm=T)


  ## collect stats from Nj
  if(0 %in% cellpickedbb[sampled_J]){
    # if the bootstrap sample did not sample all J cells, then we discount that from 
    # the one-pass bias and rmse calculations
    picked_discount[which(cellpickedbb == 0)] <- 
      picked_discount[which(cellpickedbb == 0)]+ 1
    newsampJ <- 
      sampled_J[-which(sampled_J==which(cellpickedbb == 0))]
    Njbias_wmrp[newsampJ] <<- 
      Njbias_wmrp[newsampJ] + colMeans(Njesttemp[,newsampJ], na.rm=T)-popcts[newsampJ]
    Njrmse_wmrp[newsampJ] <<- 
      Njrmse_wmrp[newsampJ] + (colMeans(Njesttemp[,newsampJ], na.rm=T)-popcts[newsampJ])^2
    CI_Nj_wmrp <- 
      apply(Njesttemp[,newsampJ], 2, quantile, c(0.025, 0.975))
    CIlength_Nj_wmrp[newsampJ] <<-
      CIlength_Nj_wmrp[newsampJ]+ apply(CI_Nj_wmrp,2,diff)
    CI_Njcov_wmrp[newsampJ] <<- 
      CI_Njcov_wmrp[newsampJ] + ((popcts[newsampJ] >= CI_Nj_wmrp[1,]) & (popcts[newsampJ]<=CI_Nj_wmrp[2,]))
  }else{
    Njbias_wmrp[sampled_J] <<- Njbias_wmrp[sampled_J] + colMeans(Njesttemp[,sampled_J], na.rm=T)-popcts[sampled_J]
    Njrmse_wmrp[sampled_J] <<- Njrmse_wmrp[sampled_J] + (colMeans(Njesttemp[,sampled_J], na.rm=T)-popcts[sampled_J])^2
    CI_Nj_wmrp <- apply(Njesttemp[,sampled_J], 2, quantile, c(0.025, 0.975), na.rm=T)
    CIlength_Nj_wmrp[sampled_J] <<- CIlength_Nj_wmrp[sampled_J] + apply(CI_Nj_wmrp,2,diff)
    CI_Njcov_wmrp[sampled_J] <<- CI_Njcov_wmrp[sampled_J] + ((popcts[sampled_J] >= CI_Nj_wmrp[1,]) & (popcts[sampled_J]<=CI_Nj_wmrp[2,]))
  }
  ## collect general stats
  # CI_Nj_wmrp <<- apply(Njesttemp[,sampled_J], 2, quantile, c(0.025, 0.975), na.rm=T)
  # CIlength_Nj_wmrp[i,sampled_J] <<- apply(CI_Nj_wmrp,2,diff)
  # CI_Njcov_wmrp[i,sampled_J] <<- ((popcts[sampled_J] >= CI_Nj_wmrp[1,]) & (popcts[sampled_J] <= CI_Nj_wmrp[2,]))
  CI_wmrp[i,] <<- quantile(mrp_posterior, c(0.025, 0.975),na.rm=T)
  CI_wmrp1[i,] <<- quantile(mrp_posterior1, c(0.025, 0.975),na.rm=T)
  CI_wmrp2[i,] <<- quantile(mrp_posterior2, c(0.025, 0.975),na.rm=T)
  CI_wmrp3[i,] <<- quantile(mrp_posterior3, c(0.025, 0.975),na.rm=T)
  CI_wmrp4[i,] <<- quantile(mrp_posterior4, c(0.025, 0.975),na.rm=T)
  CI_wmrp5[i,] <<- quantile(mrp_posterior5, c(0.025, 0.975),na.rm=T)
  
  if(is.na(between(true_mean, CI_wmrp[i,1], CI_wmrp[i,2]))){break}
  if(between(true_mean, CI_wmrp[i,1], CI_wmrp[i,2])){CI_covered_wmrp[i] <<- 1}else if(sum(is.na(CI_wmrp[i,])!=0)){CI_covered_wmrp[i] <<- NA}else{CI_covered_wmrp[i] <<- 0}
  if(between(true_mean1, CI_wmrp1[i,1], CI_wmrp1[i,2])){CI_covered_wmrp1[i] <<- 1} else{CI_covered_wmrp1[i] <<- 0}
  if(between(true_mean2, CI_wmrp2[i,1], CI_wmrp2[i,2])){CI_covered_wmrp2[i] <<- 1} else{CI_covered_wmrp2[i] <<- 0}
  if(between(true_mean3, CI_wmrp3[i,1], CI_wmrp3[i,2])){CI_covered_wmrp3[i] <<- 1} else{CI_covered_wmrp3[i] <<- 0}
  if(between(true_mean4, CI_wmrp4[i,1], CI_wmrp4[i,2])){CI_covered_wmrp4[i] <<- 1} else{CI_covered_wmrp4[i] <<- 0}
  if(between(true_mean5, CI_wmrp5[i,1], CI_wmrp5[i,2])){CI_covered_wmrp5[i] <<- 1} else{CI_covered_wmrp5[i] <<- 0}
  
}