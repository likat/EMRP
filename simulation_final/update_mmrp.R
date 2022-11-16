containerfun(basename="y_bar_mmrp", subgroups=grps, type="vector",dims = iter)
containerfun(basename="CI_mmrp", subgroups=grps, type="matrix",dims = c(iter,2))
containerfun(basename="CI_covered_mmrp",subgroups=grps, type="vector",dims = iter)
Njbias_mmrp <- rep(0,J)
Njrmse_mmrp <- rep(0,J)
CIlength_Nj_mmrp <-  rep(0,J)
CI_Njcov_mmrp <-  rep(0,J)

update.mmrp <- function(){
  #-- calculate weights
  zrelfreq <- matrix(0,nrow = Zcat, ncol=2)
  zrelfreq[sampled_z,] <- t(table(sampled$X.1, sampled$Z_categorical))*(as.vector(1/table(sampled$Z_categorical)))

  #-- extract estimates of cell means from the stan model
  wts_mrp <- matrix(0,nrow = Zcat,ncol=2,byrow=F)
  Njesttemp <- matrix(NA,nrow=sim, ncol=J)
  mrp_posterior <- mrp_posterior1<- mrp_posterior2 <- mrp_posterior3 <- mrp_posterior4 <-rep(NA, sim)
  
  for(s in 1:sim){
    
    for(m in sampled_z){
      wts_mrp[m,] <- rmultinom(1,zpopcts_lengthZ[m],prob = zrelfreq[m,])
    }
    
    current_cellmean <- cellmeans_emrp[s,]
    wts <- as.vector(t(wts_mrp))[sampled_J]
    Njesttemp[s,sampled_J] <- wts
    wts <- normalize(wts)
    wts1 <- normalize((wts*cmat1[sampled_J]))
    wts2 <- normalize((wts*cmat2[sampled_J]))
    wts3 <- normalize((wts*cmat3[sampled_J]))
    wts4 <- normalize((wts*cmat4[sampled_J]))

    mrp_posterior[s] <- crossprod(wts, current_cellmean)
    mrp_posterior1[s] <- crossprod(wts1, current_cellmean)
    mrp_posterior2[s] <- crossprod(wts2, current_cellmean)
    mrp_posterior3[s] <- crossprod(wts3, current_cellmean)
    mrp_posterior4[s] <- crossprod(wts4, current_cellmean)

  }
  y_bar_mmrp[i] <<- mean(mrp_posterior,na.rm=T)
  y_bar_mmrp1[i] <<- mean(mrp_posterior1,na.rm=T)
  y_bar_mmrp2[i] <<- mean(mrp_posterior2,na.rm=T)
  y_bar_mmrp3[i] <<- mean(mrp_posterior3,na.rm=T)
  y_bar_mmrp4[i] <<- mean(mrp_posterior4,na.rm=T)

  # Nj summary
  Njbias_mmrp[sampled_J] <<- Njbias_mmrp[sampled_J] + colMeans(Njesttemp[,sampled_J], na.rm=T)-popcts[sampled_J]
  Njrmse_mmrp[sampled_J] <<- Njrmse_mmrp[sampled_J] + (colMeans(Njesttemp[,sampled_J], na.rm=T)-popcts[sampled_J])^2
  CI_mmrp_Nj <- apply(Njesttemp[,sampled_J], 2, quantile, c(0.025, 0.975), na.rm=T)
  CIlength_Nj_mmrp[sampled_J] <<- CIlength_Nj_mmrp[sampled_J] + apply(CI_mmrp_Nj,2,diff)
  CI_Njcov_mmrp[sampled_J] <<- CI_Njcov_mmrp[sampled_J]+ ((popcts[sampled_J] >= CI_mmrp_Nj[1,]) & (popcts[sampled_J] <= CI_mmrp_Nj[2,]))
  
  CI_mmrp[i,] <<- quantile(mrp_posterior, c(0.025, 0.975),na.rm=T)
  CI_mmrp1[i,] <<- quantile(mrp_posterior1, c(0.025, 0.975),na.rm=T)
  CI_mmrp2[i,] <<- quantile(mrp_posterior2, c(0.025, 0.975),na.rm=T)
  CI_mmrp3[i,] <<- quantile(mrp_posterior3, c(0.025, 0.975),na.rm=T)
  CI_mmrp4[i,] <<- quantile(mrp_posterior4, c(0.025, 0.975),na.rm=T)

  if(is.na(between(true_mean, CI_mmrp[i,1], CI_mmrp[i,2]))){break}
  if(between(true_mean, CI_mmrp[i,1], CI_mmrp[i,2])){CI_covered_mmrp[i] <<- 1}else{CI_covered_mmrp[i] <<- 0}
  if(between(true_mean1, CI_mmrp1[i,1], CI_mmrp1[i,2])){CI_covered_mmrp1[i] <<- 1} else{CI_covered_mmrp1[i] <<- 0}
  if(between(true_mean2, CI_mmrp2[i,1], CI_mmrp2[i,2])){CI_covered_mmrp2[i] <<- 1} else{CI_covered_mmrp2[i] <<- 0}
  if(between(true_mean3, CI_mmrp3[i,1], CI_mmrp3[i,2])){CI_covered_mmrp3[i] <<- 1} else{CI_covered_mmrp3[i] <<- 0}
  if(between(true_mean4, CI_mmrp4[i,1], CI_mmrp4[i,2])){CI_covered_mmrp4[i] <<- 1} else{CI_covered_mmrp4[i] <<- 0}

}