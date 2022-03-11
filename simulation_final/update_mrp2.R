containerfun(basename="y_bar_mrp2", subgroups=grps, type="vector",dims = iter)
containerfun(basename="CI_mrp2", subgroups=grps, type="matrix",dims = c(iter,2))
containerfun(basename="CI_covered_mrp2",subgroups=grps, type="vector",dims = iter)
Njbias_mrp2 <- rep(0,J)
Njrmse_mrp2 <- rep(0,J)
CIlength_Nj_mrp2 <-  rep(0,J)
CI_Njcov_mrp2 <-  rep(0,J)


update.mrp2 <- function(){
  
  # Nj_mrp2 are the estimated cell frequencies
  mrp_posterior <- mrp_posterior1<- mrp_posterior2 <- mrp_posterior3 <- mrp_posterior4 <-mrp_posterior5 <- rep(NA, sim)
  zpopsampj <- zpopcts_lengthJ[sampled_J]
  
  for(s in 1:sim){
    wts_mrp <- normalize(Nj_mrp2[s,] * zpopsampj)
    
    current_cellmean <- cellmeans_emrp[s,]
    wts_mrp1 <- normalize(wts_mrp*cmat1[sampled_J])
    wts_mrp2 <- normalize(wts_mrp*cmat2[sampled_J])
    wts_mrp3 <- normalize(wts_mrp*cmat3[sampled_J])
    wts_mrp4 <- normalize(wts_mrp*cmat4[sampled_J]) 
    wts_mrp5 <- normalize(wts_mrp*cmat5[sampled_J])
    
    mrp_posterior[s] <- crossprod(wts_mrp, current_cellmean)
    mrp_posterior1[s] <- crossprod(wts_mrp1, current_cellmean)
    mrp_posterior2[s] <- crossprod(wts_mrp2, current_cellmean)
    mrp_posterior3[s] <- crossprod(wts_mrp3, current_cellmean)
    mrp_posterior4[s] <- crossprod(wts_mrp4, current_cellmean)
    mrp_posterior5[s] <- crossprod(wts_mrp5, current_cellmean)
  }
  y_bar_mrp2[i] <<- mean(mrp_posterior,na.rm=T)
  y_bar_mrp21[i] <<- mean(mrp_posterior1,na.rm=T)
  y_bar_mrp22[i] <<- mean(mrp_posterior2,na.rm=T)
  y_bar_mrp23[i] <<- mean(mrp_posterior3,na.rm=T)
  y_bar_mrp24[i] <<- mean(mrp_posterior4,na.rm=T)
  y_bar_mrp25[i] <<- mean(mrp_posterior5,na.rm=T)
  
  # Nj summary
  Njbias_mrp2[sampled_J] <<- Njbias_mrp2[sampled_J] + colMeans(Nj_mrp2*zpopsampj, na.rm=T)-popcts[sampled_J]
  Njrmse_mrp2[sampled_J] <<- Njrmse_mrp2[sampled_J] + (colMeans(Nj_mrp2*zpopsampj, na.rm=T)-popcts[sampled_J])^2
  CI_mrp2_Nj <- apply(Nj_mrp2*zpopsampj, 2, quantile, c(0.025, 0.975), na.rm=T)
  CIlength_Nj_mrp2[sampled_J] <<- CIlength_Nj_mrp2[sampled_J] + apply(CI_mrp2_Nj,2,diff)
  CI_Njcov_mrp2[sampled_J] <<- CI_Njcov_mrp2[sampled_J]+ ((popcts[sampled_J] >= CI_mrp2_Nj[1,]) & (popcts[sampled_J] <= CI_mrp2_Nj[2,]))
  
  CI_mrp2[i,] <<- quantile(mrp_posterior, c(0.025, 0.975),na.rm=T)
  CI_mrp21[i,] <<- quantile(mrp_posterior1, c(0.025, 0.975),na.rm=T)
  CI_mrp22[i,] <<- quantile(mrp_posterior2, c(0.025, 0.975),na.rm=T)
  CI_mrp23[i,] <<- quantile(mrp_posterior3, c(0.025, 0.975),na.rm=T)
  CI_mrp24[i,] <<- quantile(mrp_posterior4, c(0.025, 0.975),na.rm=T)
  CI_mrp25[i,] <<- quantile(mrp_posterior5, c(0.025, 0.975),na.rm=T)
  
  # sampsub <<- unique(sampled$subgroup)
  if(is.na(between(true_mean, CI_mrp2[i,1], CI_mrp2[i,2]))){break}
  if(between(true_mean, CI_mrp2[i,1], CI_mrp2[i,2])){CI_covered_mrp2[i] <<- 1}else if(sum(is.na(CI_mrp2[i,])!=0)){CI_covered_mrp2[i] <<- NA}else{CI_covered_mrp2[i] <<- 0}
  if(between(true_mean1, CI_mrp21[i,1], CI_mrp21[i,2])){CI_covered_mrp21[i] <<- 1} else{CI_covered_mrp21[i] <<- 0}
  if(between(true_mean2, CI_mrp22[i,1], CI_mrp22[i,2])){CI_covered_mrp22[i] <<- 1} else{CI_covered_mrp22[i] <<- 0}
  if(between(true_mean3, CI_mrp23[i,1], CI_mrp23[i,2])){CI_covered_mrp23[i] <<- 1} else{CI_covered_mrp23[i] <<- 0}
  if(between(true_mean4, CI_mrp24[i,1], CI_mrp24[i,2])){CI_covered_mrp24[i] <<- 1} else{CI_covered_mrp24[i] <<- 0}
  if(between(true_mean5, CI_mrp25[i,1], CI_mrp25[i,2])){CI_covered_mrp25[i] <<- 1} else{CI_covered_mrp25[i] <<- 0}
  #-- remove stanfit object to prevent segfault and conserve memory
  # rm(list=c(modelname[i], nqmodelname[i], "cellmeans_stan", "Nj_mrp2", "mrp_posterior","mrp_posterior1","mrp_posterior2","mrp_posterior3","mrp_posterior4",
  #           "mrp_posterior5"))
  
}