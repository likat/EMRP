containerfun(basename="y_bar_mrp", subgroups=grps, type="vector",dims = iter)
containerfun(basename="CI_mrp", subgroups=grps, type="matrix",dims = c(iter,2))
containerfun(basename="CI_covered_mrp",subgroups=grps, type="vector",dims = iter)

update.mrp <- function(){
  
    wts<- zpopcts_lengthZ
    mrp_posterior <- mrp_posterior1<- mrp_posterior2 <- mrp_posterior3 <- mrp_posterior4 <-mrp_posterior5 <-rep(NA, sim)
    
    wts_mrp <- normalize(wts[sampled_z])
    wts_mrp1 <- normalize((wts_mrp*cmat1z[sampled_z]))
    wts_mrp2 <- normalize((wts_mrp*cmat2z[sampled_z]))
    wts_mrp3 <- normalize((wts_mrp*cmat3z[sampled_z]))
    wts_mrp4 <- normalize((wts_mrp*cmat4z[sampled_z]))
    wts_mrp5 <- normalize((wts_mrp*cmat5z[sampled_z]))
    
    for(s in 1:sim){
      current_cellmean <- cellmeans_mrp[s,]
      mrp_posterior[s] <- crossprod(wts_mrp, current_cellmean)
      mrp_posterior1[s] <- crossprod(wts_mrp1, current_cellmean)
      mrp_posterior2[s] <- crossprod(wts_mrp2, current_cellmean)
      mrp_posterior3[s] <- crossprod(wts_mrp3, current_cellmean)
      mrp_posterior4[s] <- crossprod(wts_mrp4, current_cellmean)
      mrp_posterior5[s] <- crossprod(wts_mrp5, current_cellmean)
      
    }
    y_bar_mrp[i] <<- mean(mrp_posterior,na.rm=T)
    y_bar_mrp1[i] <<- mean(mrp_posterior1,na.rm=T)
    y_bar_mrp2[i] <<- mean(mrp_posterior2,na.rm=T)
    y_bar_mrp3[i] <<- mean(mrp_posterior3,na.rm=T)
    y_bar_mrp4[i] <<- mean(mrp_posterior4,na.rm=T)
    y_bar_mrp5[i] <<- mean(mrp_posterior5,na.rm=T)
    CI_mrp[i,] <<- quantile(mrp_posterior, c(0.025, 0.975),na.rm=T)
    CI_mrp1[i,] <<- quantile(mrp_posterior1, c(0.025, 0.975),na.rm=T)
    CI_mrp2[i,] <<- quantile(mrp_posterior2, c(0.025, 0.975),na.rm=T)
    CI_mrp3[i,] <<- quantile(mrp_posterior3, c(0.025, 0.975),na.rm=T)
    CI_mrp4[i,] <<- quantile(mrp_posterior4, c(0.025, 0.975),na.rm=T)
    CI_mrp5[i,] <<- quantile(mrp_posterior5, c(0.025, 0.975),na.rm=T)
    
    # sampsub <- unique(sampled$subgroup)
    if(is.na(between(true_mean, CI_mrp[i,1], CI_mrp[i,2]))){break}
    if(between(true_mean, CI_mrp[i,1], CI_mrp[i,2])){CI_covered_mrp[i] <<- 1}else if(sum(is.na(CI_mrp[i,])!=0)){CI_covered_mrp[i] <<- NA}else{CI_covered_mrp[i] <<- 0}
    if(between(true_mean1, CI_mrp1[i,1], CI_mrp1[i,2])){CI_covered_mrp1[i] <<- 1} else{CI_covered_mrp1[i] <<- 0}
    if(between(true_mean2, CI_mrp2[i,1], CI_mrp2[i,2])){CI_covered_mrp2[i] <<- 1} else{CI_covered_mrp2[i] <<- 0}
    if(between(true_mean3, CI_mrp3[i,1], CI_mrp3[i,2])){CI_covered_mrp3[i] <<- 1} else{CI_covered_mrp3[i] <<- 0}
    if(between(true_mean4, CI_mrp4[i,1], CI_mrp4[i,2])){CI_covered_mrp4[i] <<- 1} else{CI_covered_mrp4[i] <<- 0}
    if(between(true_mean5, CI_mrp5[i,1], CI_mrp5[i,2])){CI_covered_mrp5[i] <<- 1} else{CI_covered_mrp5[i] <<- 0}
    #-- remove stanfit object to prevent segfault and conserve memory
    # rm(list=c(modelname[i],"cellmeans_mrp",
    #           "mrp_posterior","mrp_posterior1", "mrp_posterior2","mrp_posterior3","mrp_posterior4","mrp_posterior5"))
    
}