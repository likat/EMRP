# library loading
# Libraries -------
library(dplyr)
library(stringr)
require(rstan)
require(polyapost)
require(LaplacesDemon)

# CONTAINER CREATION FUNCTION------
# purpose: create containers to store results in
# INPUTS:
#   basename <string> = what is the base name of the container?
#   subgroups <real> = number of subgroups we are estimating
#   type <string> = matrix or vector?
#   dim <real> or <vector> = <real> if type is vector, <vector>[2] if type is matrix
# OUTPUT 
#   assign subgroups+1 objects to the global environment
# USAGE
# containerfun(basename="y_bar", subgroups=grps, type="vector",dims = sim)
# containerfun(basename="CI", subgroups=grps, type="matrix",dims = c(sim,2))

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

# EXPIT ------
# purpose: expit function
expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

# DLL unloading -----
## Workaround for fitting many stan models in the cluster
# to fix the dyn.load() unable to load shared object warning message
# adapted from https://github.com/stan-dev/rstan/issues/448
unload.dll <- function(model=model){
  dso_filename = model@dso@dso_filename
  loaded_dlls = getLoadedDLLs()
  if (dso_filename %in% names(loaded_dlls)) {
    message("Unloading DLL for model dso ", dso_filename)
    model.dll = loaded_dlls[[dso_filename]][['path']]
    dyn.unload(model.dll)
  } else {
    message("No loaded DLL for model dso ", dso_filename)
  }
  
  loaded_dlls = getLoadedDLLs()
  loaded_dlls = loaded_dlls[str_detect(names(loaded_dlls), '^file')]
  if (length(loaded_dlls) > 10) {
    for (dll in head(loaded_dlls, -10)) {
      message("Unloading DLL ", dll[['name']], ": ", dll[['path']])
      dyn.unload(dll[['path']])
    }
  }
  message("DLL Count = ", length(getLoadedDLLs()), ": [", str_c(names(loaded_dlls), collapse = ","), "]")
  
}

# NORMALIZE -------
# purpose: normalize a vector to sum to 1
# input: double vector[n]
# output: normalized vector[n]
normalize <- function(vec){
  return(vec/sum(vec))
}

# NOT IN--------
# purpose: boolean operator; are elements in 'vector' not in 'othervector'?
# use: double vector[n] %notin% another double vector[n]
# output: T or F 
`%notin%` <- Negate(`%in%`)



# RESMAT FUNCTION for MRP/EMRP methods ----
resmat.emrp.mrp <- function(label, subgroups=5,type){
  labs <- rep(NA,subgroups+1)
  labelvec <- c(seq(1,subgroups))
  
  # make labels for nj containers 
  Njbiastext <- paste0("Njbias_",label)
  Njrmsetext <- paste0("Njrmse_", label)
  NjCIlengthtext <- paste0("CIlength_Nj_",label)
  Njcovtext <- paste0("CI_Njcov_", label)
  # containers for final output
  ybarmat <- matrix(nrow=iter,ncol=subgroups+1)
  truemeanmat <- c(true_mean,true_mean1,true_mean2,true_mean3,true_mean4,true_mean5)
  CImat <-vector(mode = "list", length = subgroups+1)
  CIcovmat <- matrix(nrow=iter,ncol=subgroups+1)

  labs[1] <- paste0("Overall")
  # get name from the inputs
  ybarext <- paste0("y_bar_",label)
  CIext <- paste0("CI_",label)
  CIcoveredext <- paste0("CI_covered_",label)
  
  # assign to the appropriate containers
  ybarmat[,1] <- get(ybarext)
  CImat[[1]] <- get(CIext)
  CIcovmat[,1] <- get(CIcoveredext)
  
  
  for(j in 1:(subgroups)){ # Change here if subgroups are different from 5
    # create subgroup label
    grplab <- labelvec[j]
    labs[j+1] <- paste0("group",grplab)
    # get name from the inputs
    ybarext <- paste0("y_bar_",label,grplab)
    CIext <- paste0("CI_",label,grplab)
    CIcoveredext <- paste0("CI_covered_",label,grplab)

    # assign to the appropriate containers
    ybarmat[,j+1] <- get(ybarext)
    CImat[[j+1]] <- get(CIext)
    CIcovmat[,j+1] <- get(CIcoveredext)

  }
  # labs <- c("Overall", "group1", "group2", "group3", "group4", "group5")
  errormat <- t(t(ybarmat)-truemeanmat)
  bias <- colMeans(errormat)
  rMSE <- sqrt(colMeans(errormat^2))
  SE <- sqrt(apply(ybarmat,2, var))
  CI_length <- vector(length = ncol(ybarmat))
  for(i in 1:ncol(ybarmat)){
    CItemp <- CImat[[i]]
    CI_length[i] <- mean(CItemp[,2]-CItemp[,1])
  }
  coverage_rate <- colMeans(CIcovmat)
  res <-rbind(rMSE, bias, SE, CI_length, coverage_rate)
  colnames(res) <- labs

  if(label == "mrp"){
    write.csv(res, paste0(type,"/results/","res_",label,".csv"))
  }else{
    ## Nj Summary
    Njbias <- get(Njbiastext) / cellpicked
    Njrmse <- sqrt(get(Njrmsetext)/ cellpicked)
    Njse <- sqrt(Njrmse^2 - Njbias^2)
    
    write.csv(res, paste0(type,"/results/","res_",label,".csv"))
    write.csv(get(NjCIlengthtext)/iter,paste0(type,"/results/","NjCIlength_", label,".csv"))
    write.csv(get(Njcovtext)/iter, paste0(type,"/results/","Njcov_", label,".csv"))
    write.csv(Njbias,paste0(type,"/results/","Njbias_", label,".csv"))
    write.csv(Njse, paste0(type,"/results/","Njse_", label,".csv"))
    write.csv(Njrmse, paste0(type,"/results/","Njrmse_", label,".csv"))
  }
  

}

# RESMAT FUNCTION for WFPBBY -----
resmat.wfpbby <- function(label,subgroups=5,type){
  labs <- rep(NA,subgroups+1)
  labelvec <- c(NULL,seq(1,subgroups))
  
  # containers for final output
  ybarmat <- matrix(nrow=iter,ncol=subgroups+1)
  truemeanmat <- c(true_mean,true_mean1,true_mean2,true_mean3,true_mean4,true_mean5)
  CImat <-vector(mode = "list", length = subgroups+1)
  CIcovmat <- matrix(nrow=iter,ncol=subgroups+1)
  CIestmat <- vector(mode = "list", length = subgroups+1)
  CIestcovmat <- matrix(nrow=iter,ncol=subgroups+1)
  SEestrow <- rep(NA,subgroups+1)
  bayesvarmat <- matrix(nrow=iter,ncol=subgroups+1)
  
  labs[1] <- "Overall"
  ybarext <- paste0("y_bar_",label)
  CIext <- paste0("CI_",label)
  CIcoveredext <- paste0("CI_covered_",label)
  CIestext <- paste0("CI_est_",label)
  CIestcoveredext <- paste0("CI_estcovered_",label)
  SEestext <- paste0("SE_est_",label)
  bayesvarext <- paste0("bayesvar_",label)
  
  # assign to the appropriate containers
  ybarmat[,1] <- get(ybarext)
  CImat[[1]] <- get(CIext)
  CIcovmat[,1] <- get(CIcoveredext)
  CIestmat[[1]] <- get(CIestext)
  CIestcovmat[,1] <- get(CIestcoveredext)
  SEestrow[1] <- mean(get(SEestext))
  bayesvarmat[,1] <- mean(get(bayesvarext))
  
  
  for(j in 1:(subgroups)){ # Change here if subgroups are different from 5
    # create subgroup label
    grplab <- labelvec[j]
    labs[j+1] <- paste0("group",grplab)
    # get name from the inputs
    ybarext <- paste0("y_bar_",label,grplab)
    CIext <- paste0("CI_",label,grplab)
    CIcoveredext <- paste0("CI_covered_",label,grplab)
    CIestext <- paste0("CI_est_",label,grplab)
    CIestcoveredext <- paste0("CI_estcovered_",label,grplab)
    SEestext <- paste0("SE_est_",label,grplab)
    bayesvarext <- paste0("bayesvar_",label,grplab)
    
    # assign to the appropriate containers
    ybarmat[,j+1] <- get(ybarext)
    CImat[[j+1]] <- get(CIext)
    CIcovmat[,j+1] <- get(CIcoveredext)
    CIestmat[[j+1]] <- get(CIestext)
    CIestcovmat[,j+1] <- get(CIestcoveredext)
    SEestrow[j+1] <- mean(get(SEestext))
    bayesvarmat[,j+1] <- mean(get(bayesvarext))

  }
  # labs <- c("Overall", "group1", "group2", "group3", "group4", "group5")
  errormat <- t(t(ybarmat)-truemeanmat)
  bias <- colMeans(errormat)
  rMSE <- sqrt(colMeans(errormat^2))
  bayesse <- sqrt(colMeans(bayesvarmat))
  SE <- sqrt(apply(ybarmat,2, var))
  CI_length <- vector(length = ncol(ybarmat))
  CIestlength <- vector(length =ncol(ybarmat))
  for(i in 1:ncol(ybarmat)){
    CItemp <- CImat[[i]]
    CI_length[i] <- mean(CItemp[,2]-CItemp[,1])
    CIesttemp <- CIestmat[[i]]
    CIestlength[i] <- mean(CIesttemp[,2]-CIesttemp[,1])
    
  }
  coverage_rate <- colMeans(CIcovmat)
  estcoverage_rate <- colMeans(CIestcovmat)
  res <-rbind(rMSE, bias, SE, SEestrow,bayesse, CI_length, CIestlength, coverage_rate,estcoverage_rate)
  colnames(res) <- labs
  # return(res)
  write.csv(res, paste0(type,"/results/","res_",label,".csv"))
  
}

