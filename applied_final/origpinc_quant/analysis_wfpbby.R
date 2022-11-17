#=====================================#
# WFPBB direct imputation applications
#=====================================#
set.seed(51)
#---- SETUP -----
## load packages
require(dplyr)
require(polyapost)
require(LaplacesDemon)

#-- source cleaned data and group indices
type = "origpinc_quantpinc"
source(paste0("BASELINE_", type, ".R"))
N <- sum(acs$perwt)
M <- length(unique(samp$sampled_x1_label))
n <- nrow(samp)
L <- 5000

Tfact <- 30
F_draw <- 20
fpbb_temp <- rep(0, ncol=n*Tfact)

#-- helper functions
`%notin%` <- Negate(`%in%`)
normalize <- function(vec){
  res <- vec
  res[!is.na(vec)] <- vec[!is.na(vec)]/sum(vec[!is.na(vec)])
  return(res)
}
resmatfun <- function(x, id=seq(1,ncol(x))){
  estx <- colMeans(x,na.rm=T)
  varx <- apply(x,2,var,na.rm=T)
  varxt <- varx *(1+1/L)
  CIlower <- apply(x,2,quantile,0.025,na.rm=T)
  CIupper <- apply(x,2,quantile,0.975,na.rm=T)
  CIlowerest <- estx - qt(0.975,df=L)*sqrt(varxt)
  CIupperest <-  estx + qt(0.975,df=L)*sqrt(varxt)
  resmat <- data.frame(
    Estimate = estx,
    SE = sqrt(varx),
    CIlower = CIlower,
    CIupper = CIupper,
    CIlowerest = CIlowerest,
    CIupperest = CIupperest
  )
  rownames(resmat) <- id
  return(resmat)
}

# initialize containers for results
ybar <- rep(0, L) # overall
ybar1 <- rep(0, L) # group 1
ybar2 <- rep(0, L) # group 2
ybar3 <- rep(0, L) # group 3
ybar4 <- rep(0, L) # group 4


#----- BEGIN SIMULATIONS ----
## beginning of runtime
t0 <- Sys.time()

samp$foodindsev <- as.numeric(samp$foodindsev)-1

#------ DRAW SYNTHETIC POPULATIONS ----
for(l in 1:L){
  
  ## (1) Bootstrap from parent sample
  samp$repwts <- NULL
  samp$repwts <- n*normalize(BayesianBootstrap(c(1:n),1))*samp$wts
  boottbl <- 
    samp %>%
    group_by(foodindsev,J_cell, grp1id, grp2id, grp3id, grp4id) %>%
    summarise(bootwts = sum(repwts))
  boottbl <- boottbl[boottbl$bootwts !=0,]
  wts_new <- (n*Tfact)*normalize(boottbl$bootwts)
  bootsamp_size <- nrow(boottbl)
  boot_cats <- 1:bootsamp_size
  
  ## (2) WFPBB with new weights and unique indices
  bbmean <- bbmean1 <- bbmean2 <- bbmean3 <- bbmean4  <- 0
  
  for(f in 1:F_draw){
    # fpbb_temp <- wtpolyap(ysamp = boot_cats, wts = twt, k = (n*Tfact)-bootsamp_size)
    fpbb_temp <- wtpolyap(ysamp = boot_cats, wts = wts_new, k = (n*Tfact)-bootsamp_size)
    # tempdf <- samp[fpbb_temp, c("foodindsev", "grp1id","grp2id","grp3id","grp4id")]
    synthpop <- boottbl[fpbb_temp,]
    bbmean <- bbmean + mean(synthpop$foodindsev)/F_draw
    bbmean1 <- bbmean1 + mean(synthpop$foodindsev[synthpop$grp1id ==1])/F_draw
    bbmean2 <- bbmean2 + mean(synthpop$foodindsev[synthpop$grp2id ==1])/F_draw
    bbmean3 <- bbmean3 + mean(synthpop$foodindsev[synthpop$grp3id ==1])/F_draw
    bbmean4 <- bbmean4 + mean(synthpop$foodindsev[synthpop$grp4id ==1])/F_draw
  }
  
  ybar[l] <- mean(bbmean)
  ybar1[l] <- mean(bbmean1)
  ybar2[l] <- mean(bbmean2)
  ybar3[l] <- mean(bbmean3)
  ybar4[l] <- mean(bbmean4)
}  
#-- (4) Print results
resmatid <- c("Overall","Group 1","Group 2","Group 3","Group 4")

allgrp_wfpbby <- cbind(ybar, ybar1, ybar2, ybar3, ybar4)
resmat_wfpbby <- resmatfun(allgrp_wfpbby,id = resmatid)
write.csv(resmat_wfpbby, paste0("results/BASELINE_", type,"_wfpbby.csv"))
