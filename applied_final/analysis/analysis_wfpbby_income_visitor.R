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
type = "income_visitor"
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
ybar5 <- rep(0, L) # group 1
ybar6 <- rep(0, L) # group 2
ybar7 <- rep(0, L) # group 3
ybar8 <- rep(0, L) # group 4
ybar9 <- rep(0, L) # group 3
ybar10 <- rep(0, L) # group 4


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
    group_by(foodindsev,J_cell, grp1id, grp2id, grp3id, grp4id, grp5id, grp6id, grp7id, grp8id,
             grp9id, grp10id) %>%
    summarise(bootwts = sum(repwts))
  boottbl <- boottbl[boottbl$bootwts !=0,]
  wts_new <- (n*Tfact)*normalize(boottbl$bootwts)
  bootsamp_size <- nrow(boottbl)
  boot_cats <- 1:bootsamp_size
  
  ## (2) WFPBB with new weights and unique indices
  bbmean <- bbmean1 <- bbmean2 <- bbmean3 <- bbmean4  <- bbmean5 <- bbmean6 <- bbmean7 <- bbmean8  <- bbmean9 <- bbmean10 <- 0
  
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
	bbmean5 <- bbmean5 + mean(synthpop$foodindsev[synthpop$grp5id ==1])/F_draw
    bbmean6 <- bbmean6 + mean(synthpop$foodindsev[synthpop$grp6id ==1])/F_draw
    bbmean7 <- bbmean7 + mean(synthpop$foodindsev[synthpop$grp7id ==1])/F_draw
    bbmean8 <- bbmean8 + mean(synthpop$foodindsev[synthpop$grp8id ==1])/F_draw
    bbmean9 <- bbmean9 + mean(synthpop$foodindsev[synthpop$grp9id ==1])/F_draw
    bbmean10 <- bbmean10 + mean(synthpop$foodindsev[synthpop$grp10id ==1])/F_draw
    
  }
  
  ybar[l] <- mean(bbmean)
  ybar1[l] <- mean(bbmean1)
  ybar2[l] <- mean(bbmean2)
  ybar3[l] <- mean(bbmean3)
  ybar4[l] <- mean(bbmean4)
  ybar5[l] <- mean(bbmean5)
  ybar6[l] <- mean(bbmean6)
  ybar7[l] <- mean(bbmean7)
  ybar8[l] <- mean(bbmean8)
  ybar9[l] <- mean(bbmean9)
  ybar10[l] <- mean(bbmean10)
}  
#-- (4) Print results
resmatid <- c("Overall",grplabs)

allgrp_wfpbby <- cbind(ybar, ybar1, ybar2, ybar3, ybar4,
ybar5, ybar6, ybar7, ybar8,ybar9, ybar10)
resmat_wfpbby <- resmatfun(allgrp_wfpbby,id = resmatid)
write.csv(resmat_wfpbby, paste0("results/BASELINE_", type,"_wfpbby.csv"))
