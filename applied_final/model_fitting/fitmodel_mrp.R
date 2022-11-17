#================================
# Embedded MRP with applications
#================================
#-- last updated: aug 2021 --
#-- uses public use baseline data
#-- just (Z), no X (svefreq2), use MRP to analyze
set.seed(51)
options(mc.cores = 4)
#-- load required packages
require(dplyr)
require(foreign)
require(MASS)
require(polyapost)
require(mvtnorm)
require(rstan)
require(LaplacesDemon)

#-- helper functions
`%notin%` <- Negate(`%in%`)
normalize <- function(vec){
  # does exactly what it's called
  return(vec/sum(vec))
}
resmatfun <- function(x, id=seq(1,ncol(x))){
  estx <- colMeans(x,na.rm=T)
  varx <- apply(x,2,var,na.rm=T)
  CIlower <- apply(x,2,quantile,0.025,na.rm=T)
  CIupper <- apply(x,2,quantile,0.975,na.rm=T)
  resmat <- data.frame(
    Estimate = estx,
    SE = sqrt(varx),
    CIlower = CIlower,
    CIupper = CIupper
  )
  rownames(resmat) <- id
  return(resmat)
}

#-- source cleaned data and group indices
source("BASELINE_origpinc_quantpinc.R")

# initialize containers for results
sim <- 5000
staniters = 10000
M <- length(unique(samp$sampled_x1_label))
N <- sum(acs$perwt)

#===== (1) Run multilevel regression via Stan=====
# parameters for stan model
Ma <- nlevels(samp$age3) # combinations of x1: age, sex, educat, povgap
Mb <- nlevels(samp$sex2)
Mc <- nlevels(samp$race3)
Md <- nlevels(samp$educat4)
Me <- nlevels(samp$povgap4)
Mab <- nlevels(samp$abcat)
Mbc <- nlevels(samp$bccat)
Mcd <- nlevels(samp$cdcat)
Mdx <- nlevels(samp$dfcat)
L <- nlevels(samp$svefreq2)
n <- nrow(samp)

y <- as.numeric(samp$foodindsev)-1

grptbl <- samp %>%
  group_by(age3,sex2,race3,educat4,povgap4,abcat,bccat,cdcat) %>%
  summarise(n(),zcts = mean(x1_popcts), zlab = mean(sampled_x1_label))
grptbl$sampJ <- 1:nrow(grptbl)
samp <- left_join(samp,grptbl)
J_stan <- nrow(grptbl) # if fitting with sampgroup
agroup <- as.numeric(grptbl$age3)
bgroup <- as.numeric(grptbl$sex2) -1
cgroup <- as.numeric(grptbl$race3) 
dgroup <- as.numeric(grptbl$educat4)
egroup <- as.numeric(grptbl$povgap4)
abgroup <- as.numeric(grptbl$abcat)
bcgroup <- as.numeric(grptbl$bccat)
cdgroup <- as.numeric(grptbl$cdcat)
cell_label <- samp$sampled_x1_label
zpopcts <- grptbl$zcts

stanfit <- stan(file = "outcome_model_BASELINE.stan",
                data = c("Ma", "Mb", "Mc", "Md", "Me", "L", "Mab","Mbc","Mcd",
                         "J_stan", "n", "y", 
                         "agroup", "bgroup","cgroup", "dgroup", "egroup",
                         "abgroup","bcgroup","cdgroup","cell_label"),
                iter=staniters,
                warmup=staniters-sim/4,control=list(adapt_delta=0.99),chains=4)

# extract and save estimates of cell means from the stan model
stanpars <-
  rstan::extract(object=stanfit, permuted = TRUE)
saveRDS(stanpars, "mrppars.rds")
