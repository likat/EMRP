#================================
# Embedded MRP with applications
#================================
#-- last updated: nov 2022 --
#-- uses public use LSW data

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
Tfact <- 30
F_draw = 20
M <- length(unique(samp$sampled_x1_label))
N <- sum(acs$perwt)
wmrp_Njhatmat <- matrix(0,nrow=sim, ncol=J)
mmrp_Njhatmat <- matrix(0,nrow=sim, ncol=J)
mrp2_Njhatmat <- matrix(0,nrow=sim, ncol=J)
#== WFPBB-MRP
wmrp_posterior <- rep(0,sim)
wmrp_posterior1 <- rep(0, sim)
wmrp_posterior2 <- rep(0, sim)
wmrp_posterior3 <- rep(0, sim)
wmrp_posterior4 <- rep(0, sim)

#== Multinomial-MRP
mmrp_posterior <- mmrp_posterior1 <- mmrp_posterior2 <- mmrp_posterior3 <- mmrp_posterior4 <- rep(0,sim)

#== Two-stage MRP
mrp2_posterior <- mrp2_posterior1 <- mrp2_posterior2 <- mrp2_posterior3 <- mrp2_posterior4 <- rep(0,sim)

#===== (0) Regress X~Z for two-stage MRP====
# parameters for Nj stan model
Ma <- nlevels(samp$age3) 
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
y <- as.numeric(samp$svefreq2)-1
xj <- samp$J_cell
grptbl <- samp %>%
  group_by(age3,sex2,race3,educat4,povgap4,svefreq2,abcat,bccat,cdcat,dfcat) %>%
  summarise(zpopcts = mean(x1_popcts))
J_stan <- nrow(grptbl)
agroup <- as.numeric(grptbl$age3)
bgroup <- as.numeric(grptbl$sex2) -1
cgroup <- as.numeric(grptbl$race3)
dgroup <- as.numeric(grptbl$educat4)
egroup <- as.numeric(grptbl$povgap4)
xgroup <- as.numeric(grptbl$svefreq2)-1
abgroup <- as.numeric(grptbl$abcat)
bcgroup <- as.numeric(grptbl$bccat)
cdgroup <- as.numeric(grptbl$cdcat)
dxgroup <- as.numeric(grptbl$dfcat)
zpopcts <- grptbl$zpopcts


stanfitNq <- stan(file = "mrp2_Nq_BASELINE.stan",
                  data = c("Ma", "Mb", "Mc", "Md", "Me", "L","Mab","Mbc","Mcd", "Mdx", 
                           "J_stan", "n", "y", "xj",
                           "agroup", "bgroup", "cgroup", "dgroup", "egroup", "xgroup",
                           "abgroup", "bcgroup","cdgroup","dxgroup"),
                  iter=staniters,
                  #pars = c("propmk"),
                  warmup=staniters-sim/4,control=list(adapt_delta=0.99, max_treedepth=13),chains=4)
nq_pars <- 
  rstan::extract(stanfitNq,permuted = TRUE, inc_warmup = FALSE,include = TRUE) 

# extract cell frequency estimates from MRP
nq_draws <- nq_pars$propmk 
saveRDS(nq_draws, "nq_stan.rds")
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
  group_by(age3,sex2,race3,educat4,povgap4,svefreq2,abcat,bccat,cdcat,dfcat) %>%
  summarise(n())
grptbl$sampJ <- 1:nrow(grptbl)
samp <- left_join(samp,grptbl)
J_stan <- nrow(grptbl) # if fitting with sampgroup
agroup <- as.numeric(grptbl$age3)
bgroup <- as.numeric(grptbl$sex2) -1
cgroup <- as.numeric(grptbl$race3) 
dgroup <- as.numeric(grptbl$educat4)
egroup <- as.numeric(grptbl$povgap4)
xgroup <- as.numeric(grptbl$svefreq2)-1
abgroup <- as.numeric(grptbl$abcat)
bcgroup <- as.numeric(grptbl$bccat)
cdgroup <- as.numeric(grptbl$cdcat)
dxgroup <- as.numeric(grptbl$dfcat)
cell_label <- samp$sampJ
grp1id <- grp1id[sort(unique(samp$J_cell))]
grp2id <- grp2id[sort(unique(samp$J_cell))]
grp3id <- grp3id[sort(unique(samp$J_cell))]
grp4id <- grp4id[sort(unique(samp$J_cell))]

stanfit <- stan(file = "outcome_model_BASELINE2.stan",
                data = c("Ma", "Mb", "Mc", "Md", "Me", "L", "Mab", "Mbc","Mcd","Mdx",
                         "J_stan", "n", "y", 
                         "agroup", "bgroup","cgroup", "dgroup", "egroup", "xgroup", 
                         "abgroup", "bcgroup","cdgroup","dxgroup","cell_label"),
                iter=staniters,
                pars = c("cellmean","y_sim", "p", "beta","alphaa","alphac","alphad","alphae",
                        "alphadx",
                         "sigma_a","sigma_c","sigma_d","sigma_e","sigma_dx"),#"log_lik"),
                warmup=staniters-sim/4,control=list(adapt_delta=0.99),chains=4)

# extract estimates of cell means from the stan model
stanpars <-
  rstan::extract(object=stanfit, permuted = TRUE)#, inc_warmup = FALSE,include = TRUE)
cellmeans_stan <- stanpars$cellmean
saveRDS(stanpars, "emrp_stan.rds")
