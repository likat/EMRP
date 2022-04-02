#================================
# Embedded MRP with applications
#================================
#-- last updated: aug 2021 --
#-- uses public use baseline data

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
library(loo)

#-- helper functions
`%notin%` <- Negate(`%in%`)
normalize <- function(vec){
  # does exactly what it's called
  return(vec/sum(vec))
}
#-- source cleaned data and group indices
# setwd("~/Desktop/fwbb_mrp/EMRP PAPER BACKUP/realdata-output/FINALCODE_applications/")
source("../BASELINE.R")

# initialize containers for results
sim <- 5000
staniters = 10000
# sim <- 1000
# staniters = 2000
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
# sampgroup <- as.numeric(grptbl$samptype) -1 # for testing samptype sig
cell_label <- samp$sampJ
grp1id <- grp1id[sort(unique(samp$J_cell))]
grp2id <- grp2id[sort(unique(samp$J_cell))]
grp3id <- grp3id[sort(unique(samp$J_cell))]
grp4id <- grp4id[sort(unique(samp$J_cell))]
grp5id <- grp5id[sort(unique(samp$J_cell))]
grp6id <- grp6id[sort(unique(samp$J_cell))]
grp7id <- grp7id[sort(unique(samp$J_cell))]
# xsampgroup <- as.numeric(grptbl$samptype)-1 # for testing samptype sig

# age:sex, education:visit
stanfit <- stan(file = "outcome_model_BASELINE2.stan",
                data = c("Ma", "Mb", "Mc", "Md", "Me", "L", "Mab", "Mbc","Mcd","Mdx",
                         "J_stan", "n", "y", 
                         "agroup", "bgroup","cgroup", "dgroup", "egroup", "xgroup", 
                         "abgroup", "bcgroup","cdgroup","dxgroup","cell_label"),
                iter=staniters,
                # pars = c("log_lik"),
                warmup=staniters-sim/4,control=list(adapt_delta=0.99),chains=4)
loglikemrp <- extract_log_lik(stanfit,merge_chains=F)
r_eff_emrp <- relative_eff(exp(loglikemrp))
loo_emrp <- loo(loglikemrp, r_eff = r_eff_emrp, cores = 4)
print(loo_emrp)
stanpars <-
  rstan::extract(object=stanfit, permuted = TRUE)#, inc_warmup = FALSE,include = TRUE)
write.csv(stanpars$y_sim, "ysimraw.csv")
# write.csv(stanpars$theta_rep, "predthetaraw.csv")
write.csv(stanpars$p, "probraw.csv")
write.csv(stanpars$beta[,1], "betaintraw.csv")
write.csv(stanpars$beta[,2], "betasexraw.csv")
write.csv(stanpars$beta[,3], "betavisitraw.csv")
write.csv(stanpars$alphaa, "araw.csv")
write.csv(stanpars$alphac, "craw.csv")
write.csv(stanpars$alphad, "draw.csv")
write.csv(stanpars$alphae, "aeraw.csv")
# write.csv(stanpars$alphaab, "abraw.csv")
# write.csv(stanpars$alphabc, "bcraw.csv")
# write.csv(stanpars$alphacd, "cdraw.csv")
write.csv(stanpars$alphadx, "dxraw.csv")
write.csv(stanpars$sigma_a, "sigaraw.csv")
write.csv(stanpars$sigma_c, "sigcraw.csv")
write.csv(stanpars$sigma_d, "sigdraw.csv")
write.csv(stanpars$sigma_e, "sigeraw.csv")
# write.csv(stanpars$sigma_ab, "sigabraw.csv")
# write.csv(stanpars$sigma_bc, "sigbcraw.csv")
# write.csv(stanpars$sigma_cd, "sigcdraw.csv")
write.csv(stanpars$sigma_dx, "sigdxraw.csv")

#==== MRP

M <- length(unique(samp$sampled_x1_label))
#== redefine groupids in terms of Z
x1tbl <- samp %>% group_by(J_cell) %>% summarise(z = mean(sampled_x1_label))

grp1z <-grp2z <- grp3z <- grp4z <- grp5z <- grp6z<- grp7z <- rep(0,M)
grp1z[unique(x1tbl$z[grp1id==1])] <- 1
grp2z[unique(x1tbl$z[grp2id==1])] <- 1
grp3z[unique(x1tbl$z[grp3id==1])] <- 1
grp4z[unique(x1tbl$z[grp4id==1])] <- 1
grp5z[unique(x1tbl$z[grp5id==1])] <- 1
grp6z[unique(x1tbl$z[grp6id==1])] <- 1
grp7z[unique(x1tbl$z[grp7id==1])] <- 1



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
grp1id <- grp1id[sort(unique(samp$J_cell))]
grp2id <- grp2id[sort(unique(samp$J_cell))]
grp3id <- grp3id[sort(unique(samp$J_cell))]
grp4id <- grp4id[sort(unique(samp$J_cell))]
grp5id <- grp5id[sort(unique(samp$J_cell))]
grp6id <- grp6id[sort(unique(samp$J_cell))]
grp7id <- grp7id[sort(unique(samp$J_cell))]
zpopcts <- grptbl$zcts

stanfit <- stan(file = "outcome_model_BASELINE.stan",
                data = c("Ma", "Mb", "Mc", "Md", "Me", "L", "Mab","Mbc","Mcd",
                         "J_stan", "n", "y", 
                         "agroup", "bgroup","cgroup", "dgroup", "egroup",
                         "abgroup","bcgroup","cdgroup","cell_label"), 
                pars = c("log_lik"),
                iter=staniters,
                warmup=staniters-sim/4,control=list(adapt_delta=0.99),chains=4)

loglikmrp <- extract_log_lik(stanfit,merge_chains=F)
r_eff_mrp <- relative_eff(exp(loglikmrp))
loo_mrp <- loo(loglikmrp, r_eff = r_eff_mrp, cores = 4)
print(loo_mrp)

comp <- loo_compare(loo_emrp, loo_mrp)
print(comp)
stanpars <-
  rstan::extract(object=stanfit, permuted = TRUE)
write.csv(stanpars, "mrppars.csv")



