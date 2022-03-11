##################
# MAIN CODE BODY #
##################
## Interaction terms in the X-Z relationship
type = "int"
options(mc.cores =2)
# Libraries -------
library(dplyr)
library(stringr)
require(rstan)
require(polyapost)
require(LaplacesDemon)
# Simulation parameters ------
set.seed(50)
iter = 200      # how many parent samples to draw
sim = 1000       # how many iterations to keep from multilevel models
staniters = 2000 # how many iterations to run rstan for
L <- 1000      # how many synthetic populations to create
F_draw <- 20     # how many subpopulations for each L
grps=5

# Generate the population----
source("helpers.R")
source(paste0("genPop_",type,".R"))
N <- nrow(df)
true_mean <- mean(df$Y)
true_mean1 <- mean(df$Y[df$J_cell%in% group1])
true_mean2 <- mean(df$Y[df$J_cell%in% group2])
true_mean3 <- mean(df$Y[df$J_cell%in% group3])
true_mean4 <- mean(df$Y[df$J_cell%in% group4])
true_mean5 <- mean(df$Y[df$J_cell%in% group5])
true_wts <-
  df %>% 
  group_by(J_cell) %>% 
  summarise(popcts = mean(Nj),
            zpopcts=mean(zcts),
            wts = mean(Nj)/N, 
            cellmeans = mean(Y))#, 
# pxz=mean(X.1*(1-pxz)+(1-X.1)*(pxz)))
# true_nq <- true_wts$pxz
popcts <- true_wts$popcts
zpopcts_lengthJ <- true_wts$zpopcts
true_wtsz <-
  df %>% 
  group_by(Z_categorical) %>% 
  summarise(
    zpopcts=mean(zcts))
zpopcts_lengthZ <- true_wtsz$zpopcts
# Source code for methods -----
source("rstan_mrp.R")
source("rstan_emrp.R")
source("rstan_mrp2.R")

## update functions will declare empty containers on their own
source("update_mrp.R")
source("update_mmrp.R")
source("update_mrp2.R")
source("update_wfpbby.R")
source("update_wfpbbmrp.R")

cellpicked <- rep(0,J) # for calculating Nj stats using one-pass algorithms

  set.seed(50)
# Simulation body -----
for(i in 1:iter){
  
  #-- draw sample
  df$I <- rbinom(n=N, size=1, prob = df$p_include)
  sampled <- df[df$I==1,]
  n <- nrow(sampled)
  # relabeling Z-categories for mrp 
    sampled_z <- sort(unique(sampled$Z_categorical))
    sampled_Zlength <- length(sampled_z)
    # create weights for wfpbb
    temptbl <- sampled %>% group_by(Z_categorical) %>% summarise(wts = mean(zcts)/n())
    temptbl$sampledzlab <- seq(1, sampled_Zlength)
    sampled <- left_join(sampled,temptbl)
    
  # relabeling j-categories for emrp
    temptbl2 <- sampled %>% group_by(J_cell) %>% summarise(sampled_J=mean(J_cell))
    sampled_J <- temptbl2$J_cell
    sampled_Jlength <- length(sampled_J)
    temptbl2$sampledJlab <- seq(1, sampled_Jlength)
    sampled <- left_join(sampled,temptbl2)

  cellpicked[sampled_J] <- cellpicked[sampled_J] + 1
    
  #-- run rstan multilevel regression for MRP
  cellmeans_mrp <- rstan.mrp(samp=sampled)
  
  #-- run rstan multilevel regression for EMRP 
  cellmeans_emrp <- rstan.emrp()
  
  #-- find Nj estimates for 2 stage MRP
  Nj_mrp2 <- nj.mrp2()

  #-- use WFPBBY to analyze
  update.wfpbby()
  
  #-- use classical MRP to analyze
  update.mrp()
  
  #-- use multinomial-mrp to analyze the sample
  update.mmrp()
  
  #-- use 2stage MRP to analyze
  update.mrp2()
  
  #-- use WFPBB-MRP to analyze
  update.wfpbbmrp()
  
}

# Compile results and write
  resmat.emrp.mrp(label="mrp", subgroups=grps, type = type)
  resmat.emrp.mrp(label="mmrp", subgroups=grps,type = type)
  resmat.emrp.mrp(label="mrp2", subgroups=grps,type = type)
  resmat.emrp.mrp(label="wmrp", subgroups=grps,type = type)
  resmat.wfpbby(label="wfpbby", subgroups=grps,type = type)
  