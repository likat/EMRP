# classic MRP
#-- load required packages
require(dplyr)
require(rstan)
options(mc.cores = 2)

#-- load population generating function:
source("genPop_hixz_hixy.R")

# weights from population
poptblwt <- 
  df %>% 
  group_by(Z_categorical) %>% 
  summarise(zwts= n()/N)
wts <- poptblwt$zwts
#-- parameters for testing the function
set.seed(50)
iter = 200
sim = 1000
staniters = 2000

#-- random integers for naming the stanfit to avoid segfault issues
stanfitrand <- sample(1:iter, iter, replace=F)
modelname <- paste("stanfit",stanfitrand, sep="")

#-- helper function to normalize
normalize <- function(vec){
  return(vec/sum(vec))
}

`%notin%` <- Negate(`%in%`)

# population size, number cells
N <- nrow(df)

# what is mean(Y) for the subgroup?
true_mean <- mean(df$Y)
true_mean1 <- mean(df$Y[df$J_cell%in% group1])
true_mean2 <- mean(df$Y[df$J_cell%in% group2])
true_mean3 <- mean(df$Y[df$J_cell%in% group3])
true_mean4 <- mean(df$Y[df$J_cell%in% group4])
true_mean5 <- mean(df$Y[df$J_cell%in% group5])

# initialize containers for results
y_bar <- y_bar1 <- y_bar2 <- y_bar3 <- y_bar4 <- y_bar5 <-rep(NA, iter)
mrp_posterior <- mrp_posterior1<- mrp_posterior2 <- mrp_posterior3 <- mrp_posterior4 <-mrp_posterior5<-rep(NA, sim)


CI       <- matrix(nrow = iter, ncol = 2)   # credible interval
CI1       <- matrix(nrow = iter, ncol = 2)   # credible interval
CI2       <- matrix(nrow = iter, ncol = 2)   # credible interval
CI3       <- matrix(nrow = iter, ncol = 2)   # credible interval
CI4       <- matrix(nrow = iter, ncol = 2)   # credible interval
CI5       <- matrix(nrow = iter, ncol = 2)   # credible interval

CI_covered    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered1    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered2    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered3    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered4    <- rep(NA,iter)    # 1 if covered by CI, NA else
CI_covered5    <- rep(NA,iter)    # 1 if covered by CI, NA else

bias <- 0
rMSE <- 0
SE <- 0


for(i in 1:iter){
  #-- draw sample
  df$I <- rbinom(n=N, size=1, prob = df$p_include)
  sampled <- df[df$I==1,-which(colnames(df)%in% c("Nj","p_include","pxz"))]
  
  sampled_z <- sort(unique(sampled$Z_categorical))
  sampled_Zlength <- length(sampled_z)
  temptbl2 <- sampled %>% group_by(Z_categorical) %>% summarise(sampled_z=mean(Z_categorical))
  temptbl2$sampledzlab <- seq(1, sampled_Zlength)
  sampled <- left_join(sampled,temptbl2)
  
  #-- Parameters for stan model
  Ma <- nlevels(sampled$Za_categorical)
  Mb <- nlevels(sampled$Zb_categorical)
  Mc <- nlevels(sampled$Zc_categorical)
  J_stan <- sampled_Zlength
  n <- nrow(sampled)
  y <- as.vector(sampled$Y)
  cell_label<-sampled$sampledzlab
  
  grptbl <- sampled %>% 
    group_by(Z_categorical) %>% 
    summarise(zalab = mean(as.numeric(Za_categorical)),
              zblab = mean(as.numeric(Zb_categorical)),
              zclab = mean(as.numeric(Zc_categorical)-1),
    )
  agroup <- grptbl$zalab
  bgroup <- grptbl$zblab
  cgroup <- grptbl$zclab

  
  #-- run the model
  assign(modelname[i],stan(file= "mrp_sim_Zonly.stan",
                           data = c("Ma", "Mb", "Mc","J_stan","n","y","agroup","bgroup","cgroup", "cell_label"),
                           iter=staniters, pars=c("cellmean"),warmup = staniters-sim/2, control=list(adapt_delta=0.99, max_treedepth=13),chains=2))

  cellmeans_stan <- extract(get(modelname[i]),permuted = TRUE, inc_warmup = FALSE,include = TRUE)$cellmean
  wts_mrp <- matrix(0,nrow = Zcat,ncol=2,byrow=F)
  mrp_posterior <- mrp_posterior1<- mrp_posterior2 <- mrp_posterior3 <- mrp_posterior4 <-mrp_posterior5 <-rep(NA, sim)
  
    wts_mrp <- normalize(wts[sampled_z])
    wts_mrp1 <- normalize((wts_mrp*cmat1z[sampled_z]))
    wts_mrp2 <- normalize((wts_mrp*cmat2z[sampled_z]))
    wts_mrp3 <- normalize((wts_mrp*cmat3z[sampled_z]))
    wts_mrp4 <- normalize((wts_mrp*cmat4z[sampled_z]))
    wts_mrp5 <- normalize((wts_mrp*cmat5z[sampled_z]))
  for(s in 1:sim){
    current_cellmean <- cellmeans_stan[s,]
    
    mrp_posterior[s] <- crossprod(wts_mrp, current_cellmean)
    mrp_posterior1[s] <- crossprod(wts_mrp1, current_cellmean)
    mrp_posterior2[s] <- crossprod(wts_mrp2, current_cellmean)
    mrp_posterior3[s] <- crossprod(wts_mrp3, current_cellmean)
    mrp_posterior4[s] <- crossprod(wts_mrp4, current_cellmean)
    mrp_posterior5[s] <- crossprod(wts_mrp5, current_cellmean)

  }
  y_bar[i] <- mean(mrp_posterior,na.rm=T)
  y_bar1[i] <- mean(mrp_posterior1,na.rm=T)
  y_bar2[i] <- mean(mrp_posterior2,na.rm=T)
  y_bar3[i] <- mean(mrp_posterior3,na.rm=T)
  y_bar4[i] <- mean(mrp_posterior4,na.rm=T)
  y_bar5[i] <- mean(mrp_posterior5,na.rm=T)
  CI[i,] <- quantile(mrp_posterior, c(0.025, 0.975),na.rm=T)
  CI1[i,] <- quantile(mrp_posterior1, c(0.025, 0.975),na.rm=T)
  CI2[i,] <- quantile(mrp_posterior2, c(0.025, 0.975),na.rm=T)
  CI3[i,] <- quantile(mrp_posterior3, c(0.025, 0.975),na.rm=T)
  CI4[i,] <- quantile(mrp_posterior4, c(0.025, 0.975),na.rm=T)
  CI5[i,] <- quantile(mrp_posterior5, c(0.025, 0.975),na.rm=T)

  # sampsub <- unique(sampled$subgroup)
  if(is.na(between(true_mean, CI[i,1], CI[i,2]))){break}
  if(between(true_mean, CI[i,1], CI[i,2])){CI_covered[i] <- 1}else if(sum(is.na(CI[i,])!=0)){CI_covered[i] <- NA}else{CI_covered[i] <- 0}
  if(between(true_mean1, CI1[i,1], CI1[i,2])){CI_covered1[i] <- 1} else{CI_covered1[i] <- 0}
  if(between(true_mean2, CI2[i,1], CI2[i,2])){CI_covered2[i] <- 1} else{CI_covered2[i] <- 0}
  if(between(true_mean3, CI3[i,1], CI3[i,2])){CI_covered3[i] <- 1} else{CI_covered3[i] <- 0}
  if(between(true_mean4, CI4[i,1], CI4[i,2])){CI_covered4[i] <- 1} else{CI_covered4[i] <- 0}
  if(between(true_mean5, CI5[i,1], CI5[i,2])){CI_covered5[i] <- 1} else{CI_covered5[i] <- 0}
  #-- remove stanfit object to prevent segfault and conserve memory
  rm(list=c(modelname[i],"cellmeans_stan",
            "mrp_posterior","mrp_posterior1", "mrp_posterior2","mrp_posterior3","mrp_posterior4","mrp_posterior5"))
  
}    

# CI is a list
#ybarmat, CIcovered is a matrix
ybarmat <- cbind(y_bar, y_bar1, y_bar2, y_bar3, y_bar4, y_bar5)
truemeanmat <- c(true_mean,true_mean1,true_mean2,true_mean3,true_mean4,true_mean5)
labs <- c("Overall", "group 1", "group2", "group3", "group4", "group5")
CImat <- list(CI,CI1,CI2,CI3,CI4,CI5)
CIcovmat <- cbind(CI_covered,CI_covered1,CI_covered2,CI_covered3,CI_covered4,CI_covered5)
resmatfun <- function(ybarmat, truemean, CI, CIcovered,labs){
  errormat <- t(t(ybarmat)-truemean)
  bias <- colMeans(errormat)
  rMSE <- sqrt(colMeans(errormat^2))
  SE <- t(t(ybarmat)-colMeans(ybarmat))
  SE <- sqrt(colMeans(SE^2))
  CI_length <- vector(length = ncol(ybarmat))
  for(i in 1:ncol(ybarmat)){
    CItemp <- CI[[i]]
    CI_length[i] <- mean(CItemp[,2]-CItemp[,1])
  }
  coverage_rate <- colMeans(CIcovered)
  res <-rbind(rMSE, bias, SE, CI_length, coverage_rate)
  colnames(res) <- labs
  return(res)
}
resultmat <- resmatfun(ybarmat = ybarmat, truemean = truemeanmat, CI = CImat, CIcovered = CIcovmat, labs = labs)

# PRINT RESULTS
print("Result Summary")
resultmat  

warnings()

write.csv(resultmat, "res_mrp.csv")

