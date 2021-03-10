#-- RESULTS for case 2 ------
library(ggplot2)
setwd("~/Desktop/fwbb_mrp/EMRP PAPER BACKUP/run2v2/res/case2")

#====== Result tables =========
wmrp <- read.csv("res_wmrp.csv",row.names=1)
mrp2 <- read.csv("res_mrp2.csv",row.names=1)
  names(mrp2) <- names(wmrp)
mmrp <- read.csv("res_mmrp.csv",row.names=1)
  names(mmrp) <- names(wmrp)
wfpbby <- read.csv("res_wfpbby.csv",row.names=1)

rmse <- rbind(wfpbby[1,], wmrp[1,], mmrp[1,], mrp2[1,])
rmse <- rmse[,c(1,6,5,4,3,2)]
rownames(rmse) <- c("WFPBB", "WFPBB-MRP", "Multinomial-MRP", "Two-Stage MRP")
bias <- rbind(wfpbby[2,], wmrp[2,], mmrp[2,], mrp2[2,])
bias <- bias[,c(1,6,5,4,3,2)]
rownames(bias) <- c("WFPBB", "WFPBB-MRP", "Multinomial-MRP", "Two-Stage MRP")
SE <- rbind(wfpbby[3,], wmrp[3,], mmrp[3,], mrp2[3,])
SE <- SE[,c(1,6,5,4,3,2)]
rownames(SE) <- c("WFPBB", "WFPBB-MRP", "Multinomial-MRP", "Two-Stage MRP")
CI_length <- rbind(wfpbby[4,], wmrp[4,], mmrp[4,], mrp2[4,])
CI_length <- CI_length[,c(1,6,5,4,3,2)]
rownames(CI_length)<- c("WFPBB", "WFPBB-MRP", "Multinomial-MRP", "Two-Stage MRP")
coverage_rate <- rbind(wfpbby[5,], wmrp[5,], mmrp[5,], mrp2[5,])
coverage_rate <- coverage_rate[,c(1,6,5,4,3,2)]
rownames(coverage_rate) <- c("WFPBB", "WFPBB-MRP", "Multinomial-MRP", "Two-Stage MRP")

# Source genPop to get ordered cell count indices
# "ordered_cellcts_hihi" := orders cells from small theoretical sample size to large
source("../../case2-xz-main/genPop_hixz_hixy.R")

# labeling heatmaps
method <- c( "WFPBB", "WFPBB-MRP", "Multinomial-MRP","Two-Stage-MRP")
inference <- c("Overall","Group 5", "Group 4", "Group 3", "Group 2",  "Group 1")
data <- expand.grid(Method=method, Inference=inference)

#====== Overall summaries =======
## Use Bias, rMSE, coverage rate for loxz/hixy (most difference observed)
#- rMSE
# start from overall WFPBB -> overall W-MRP -> overall M-MRP... -> group 5 WFPBB -> group 5 W-MRP...
data$Z <- as.numeric(as.matrix(rmse))
ggplot(data, aes(Method, Inference, fill= Z)) + 
  geom_tile() + labs(title = "rMSE for main XZ case", fill = "rMSE")+scale_fill_continuous(high = "#28629c", low = "#56B1F7")

#- bias
## hi xz / hi xy
data$Z <- as.numeric(as.matrix(bias))
ggplot(data, aes(Method, Inference, fill= abs(Z))) + 
  geom_tile() + 
  labs(title = "Absolute bias for main XZ case", fill = "bias")+scale_fill_continuous(high = "#28629c", low = "#56B1F7")

#- coverage
## hi xz / hi xy
data$Z <- as.numeric(as.matrix(coverage_rate))
ggplot(data, aes(Method, Inference, fill= Z)) + 
  geom_tile()+ labs(title = "95% Interval coverage for main XZ case", fill = "coverage")+scale_fill_continuous(high = "#28629c", low = "#56B1F7")

#====== read in Nj ======
wmrp_njbias <- read.csv("Njbias_wmrp.csv", header = T, row.names = 1) %>% unlist %>% as.vector
mmrp_njbias <- read.csv("Njbias_mmrp.csv", header = T, row.names = 1)%>% unlist %>% as.vector
mrp2_njbias <- read.csv("Njbias_mrp2.csv", header = T, row.names = 1)%>% unlist %>% as.vector
wmrp_njse <- read.csv("Njse_wmrp.csv", header = T, row.names = 1)%>% unlist %>% as.vector
mmrp_njse <- read.csv("Njse_mmrp.csv", header = T, row.names = 1)%>% unlist %>% as.vector
mrp2_njse <- read.csv("Njse_mrp2.csv", header = T, row.names = 1)%>% unlist %>% as.vector
wmrp_njrmse <- read.csv("Njrmse_wmrp.csv", header = T, row.names = 1)%>% unlist %>% as.vector
mmrp_njrmse <- read.csv("Njrmse_mmrp.csv", header = T, row.names = 1)%>% unlist %>% as.vector
mrp2_njrmse <- read.csv("Njrmse_mrp2.csv", header = T, row.names = 1)%>% unlist %>% as.vector
wmrp_njcov <- read.csv("Njcov_wmrp.csv", header = T, row.names = 1)%>% unlist %>% as.vector
mmrp_njcov <- read.csv("Njcov_mmrp.csv", header = T, row.names = 1)%>% unlist %>% as.vector
mrp2_njcov <- read.csv("Njcov_mrp2.csv", header = T, row.names = 1)%>% unlist %>% as.vector

##======= Nj Bias =====
#--- faceted plots
## note: J=99 for this population, so append NA to end of each result vector
x <- y <- factor(1:10)
method <- c("WFPBB-MRP", "Multinomial-MRP", "Two-Stage-MRP")
data <- expand.grid(X=x,Y=y, Method = method)
wmdat <- wmrp_njbias %>% abs()
m2dat <- mrp2_njbias%>% abs()
mmdat <- mmrp_njbias%>% abs()
wmdat <- c(wmdat[ordered_cellcts_hihi])#,NA)
m2dat <- c(m2dat[ordered_cellcts_hihi])#,NA)
mmdat <- c(mmdat[ordered_cellcts_hihi])#,NA)
data$Z <- c(wmdat, mmdat, m2dat)
ggplot(data, aes(x=X,y=Y, fill= Z)) + 
  geom_tile()+
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "Absolute bias for cell counts in main XZ case",fill="bias")+scale_fill_continuous(high = "#28629c", low = "#56B1F7")+facet_grid(. ~ Method)

##===== Nj rMSE ======

#-- faceted plots  
x <- y <- factor(1:10)
method <- c("WFPBB-MRP", "Multinomial-MRP", "Two-Stage-MRP")
data <- expand.grid(X=x,Y=y, Method = method)
wmdat <- wmrp_njrmse
mmdat <- mmrp_njrmse
m2dat <- mrp2_njrmse
wmdat <- c(wmdat[ordered_cellcts_hihi])#,NA)
m2dat <- c(m2dat[ordered_cellcts_hihi])#,NA)
mmdat <- c(mmdat[ordered_cellcts_hihi])#,NA)
data$Z <- c(wmdat, mmdat, m2dat)
ggplot(data, aes(x=X,y=Y, fill= Z)) + 
  geom_tile()+
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "rMSE of estimated population cell count, main XZ case", fill = "rMSE")+
  scale_fill_continuous(high = "#28629c", low = "#56B1F7")+
  facet_grid(. ~ Method)

##===== Nj SE ====
#-- faceted plots  
x <- y <- factor(1:10)
method <- c("WFPBB-MRP", "Multinomial-MRP", "Two-Stage-MRP")
data <- expand.grid(X=x,Y=y, Method = method)
wmdat <- wmrp_njse
mmdat <- mmrp_njse
m2dat <- mrp2_njse
wmdat <- c(wmdat[ordered_cellcts_hihi])#,NA)
m2dat <- c(m2dat[ordered_cellcts_hihi])#,NA)
mmdat <- c(mmdat[ordered_cellcts_hihi])#,NA)
data$Z <- c(wmdat, mmdat, m2dat)
ggplot(data, aes(x=X,y=Y, fill= Z)) + 
  geom_tile()+
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "Standard error of estimated population cell count, main XZ case", fill = "coverage")+
  scale_fill_continuous(high = "#28629c", low = "#56B1F7")+
  facet_grid(. ~ Method)

#==== Nj coverage rates======

#-- faceted plots  
x <- y <- factor(1:10)
method <- c("WFPBB-MRP", "Multinomial-MRP", "Two-Stage-MRP")
data <- expand.grid(X=x,Y=y, Method = method)
wmdat <- wmrp_njcov
mmdat <- mmrp_njcov
m2dat <- mrp2_njcov
wmdat <- c(wmdat[ordered_cellcts_hihi])#,NA)
m2dat <- c(m2dat[ordered_cellcts_hihi])#,NA)
mmdat <- c(mmdat[ordered_cellcts_hihi])#,NA)
data$Z <- c(wmdat, mmdat, m2dat)
ggplot(data, aes(x=X,y=Y, fill= Z)) + 
  geom_tile()+
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "95% interval coverage rate of estimated population cell count, main XZ case", fill = "coverage")+
  scale_fill_continuous(high = "#28629c", low = "#56B1F7")+
  facet_grid(. ~ Method)

#===== Individual plots ======

#-- Bias
#- wmrp
x <- y <- factor(1:10)
data <- expand.grid(X=x,Y=y)
dat <- wmrp_njbias
dat <- abs(dat[ordered_cellcts_hihi])
data$Z <- dat
ggplot(data, aes(x=X,y=Y, fill= Z)) + 
  geom_tile()+
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "Cell count bias for High XZ case, W-MRP")+scale_fill_continuous(high = "#28629c", low = "#56B1F7")

#- m-mrp
x <- y <- factor(1:10)
data <- expand.grid(X=x,Y=y)
dat <- mmrp_njbias
dat <- abs(dat[ordered_cellcts_hihi])
data$Z <- dat
ggplot(data, aes(x=X,y=Y, fill= Z)) +
  geom_tile()+
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "Cell count bias for High XZ case, M-MRP")+scale_fill_continuous(high = "#28629c", low = "#56B1F7")

#- 2-mrp
x <- y <- factor(1:10)
data <- expand.grid(X=x,Y=y)
dat <- mrp2_njbias
dat <- abs(dat[ordered_cellcts_hihi])
data$Z <- dat
ggplot(data, aes(x=X,y=Y, fill= Z)) + 
  geom_tile()+
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "Cell count bias for High XZ case, 2-MRP")+scale_fill_continuous(high = "#28629c", low = "#56B1F7")




#---- rMSE
#- wmrp
x <- y <- factor(1:10)
data <- expand.grid(X=x,Y=y)
dat <- wmrp_njrmse
dat <- abs(dat[ordered_cellcts_hihi])
data$Z <- dat
ggplot(data, aes(x=X,y=Y, fill= Z)) + 
  geom_tile()+
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "rMSE for cell counts in main XZ case, WFPBB-MRP")+scale_fill_continuous(high = "#28629c", low = "#56B1F7")

#- m-mrp
x <- y <- factor(1:10)
data <- expand.grid(X=x,Y=y)
dat <- mmrp_njrmse
dat <- abs(dat[ordered_cellcts_hihi])
data$Z <- dat
ggplot(data, aes(x=X,y=Y, fill= Z)) + 
  geom_tile()+
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "rMSE for cell counts in main XZ case, Multinomial-MRP")+scale_fill_continuous(high = "#28629c", low = "#56B1F7")

#- 2-mrp
x <- y <- factor(1:10)
data <- expand.grid(X=x,Y=y)
dat <- mrp2_njrmse
dat <- abs(dat[ordered_cellcts_hihi])
data$Z <- dat
ggplot(data, aes(x=X,y=Y, fill= Z)) + 
  geom_tile()+
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "rMSE for cell counts in main XZ case, two-stage MRP")+scale_fill_continuous(high = "#28629c", low = "#56B1F7")

#----- SE
#- wmrp
x <- y <- factor(1:10)
data <- expand.grid(X=x,Y=y)
dat <- wmrp_njse
dat <- dat[ordered_cellcts_hihi]
data$Z <- dat
ggplot(data, aes(x=X,y=Y, fill= Z)) + 
  geom_tile()+
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "Cell count SE for main XZ case, WFPBB-MRP")+scale_fill_continuous(high = "#28629c", low = "#56B1F7")

#- m-mrp
x <- y <- factor(1:10)
data <- expand.grid(X=x,Y=y)
dat <- mmrp_njse
dat <- dat[ordered_cellcts_hihi]
data$Z <- dat
ggplot(data, aes(x=X,y=Y, fill= Z)) +
  geom_tile()+
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "Cell count SE for main XZ case, Multinomial-MRP")+scale_fill_continuous(high = "#28629c", low = "#56B1F7")

#- 2-mrp
x <- y <- factor(1:10)
data <- expand.grid(X=x,Y=y)
dat <- mrp2_njse
dat <- dat[ordered_cellcts_hihi]
data$Z <- dat
ggplot(data, aes(x=X,y=Y, fill= Z)) + 
  geom_tile()+
  theme(axis.ticks.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  labs(title = "Cell count SE for main XZ case, Two-stage MRP")+scale_fill_continuous(high = "#28629c", low = "#56B1F7")
