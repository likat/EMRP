####### Code for population generation MAIN-BIN ########
### Y is continuous
### Z split into three factors, X binary
### Data are *MAR*
library(MASS)
library(dplyr)
set.seed(42)

## expit helper function
expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

N = 1e4 # population size

## [Z] Poststratifiers known in the population
# 3 separate variables: A, B, C
acats <- 5; bcats <- 5; ccats <- 2
Zcat = acats*bcats*ccats # should be 50

## [X] Poststratifier known in sample only    
Xcat = 2

#-- [Y] DGP parameters
sig_y=1 # sd for data
beta0 = 0
alpha_a <- rnorm(5)
alpha_b <- rnorm(5)
# alpha_a <- c(1.35, -0.5, 0.35, 0.6, 0.4)
# alpha_b <- c(-0.1, 1.5, -0.1, 2, -0.05)

alpha_c <-c(0,0.24)
alpha_x <- c(0,-1.3)

#-- P(I=1|Z) Inclusion probabilities
pinc_tbl = data.frame(Z_categorical=c(1:50), p_include = c(
  sample(seq(0.01,0.1,length.out=10),size = 5,replace = T), #1:5
  sample(seq(0.11,0.4,length.out=10),size = 15,replace = T),#6:20
  sample(seq(0.21,0.6,length.out=10),size = 20,replace = T),#21:40
  sample(seq(0.51,0.8,length.out=10),size = 5,replace = T),#41:45
  sample(seq(0.8,0.99,length.out=10),size = 5,replace = T)))#46:50

#-- Parameters to generate Z
p_a = runif(acats, 0.25,1)
p_b = runif(bcats, 0.25,1)
p_c = runif(ccats, 0.25,1)

Za <- t(rmultinom(n=N, size=1, p=p_a)) %>% as.data.frame()
Zb <- t(rmultinom(n=N, size=1, p=p_b)) %>% as.data.frame()
Zc <- t(rmultinom(n=N, size=1, p=p_c)) %>% as.data.frame()

for(a in 1:acats){
  colnames(Za)[a] = paste("A", a, sep="")
}
for(b in 1:bcats){
  colnames(Zb)[b] = paste("B", b, sep="")
}
for(c in 1:ccats){
  colnames(Zc)[c] = paste("C", c, sep="")
}

Z = cbind(Za, Zb, Zc)

# Create categorical levels for easier reference
Za_categorical <- Za
Zb_categorical <- Zb
Zc_categorical <- Zc
for(Z1 in 1:acats){
  Za_categorical[,Z1] <- Z1*Za_categorical[,Z1]
}
for(Z2 in 1:bcats){
  Zb_categorical[,Z2] <- Z2*Zb_categorical[,Z2]
}
for(Z3 in 1:ccats){
  Zc_categorical[,Z3] <- Z3*Zc_categorical[,Z3]
}
Za_categorical <- rowSums(Za_categorical)
Zb_categorical <- rowSums(Zb_categorical)
Zc_categorical <- rowSums(Zc_categorical)

df <- 
  as.data.frame(cbind(Z,Za_categorical, Zb_categorical, Zc_categorical)) %>% 
  mutate(Za_categorical = factor(Za_categorical),
         Zb_categorical = factor(Zb_categorical),
         Zc_categorical = factor(Zc_categorical)) %>% 
  arrange(Za_categorical, Zb_categorical, Zc_categorical) %>% as.data.frame()

# get counts for each category of Z
tempdf <- 
  df %>% group_by(Za_categorical, Zb_categorical, Zc_categorical) %>% 
  summarise(zcts=n())  %>%  as.data.frame()
tempdf$Z_categorical <- 1:nrow(tempdf)

#-- [X|Z] generation
bxz0 <- -0.5
# bxza <- rnorm(acats+1)-0.5
bxza <- c(1.7, 0.25, 0.2, -0.75, -1.7)

# bxza <- bxza[-1] # get rid of strange value (0.0025)
# bxzb <- rnorm(bcats) + 0.5
bxzb <- c(2.3, 1.5, 0.15, 0.2, 0.9)
# bxzc <- c(0,rnorm(1)+0.1)
bxzc <- c(0,-1)
# bxzac = c(0, rnorm(acats-1) - 0.2) 
# bxzbc = c(0, rnorm(bcats-1) + 0.2) 
# bxzac = c(0, -0.6, 0.5, 0.35, -0.4) 
# bxzbc = c(0, 1.7, 0.1, 2, -0.75)

bxzac <- rep(0, acats) #MAIN
bxzbc <- rep(0,bcats)  #MAIN

df <-
  df %>% 
  left_join(tempdf)
df$pxz <-  expit(bxz0+bxza[df$Za_categorical] +
                   bxzb[df$Zb_categorical] +
                   bxzc[df$Zc_categorical])
df$X.0 = rbinom(N, 1, 1-df$pxz)
df$X.1 = as.numeric(df$X.0 != 1)

# join tables to get zcts, p_include
df <- left_join(df,pinc_tbl)

#-- [Y] generation
df$mu <- beta0 +
  alpha_a[df$Za_categorical] +
  alpha_b[df$Zb_categorical] +
  alpha_c[df$Zc_categorical] +
  alpha_x[df$X.1+1]
# df$Y <-  df$mu + rnorm(N, 0, sig_y) # CTS case
df$Y <-  rbinom(N,1,expit(df$mu))  # BIN case

#---- Specification of cell levels between XxZ
popcts_tl = df %>% 
  group_by(Za_categorical, Zb_categorical, Zc_categorical, X.1) %>%  
  summarise(Nj = n(), zcategorical=mean(as.numeric(Z_categorical))) %>% 
  as.data.frame()
J <- nrow(popcts_tl)
popcts_tl$J_cell = 1:J
df = left_join(df, popcts_tl) %>% arrange(J_cell)

## order by the theoretical sampled counts E[nj]; subgroup 1 is 
# composed of 20 of cells with least sampled, subgroup 5
# composed of 20 of cells with most sampled
poptbl <- df %>% 
  group_by(J_cell) %>% 
  summarise(pinc = mean(p_include), cts = n(), cellmeans = mean(Y), cellvars = var(Y))

# order subgroups by theoretical sampled cell count
temp <- poptbl$pinc 
ordered_groupcts <- order(temp)
ordered_cellcts_hihi <- ordered_groupcts
poptbl <- df %>%
  group_by(J_cell) %>%
  summarise(
    pinc = mean(p_include),
    pxz = mean(pxz),
    cts = n(),
    sampct = pinc*cts,
    wts = cts/sampct,
    Za = mean(as.numeric(Za_categorical)),
    Zb =  mean(as.numeric(Zb_categorical)),
    Zc =  mean(as.numeric(Zc_categorical)),
    zcategorical=mean(as.numeric(Z_categorical)),
    cellmeans = mean(Y),
    cellvars=var(Y),
    avgX1=mean(X.1))

grpincr <- ceiling(J/5)
# Zcdraws <- ceiling((1+5+10+15+19)/20 * grpincr)
# Zccells <- sample(poptbl$J_cell[poptbl$Zc==2], Zcdraws)
# nonZccells <- sample(poptbl$J_cell[poptbl$Zc==1], J-Zcdraws)
# group1 <- c(nonZccells[1:19],Zccells[1])
# group2 <- c(nonZccells[20:34],Zccells[2:6])
# group3 <- c(nonZccells[35:44],Zccells[7:16])
# group4 <- c(nonZccells[45:49],Zccells[17:31])
# group5 <- c(nonZccells[50],Zccells[32:50])

# ordered by pinclude, X = 0
pool10 <- ordered_groupcts[1:(2*grpincr)][ordered_groupcts[1:(2*grpincr)] %% 2 ==0 ]
pool20 <- ordered_groupcts[(grpincr+1):(3*grpincr)][ordered_groupcts[(grpincr+1):(3*grpincr)] %% 2 ==0 ]
pool30 <- ordered_groupcts[(2*grpincr+1):(4*grpincr)][ordered_groupcts[(2*grpincr+1):(4*grpincr)] %% 2 ==0 ]
pool40 <- ordered_groupcts[(3*grpincr+1):(5*grpincr)][ordered_groupcts[(3*grpincr+1):(5*grpincr)] %% 2 ==0]
# pool50 <- ordered_groupcts[(4*grpincr+1):J][ordered_groupcts[(4*grpincr+1):J] %% 2 ==0]

# ordered by pinclude, X = 1
pool11 <- ordered_groupcts[1:(2*grpincr)][ordered_groupcts[1:(2*grpincr)] %% 2 ==1 ]
pool21 <- ordered_groupcts[(grpincr+1):(3*grpincr)][ordered_groupcts[(grpincr+1):(3*grpincr)] %% 2 ==1 ]
pool31 <- ordered_groupcts[(2*grpincr+1):(4*grpincr)][ordered_groupcts[(2*grpincr+1):(4*grpincr)] %% 2 ==1 ]
pool41 <- ordered_groupcts[(3*grpincr+1):(5*grpincr)][ordered_groupcts[(3*grpincr+1):(5*grpincr)] %% 2 ==1]
# pool51 <- ordered_groupcts[(4*grpincr+1):J][ordered_groupcts[(4*grpincr+1):J] %% 2 ==1]

group1 <- c(sample(pool10,5), sample(pool11,15))
group2 <- c(sample(pool20,5), sample(pool21,15))
group3 <- c(sample(pool30,15), sample(pool31,5))
group4 <- c(sample(pool40,5), sample(pool41,15))
group5 <- c(sample(pool40,15), sample(pool41,5))

cmat1 <-cmat2 <- cmat3 <- cmat4 <- cmat5 <- rep(0,J)
cmat1[group1] <- 1
cmat2[group2] <- 1
cmat3[group3] <- 1
cmat4[group4] <- 1
cmat5[group5] <- 1

 cmat1z <-cmat2z <- cmat3z <- cmat4z <- cmat5z <- rep(0,Zcat)
 cmat1z[unique(popcts_tl$zcategorical[group1])] <- 1
 cmat2z[unique(popcts_tl$zcategorical[group2])] <- 1
 cmat3z[unique(popcts_tl$zcategorical[group3])] <- 1
 cmat4z[unique(popcts_tl$zcategorical[group4])] <- 1
 cmat5z[unique(popcts_tl$zcategorical[group5])] <- 1
# 
# sampcat_groupings <- list(group1, group2, group3, group4, group5)
# df$subgroup <- 0
# for(i in 1:5){
#   df$subgroup[df$J_cell %in% sampcat_groupings[[i]]] <- i
# }


# poptbl <- df %>%
#   group_by(J_cell,subgroup) %>%
#   summarise(
#     pinc = mean(p_include),
#     pxz = mean(pxz),
#     cts = n(),
#     sampct = pinc*cts,
#     wts = cts/sampct,
#     Za = mean(as.numeric(Za_categorical)),
#     Zb =  mean(as.numeric(Zb_categorical)),
#     Zc =  mean(as.numeric(Zc_categorical)),
#     cellmeans = mean(Y),
#     cellvars=var(Y),
#     avgX1=mean(X.1))

## order subgroups by proportion of X=1 
# group1 <- which(poptbl$Za==1 & poptbl$Zb==3)
# group1 <- which(poptbl$Zb==2 & poptbl$avgX1==1)
# group2 <- which(poptbl$Zb==4)
# group2 <- which(poptbl$Zb==4)
# group3 <- which(poptbl$Za==2)
# group4 <- which(poptbl$Za == 2 & poptbl$Zc == 2)
# group5 <- which(poptbl$Zb==4 & poptbl$Zc==2)
# 
# cmat1 <-cmat2 <- cmat3 <- cmat4 <- cmat5 <- rep(0,J)
# cmat1[group1] <- 1
# cmat2[group2] <- 1
# cmat3[group3] <- 1
# cmat4[group4] <- 1
# cmat5[group5] <- 1
# cmat1z <-cmat2z <- cmat3z <- cmat4z <- cmat5z <- rep(0,Zcat)
# cmat1z[unique(popcts_tl$zcategorical[group1])] <- 1
# cmat2z[unique(popcts_tl$zcategorical[group2])] <- 1
# cmat3z[unique(popcts_tl$zcategorical[group3])] <- 1
# cmat4z[unique(popcts_tl$zcategorical[group4])] <- 1
# cmat5z[unique(popcts_tl$zcategorical[group5])] <- 1


#-- Finally, deselect all irrelevant columns from df to save memory
df <- df %>%
  dplyr::select(
    Z_categorical,Za_categorical,Zb_categorical,Zc_categorical,pxz,
    p_include, X.1, Y, J_cell, zcts, Nj
  )
# summary(glm(X.1 ~ Zb_categorical*Zc_categorical + Za_categorical*Zc_categorical, data=df, family = "binomial"))
# mod1 <- glm(X.1 ~ Zb_categorical*Zc_categorical + Za_categorical*Zc_categorical, data=samp, family = "binomial")
# pred1 <- predict.glm(mod1,newdata=df,type="response")
# mod2 <- glm(X.1 ~ Za_categorical + Zb_categorical + Zc_categorical, data=samp, family = "binomial")
# pred2 <- predict.glm(mod2,newdata=df, type="response")
# pROC::roc(df$X.1~pred1, data=df)
# pROC::roc(df$X.1~pred2, data=df)
# 
# 
# mod1 <- glm(X.1 ~ Zb_categorical*Zc_categorical + Za_categorical*Zc_categorical, data=df, family = "binomial")
# pred1 <- predict.glm(mod1,newdata=df,type="response")
# mod2 <- glm(X.1 ~ Za_categorical + Zb_categorical + Zc_categorical, data=df, family = "binomial")
# pred2 <- predict.glm(mod2,newdata=df, type="response")
# pROC::roc(df$X.1~pred1, data=df)
# pROC::roc(df$X.1~pred2, data=df)
