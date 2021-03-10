####### Code for population generation with multiple categories in Z ########
### Z split into three factors, X binary
### Y_i|j = beta_1j*Z_i + beta_2j*X_i
### p is a vector of corresponding probabilities for each Z category (length proportional to Zcat)
### first 1:Z_cat columns are dummy variables for the categories of Z
## returns df as sorted by Zcat, X
### Data are *MAR*

# issues:
# need to adjust parameters so that in the case where there is low cor(x,z)
# and high cor(Y,X), wfpbb-mrp outperforms wfpbb because we actually account for X

require(MASS)
require(dplyr)
#set.seed(50)
set.seed(42)

## expit function
expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

N = 1e4 # population size

#-- Enumerate the categories of population covariates
## [Z] Poststratifiers known in the population
# 3 separate variables: A, B, C
acats <- 5; bcats <- 5; ccats <- 2
Zcat = acats*bcats*ccats # should be 50

## [X] Poststratifier known in sample only    
Xcat = 2; 
## [J] Number of Z:X combinations (cells)
# J = Zcat*Xcat

#-- DGP parameters
# dispersion
sig_a = 1 # sd for random intercepts given Z category
sig_b = 1
sig_c = 1
sig_bc =1
sig_y=1 # sd for data

# global mean
beta0 = 0

# main effects
alpha_a <- c(-1,-0.5,0,0.5,1)
alpha_b <- c(-1,-0.5,0,0.5,1)
alpha_c <-c(0,1)
#alpha_c <-c(0,10)
# alpha_x <- c(0,3) ## LOW cor(X,Y) CASE
alpha_x <- c(0,3) # HIGH cor(X,Y) CASE
#alpha_x <- c(0,20)
# interactions
# alpha_bc <- rnorm(bcats*ccats,0,sig_bc)
alpha_bc <- rep(0,bcats*ccats) ## deprecated, keep at 0
alpha_cx <- c(-1.5, -0.5,0.5,1.5)
# collect effects
alphavec = c(alpha_a, alpha_b, alpha_c, alpha_x,  alpha_bc, alpha_cx)

#-- Parameters relevant to selection
## for p(x|z=m) = expit(bxz*m)
## HIGH COR(X,Z)
## this increases the difference between pop and sampled means
# bxz0 <- -1.1
# bxza <- -0.5
# bxzb <- 0.5
# bxzc <- 0.1
# ## LOW COR(X,Z)
## this reduces the difference in pop and sampled means
# bxz0 <- -1.1
# bxza <- -0.05
# bxzb <- 0.05
# bxzc <- 0.05
## INDEP X,Z
# bxz0 <- -1.10 # pxz =  0.25
# bxz0 <- 0 # pxz= 0.5
# bxz0 <- 1.10 #pxz = 0.75
# bxza <- 0
# bxzb <- 0
# bxzc <- 0


## MAR CASEfor p(I=1|Z=m)
## p-include parameters
#pinc <- sample(seq(0.01,0.9,length.out=Zcat), size=Zcat, replace=F)

pinc_tbl = data.frame(Z_categorical=1:50, p_include = c(
sample(seq(0.01,0.1,length.out=10),size = 5,replace = T),
sample(seq(0.11,0.4,length.out=10),size = 10,replace = T),
sample(seq(0.21,0.6,length.out=10),size = 20,replace = T),
sample(seq(0.51,0.8,length.out=10),size = 10,replace = T),
sample(seq(0.8,0.99,length.out=10),size = 5,replace = T)))
# pinc_tbl = data.frame(Z_categorical=1:50, p_include = c(
#   sample(seq(0.01,0.1,10),size = 5,replace = T),
#   sample(seq(0.11,0.4,10),size = 10,replace = T),
#   sample(seq(0.21,0.6,10),size = 20,replace = T),
#   sample(seq(0.51,0.8,10),size = 10,replace = T),
#   sample(seq(0.8,0.99,10),size = 5,replace = T)))

# pinc <- expit(-3.5+0.1*seq(1,50))

# pinc_tbl <- data.frame(Z_categorical=1:50, p_include = c(rep(0.05, 20),
#                                                          rep(0.5, 10),
#                                                          rep(0.9, 15),
#                                                          rep(0.2, 5)))
# pinc_tbl <- data.frame(Z_categorical=1:50, p_include = c(rep(0.05, 15),
#                                                          rep(0.1,5),
#                                                          rep(0.3, 10),
#                                                          rep(0.8,10),
#                                                          rep(0.9,10)))

# binc0 <- -4
# bincz <- 0.12
# bincza <- 0.4
# binczb <- 0.6
# binczc <- 0.4
# bincx <- 0 # !=0 if MNAR
# binc0 <- -0.5
# bincza <- 0.5
# binczb <- -0.6
# binczc <- -0.2
# bincx <- 0 # !=0 if MNAR

#-- Parameters to generate Z
p_a = runif(acats, 0.5,1)
p_b = runif(bcats, 0.5,1)
p_c = runif(ccats, 0.5,1)

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
designmat <- df[,-((ncol(df)-Xcat):ncol(df))]

# get counts for each category of Z
tempdf <- 
  df %>% group_by(Za_categorical, Zb_categorical, Zc_categorical) %>% 
  summarise(zcts=n())  %>%  as.data.frame()
tempdf$Z_categorical <- 1:nrow(tempdf)

bxz0 <- -1.1
# bxza <- rnorm(acats)-0.5
bxza <- c(-1.59,-0.56,1.27,1.07, -2.41)
# bxzb <- rnorm(bcats) + 0.5
bxzb <- c(0.74, 0.33, -1.04,-1.37,2.27)
# bxzc <- rnorm(ccats) + 0.1
bxzc <- c(0.98, 2.47)

# bxzac <- c(0, -0.46,-0.61,0.75,-0.81)
# bxzbc <- c(0, 0.47, 0.55, 0.98, -0.34)
# join tables to get zcts, p_include

df <-
  df %>% 
  left_join(tempdf) %>% 
  mutate(pxz=
           expit(bxz0+bxza[Za_categorical] +
                   bxzb[Zb_categorical] +
                   bxzc[Zc_categorical]) 
 #                  bxzac[Za_categorical]*as.numeric(Zc_categorical==2)+
#                   bxzbc[Zb_categorical]*as.numeric(Zc_categorical==2))
  )
# mutate(pxz= expit(bxz0+bxz * Z_categorical))
designmat$X.0 = rbinom(N, 1, 1-df$pxz)
designmat$X.1 = as.numeric(designmat$X.0 != 1)

df$X.0 <- designmat$X.0
df$X.1 <-designmat$X.1
# create interaction terms in the design matrix
designmat <- model.matrix(~.^2-1, designmat)
selecta <- which(grepl("^A[[:digit:]]$",colnames(designmat)))
selectb <- which(grepl("^B[[:digit:]]$",colnames(designmat)))
selectc <- which(grepl("^C[[:digit:]]$",colnames(designmat)))
selectx <- which(grepl("^X.[[:digit:]]$",colnames(designmat)))
selectbc <- which(grepl("B[[:digit:]]:C[[:digit:]]",colnames(designmat)))
selectcx <- which(grepl("C[[:digit:]]:X.[[:digit:]]",colnames(designmat)))

# alphavec = c(alpha_a, alpha_b, alpha_c, alpha_x,  alpha_bc, alpha_cx)

designmat <- designmat[,c(selecta, selectb, selectc, selectx, selectbc, selectcx)] %>% as.data.frame()
df$mu <- beta0 + as.matrix(designmat) %*% alphavec

#-- [Y] Generation
df$Y <-  df$mu + rnorm(N, 0, sig_y)

#-- MNAR modification if needed
df <- left_join(df,pinc_tbl)
# df <- df %>% mutate(
#   p_include = expit(
#     binc0 + bincz*as.numeric(Z_categorical)
# bincza*as.numeric(Za_categorical) +
# binczb * as.numeric(Zb_categorical) +
# binczc * as.numeric(Zc_categorical) + 
# bincx * X.1
#   )
# ) 

# data.frame(Z_categorical = 1:50) %>% 
# mutate(p_include = expit(1-Z_categorical*bincz))# if MNAR, X effect added later

#-- Identifiers for interaction terms
bc_tbl <- 
  df %>% group_by(Zb_categorical, Zc_categorical) %>% summarise(cts=n()) %>% dplyr::select(-cts)
cx_tbl <- 
  df %>% group_by(Zc_categorical, X.1) %>% summarise(cts=n()) %>% dplyr::select(-cts)
bc_tbl$bc <- factor(1:nrow(bc_tbl))
cx_tbl$cx <- factor(1:nrow(cx_tbl))
df <- df %>%  left_join(bc_tbl)%>% left_join(cx_tbl)

#---- NOTE WE NEED TO FIX J < 100
popcts_tl = df %>% 
  group_by(Za_categorical, Zb_categorical, Zc_categorical, X.1) %>%  
  summarise(Nj = n()) %>% 
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
temp <- poptbl$pinc * poptbl$cts
ordered_groupcts <- order(temp)

grpincr <- ceiling(J/5)
group1 <- ordered_groupcts[1:grpincr]
group2 <- ordered_groupcts[(grpincr+1):(2*grpincr)]
group3 <- ordered_groupcts[(2*grpincr+1):(3*grpincr)]
group4 <- ordered_groupcts[(3*grpincr+1):(4*grpincr)]
group5 <- ordered_groupcts[(4*grpincr+1):J]
df$subgroup6 <- df$Za_categorical==1 & df$X.1 == 1
df$subgroup7 <- df$Zc_categorical==1 & df$X.1 == 1
grptbl <- df %>% 
  group_by(J_cell, Za_categorical, Zc_categorical,X.1) %>% 
  summarise(group6=mean(Za_categorical==1 & X.1 ==1),
            group7=mean(Zc_categorical==1 & X.1==1))
cmat6 <- grptbl$group6
cmat7 <- grptbl$group7

cmat1 <-cmat2 <- cmat3 <- cmat4 <- cmat5 <- rep(0,J)
cmat1[group1] <- 1
cmat2[group2] <- 1
cmat3[group3] <- 1
cmat4[group4] <- 1
cmat5[group5] <- 1

# cmat1z <-cmat2z <- cmat3z <- cmat4z <- cmat5z <- rep(0,Zcat)
# cmat1z[unique(popcts_tl$zcategorical[group1])] <- 1
# cmat2z[unique(popcts_tl$zcategorical[group2])] <- 1
# cmat3z[unique(popcts_tl$zcategorical[group3])] <- 1
# cmat4z[unique(popcts_tl$zcategorical[group4])] <- 1
# cmat5z[unique(popcts_tl$zcategorical[group5])] <- 1
sampcat_groupings <- list(group1, group2, group3, group4, group5)
df$subgroup <- 0
for(i in 1:5){
  df$subgroup[df$J_cell %in% sampcat_groupings[[i]]] <- i
}

poptbl <- df %>%
  group_by(J_cell,subgroup) %>%
  summarise(
    pinc = mean(p_include),
    pxz = mean(pxz),
    cts = n(),
    sampct = pinc*cts,
    wts = cts/sampct,
    Za = mean(as.numeric(Za_categorical)),
    Zb =  mean(as.numeric(Zb_categorical)),
    Zc =  mean(as.numeric(Zc_categorical)),
    cellmeans = mean(Y),
    cellvars=var(Y),
    avgX1=mean(X.1))

ordered_cellcts_hihi <- order(poptbl$sampct)

poptbl %>%
  group_by(subgroup) %>%
  summarise(
    varthj = var(cellmeans),
    sdpinc = sd(pinc),
    sdsampct = sd(sampct),
    meanpinc = mean(pinc),
    sdcts = sd(cts),
    avgvar=mean(cellvars),
    avgcellsamp=mean(pinc * cts),
    avgsamp= sum(pinc*cts),
    avgX1=mean(avgX1))
#  for(i in 1:5){
#    print(table(poptbl$Za[poptbl$subgroup==i]))
#  }
# for(i in 1:5){
#   print(table(poptbl$Zb[poptbl$subgroup==i]))
# }
# for(i in 1:5){
#   print(table(poptbl$Zc[poptbl$subgroup==i]))
# }
# for(i in 1:5){
#   print(table(poptbl$avgX1[poptbl$subgroup==i]))
# }
#for(i in 1:5){
#  print(sd(poptbl$pxz[poptbl$subgroup==i]))
#}
#for(i in 1:5){
#  print(mean(poptbl$pxz[poptbl$subgroup==i]))
#}


#-- Finally, deselect all irrelevant columns from df to save memory
df <- df %>%
  dplyr::select(
    Z_categorical,Za_categorical,Zb_categorical,Zc_categorical,
    bc,cx, pxz,
    p_include, X.1, Y, J_cell, zcts, subgroup,subgroup6, subgroup7,Nj
  )
# rm(list=c("bc_tbl", "cx_tbl","tempdf", "designmat", "Z", "Za", "Zb", "Zc", "Za_categorical", "Zb_categorical", "Zc_categorical",
#           "bxza", "bxzb", "bxzc", "popcts_tl","pinc_tbl", "poptbl", "temp"))
#summary(glm(X.1 ~ Za_categorical+Zb_categorical*Zc_categorical, data=df, family = "binomial"))
