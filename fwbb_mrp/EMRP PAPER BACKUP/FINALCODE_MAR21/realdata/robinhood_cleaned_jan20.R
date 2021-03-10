# ==============================================
# Code for wrangling Robin Hood and ACS datasets
# ==============================================
#- data used for cases 2 and 3 in applications
#- "sample" used to identify SRBI and Agency subsets

require(foreign)
require(dplyr)
acs   <- read.dta("acs_nyc_2011_wpov1.dta")
#acs = readRDS("~/Box/Kli/EMRP/robinhood_jan2020/acs.rds")
#lsw = readRDS("~/Box/Kli/EMRP/robinhood_jan2020/lsw.rds")
srbi  <- read.dta("SRBIandAGENCYbaselinew1w2w3datawithwghts072114.dta")
#-- checking data characteristics
# DIMENSIONS
# dim(acs) # nrow = 66246, ncol = 75
# dim(srbi) # nrow = 2228, ncol = 2579
srbi_nowt <- srbi[ ,-grep("wt", colnames(srbi))] # exclude weight variables

# agecat4: 4 levels of age
# sex2: 1= Male, 2=Female
# race5: 1=White, 2=Black, 3=Asian or Native Hawaiian/Pacific Islander, 4=American Indian/Alaska Native, 5=Other
# educat4: factored educat variable
# povgap4: factored povgap variable (levels for <50, 50-100, 100-200, 200-300, 300+)

#-- X2: svefreq: 0 = once or twice in past year, 1 =  more frequently than once or twice yearly
# 0. No service use
# 1. Every day
# 2. Every week
# 3. Every month
# 4. Several/few times a year
# 5. Once or twice in the past year
#-- Y: bestlife_num^3 (bestlife as continuous)

# Cleaning and variable manipulation for SAMPLE
x1_poststratifiers <- c("age3", "sex2",  "race3","educat4", "povgap4")
all_poststratifiers <- c(x1_poststratifiers, "svefreq2") # combinations = J_cell

pop <- 
  acs %>% 
  dplyr::select(age, sex, race, educat, poverty, perwt) %>% 
  as.data.frame()

pop$age <- as.character(pop$age)
pop[pop$age == "Less than 1 year old", "age"] = 0
pop[pop$age == "90 (90+ in 1980 and 1990)", "age"] = 90
pop$age <- as.numeric(pop$age)
# popage_quartiles <- quantile(pop$age, c(0.25,0.5,0.75))
#popage_tertile <- quantile(pop$age, c(0.33, 0.66))
pop <- 
  pop %>% 
  transmute(
    age3 = factor(case_when(
      age <= 35 ~1,
      age <= 50 ~2,
      TRUE ~ 3
    )),
    sex2 = factor(case_when(
      sex == "Male" ~ 1,
      sex == "Female" ~2
    )),
    race3 = factor(case_when(
      race == "White" ~ 1,
      race == "Black/Negro" ~ 2,
      TRUE ~ 3
    )),
    educat4 = factor(
      educat
    ),
    povgap4 = factor(case_when(
      poverty < 100 ~ 1,
      poverty < 200 ~ 2,
      poverty < 300 ~ 3,
      poverty >= 300 ~ 4
    )),
    perwt =perwt
  ) %>%
  filter(complete.cases(.)) %>% 
  as.data.frame()

srbi_nowt = srbi_nowt %>%
 mutate(
   foodindsev=ifelse(foodindsev=="NA",sample(foodindsev[! foodindsev=="NA"], sum(foodindsev=="NA")),foodindsev)
   #bestlife=ifelse(bestlife=="NA",sample(bestlife[! bestlife=="NA"], sum(bestlife=="NA")),bestlife)
)

samp <-   
  srbi_nowt %>% 
  dplyr::select(age,r_gender, race_dc, educat, povgap,sample, svefreq, foodindsev) %>% 
  transmute(
    age3 = factor(case_when(
      age <= 35 ~1,
      age <= 50 ~2,
      TRUE ~ 3
    )),
    sex2 = factor(case_when(
      r_gender == "1. Male" ~ 1,
      r_gender == "2.Female" ~ 2)),
    race3 = factor(case_when(
      race_dc == 1 ~ 1,
      race_dc == 2 ~ 2,
      race_dc > 2 ~ 3
    )),
    educat4 = factor(
      educat),
    povgap4 = factor(
      case_when(
        povgap == "1 Under 50%" ~1,
        povgap == "2 50-100%" ~1,
        povgap == "3 100-200%" ~2,
        povgap == "4 200-300%" ~3,
        povgap == "5 300%+" ~4
      )
    ),
    samptype = factor(sample),
    # svefreq2 = factor(case_when(
    #   svefreq == 5 ~ 0,
    #   TRUE ~ 1)),
    svefreq2 = factor(case_when(
      svefreq <= 2 ~ 1,
      TRUE ~ 2
    )),
    foodindsev = as.numeric(as.character(foodindsev))
    #bestlife_num = as.numeric(as.character(bestlife))
  ) %>% 
  filter(complete.cases(.)) %>% 
  as.data.frame()
sampt1 <-   
  srbi_nowt %>% 
  dplyr::select(age,r_gender, race_dc, educat, povgap,sample, svefreq, foodindsev) %>% 
  transmute(
    age3 = factor(case_when(
      age <= 35 ~1,
      age <= 50 ~2,
      TRUE ~ 3
    )),
    sex2 = factor(case_when(
      r_gender == "1. Male" ~ 1,
      r_gender == "2.Female" ~ 2)),
    race3 = factor(case_when(
      race_dc == 1 ~ 1,
      race_dc == 2 ~ 2,
      race_dc > 2 ~ 3
    )),
    educat4 = factor(
      educat),
    povgap4 = factor(
      case_when(
        povgap == "1 Under 50%" ~1,
        povgap == "2 50-100%" ~1,
        povgap == "3 100-200%" ~2,
        povgap == "4 200-300%" ~3,
        povgap == "5 300%+" ~4
      )
    ),
    samptype = factor(sample),
    # svefreq2 = factor(case_when(
    #   svefreq == 5 ~ 0,
    #   TRUE ~ 1)),
    svefreq2 = factor(case_when(
      svefreq <= 2 ~ 1,
      TRUE ~ 2
    )),
    foodindsev = as.numeric(as.character(foodindsev))
    #bestlife_num = as.numeric(as.character(bestlife))
  ) %>% 
  # filter(complete.cases(.)) %>% 
  as.data.frame()

#-- SRBI only X1 table to create weights
samp_X1_srbi <- 
  samp %>% filter(samptype=="SRBI") %>% 
  dplyr::select(c("age3", "sex2",  "race3","educat4", "povgap4")) %>% 
  group_by_all() %>% 
  summarise(x1_sampcts = n()) 

samp_aetbl <- 
  samp %>% dplyr::group_by(educat4, povgap4) %>% summarise(n()) %>% dplyr::select(-`n()`)
samp_aetbl$aecat <- factor(1:nrow(samp_aetbl))


pop_X1 <- 
  pop %>% dplyr::group_by(age3, sex2, race3,educat4,povgap4) %>%   
  summarise(x1_popcts = sum(perwt))

wt_tbl <- 
  left_join(samp_X1_srbi, pop_X1) %>% 
  mutate(wts = x1_popcts/x1_sampcts)

wt_tbl$sampled_x1_label <- seq(1:nrow(wt_tbl))

#-- create easy way of referencing the cells created by {X1} x {X2}  
#-- assemble new sampled dataset with J_cell membership and sampled X1 counts
# samp <- left_join(samp, wt_tbl) %>% 
#   full_join(samp_allX) %>% arrange(age3, sex2, race5, educat4, povgap4, svefreq2)
samp <- left_join(samp, wt_tbl) %>%  left_join(samp_aetbl) 
samp <- samp[-which(is.na(samp$wts)),] # 11 Agency observations excluded
samp_allX <- samp %>%
  dplyr::select(all_poststratifiers)  %>%
  group_by_all() %>% 
  summarise(nj = n())

# number of cells
J = nrow(samp_allX) # 487
#N <- nrow(pop)
N = sum(pop$perwt)
samp_allX$J_cell = seq(1,J)

samp <- samp %>% 
  full_join(samp_allX)
group_x1_label <- samp %>% group_by(J_cell) %>% summarise(x1lab = mean(sampled_x1_label))
group_x1_label <- group_x1_label$x1lab
samp <- samp[!is.na(samp$J_cell),] # removes 11 observations who have Z cells unique to Agency (and not SRBI)

agencysamp <- samp %>% filter(samptype == "Agency")
srbisamp <- samp %>% filter(samptype == "SRBI")


grptbl <- samp %>% group_by(age3,sex2,race3,educat4,povgap4,svefreq2) %>% summarise(
  grp1id = mean(age3==1),#NEW
  grp2id = mean(educat4 %in% c(1,2)), #NEW
  grp3id = mean(povgap4==1),
  grp4id = mean(race3==2),
  grp5id = mean(povgap4 %in% c(1,2) & sex2==1),
  grp6id = mean(svefreq2==1),
  grp7id = mean(sex2==2 & svefreq2==1))
grp1id <- grptbl$grp1id
grp2id <- grptbl$grp2id
grp3id <- grptbl$grp3id
grp4id <- grptbl$grp4id
grp5id <- grptbl$grp5id
grp6id <- grptbl$grp6id
grp7id <- grptbl$grp7id
samp <- left_join(samp, grptbl)

table(samp$grp1id, samp$samptype)[2,]
table(samp$grp2id, samp$samptype)[2,]
table(samp$grp3id, samp$samptype)[2,]
table(samp$grp4id, samp$samptype)[2,]
table(samp$grp5id, samp$samptype)[2,]
table(samp$grp6id, samp$samptype)[2,]
table(samp$grp7id, samp$samptype)[2,]
# samp %>% group_by(age3,povgap4, samptype) %>% summarise(mean(bestlife_num))
# samp %>% group_by(age3,educat4, samptype) %>% summarise(mean(bestlife_num))
# View(samp %>% group_by(educat4, povgap4,samptype) %>% summarise(cts=n(),mean(bestlife_num)))
# View(samp %>% group_by(povgap4,educat4, hispan2,samptype) %>% summarise(mean(bestlife_num)))
# View(samp %>% group_by(povgap4,sex2, hispan2,samptype) %>% summarise(mean(bestlife_num)))

# temp <- left_join(samp, grptbl)

# View(samp %>% group_by(age3, povgap4,samptype) %>% summarise(mean(bestlife_num), cts=n()))
# View(temp %>% group_by(age3,educat4,samptype) %>% summarise(mean(nj), mean(bestlife_num)))
# temp %>% group_by(grp1id) %>% summarise(mean(nj),var(bestlife_num), mean(bestlife_num), n())
# temp %>% group_by(grp2id) %>% summarise(mean(nj),var(bestlife_num), mean(bestlife_num), n())
# temp %>% group_by(grp3id) %>% summarise(mean(nj),var(bestlife_num), mean(bestlife_num), n())
# temp %>% group_by(grp4id) %>% summarise(mean(nj),var(bestlife_num), mean(bestlife_num), n())

# temp %>% group_by(grp1id, samptype) %>% summarise(mean(nj),var(bestlife_num), mean(bestlife_num),n())
# temp %>% group_by(grp2id, samptype) %>% summarise(mean(nj),var(bestlife_num), mean(bestlife_num),n())
# temp %>% group_by(grp3id, samptype) %>% summarise(mean(nj),var(bestlife_num), mean(bestlife_num),n())
# temp %>% group_by(grp4id, samptype) %>% summarise(mean(nj),var(bestlife_num), mean(bestlife_num),n())


# #--- LASSO -----something is wrong here
# library(glmnet)
# 
# # Loading the data
# x_vars <- model.matrix(foodindsev~ (age3 + sex2 + race3 + educat4 + povgap4 + svefreq2)^2-1, data = samp)
# y_var <- samp$foodindsev #^3
# lambda_seq <- seq(0.001, 1, by=1)
# #
# set.seed(50)
# pterm <- c(rep(0,17),rep(1,ncol(x_vars)-17))
# cv_output <- cv.glmnet(x_vars, y_var,family = "binomial", type.measure = "class",
#                        alpha = 1, lambda = lambda_seq,
#                        nfolds = 10, penalty.factor=pterm)
# # identifying best lamda
# best_lam <- cv_output$lambda.min
# lasso_best <- glmnet(x_vars, y_var, family = "binomial", alpha = 1, lambda = best_lam,penalty.factor=pterm)
# Math.cbrt <- function(x) {
#   sign(x) * abs(x)^(1/3)
# }
# Math.cbrt(lasso_best$beta) # age:povgap

#
#
# #--- STEPWISE SELECTION -----
# library(MASS)
# # Fit the full model
# full.model <- glm(foodindsev~ (age3 + sex2 + race3 + educat4 + povgap4 + svefreq2)^2, data = samp, family = "binomial")
# # Stepwise regression model
# step.model <- stepAIC(full.model, scope = list(
#   upper= ~(age3 + sex2 + race3 + educat4 + povgap4 + svefreq2)^2,
#   lower= ~age3 + sex2 + race3 + educat4 + povgap4 + svefreq2),
#   direction = "both", trace = FALSE)
# summary(step.model) 
