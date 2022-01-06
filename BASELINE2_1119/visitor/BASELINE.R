#---
#====== CODE FOR BASELINE DATA ANALYSIS ========
#---
library(tidyverse)
library(dplyr)
require(foreign)
acs   <- read.dta("acs_nyc_2011_wpov1.dta")
baseline <- read.csv("Baseline+Public+Use.csv")

## 


set.seed(51)
# Impute missing outcome values and education by random draw
baseline <- baseline %>% mutate(
    # ifelse(test, yes, no)
  foodindsev = ifelse(qf2>3, sample(qf2[qf2 <=3], sum(qf2>3)), qf2),
  qi5_imp    = ifelse(qi5>8, sample(qi5[qi5<=8], sum(qi5>8)), qi5) #education
  # famsize  = adultx + childx,
  # opmfam = opmres_tc/famsize
  #bestlife=ifelse(bestlife=="NA",sample(bestlife[! bestlife=="NA"], sum(bestlife=="NA")),bestlife)
)

## Define opmres cutoffs for variable
# opmresquant <- quantile(baseline$opmres_tc)
# opmresquant <- c(0,20000, 50000, 100000, 200000)
opmresquant <- c(0,35000, 55000, 100000)
x1_poststratifiers <- c("age3", "sex2",  "race3","educat4", "povgap4")
all_poststratifiers <- c(x1_poststratifiers, "svefreq2") # combinations = J_cell

pop <- 
  acs %>% 
  dplyr::select(age, sex, race, educat,famsize, opmres, perwt) %>% 
  as.data.frame()

pop$age <- as.character(pop$age)
pop[pop$age == "Less than 1 year old", "age"] = 0
pop[pop$age == "90 (90+ in 1980 and 1990)", "age"] = 90
pop$age <- as.numeric(pop$age)
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
    # famsize = case_when(
    #   famsize == "1 family member present" ~ 1,
    #   famsize == "2 family members present" ~ 2,
    #   famsize == 3 ~ 3,
    #   famsize == 4 ~ 4,
    #   famsize == 5 ~ 5,
    #   famsize == 6 ~ 6,
    #   famsize == 7 ~ 7,
    #   famsize == 8 ~ 8,
    #   famsize == 9 ~ 9,
    #   famsize == 10 ~ 10,
    #   famsize == 11 ~ 11,
    #   famsize == 12 ~ 12,
    #   famsize == 13 ~ 13,
    #   famsize == 15 ~ 15,
    #   famsize == 16 ~ 16
    # ),
    # opmfam = opmres/famsize,
    povgap4 = factor(case_when(
      opmres <= opmresquant[2] ~1,
      opmres <= opmresquant[3] ~2,
      opmres <= opmresquant[4] ~3,
      TRUE ~4
    )),
    perwt =perwt
  ) %>%
  filter(complete.cases(.)) %>% 
  as.data.frame()


# Proposal: 4-category education variable
# 4 = Less than HS (<=2)
# 3 = HS or equivalent ( 3 or 4)
# 2 = Some college or associates (5,6)
# 1 = Bachelors or higher (7,8)
samp <- 
  baseline %>% 
  transmute(
    age3 = factor(case_when(
        qa3_age_tc <= 35 ~1,
        qa3_age_tc <= 50 ~2,
        TRUE ~ 3
    )),
    sex2 = factor(case_when(
      imp_female== 0 ~ 1,
      imp_female== 1 ~ 2
      )),
    race3 = factor(case_when(
      imp_race==1 ~ 1,
      imp_race == 2 ~ 2,
      imp_race > 2 ~ 3
    )),
    educat4 = factor(case_when(
      qi5_imp <= 2 ~ 1,
      qi5_imp <= 4 ~ 2,
      qi5_imp <= 6 ~ 3,
      qi5_imp <= 8 ~ 4
    )),
    povgap4 = factor(
      case_when(
        opmres_tc <= opmresquant[2] ~1,
        opmres_tc <= opmresquant[3] ~2,
        opmres_tc <= opmresquant[4] ~3,
        TRUE ~4 # no missing values in original dataset
      )
    ),
    pov2 = factor(opmpov),
    # famsize = factor(famsize),
    svefreq2 = factor(case_when(
      qc6 < 98 ~ 2,
      TRUE ~ 1
    )),
    foodindsev = factor(case_when(
      foodindsev >=2 ~ 0,
      foodindsev ==1 ~ 1
    )),
    qweight_pu = qweight_pu
  )

#--- STEPWISE SELECTION -----
# library(MASS)
# # Fit the full model
# full.model <- glm(foodindsev~ (age3 + sex2 + race3 + educat4 + povgap4 + svefreq2)^2, data = samp, family = "binomial")
# # Stepwise regression model
# step.model <- stepAIC(full.model, scope = list(
#   upper= ~(age3 + sex2 + race3 + educat4 + povgap4 + svefreq2)^2,
#   lower= ~age3 + sex2 + race3 + educat4 + povgap4 + svefreq2),
#   direction = "both", trace = FALSE)
# summary(step.model)
# step.modelF <- stepAIC(full.model, scope = list(
#   upper= ~(age3 + sex2 + race3 + educat4 + povgap4 + svefreq2)^2,
#   lower= ~age3 + sex2 + race3 + educat4 + povgap4 + svefreq2),
#   direction = "forward", trace = FALSE)
# summary(step.modelF)
# step.modelB <- stepAIC(full.model, scope = list(
#   upper= ~(age3 + sex2 + race3 + educat4 + povgap4 + svefreq2)^2,
#   lower= ~age3 + sex2 + race3 + educat4 + povgap4 + svefreq2),
#   direction = "backward", trace = FALSE)
# summary(step.modelB)
# 
# # glmchosen <- glm(foodindsev~(age3 + sex2 + race3 + educat4 + povgap4 + svefreq2)^2,family="binomial", data = samp)
# glmchosen <- glm(step.modelB,family="binomial", data = samp)
# 
# glmchosenpred <- predict.glm(glmchosen,type = "response")
# test_roc = roc(samp$foodindsev ~ glmchosenpred, plot = TRUE, print.auc = TRUE)
# library(pROC)
# pROC::roc(as.numeric(samp$foodindsev)-1,glmchosenpred)
#age:sex
#sex:race
#race:educat
#povgap:svefreq

#---- CREATE WEIGHTS ------
samp_X1 <- 
  samp %>% 
  dplyr::select(c("age3", "sex2",  "race3","educat4", "povgap4")) %>% 
  group_by_all() %>% 
  summarise(x1_sampcts = n()) 
pop_X1 <- 
  pop %>% dplyr::group_by(age3, sex2, race3,educat4,povgap4) %>%   
  summarise(x1_popcts = sum(perwt))

wt_tbl <- 
  left_join(samp_X1, pop_X1) %>% 
  mutate(wts = x1_popcts/x1_sampcts)

wt_tbl$sampled_x1_label <- seq(1:nrow(wt_tbl))

#---- CREATE INTERACTION TERMS ----
#age:sex
#sex:race
#race:educat
#povgap:svefreq

samp_abtbl <-
  samp %>% dplyr::group_by(age3, sex2) %>% summarise(n()) %>% dplyr::select(-`n()`)
  samp_abtbl$abcat <- 1
  samp_abtbl$abcat[samp_abtbl$sex2==2] <- 1:3
  # samp_abtbl$abcat <- 1:nrow(samp_abtbl)
samp_bctbl <-
  samp %>% dplyr::group_by(sex2,race3) %>% summarise(n()) %>% dplyr::select(-`n()`)
  samp_bctbl$bccat <- 1
  samp_bctbl$bccat[samp_bctbl$sex2==2] <- 1:3
  # samp_bctbl$bccat <- 1:nrow(samp_bctbl)
samp_cdtbl <-
  samp %>% dplyr::group_by(race3,educat4) %>% summarise(n()) %>% dplyr::select(-`n()`)
  samp_cdtbl$cdcat <- 1:nrow(samp_cdtbl)
samp_dftbl <-
  samp %>% dplyr::group_by(povgap4, svefreq2) %>% summarise(n()) %>% dplyr::select(-`n()`)
  samp_dftbl$dfcat <- 1
  samp_dftbl$dfcat[samp_dftbl$svefreq2==2] <- 1:4
  # samp_dftbl$dfcat <- 1:nrow(samp_dftbl)

# samp_dftbl <- samp_dftbl[samp_dftbl$svefreq2==2,]
# samp_dftbl$dfcat <- 1:nrow(samp_dftbl)

samp <- left_join(samp, wt_tbl) %>%  
  left_join(samp_abtbl) %>% 
  left_join(samp_dftbl) %>% 
  left_join(samp_cdtbl) %>% 
  left_join(samp_bctbl)
  
samp$abcat[is.na(samp$abcat)] <- nrow(samp_abtbl)+1
samp$bccat[is.na(samp$bccat)] <- nrow(samp_bctbl)+1
samp$cdcat[is.na(samp$cdcat)] <- nrow(samp_cdtbl)+1
samp$dfcat[is.na(samp$dfcat)] <- nrow(samp_dftbl)+1

samp$abcat <- factor(samp$abcat)
samp$bccat <- factor(samp$bccat)
samp$cdcat <- factor(samp$cdcat)
samp$dfcat <- factor(samp$dfcat)

samp_allX <- samp %>%
  dplyr::select(c("age3", "sex2",  "race3","educat4", "povgap4","svefreq2"))  %>%
  group_by_all() %>% 
  summarise(nj = n())

# number of cells
J = nrow(samp_allX) # 468
#N <- nrow(pop)
N = sum(pop$perwt)
samp_allX$J_cell = seq(1,J)

samp <- samp %>% 
  full_join(samp_allX)
group_x1_label <- samp %>% group_by(J_cell) %>% summarise(x1lab = mean(sampled_x1_label))
group_x1_label <- group_x1_label$x1lab


#---- GROUP INDICATORS -----
# suggested groups: 
## choose a variety of small (larger diff from WFPBB) and large groups (smaller diff from WFPBB)
#  BLACK, and LOW INCOME and LOW EDUCATION (n=202) GOOD
#  BLACK OR OTHER, YOUNG, LOW INCOME (n=131) NOT GREAT
#  LOW INCOME, LOW EDUCATION, YOUNG (n=31)
#  FEMALE, BLACK, LOW INCOME (n=202)
#  FEMALE, LOW INCOME, AGE = 2 (n=77)
#  FEMALE, povgap = 2, age=1 (n=69 + 51 = 120)
#age:sex
#sex:race
#race:educat
#povgap:svefreq

# vary n (n < 100, n > 100)

# grp1: svefreq2 == 2
# grp2: svefreq2 == 2 & povgap4 == 4
# grp3: svefreq2 == 1 & povgap4 == 1
# grp4: race3 == 2 & povgap4 == 1
# grp5: sex2 == 1 & svefreq2 == 2
# grp6: age3 == 1 & svefreq2 == 2
# grp7: educat4 == 4 & svefreq2 == 2

grptbl <- samp %>% group_by(age3,sex2,race3,educat4,povgap4,svefreq2) %>% summarise(
  grp1id = mean(svefreq2 == 2), #visitor n=626
  grp2id = mean(svefreq2 == 2 & povgap4 == 4), #n=79
  grp3id = mean(svefreq2 == 1 & povgap4 == 1), #poor nonvis n=324
  grp4id = mean(race3 == 2 & povgap4 == 1), # n=202
  grp5id = mean(sex2 == 1 & svefreq2 == 2), # male vis n=218
  grp6id = mean(age3 == 1 & svefreq2 == 2), # n=212
  grp7id = mean(educat4 == 4 & svefreq2 == 2)) # highed vis n=123
grplabs <- c("visitor", "rich vis", "poor nonvis", "poor black",
             "male vis", "young vis", "highed vis")
# grptbl <- samp %>% group_by(age3,sex2,race3,educat4,povgap4,svefreq2) %>% summarise(
#   grp1id = mean(age3==1),#used
#   grp2id = mean(educat4 %in% c(1,2)), #used
#   grp3id = mean(povgap4==1),# used
#   grp4id = mean(race3==2), # used
#   grp5id = mean(povgap4 %in% c(1,2) & sex2==1),
#   grp6id = mean(svefreq2==1),
#   grp7id = mean(sex2==2 & svefreq2==1)) # used

grp1id <- grptbl$grp1id
grp2id <- grptbl$grp2id
grp3id <- grptbl$grp3id
grp4id <- grptbl$grp4id
grp5id <- grptbl$grp5id
grp6id <- grptbl$grp6id
grp7id <- grptbl$grp7id
samp <- left_join(samp, grptbl)

library(survey)
ybar <- rep(0,5)
cilower <- rep(0,5)
ciupper <- rep(0,5)

wtdes <- svydesign(ids=~0, weights =~qweight_pu, data = samp)
bwtdes <- as.svrepdesign(wtdes, type="bootstrap", replicates=50)
ybar[1] <- svymean(~foodindsev,bwtdes)[2]
cilower[1] <- confint(svymean(~foodindsev,bwtdes))[2]
ciupper[1] <- confint(svymean(~foodindsev,bwtdes))[4]

#             mean     SE
# foodindsev1 0.093374 0.0091
tempdf <- samp[samp$grp1id==1,]
wtdes <- svydesign(ids=~0, weights =~qweight_pu, data = tempdf)
bwtdes <- as.svrepdesign(wtdes, type="bootstrap", replicates=50)
ybar[2] <- svymean(~foodindsev,bwtdes)[2]
cilower[2] <- confint(svymean(~foodindsev,bwtdes))[2]
ciupper[2] <- confint(svymean(~foodindsev,bwtdes))[4]


# foodindsev1 0.20136 0.0226
tempdf <- samp[samp$grp3id==1,]
wtdes <- svydesign(ids=~0, weights =~qweight_pu, data = tempdf)
bwtdes <- as.svrepdesign(wtdes, type="bootstrap", replicates=50)
ybar[3] <- svymean(~foodindsev,bwtdes)[2]
cilower[3] <- confint(svymean(~foodindsev,bwtdes))[2]
ciupper[3] <- confint(svymean(~foodindsev,bwtdes))[4]

# foodindsev1 0.13168 0.0178
tempdf <- samp[samp$grp5id==1,]
wtdes <- svydesign(ids=~0, weights =~qweight_pu, data = tempdf)
bwtdes <- as.svrepdesign(wtdes, type="bootstrap", replicates=50)
ybar[4] <- svymean(~foodindsev,bwtdes)[2]
cilower[4] <- confint(svymean(~foodindsev,bwtdes))[2]
ciupper[4] <- confint(svymean(~foodindsev,bwtdes))[4]

# foodindsev1 0.18062 0.0415
tempdf <- samp[samp$grp7id==1,]
wtdes <- svydesign(ids=~0, weights =~qweight_pu, data = tempdf)
bwtdes <- as.svrepdesign(wtdes, type="bootstrap", replicates=50)
ybar[5] <- svymean(~foodindsev,bwtdes)[2]
cilower[5] <- confint(svymean(~foodindsev,bwtdes))[2]
ciupper[5] <- confint(svymean(~foodindsev,bwtdes))[4]
# foodindsev1 0.092462 0.0373

svywtdres <- data.frame(
  group = c("overall",grplabs[c(1,3,5,7)]),
  mean = ybar,
  CIlower = cilower,
  CIupper = ciupper
) %>% mutate(CIlength = CIupper-CIlower)

write.csv(svywtdres, "svywtdres.csv")