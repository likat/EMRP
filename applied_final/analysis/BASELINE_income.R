#---
#====== CODE FOR BASELINE DATA ANALYSIS ========
#---
library(tidyverse)
library(dplyr)
require(foreign)
library(survey)
acs   <- read.dta("../data/acs_nyc_2011_wpov1.dta")
baseline <- read.csv("../data/Baseline+Public+Use.csv")
type = "income"

set.seed(51)
# Impute missing outcome values and education by random draw
baseline <- baseline %>% mutate(
  # ifelse(test, yes, no)
  foodindsev = ifelse(qf2>3, sample(qf2[qf2 <=3], sum(qf2>3)), qf2),
  qi5_imp    = ifelse(qi5>8, sample(qi5[qi5<=8], sum(qi5>8)), qi5) #education
)

## Define opmres cutoffs for variable

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

samp_abtbl <-
  samp %>% dplyr::group_by(age3, sex2) %>% summarise(n()) %>% dplyr::select(-`n()`)
samp_abtbl$abcat <- 1
samp_abtbl$abcat[samp_abtbl$sex2==2] <- 1:3
samp_bctbl <-
  samp %>% dplyr::group_by(sex2,race3) %>% summarise(n()) %>% dplyr::select(-`n()`)
samp_bctbl$bccat <- 1
samp_bctbl$bccat[samp_bctbl$sex2==2] <- 1:3
samp_cdtbl <-
  samp %>% dplyr::group_by(race3,educat4) %>% summarise(n()) %>% dplyr::select(-`n()`)
samp_cdtbl$cdcat <- 1:nrow(samp_cdtbl)
samp_dftbl <-
  samp %>% dplyr::group_by(povgap4, svefreq2) %>% summarise(n()) %>% dplyr::select(-`n()`)
samp_dftbl$dfcat <- 1
samp_dftbl$dfcat[samp_dftbl$svefreq2==2] <- 1:4

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

grptbl <- 
  samp %>% 
  group_by(age3, sex2, race3,educat4,povgap4, svefreq2,J_cell) %>% 
  summarise(nj = mean(nj),
            # pinc_est = nj / sum(qweight_pu))
            pinc_est = mean(1/wts))

# quantsampcts <- quantile(grptbl$nj)
quantsampcts <- quantile(unique(grptbl$nj))
quantpinc <- quantile(unique(grptbl$pinc_est))
# quantpinc <-quantile(unique(grptbl$pinc_est), c(0.1,0.2, 0.3, 0.4, 0.5, 0.6,0.7,0.8,0.9))
# educat43          -2.6721490  0.001113 **
# educat44          -2.0533153  0.031475 *
# age33:educat43     1.3980926 0.005927 **
# age33:povgap42    -1.8676514 0.000925 ***
# age32:povgap44    -1.4490616  *
# age33:povgap44    -1.4794415 *
# educat44:povgap43 -1.5538706 0.020421 *
# educat43:povgap44  2.4128150 0.023203 *
sampdes <- svydesign(ids=~0, weights=~qweight_pu, data=samp)

modxz <- svyglm(svefreq2 ~( age3 + sex2 + race3 + educat4 + povgap4)^2, design = sampdes, family="binomial")
summary(modxz)
# mody <- svyglm(foodindsev ~svefreq2+ svefreq2*povgap4+age3 + sex2 + race3 + educat4 + povgap4, design = sampdes, family="binomial")
mody <- glm(foodindsev ~svefreq2+ svefreq2*povgap4+age3 + sex2 + race3 + educat4 + povgap4, data=samp,  family="binomial")

summary(mody)


## 
# modymrp <- svyglm(foodindsev ~age3 + sex2 + race3 + educat4 + povgap4, design = sampdes, family="binomial")
modymrp <- glm(foodindsev ~age3 + sex2 + race3 + educat4 + povgap4, data=samp, family="binomial")
summary(modymrp)

pred_mody <- predict.glm(mody,type = "response")
pred_modymrp <- predict.glm(modymrp,type = "response")
pROC::auc(samp$foodindsev, pred_mody)
pROC::auc(samp$foodindsev, pred_modymrp)
# 0,35000, 55000, 100000

grptbl <- samp %>% group_by(age3,sex2,race3,educat4,povgap4,svefreq2) %>% summarise(
  grp1id = mean(povgap4==1), # 0-35k
  grp2id = mean(povgap4==2), # 35-55 k
  grp3id = mean(povgap4==3), # 55-100k
  grp4id = mean(povgap4==4)) # >100k
grplabs <- c("<$35k", "$35-55k", "$55-100k", ">$100k")
grp1id <- grptbl$grp1id
grp2id <- grptbl$grp2id
grp3id <- grptbl$grp3id
grp4id <- grptbl$grp4id
samp <- left_join(samp, grptbl)

# write.csv(samp,"samp.csv")
## using the person-level weights from BASELINE code, which 
# calibrate the SRBI / Agency samples to the ACS data
res.table.foodindsev <- function(des){
  sub1 <- subset(des,grp1id==1)
  sub2 <- subset(des,grp2id==1)
  sub3 <- subset(des,grp3id==1)
  sub4 <- subset(des,grp4id==1)
  resall <- svymean(~foodindsev2, des)
  res1 <- svymean(~foodindsev2, sub1)
  res2 <- svymean(~foodindsev2, sub2)
  res3 <- svymean(~foodindsev2, sub3)
  res4 <- svymean(~foodindsev2, sub4)
  
  ciall <- svyciprop(~foodindsev2, des,method="logit")
  ci1 <- svyciprop(~foodindsev2, sub1, method = "logit")
  ci2 <- svyciprop(~foodindsev2, sub2, method = "logit")
  ci3 <- svyciprop(~foodindsev2, sub3, method = "logit")
  ci4 <- svyciprop(~foodindsev2, sub4, method = "logit")
  # ciall <- attributes(svyciprop(~I(foodindsev2==1), des, method = "logit"))$ci
  # ci1 <- attributes(svyciprop(~I(foodindsev2==1), sub1, method = "logit"))$ci
  # ci2 <- attributes(svyciprop(~I(foodindsev2==1), sub2, method = "logit"))$ci
  # ci3 <- attributes(svyciprop(~I(foodindsev2==1), sub3, method = "logit"))$ci
  # ci4 <- attributes(svyciprop(~I(foodindsev2==1), sub4, method = "logit"))$ci
  resmat <- data.frame(mean = c(resall,res1,res2,res3,res4),
                       SE = c(SE(resall),SE(res1), SE(res2), SE(res3), SE(res4)),
                       row.names=c("overall", "1","2","3","4"))
  CImat <- matrix(nrow = 5,ncol=2)
  CImat[1,] <- confint(ciall)
  CImat[2,] <- confint(ci1)
  CImat[3,] <- confint(ci2)
  CImat[4,] <- confint(ci3)
  CImat[5,] <- confint(ci4)
  
  return(list(resmat*100, CImat))
}
samp$foodindsev2 <- as.numeric(samp$foodindsev)-1
unwtdes<- svydesign(ids=~0,  data = samp)
# samp$qweight_pu <- samp$qweight_pu/sum(samp$qweight_pu) * 2228
wtdes_person <- svydesign(ids=~0, weights =~qweight_pu, data = samp)
# wtdes_person <- svydesign(ids=~0, weights =~wts, data = samp)
unwtdres <- res.table.foodindsev(unwtdes)
wtdres_person <- res.table.foodindsev(wtdes_person)
sum(samp$grp1id)
sum(samp$grp2id)
sum(samp$grp3id)
sum(samp$grp4id)

unwtdres;wtdres_person
saveRDS(unwtdres, paste0("results/res_", type, "_unwtd.rds"))
saveRDS(wtdres_person, paste0("results/res_", type, "_wtd.rds"))
# table(grptbl$svefreq2[grptbl$grp1id==1])
# table(grptbl$svefreq2[grptbl$grp2id==1])
# table(grptbl$svefreq2[grptbl$grp3id==1])
# table(grptbl$svefreq2[grptbl$grp4id==1])
# 
cattbl <- tableone::CreateCatTable(data=samp, vars = c("age3", "sex2","race3","educat4","povgap4","svefreq2"))
print(cattbl, format = 'p', quote=T, cramVars=c("sex2", "svefreq2"))

cattbl1 <- tableone::CreateCatTable(data=samp[samp$grp1id==1,], vars = c("age3", "sex2","race3","educat4","povgap4","svefreq2"))
print(cattbl1, format = 'p', quote=T, cramVars=c("sex2", "svefreq2"))

cattbl2 <- tableone::CreateCatTable(data=samp[samp$grp2id==1,], vars = c("age3", "sex2","race3","educat4","povgap4","svefreq2"))
print(cattbl2, format = 'p', quote=T, cramVars=c("sex2", "svefreq2"))

cattbl3 <- tableone::CreateCatTable(data=samp[samp$grp3id==1,], vars = c("age3", "sex2","race3","educat4","povgap4","svefreq2"))
print(cattbl3, format = 'p', quote=T, cramVars=c("sex2", "svefreq2"))
cattbl4 <- tableone::CreateCatTable(data=samp[samp$grp4id==1,], vars = c("age3", "sex2","race3","educat4","povgap4","svefreq2"))
print(cattbl4, format = 'p', quote=T, cramVars=c("sex2", "svefreq2"))
# 
