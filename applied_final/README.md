# EMRP Application: Longitudinal Study for Wellbeing

Code for EMRP analyses applied to the Longitudinal Study for Wellbeing (LSW), organized into different folders:

 * "data/" contains data from the ACS-NYC 2011 and the baseline survey
 * "income/" contains code for estimating the food insecurity prevalence among different income groups
  	- "analysis_emrp_income.R" is the code for running EMRP methods (WFPBB-MRP, multinomial-MRP, and two-stage-MRP)
	- "analysis_wfpbby.R" estimates subgroup food insecurity using direct WFPBB
 	- "BASELINE_income.R" is the code that prepares both LSW and ACS-NYC datasets for analysis and defines the subgroups for analysis. 
	- "samp.csv" is the cleaned baseline sample output from "BASELINE_income.R"
	- "mrp2_Nq_BASELINE.stan" and "outcome_model_BASELINE2.stan" are scripts for running the MCMC sampler through RStan for the X|Z model in two-stage MRP and the Y|X,Z model in EMRP methods, respectively 
 * "income_nosvefreq/" contains code for estimating the food insecurity prevalence among different income groups using classical MRP (no incorporation of agency visit frequency)
 * "visitor/" contains code for estimating the food insecurity prevalence among different sociodemographic-visitor groups. Naming conventions and definitions are similar to those of "income/"
 