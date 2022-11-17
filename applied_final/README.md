# EMRP Application: Longitudinal Study for Wellbeing

Code for EMRP analyses applied to the Longitudinal Study for Wellbeing (LSW). If you are reproducing the results from the paper, please create a "results/" folder after pulling the repo and run the code in the following order:

1. Run "model_fitting/fitmodel_emrp.R" and "model_fitting/firmodel_mrp.R" code to obtain ".rds" files from the model fit. These files will be used in the analyses later.
2. For analysis set one (using estimated inclusion probabilities), run all "analysis_*" codes in "origpinc_quant/".
3. For analysis set two (using sociodemographic-visitor subgroups), run all "analysis_*" codes in "visitor/".

**A comprehensive description of this repository**:

 * "data/" contains data from the ACS-NYC 2011 and the baseline survey
 	- "Baseline+Public+Use.csv" is the LSW public-use code.
	- "acs_nyc_2011_wpov1.dta" is the ACS data that we use to derive population characteristics (Z).
 * "loo/" contains code for doing leave-one-out diagnostics.
 * "model_fitting/" includes code used for fitting the EMRP outcome model (Y|Z,X), MRP outcome model (Y|Z), and the (X|Z) model for two-stage MRP.
 	- "BASELINE_origpinc_quantpinc.R" will clean the LSW data and create weights. You do not have to run it by itself; it is a script sourced by "fitmodel*" codes.
	- "fitmodel_emrp.R" will fit the (Y|Z,X) and (X|Z) models using RStan ("outcome_model_BASELINE2.stan"; "mrp2_Nq_BASELINE.stan", respectively) and write the stanfits to "emrp_stan.rds" and "nq_stan.rds". These RDS files will be used in EMRP analyses.
	- "fitmodel_mrp" fits the (Y|Z) model using RStan ("outcome_model_BASELINE.stan") and write the stanfit to "mrppars.rds" to be used in the classical MRP analysis.
 * "origpinc_quant/" conducts analysis on the subgroups formed by inclusion probabilities (given that the user has already run the "fitmodel") codes. 
 	- "BASELINE_origpinc_quantpinc.R" is identical to the code in "model_fitting" and again will not need to be explicitly called. It is sourced in all analysis codes.
	- "analysis_MRP_rstanavail.R" will conduct classical MRP analysis with model estimates stored in "mrppars.rds"
	- "analysis_emrp_rstanavail.R" will conduct EMRP (WFPBB-MRP, multinomial-MRP, and two-stage-MRP) analysis with model estimates stored in "emrp_stan.rds" and "nq_stan.rds".
	- "analysis_wfpbby.R" will conduct direct imputation of the outcome using the WFPBB.
* "visitor/" codes are direct analogs of those in "origpinc_quant/"
