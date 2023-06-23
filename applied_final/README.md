# EMRP Application: Longitudinal Study for Wellbeing

This is code for the EMRP analyses applied to the Longitudinal Study for Wellbeing (LSW). If you are reproducing the results from the paper, please create a "results/" folder after pulling the repo and run the code in the following order:

1. Run "model_fitting/fitmodel_emrp.R" and "model_fitting/fitmodel_mrp.R" code to obtain ".rds" files from the model fit. These files will be used in the analyses later.
2. For analysis set one (using income subdomains), run all "analysis_\*" codes in "analysis/\*_income".
3. For analysis set two (using income-visitor subdomains), run all "analysis_\*" codes in "analysis/\*_income_visitor".

## Documentation

 * **"data/"** contains data from the ACS-NYC 2011 and the baseline survey
 	- "Baseline+Public+Use.csv" is the LSW public-use code.
	- "acs_nyc_2011_wpov1.dta" is the ACS data that we use to derive population characteristics (Z).
 * **"loo/"** contains code for doing leave-one-out diagnostics.
 * **"model_fitting/"** includes code used for fitting the EMRP outcome model (Y|Z,X), MRP outcome model (Y|Z), and the (X|Z) model for two-stage MRP.
 	- "BASELINE_income_visitor.R" will clean the LSW data and create weights. You do not have to run it by itself; it is a script sourced by "fitmodel*" codes.
	- "fitmodel_emrp.R" will fit the (Y|Z,X) and (X|Z) models using RStan ("outcome_model_BASELINE2.stan"; "mrp2_Nq_BASELINE.stan", respectively) and write the stanfits to "emrp_stan.rds" and "nq_stan.rds". These RDS files will be used in EMRP analyses.
	- "fitmodel_mrp" fits the (Y|Z) model using RStan ("outcome_model_BASELINE.stan") and write the stanfit to "mrppars.rds" to be used in the classical MRP analysis.
 * **"analysis/"** conducts analysis on the subdomains formed by income brackets (given that the user has already run the "fitmodel") codes. 
 	- "BASELINE_income.R" and "BASELINE_income_visitor.R" are used to clean and define subdomain indices and will not need to be explicitly called as it is sourced in all analysis codes.
	- "analysis_MRP_rstanavail_income.R" and "analysis_MRP_rstanavail_income_visitor.R" will conduct classical MRP analysis with model estimates stored in "mrppars.rds" for the income and income-visitor subdomains, respectively.
	- "analysis_emrp_rstanavail_income.R" and "analysis_emrp_rstanavail_income_visitor.R" will conduct EMRP (WFPBB-MRP, multinomial-MRP, and two-stage-MRP) analysis with model estimates stored in "emrp_stan.rds" and "nq_stan.rds" for the income and income-visitor subdomains, respectively.
	- "analysis_wfpbby_income.R" and "analysis_wfpbby_income_visitor.R" will conduct direct imputation of the outcome using the WFPBB for the income and income-visitor subdomains, respectively.
