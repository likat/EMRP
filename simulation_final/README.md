# Embedded MRP Simulation Code

Multilevel Regression and Poststratification (MRP) is a powerful method used in small area estimation known for its stabilizing and bias reduction properties. However, the structure of its multilevel model is constrained by availability of the population distribution of poststratifying variables, which are often derived from auxiliary data such as a census. Embedded-MRP (EMRP) is a class of methods that "embed" the estimation of the population distribution of such unavailable predictive variables in the MRP workflow. 

We see the greatest improvement in bias from EMRP relative to Classical MRP when the distribution of $X|Z$ (following notation in our paper; the distribution of the unobserved poststratifier given the observed poststratifiers) is different in the sample relative to the population. In the simulation, we explore different cases of subdomain construction: one where the subdomain is defined by levels of $Z$ variables only and one where the subdomain is defined by the joint $Z,X$ distribution.

The simulations compare the following Embedded MRP methods:

 - WFPBB-MRP
 - Multinomial-MRP
 - Two-stage MRP

To replicate the simulation results, create a "results/" folder after forking the repo and run the "main" R scripts. "Main" scripts will call on the other R scripts listed in the documentation.

## Documentation

 - **"main_\*.R"**: Main simulation body. Calls on all functions listed in the documentation. 
 - **"genPop\*.R"**: Function for simulating the population.
 - **"rstan_emrp.R"**: Defines the function "rstan.emrp," which fits the multilevel model for Y|Z,X using "rstan_emrp.R" and returns estimated cell means for use with Embedded MRP methods.
 - **"rstan_mrp.R"**: Defines the function "rstan.mrp," which fits the multilevel model for Y|Z using "rstan_mrp.R" and returns estimated cell means for use with Classical MRP.
 - **"rstan_mrp2.R"**: Defines the function "nj.mrp2," which fits the multilevel model for X|Z using "rstan_mrp2.R" and returns estimated cell counts for use with Two-Stage MRP.
 - **"update_\*.R"**: Functions for updating containers for a given method. 
 - **"helpers.R"**: Various helper functions used in the simulation. 