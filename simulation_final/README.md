# Embedded MRP Simulation Code

Multilevel Regression and Poststratification (MRP) is a powerful method used in small area estimation known for its stabilizing and bias reduction properties. However, the structure of its multilevel model is constrained by availability of the population distribution of poststratifying variables, which are often derived from auxiliary data such as a census. Embedded-MRP is a class of methods that "embed" the estimation of the population distribution of such unavailable predictive variables in the MRP workflow. We compare the following Embedded MRP methods in a simulation:

 - WFPBB-MRP
 - Multinomial-MRP
 - Two-stage MRP

## Documentation

 - **"main_\*.R"**: Main simulation body. Calls on all functions listed in the documentation. 
 - **"genPop\*.R"**: Function for simulating the population.
 - **"rstan_emrp.R"**: Defines the function "rstan.emrp," which fits the multilevel model for Y|Z,X and returns estimated cell means for use with Embedded MRP methods.
 - **"rstan_mrp.R"**: Defines the function "rstan.mrp," which fits the multilevel model for Y|Z and returns estimated cell means for use with Classical MRP.
 - **"rstan_mrp2.R"**: Defines the function "nj.mrp2," which fits the multilevel model for X|Z and returns estimated cell counts for use with Two-Stage MRP.
 - **"update_\*.R"**: Functions for updating containers for a given method. 
 - **"helpers.R"**: Various helper functions used in the simulation. 