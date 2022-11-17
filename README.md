# Embedded Multilevel Regression and Poststratification

Multilevel Regression and Poststratification (MRP) is a powerful method used in small area estimation known for its stabilizing and bias reduction properties. However, the structure of its multilevel model is constrained by availability of the population distribution of poststratifying variables, which are often derived from auxiliary data such as a census. Embedded-MRP is a class of methods that "embed" the estimation of the population distribution of unavailable predictive variables (which we term $X$) in the MRP workflow using information from fully observed poststratifying variables ($Z$). We have a total of $J$ unique cells constructed by combinations of $(Z,X)$.

EMRP estimators fit the same outcome model (Y|Z,X) to derive $\hat{\theta}_j$, but use different methods to estimate $\hat{N}_j$. We consider these three:

 - WFPBB-MRP: use the weighted finite population Bayesian bootstrap (WFPBB) with weights constructed from cross-tabs of $Z$
 - Multinomial-MRP: use draws from the multinomial distribution using observed cell frequencies
 - Two-stage MRP: fit a $X|Z$ model and apply the estimated probabilities to observed $Z$ cross-tabs.
 
We publish the code used to generate the results in our paper, ``Embedded Multilevel Regression and Poststratification: Model-based Inference with Incomplete Auxiliary Information" (under review) in this repository. 