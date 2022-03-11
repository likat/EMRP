# MRP RSTAN -------
# purpose: run multilevel regression model with (Z) for classical MRP
# call: "mrp_sim.stan"
# input: 
#     sample
# output: <double> estimated cell means[M] (number of Z categories)

rstan.mrp <- function(samp=sampled){
  Ma <- nlevels(samp$Za_categorical)
  Mb <- nlevels(samp$Zb_categorical)
  Mc <- nlevels(samp$Zc_categorical)
  J_stan <- sampled_Zlength
  n <- nrow(samp)
  y <- as.vector(samp$Y)
  cell_label<-samp$sampledzlab
  
  grptbl <- samp %>% 
    group_by(Z_categorical) %>% 
    summarise(zalab = mean(as.numeric(Za_categorical)),
              zblab = mean(as.numeric(Zb_categorical)),
              zclab = mean(as.numeric(Zc_categorical)-1),
    )
  agroup <- grptbl$zalab
  bgroup <- grptbl$zblab
  cgroup <- grptbl$zclab
  #-- run the model
  
  # mod <- stan(file= "mrp_sim_Zonly.stan",
  #             data = c("Ma", "Mb", "Mc","J_stan","n","y","agroup","bgroup","cgroup", "cell_label"),
  #             iter=staniters, 
  #             pars=c("cellmean"),
  #             warmup = staniters-sim/2, control=list(adapt_delta=0.99, max_treedepth=13),chains=2)

  model <- stan_model(file= "mrp_sim_Zonly.stan", model_name = "mod")
  sfit <- sampling(model,
                   data = c("Ma", "Mb", "Mc","J_stan","n","y","agroup","bgroup","cgroup", "cell_label"),
                   iter=staniters, 
                   pars=c("cellmean"),
                   warmup = staniters-sim/2, control=list(adapt_delta=0.99, max_treedepth=13),chains=2)
  unload.dll(model=model)
  return(extract(sfit,permuted = TRUE, inc_warmup = FALSE,include = TRUE)$cellmean)
}
