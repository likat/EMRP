# EMRP RSTAN -------
# purpose: run multilevel regression model with (Z,X) for EMRP functions
# call: "mrp_sim.stan"
# input: 
  # sample
# output: <double> estimated cell means[sampledJ]
rstan.emrp <- function(){
  Ma <- nlevels(sampled$Za_categorical)
  Mb <- nlevels(sampled$Zb_categorical)
  J_stan <- sampled_Jlength
  n <- nrow(sampled)
  y <- as.vector(sampled$Y)
  cell_label <- sampled$sampledJlab
  grptbl <- sampled %>% 
    group_by(J_cell) %>% 
    summarise(zalab = mean(as.numeric(Za_categorical)),
              zblab = mean(as.numeric(Zb_categorical)),
              zclab = mean(as.numeric(Zc_categorical)-1),
              xlab = mean(X.1))
  agroup <- grptbl$zalab
  bgroup <- grptbl$zblab
  cgroup <- grptbl$zclab
  xgroup <- grptbl$xlab 
  # mod <- stan(file= "mrp_sim.stan",
  #                          data = c("Ma", "Mb", "J_stan","n","y","agroup","bgroup","cgroup","xgroup", "cell_label"),
  #                          iter=staniters, 
  #                          pars = c("cellmean"),
  #                          warmup = staniters-sim/2, control=list(adapt_delta=0.99, max_treedepth=13),chains=2)
  # 
  
  model <- stan_model(file= "mrp_sim.stan", model_name = "mod")
  sfit <- sampling(model,
                   data = c("Ma", "Mb", "J_stan","n","y","agroup","bgroup","cgroup","xgroup", "cell_label"),
                   iter=staniters, 
                   pars = c("cellmean"),
                   warmup = staniters-sim/2, control=list(adapt_delta=0.99, max_treedepth=13),chains=2)
  unload.dll(model=model)
  
  cellmeans <- extract(sfit,permuted = TRUE, inc_warmup = FALSE,include = TRUE)$cellmean
  return(cellmeans)
}
