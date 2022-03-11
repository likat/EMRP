

nj.mrp2 <- function(){
  #-- Step 1: Estimate cell frequency given X (Nq) via MRP
  Ma <- nlevels(sampled$Za_categorical)
  Mb <- nlevels(sampled$Zb_categorical)
  Mc <- nlevels(sampled$Zc_categorical)
  J_stan <- length(unique(sampled$sampledJlab))
  n <- nrow(sampled)
  y <- sampled$X.1
  xj <- sampled$sampledJlab
  grptbl = sampled %>% 
    group_by(sampledJlab) %>% 
    summarise(
      alab = mean(as.numeric(Za_categorical)),
      blab = mean(as.numeric(Zb_categorical)),
      clab = mean(as.numeric(Zc_categorical)-1), 
      xlab = mean(X.1),
      jlab= mean(sampledJlab)
    )
  jgroup = grptbl$jlab
  agroup = grptbl$alab
  bgroup = grptbl$blab
  cgroup=grptbl$clab
  xgroup=grptbl$xlab
  
  # fit <- stan(file = "mrp2_Nq_sim.stan",
                              # data = c("Ma","Mb","Mc","J_stan","n","y","xj","agroup","bgroup","cgroup","xgroup"),
                              # iter=staniters,pars = c("propmk"),
                              # warmup = staniters-sim/2,control=list(adapt_delta=0.99, max_treedepth=13),chains=2)
  model <- stan_model(file= "mrp2_Nq_sim.stan", model_name = "mod")
  sfit <- sampling(model,
                   data = c("Ma","Mb","Mc","J_stan","n","y","xj","agroup","bgroup","cgroup","xgroup"),
                   iter=staniters,pars = c("propmk"),
                   warmup = staniters-sim/2,control=list(adapt_delta=0.99, max_treedepth=13),chains=2)
  unload.dll(model=model)
  #-- extract estimates of cell means from the stan model
  return(extract(sfit,permuted = TRUE, inc_warmup = FALSE,include = TRUE)$propmk)
  
}
