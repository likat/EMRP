// part 1 of mrp^2: estimation of cell counts
// effectively modeling X ~ Z
// X[i] ~ bernoulli(Xbeta + alpha_x1)
data {
  int<lower=0> Ma; // number of Z_a groups (5 levels)
  int<lower=0> Mb; // number of Z_b groups (5 levels)
  int<lower=0> Mc; // number of Z_c groups (2 levels: 0,1)
  int<lower=0> J_stan; // number of ZxX categories
  int<lower=0> n; // total sample size
  int<lower=0, upper=1> y[n];      // X (binary)
  int<lower=1> xj[n]; //jgroup index at individual level
  int<lower=1, upper=Ma> agroup[J_stan]; // A index at cell level
  int<lower=1, upper=Mb> bgroup[J_stan]; // B index at cell level
  vector<lower=0, upper=1>[J_stan] cgroup; // C index at cell level 
   vector<lower=0, upper=1>[J_stan] xgroup; // X index at cell level 
}
parameters {
  vector[2] beta; //global intercept, C
  vector[Ma] az;
  vector[Mb] bz;
  real<lower=0> sigma_a; //sd of A main effect
  real<lower=0> sigma_b; //sd of B main effect
}

transformed parameters{
  vector[Ma] alphaa = az*sigma_a;
  vector[Mb] alphab= bz*sigma_b; //B main effect
  vector[J_stan] pxz = beta[1] + beta[2]*cgroup + alphaa[agroup] + alphab[bgroup]; // p(X=1|Z)
}

model {
  az ~ normal(0,1);
  bz~ normal(0,1);
  sigma_a ~ cauchy(0,1);
  sigma_b ~ cauchy(0,1);
  y ~ bernoulli_logit(pxz[xj]);
}

generated quantities {
  vector<lower=0, upper=1>[J_stan] propmk;
  for(i in 1:J_stan){
    if(xgroup[i] == 1){
      propmk[i] = inv_logit(pxz[i]);
    }
    else{
      propmk[i] = 1-inv_logit(pxz[i]);
    }
  }
}
