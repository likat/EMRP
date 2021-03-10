data{ 
  int<lower=0> Ma; // number of Z_a groups (5 levels)
  int<lower=0> Mb; // number of Z_b groups (5 levels)
  int<lower=0> Mc; // number of Z_c groups (2 levels: 0,1)
  int<lower=0> L; // number of X categories (2 levels: 0,1)
  int<lower=0> Mbc; // number of B:C interaction groups (5x2=10)
  int<lower=0> Mcx; // number of C:X groups (2x2 = 4)
  int<lower=0> J_stan; // number of cellmeans (5x5x2x2=100)
  int<lower=0> n; // total sample size
  real y[n];      // continuous outcome
  int<lower=1, upper=Ma> agroup[J_stan]; // A index at cell level
  int<lower=1, upper=Mb> bgroup[J_stan]; // B index at cell level
  vector<lower=0, upper=1>[J_stan] cgroup; // C index at cell level 
  vector<lower=0, upper=1>[J_stan] xgroup; // X grouping at cell level
  int<lower=1, upper=Mbc> bcgroup[J_stan]; //B:C index at cell level
  vector<lower=0, upper=1>[J_stan]cxgroup; //C:X index at cell level
  int<lower=1, upper=J_stan> cell_label[n];
} 

parameters{
  // vector[3] beta; //  X THIS WORKS
  vector[4] beta; 
  vector[Ma] az;
  vector[Mb] bz;
  real<lower=0> sigmay; // sd of the observations
  real<lower=0> sigma_a; //sd of A main effect
  real<lower=0> sigma_b; //sd of B main effect
}

transformed parameters{
  vector[Ma] alphaa = az*sigma_a;
  vector[Mb] alphab= bz*sigma_b; //B main effect
vector[J_stan] cellmean =  beta[1] + beta[2]*cgroup + beta[3]*xgroup + beta[4]*cxgroup + alphaa[agroup] + alphab[bgroup]; // THIS WORKS
}
 
model{
  az ~ normal(0,5);
  bz~ normal(0,5);
  sigma_a ~ normal(0,1);
  sigma_b ~ normal(0,1);
  sigmay ~ normal(0,1);
  y~ normal(cellmean[cell_label], sigmay);
}
