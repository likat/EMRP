data{ 
  int<lower=1> Ma; // number of (age) levels
  int<lower=1> Mb; // number of (sex) levels
  int<lower=1> Mc; // number of (race) levels
  int<lower=1> Md; // number of (educat) levels
  int<lower=1> Me; // number of (income) levels
  int<lower=1> L; // number of (svefreq) levels
  int<lower=1> Mab; // number of age:sex interaction groups (4x2=8)
  int<lower=1> Mbc;
  int<lower=1> Mcd;
  int<lower=1> Mdx; // number of education:visit interaction groups (4x2=8)
  int<lower=1> J_stan; // number of cellmeans
  int<lower=1> n; // total sample size
  int<lower = 0> y[n];      // binary outcome
  int<lower=1, upper=Ma> agroup[J_stan]; // A index at cell level
  vector<lower=0, upper=1>[J_stan] bgroup; // B index at cell level
  int<lower=1,upper=Mc> cgroup[J_stan]; // C index at cell level 
  int<lower=1, upper=Md> dgroup[J_stan]; // A index at cell level
  int<lower=1, upper=Me> egroup[J_stan]; // A index at cell level
 vector<lower=0, upper=1>[J_stan] xgroup; // X grouping at cell level
  int<lower=1, upper=Mab> abgroup[J_stan]; //age:sex interaction label at cell level
  int<lower=1, upper=Mbc> bcgroup[J_stan];
  int<lower=1, upper=Mcd> cdgroup[J_stan];
  int<lower=1, upper=Mdx> dxgroup[J_stan]; //education:visit interaction label at cell level
  int<lower=1, upper=J_stan> cell_label[n];
} 

parameters{
  // vector[3] beta; //int, sex,svefreq
    real b0;
  real bb;
  real bx;
  vector[Ma] az;
  vector[Mc] cz;
  vector[Md] dz;
  vector[Me] ez;
  vector[Mdx] dxz;
  real<lower=0> sigma_a; //sd of A main effect
  real<lower=0> sigma_c; //sd of C main effect
  real<lower=0> sigma_d; //sd of D main effect
  real<lower=0> sigma_e; //sd of E main effect
  real<lower=0> sigma_dx; //sd of education:visit interaction
  real<lower=0> sigma_beta0; 
  real<lower=0> sigma_betab; 
  real<lower=0> sigma_betax;
}

transformed parameters{
  vector[Ma] alphaa = az*sigma_a;
  vector[Mc] alphac= cz*sigma_c; 
  vector[Md] alphad= dz*sigma_d; 
  vector[Me] alphae= ez*sigma_e;
  vector[Mdx] alphadx= dxz*sigma_dx;
    real beta0 = b0*sigma_beta0;
  real betab = bb*sigma_betab;
  real betax = bx*sigma_betax;
 // vector[J_stan] cellmean =  inv_logit(beta[1] + beta[2]*bgroup + beta[3]*xgroup+ alphaa[agroup] + alphac[cgroup] + alphad[dgroup]+ alphae[egroup] + alphadx[dxgroup]);
 vector[J_stan] cellmean =  inv_logit(beta0 + betab*bgroup + betax*xgroup+ alphaa[agroup] + alphac[cgroup] + alphad[dgroup]+ alphae[egroup] + alphadx[dxgroup]);
}

 
model{
    b0 ~ std_normal();
  bb ~ std_normal();
  bx ~ std_normal();
  az ~ std_normal();
  cz~ std_normal();
  dz~ std_normal();
  ez~ std_normal();
  dxz ~ std_normal();
  sigma_a ~ cauchy(0,1);
  sigma_c ~ cauchy(0,1);
  sigma_d ~ cauchy(0,1);
  sigma_e ~ cauchy(0,1);
  sigma_dx ~ cauchy(0,1);
    sigma_beta0 ~ cauchy(0,5);
  sigma_betab ~ cauchy(0,5);
  sigma_betax ~ cauchy(0,5);
  y ~ bernoulli(cellmean[cell_label]);
}

generated quantities {
  vector[n] p;
  int y_sim[n];
  vector[n] log_lik;
  for (j in 1:n){
    y_sim[j] = bernoulli_rng(cellmean[cell_label[j]]);
    }
  for(l in 1:n){
    p[l] = cellmean[cell_label[l]];
  }
  for(g in 1:n){
    log_lik[g] = bernoulli_lpmf(y[g]|p[g]);
  }
}
