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
  int<lower=1> J_stan; // number of cellmeans
  int<lower=1> n; // total sample size
  int<lower = 0> y[n];      // binary outcome
  int<lower=1, upper=Ma> agroup[J_stan]; // A index at cell level
  vector<lower=0, upper=1>[J_stan] bgroup; // B index at cell level
  int<lower=1,upper=Mc> cgroup[J_stan]; // C index at cell level 
  int<lower=1, upper=Md> dgroup[J_stan]; // A index at cell level
  int<lower=1, upper=Me> egroup[J_stan]; // A index at cell level
  int<lower=1, upper=Mab> abgroup[J_stan]; //age:sex interaction label at cell level
    int<lower=1, upper=Mbc> bcgroup[J_stan];
  int<lower=1, upper=Mcd> cdgroup[J_stan];
  int<lower=1, upper=J_stan> cell_label[n];
} 

parameters{
  vector[2] beta; //int, sex
  vector[Ma] az;
  vector[Mc] cz;
  vector[Md] dz;
  vector[Me] ez;
  // vector[Mab] abz;
  // vector[Mbc] bcz;
  // vector[Mcd] cdz;
  real<lower=0> sigma_a; //sd of A main effect
  real<lower=0> sigma_c; //sd of C main effect
  real<lower=0> sigma_d; //sd of D main effect
  real<lower=0> sigma_e; //sd of E main effect
  // real<lower=0> sigma_ab; //sd of age:sex interaction
  // real<lower=0> sigma_bc;
  // real<lower=0> sigma_cd;
}

transformed parameters{
  vector[Ma] alphaa = az*sigma_a;
  vector[Mc] alphac= cz*sigma_c; 
  vector[Md] alphad= dz*sigma_d; 
  vector[Me] alphae= ez*sigma_e;
  // vector[Mab] alphaab= abz*sigma_ab; 
  // vector[Mbc] alphabc= bcz*sigma_bc; 
  // vector[Mcd] alphacd= cdz*sigma_cd; 
  vector[J_stan] cellmean =  inv_logit(beta[1] + beta[2]*bgroup + alphaa[agroup] + alphac[cgroup] + alphad[dgroup]+ alphae[egroup] );
}
 
model{
  az ~ std_normal();
  cz~ std_normal();
  dz~ std_normal();
  ez~ std_normal();
  // abz ~ std_normal();
  // bcz ~ std_normal();
  // cdz ~ std_normal();
  sigma_a ~ cauchy(0,1);
  sigma_c ~ cauchy(0,1);
  sigma_d ~ cauchy(0,1);
  sigma_e ~ cauchy(0,1);
  // sigma_ab ~ cauchy(0,1);
  // sigma_bc ~ cauchy(0,1);
  // sigma_cd ~ cauchy(0,1);
  // y ~ bernoulli(p);
    y ~ bernoulli(cellmean[cell_label]);
}

generated quantities {
  // real<lower=0,upper=1> theta_rep;
  vector[n] p;
  int y_sim[n];
  vector[n] log_lik;
  // use current estimate of theta to generate new sample
   for(l in 1:n){
    p[l] = cellmean[cell_label[l]];
  }
  for (k in 1:n)
    y_sim[k] = bernoulli_rng(cellmean[cell_label[k]]);
  // estimate theta_rep from new sample
    for(g in 1:n){
    log_lik[g] = bernoulli_lpmf(y[g]|p[g]);
  }
  // theta_rep = sum(y_sim) * 1.0 / n;
}
