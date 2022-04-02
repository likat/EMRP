// Estimation of Nj for two-stage MRP
// modeling X~Z: ("age3", "sex2",  "race3","educat4", "povgap4")

data {
  int<lower=1> Ma; // number of (age) levels
  int<lower=1> Mb; // number of (sex) levels
  int<lower=1> Mc; // number of (race) levels
  int<lower=1> Md; // number of (educat) levels
  int<lower=1> Me; // number of (povgap) levels
  int<lower=1> L; // number of (svefreq) levels
  int<lower=1> Mab; // number of age:sex interaction groups (4x2=8)
  int<lower=1> Mbc;
  int<lower=1> Mcd;
  int<lower=1> Mdx; // number of education:visit interaction groups (4x2=8)
  int<lower=0> J_stan; // number of XxZ categories
  int<lower=0> n; // total sample size
  int<lower=0, upper=1> y[n];      // svefreq(binary)
  int<lower=1> xj[n]; //jgroup index at individual level
  int<lower=1, upper=Ma> agroup[J_stan]; // A index at cell level
  vector<lower=0, upper=1>[J_stan] bgroup; // B index at cell level
  int<lower=1, upper=Mc> cgroup[J_stan]; // C index at cell level 
  int<lower=1, upper=Md> dgroup[J_stan]; // D index at cell level
  int<lower=1, upper=Me> egroup[J_stan]; // E index at cell level
  int<lower=0, upper=1> xgroup[J_stan]; // X2 index at cell level 
  int<lower=1, upper=Mab> abgroup[J_stan]; //age:sex interaction label at cell level
  int<lower=1, upper=Mbc> bcgroup[J_stan];
  int<lower=1, upper=Mcd> cdgroup[J_stan];
  int<lower=1, upper=Mdx> dxgroup[J_stan]; //education:visit interaction label at cell level
}

parameters {
  vector[2] beta; //global intercept, sex effect
  vector[Ma] az;
  vector[Mc] cz;
  vector[Md] dz;
  vector[Me] ez;
  real<lower=0> sigma_a; //sd of A main effect
  real<lower=0> sigma_c; //sd of C main effect
  real<lower=0> sigma_d; //sd of D main effect
  real<lower=0> sigma_e; //sd of E main effect
}

transformed parameters{
  vector[Ma] alphaa = az*sigma_a;
  vector[Mc] alphac= cz*sigma_c; //B main effect
  vector[Md] alphad= dz*sigma_d; //B main effect
  vector[Me] alphae= ez*sigma_e;
 vector[J_stan] pxz = beta[1] + beta[2]*bgroup + alphaa[agroup] + alphac[cgroup] + alphad[dgroup]+ alphae[egroup]; // p(X=1|Z)
}

model {
  az ~ std_normal();
  cz~ std_normal();
  dz~ std_normal();
  ez~ std_normal();
  sigma_a ~ cauchy(0,1);
  sigma_c ~ cauchy(0,1);
  sigma_d ~ cauchy(0,1);
  sigma_e ~ cauchy(0,1);
  y ~ bernoulli_logit(pxz[xj]);
}

generated quantities {
  vector<lower=0>[J_stan] propmk;
  for(i in 1:J_stan){
    if(xgroup[i] == 1){
      propmk[i] = inv_logit(pxz[i]);
    }
    else{
      propmk[i] = 1-inv_logit(pxz[i]);
    }
  }
}
