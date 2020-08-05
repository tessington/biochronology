// this version has complicated indexing so that initial size can be estimated (first assumed identical across fish).  It also allows individual variation in q separate from that expected based on k_i

data {
  int ndata;
  int nfish;
  int nyears;
  int maxfirstage;
  int minfirstage;
  real winc[ndata];
  int years[ndata];
  int startyears[nfish]; // these are years indices for birthyear
  int fishids[ndata];
  int startage[nfish];
  real cinc[ndata];
  int firstyear[nfish];
  //real w_start; // send the starting size, for testing where w_start is not estimated
  
}


parameters {
real<lower = 0> k_base;
real<lower = 0> q_base;  
real<lower = 0, upper = 1> rho;
real<upper = 1.0 > psi;
real <lower = 0> k_sigma;
real <lower = 0> beta_t;
vector[nyears] omega;
real<lower = 0> winf_sigma;
vector[nfish] eps_k_raw;
vector[nfish] eps_qi_raw;
real<lower = 0> inc_sigma;
// need to get the start size for average growth 
real<lower = 0, upper = 0.5> w_start_p; // w_start as a proportion of winf
real <lower = 0> qi_sigma;
}

transformed parameters {
  vector[nyears] eps_q;
  vector[nyears] eps_q_raw;
  vector[nyears] winf_t;
  vector[nfish] eps_k;
  vector[nfish] eps_qi;
  real w_start;
  real winf_base;
  real sigma_year;
  
  
   winf_base = q_base / k_base;
 
  
  eps_q_raw[1] = omega[1];
  for (t in 2:nyears) eps_q_raw[t] = rho * eps_q_raw[t-1] + sqrt(1- rho^2) * omega[t];

  for (t in 1:nyears) eps_q[t] = eps_q_raw[t] - mean(eps_q_raw); //center on 0
  
  

  winf_t = q_base * exp(beta_t * eps_q) / k_base ;
  
  for (i in 1:nfish) eps_k[i] = (eps_k_raw[i] - mean(eps_k_raw)) * k_sigma; // center on 0 and scale by k sigma
  
  for (i in 1:nfish) eps_qi[i] = (eps_qi_raw[i] - mean(eps_qi_raw)) * qi_sigma; // center on 0 and scale by qi_sigma
  // w_start is a proportion of w_infinity
  w_start = q_base / k_base * w_start_p;
  
}

model {
  real inc_hat[ndata];
  real k;
  real winf;
  real winf_tmp;
  real k_tmp;
  real wt[ndata];
  real firstsize[maxfirstage,nfish];
  
 // need to center the eps_q.  create eps_q raw and go from there
 
  // specify priors
  k_base ~ cauchy(0, 0.3);
  q_base ~ cauchy(0.5, 2);
  winf_sigma ~ cauchy(0, 0.5);
  beta_t ~ cauchy(0, 0.25);
  k_sigma ~ cauchy(0,0.1);
  qi_sigma ~ cauchy(0,0.1);
  psi ~ beta(2.5, 2.5);
  inc_sigma ~ cauchy(0,0.5);
  w_start_p ~ beta(1,1);
  
  
  for (t in 1:nyears) omega[t] ~ normal(0, 1);
  for (i in 1:nfish) eps_k_raw[i] ~ normal(0,1);
  for (i in 1:nfish) eps_qi_raw[i] ~ normal(0,1);
  // hard indexing part. 
 for (i in 1:nfish) {
  k_tmp = k_base* exp(eps_k[i]);
  firstsize[minfirstage,i] = w_start;
  if (startage[i] > minfirstage) {  
    for (j in minfirstage:(startage[i]-1)) {
      winf_tmp = q_base * exp(psi * eps_k[fishids[i]]+beta_t * eps_q[firstyear[i]+j] + eps_qi[fishids[i]]) / k_tmp;
      firstsize[j + 1,i] = firstsize[j,i] * exp(-k_tmp) + winf_tmp * (1 - exp(-k_tmp));
    }
  }
}

// now make a loop to get wt based on firstsize and cumulative increment, use this for wt.

for (i in 1: ndata) wt[i] = firstsize[startage[fishids[i]],fishids[i]] + cinc[i];
  // get likelihood
  
  for (i in 1:ndata) {
    k = k_base* exp(eps_k[fishids[i]]);
    winf = q_base * exp(psi * eps_k[fishids[i]]+beta_t * eps_q[years[i]] + eps_qi[fishids[i]]) / k;
    inc_hat[i] =  (1 - exp(- k)) * (winf - wt[i]);
    winc[i] ~ normal(inc_hat[i], inc_sigma);
  }
  

}
