data {
  int ndata;
  int nbeta; // number of betas that are not time related
  int ngamma; // number if unique fish
  int nyear;
  matrix[ndata, nbeta+nyear] x; // these are the intercept, age and year effects
  matrix[ndata, ngamma] z; // these are the  fish effects
  real y[ndata];


}

parameters {
  vector[nbeta] beta; //the regression parameters depth and length
  vector[ngamma] gammaraw; //the regression parameters on fishID
  vector[nyear] year_epsilon_raw;

  real <lower = -10, upper = 100 > logit_rho; // autocorrelation term in logit space
  real <lower = 0 > sigma_eta; // standard deviation of time series
  real <lower = 0 > sigma; //the standard deviation of the increments
  real <lower = 0 > sigma_fish;

}

transformed parameters {

  vector[nyear] year_effect;
  vector[nbeta+nyear] newbeta; //
  vector[nyear] yearbeta;
  vector[ndata] yhat;
  vector[ngamma] gamma;
  real sigma_year;
  vector[nyear] yearbeta_raw;
  real rho;
  
  rho = 1/(1+exp(-logit_rho));


  
  sigma_year = sigma_eta*(1-rho^2)^(0.5);
// non centered random effects of year
   year_effect = year_epsilon_raw * sigma_year;
   

   // create year vector for easier post hoc evaluation
// now create annual coefficients using random walk. 
   yearbeta_raw[1] = year_effect[1];

   for (i in 2:nyear) {
     yearbeta_raw[i] = rho * yearbeta_raw[i - 1] + year_effect[i];
   }
   
   // center so that mean is zero
   for (i in 1:nyear) {
     yearbeta[i] = yearbeta_raw[i]  - mean(yearbeta_raw);
   }

// create vector newebata that has the intercept and depth terms plus year terms
newbeta[1:nbeta]= beta;

newbeta[(nbeta+1):(nbeta+nyear)] = yearbeta[1:nyear];



  // non centered random effects.  have mean = 0
  gamma = gammaraw * sigma_fish;
  


//  linear predictors
  yhat = x * newbeta + z * gamma;
  
}

model {

// priors on variance terms
sigma ~ cauchy(0,1.0);
sigma_eta ~ cauchy(0, 1.0);
sigma_fish ~ cauchy(0, 1.0);

// prior on autocorrelation
logit_rho ~ normal(0, 1.75);

// loop through betas to assign priors
for (i in 1:nbeta){
  beta[i] ~ cauchy(0,2.5);//prior for the slopes following Gelman 2008
}

for (i in 1:nyear){
  year_epsilon_raw[i] ~ normal(0, 1);
}

// loop through basin gammas and cauchy_scales
for (i in 1:ngamma){
  gammaraw[i] ~ normal(0,1);// random effects for basin
}

// End of Prior Specification

  // calculate likelihood 

  for (i in 1:ndata) y[i] ~ normal(yhat[i], sigma);

}

generated quantities {
  
  vector[ndata] y_pred_check;
  for(i in 1:ndata) y_pred_check[i] = normal_rng(yhat[i], sigma);
  
}
