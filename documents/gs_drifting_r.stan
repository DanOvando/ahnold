data{
  
  int n_yr; // number of years
  
  vector <lower=0> [n_yr] index; //cpue index
  
  vector <lower=0> [n_yr] catches; // observed harvest
}

transformed data {
  
  vector [n_yr] log_index; // index log scale
  
  vector [n_yr] log_catches; // harvest log scale
  
  log_index = log(index);
  
  log_catches = log(catches);
  
}

parameters {

// modeling non stationary r
  real<lower=0> rinit; //initial r (year 1)
  
  real delta; // year to year changes in r
  
  real<lower=0> sigma_r; // stochastic productivity param
  
// other params for growth model

  vector [n_yr] log_bio_devs; // population deviations due to process error
  
  real log_k; // carrying capacity log  
  
  real<lower=0> sigma_p; // process error
  
// catchability, fishing mortality, and observation model for cpue index

  real<upper=0> log_q; // catchability log

  vector <lower=0, upper=0.8> [n_yr] u; // fishing mortality
  
  real<lower=0> sigma_o; // observation error in the cpue index
}

transformed parameters{
  
  real<lower=0> r[n_yr]; //drifting r
  
  real k;
  
  real q;
  
  vector [n_yr] bio; // estimated biomass vector
  
  vector [n_yr] index_hat;  //estimated index of cpue
  
  vector [n_yr] log_index_hat; // our estimate of the cpue index
  
  vector [n_yr] catches_hat; // estimated catches
  
  vector [n_yr] log_catches_hat; // estimated catches log
  
  r[1] = rinit;
  
  k = exp(log_k);
  
  q = exp(log_q);

  bio[1] = k;
  
  catches_hat[1] = bio[1]*u[1];
  
  //process
  
  for(t in 2:n_yr) {
    
    bio[t] = (bio[t-1]+r[t-1]*bio[t-1]*(1-bio[t-1]/k)-catches_hat[t-1])*exp(log_bio_devs[t]);
    
    catches_hat[t] = bio[t]*u[t];
    
    r[t] = r[t-1]+delta;

  }
  
  log_catches_hat = log(catches_hat);
  
  //observation
  
  index_hat = q*bio;
  
  log_index_hat = log(index_hat);
}

model {
  
  // process model
  
  rinit ~ cauchy(0,2.5);
  
  r[2:n_yr] ~ normal(r[1:n_yr-1], sigma_r);
  
  delta ~ normal(0,0.25);
  
  sigma_r ~ cauchy(0,2.5);
  
  log_k ~ uniform(log(4000),log(15000));
  
  log_bio_devs ~ normal(0,sigma_p);
  
  sigma_p ~ cauchy(0,2.5);
  
  log_catches ~ normal(log_catches_hat, 1e-2); //tune catches to match observed
  
  // observation model
  
  log_q ~ uniform(-12,0);
  
  log_index ~ normal(log_index_hat,sigma_o);
  
  sigma_o ~ cauchy(0,2.5);

}

generated quantities{ // draw from log_cpue_hat posterior based on successful parameter draws
  
  vector [n_yr] pp_log_index_hat;
  
  for (i in 1:n_yr){
    
    pp_log_index_hat[i] = normal_rng(log_index_hat[i], sigma_o);
    
  }
  
}


