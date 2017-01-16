/*
Bayesian Hierarchichal Regression for Project Ahnold

This script runs the STAN regression model
*/

data{
//Data section
/*
int<lower = 1> num_betas;

int<lower = 1> num_sigmas;
*/

int<lower = 1> num_pars;

int<lower = 1> num_obs;

matrix[num_obs, num_pars] x; //covariates

vector[num_obs] y;

}

parameters{
// Parameters to be estimated
vector[num_pars] betas;

real<lower = 0> sigma;


}

transformed parameters{

  /*
  real year_beta;

  year_beta = betas[year_beta_index];

  This could be a good place for creating new variables that are functions
  of the original thing, so for example the subset of betas that correspond to years


  year_beta can then be used in a likelihood in the model
  section

  */


}

model{

//print("betasize = ",num_elements(betas));

betas ~ normal(0,5);

sigma ~ cauchy(0,2.5);

y ~ normal(x*betas,sigma);

// year_beta ~ normal(0, sigma_year);

}
