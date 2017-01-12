/*
Bayesian Hierarchichal Regression for Project Ahnold

This script runs the STAN regression model
*/

data{
//Data section
int<lower = 1> num_betas;

int<lower = 1> num_sigmas;

int<lower = 1> num_pars;

int<lower = 1> num_obs;

matrix[num_obs, num_pars] X; //covariates

vector[num_obs] y;

}

parameters{
// Parameters to be estimated
vector[num_pars] betas;

vector<lower = 0>[num_sigmas] sigmas;

}

model{

betas ~ normal(0,5);

sigma ~ cauchy(0,2.5);

y ~ normal(X*betas,sigma);

}