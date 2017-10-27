data{

  int<lower = 0> n_parameters; //number of parameters

  int<lower = 0> n_observations; //number of observations

  int<lower = 0> n_observations_obs; //number of observations

  int<lower = 0> n_years; // number of years

  int<lower = 0> n_regions; //number of regions

  int year_positions[n_years]; // position of year terms

  int region_positions[n_regions]; // position of region terms

  int<lower = 1> vis_position; // position of visibility term

  int<lower = 1> kelp_position; // position of kelp terms

  int<lower = 1> intercept_position; // position of kelp terms

  matrix[n_observations,n_parameters] x; // covariates

  matrix[n_observations_obs,n_parameters] x_obs; // covariates

  vector[n_observations] log_density; // observed log densities

  int<lower=0, upper = 1> observed[n_observations_obs]; // observed or not observed


} // close data section

parameters{

// observed parameters

real intercept; // intercept for the model

real kelp_effect;

real vis_effect;

vector[n_years] year_effects; // year effects

real<lower = 0> sigma_year; // standard deviation of the year effects

vector[n_regions] region_effects; // regional effects

real<lower = 0> sigma_region; // standard eviation of the regional effects

real<lower = 0> sigma_density; // standard deviation of the observed densities

// observation parameters

real intercept_obs; // intercept for the model

real kelp_effect_obs;

real vis_effect_obs;

vector[n_years] year_effects_obs; // year effects

real<lower = 0> sigma_year_obs; // standard deviation of the year effects

vector[n_regions] region_effects_obs; // regional effects

real<lower = 0> sigma_region_obs; // standard eviation of the regional effects


} // close parameters

transformed parameters{


// real<lower = 0> variance_density;
//
// real<lower = 0> variance_year;

// real[n_years]  test;


 //
// print(betas)
 // variance_density = sigma_density^2;
 //
 // variance_year = sigma_year^2;
 //
 // variance_density = sigma_density^2;

//print("y=", y, " z=", z);


} // closed transformed paramters

model {

vector[n_parameters] betas;

vector[n_observations] log_density_hat;

vector[n_parameters] betas_obs;

vector[n_observations_obs] p_obs;

vector[n_observations_obs] p_raw;

vector[n_years] standardized_prob_seen;

vector[n_years] standardized_abundance; // year effects

matrix[n_years, n_parameters] standard_x_obs;


// observed model

betas[year_positions] = year_effects;

 betas[region_positions] = region_effects;
 //
 betas[kelp_position] = kelp_effect;
 //
 betas[vis_position] = vis_effect;
 //
 betas[intercept_position] = intercept;

// observation model

 betas_obs[year_positions] = year_effects_obs;

 betas_obs[region_positions] = region_effects_obs;
 //
 betas_obs[kelp_position] = kelp_effect_obs;
 //
 betas_obs[vis_position] = vis_effect_obs;
 //
 betas_obs[intercept_position] = intercept_obs;

 // deal with observation model
p_raw = x_obs*betas_obs;

p_obs = 1 ./ (1 + exp(-1 * (p_raw)));


for (i in 1:n_years){

 for (j in 1:n_parameters){

   standard_x_obs[i,j] = 0;
 }

}


for (i in 1:n_years){

  standard_x_obs[i,intercept_position] = 1;

  standard_x_obs[i,kelp_position] = mean(col(x_obs, kelp_position));

  standard_x_obs[i,vis_position] = mean(col(x_obs, vis_position));

  standard_x_obs[i,year_positions[i]] = 1;
  //
  // standard_x_obs[i,region_positions] = 1;

  // standard_x_obs[i,region_positions[1]] = 1;

}

standardized_prob_seen = 1 ./ (1 + exp(-standard_x_obs * betas_obs));

standardized_abundance = standardized_prob_seen .* exp(year_effects);

// print(standardized_abundance)

// deal with observed data
log_density_hat = x*betas;

target += cauchy_lpdf( sigma_density |0, 5);

target += cauchy_lpdf(sigma_year | 0, 2.5);
//
target += cauchy_lpdf(sigma_region | 0, 2.5);
//
 target += cauchy_lpdf(sigma_year_obs | 0,5);
//
 target += cauchy_lpdf(sigma_region_obs | 0,5);
//
target += normal_lpdf(year_effects | 0, sigma_year);
//
target += normal_lpdf(region_effects | 0, sigma_region);
//
target += normal_lpdf(year_effects_obs | 0, sigma_year_obs);
//
target += normal_lpdf(region_effects_obs | 0, sigma_region_obs);
//
target += normal_lpdf(log_density | log_density_hat, sigma_density);

target += bernoulli_lpmf(observed | p_obs);


//   sigma_density ~ cauchy(0, 5);
//
//   sigma_year ~ cauchy(0, 2.5);
//
//   sigma_region ~ cauchy(0, 2.5);
//
//   sigma_year_obs ~ cauchy(0, 5);
//
//   sigma_region_obs ~ cauchy(0, 5);
//
// year_effects ~ normal(0, sigma_year);
//
// region_effects ~ normal(0, sigma_region);
//
// year_effects_obs ~ normal(0, sigma_year_obs);
//
// region_effects_obs ~ normal(0, sigma_region_obs);
//
// log_density ~ normal(log_density_hat, sigma_density);

// observed ~ bernoulli_logit(p_raw);

// observed ~ bernoulli(p_obs);
// I think you can use the target += notation to make the likelihood out of an arbitrary number of betas and sigmas, so long as the type of data is known. i.e. loop over the sigmas and betas, adding the the likelihood as needed

}


