data{

  int<lower = 0> n_parameters; //number of parameters

  int<lower = 0> n_observations_obs; //number of observations

  int<lower = 0> n_years; // number of years

  int<lower = 0> n_regions; //number of regions

  int year_positions[n_years]; // position of year terms

  int region_positions[n_regions]; // position of region terms

  int<lower = 1> vis_position; // position of visibility term

  int<lower = 1> kelp_position; // position of kelp terms

  int<lower = 1> intercept_position; // position of kelp terms

  matrix[n_observations_obs,n_parameters] x_obs; // covariates

  int<lower=0, upper = 1> observed[n_observations_obs]; // observed or not observed

} // close data section

parameters{

// observation section

real intercept_obs; // intercept for the model

real kelp_effect_obs;

real vis_effect_obs;

vector[n_years] year_effects_obs; // year effects

real<lower = 0> sigma_year_obs; // standard deviation of the year effects

vector[n_regions] region_effects_obs; // regional effects

real<lower = 0> sigma_region_obs; // standard eviation of the regional effects

} // close parameters

transformed parameters{


// print(betas)
 // variance_density = sigma_density^2;
 //
 // variance_year = sigma_year^2;
 //
 // variance_density = sigma_density^2;

//print("y=", y, " z=", z);


} // closed transformed paramters

model {

vector[n_parameters] betas_obs;

// vector[n_observations_obs] p_obs;

vector[n_observations_obs] p_raw;

 betas_obs[year_positions] = year_effects_obs;

 betas_obs[region_positions] = region_effects_obs;
 //
 betas_obs[kelp_position] = kelp_effect_obs;
 //
 betas_obs[vis_position] = vis_effect_obs;
 //
 betas_obs[intercept_position] = intercept_obs;
//
p_raw = x_obs*betas_obs;
// print(max(p_raw))
// p_obs = 1 ./ (1 + exp(-1 * (p_raw)));

// for (i in 1:n_observations_obs){
// p_obs[i] = 1/(1 + p_obs[i]);
// }

  // sigma_year_obs ~ cauchy(0, 2);
  //
  // sigma_region_obs ~ cauchy(0, 2);
  //
  sigma_year_obs ~ cauchy(0, 5);
  //
  sigma_region_obs ~ cauchy(0, 5);

year_effects_obs ~ normal(0, sigma_year_obs);
//
region_effects_obs ~ normal(0, sigma_region_obs);

observed ~ bernoulli_logit(p_raw);

}
