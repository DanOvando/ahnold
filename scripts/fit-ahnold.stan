data{

  int<lower = 0> n_parameters; //number of parameters

  int<lower = 0> n_observations; //number of observations

  int<lower = 0> n_years; // number of years

  int<lower = 0> n_regions; //number of regions

  int year_positions[n_years]; // position of year terms

  int region_positions[n_regions]; // position of region terms

  int<lower = 1> vis_position; // position of visibility term

  int<lower = 1> kelp_position; // position of kelp terms

  int<lower = 1> intercept_position; // position of kelp terms

  matrix[n_observations,n_parameters] x; // covariates

  vector[n_observations] log_density; // observed log densities


} // close data section

parameters{

real intercept; // intercept for the model

real kelp_effect;

real vis_effect;

vector[n_years] year_effects; // year effects

real<lower = 0> sigma_year; // standard deviation of the year effects

vector[n_regions] region_effects; // regional effects

real<lower = 0> sigma_region; // standard eviation of the regional effects

real<lower = 0> sigma_density; // standard deviation of the observed densities

} // close parameters

transformed parameters{

vector[n_parameters] betas;

real<lower = 0> variance_density;
//
real<lower = 0> variance_year;

// real[n_years]  test;

 betas[year_positions] = year_effects;

 betas[region_positions] = region_effects;
 //
 betas[kelp_position] = kelp_effect;
 //
 betas[vis_position] = vis_effect;
 //
 betas[intercept_position] = intercept;
 //
// print(betas)
 variance_density = sigma_density^2;
 //
 variance_year = sigma_year^2;
 //
 variance_density = sigma_density^2;

//print("y=", y, " z=", z);


} // closed transformed paramters

model {

vector[n_observations] log_density_hat;

log_density_hat = x*betas;

  sigma_density ~ cauchy(0, 5);

  sigma_year ~ cauchy(0, 5);

  sigma_region ~ cauchy(0, 5);

  // variance_density ~ inv_gamma(2,2);
//
// variance_year ~ inv_gamma(alpha = 2, beta = 2);
//
// variance_region ~ inv_gamma(alpha = 2, beta = 2);
//
year_effects ~ normal(0, sigma_year);

region_effects ~ normal(0, sigma_region);

log_density ~ normal(log_density_hat, sigma_density);

}
