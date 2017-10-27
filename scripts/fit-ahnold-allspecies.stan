data{

  int<lower = 0> n_parameters; //number of parameters

  int<lower = 0> n_observations; //number of observations

  int<lower = 0> n_year_species; //number of observations

  int<lower = 0> n_species; //number of observations

  int years_per_species[n_species]; // number of years for each species

  int year_species_positions[n_year_species];

  matrix[n_observations,n_parameters] x; // covariates

  vector[n_observations] log_density; // observed log densities

} // close data section

parameters{

// observed parameters

vector[n_parameters] betas;

vector<lower = 0>[n_species] sigma_year_species; // standard deviation of the year effects

real<lower = 0> sigma_density; // standard deviation of the observed densities


} // close parameters

transformed parameters{

  vector[n_year_species] year_species_betas;

  year_species_betas = betas[year_species_positions];

}

model {

int counter;

vector[n_observations] log_density_hat;

// observed model- hierarchichal structure of year effects
counter = 1;


for (i in 1:n_species){ // hierarchical priors on year effects by species

target += normal_lpdf(year_species_betas[counter:(counter + years_per_species[i] - 1)] | 0, sigma_year_species[i]);

counter = counter + years_per_species[i];

}

log_density_hat = x*betas;
//
target += cauchy_lpdf( sigma_density |0, 5);
//
target += cauchy_lpdf(sigma_year_species | 0, 2.5);
//
target += normal_lpdf(log_density | log_density_hat, sigma_density);

 // target += bernoulli_lpmf(observed | p_obs);


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


