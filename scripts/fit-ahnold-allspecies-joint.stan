// STAN model for an integrated fit to the ahnodl model
// Goes from delta model of abundances to DiD across all species
// seen means the model fit to observed fish
// seeing means the model of the probability of observing anything
data{

  int<lower = 0> n_parameters; //number of parameters

  int<lower = 0> n_observations_seen; //number of observations

  int<lower = 0> n_observations_seeing; //number of observations

  int<lower = 0> n_species; //number of observations

  int<lower = 0> n_year_species; //number of year species combinations

  int<lower = 0> n_region_species; // number of region species combinations

  int<lower = 0> n_standard; // number of region species combinations

  int years_per_species[n_species]; // number of years for each species

  int regions_per_species[n_species]; // number of years for each species

  int year_species_positions[n_year_species]; // the location of year-species effects

  int region_species_positions[n_region_species]; // the location of region-species effects

  int species_intercepts_positions[n_species]; // the location of region-species effects

  matrix[n_observations_seen,n_parameters] x_seen; // covariates

  matrix[n_observations_seeing,n_parameters] x_seeing; // covariates

  matrix[n_standard,n_parameters] standard_matrix; // covariates

  vector[n_observations_seen] log_density; // observed log densities

  int<lower=0, upper = 1> observed[n_observations_seeing]; // observed or not observed

} // close data section

parameters{

// seen fish parameters

vector[n_parameters] betas;

vector<lower = 0>[n_species] sigma_year_species; // standard deviation of the year effects

vector<lower = 0>[n_species] sigma_region_species; // standard deviation of the year effects

real<lower = 0> sigma_density; // standard deviation of the observed densities

// seeing fish parameters

vector[n_parameters] seeing_betas;

vector<lower = 0>[n_species] seeing_sigma_year_species; // standard deviation of the year effects

vector<lower = 0>[n_species] seeing_sigma_region_species; // standard deviation of the year effects


} // close parameters

transformed parameters{

// seen fish

  vector[n_year_species] year_species_betas;

  vector[n_region_species] region_species_betas;

  vector[n_species] species_intercepts;

  vector[n_year_species] seeing_year_species_betas;

  vector[n_region_species] seeing_region_species_betas;

  vector[n_species] seeing_species_intercepts;

  year_species_betas = betas[year_species_positions];

  region_species_betas = betas[region_species_positions];

  species_intercepts = betas[species_intercepts_positions];

//seeing fish


  seeing_year_species_betas = seeing_betas[year_species_positions];

  seeing_region_species_betas = seeing_betas[region_species_positions];

  seeing_species_intercepts = seeing_betas[species_intercepts_positions];


}

model {
// seen model //////////////////////////////////////////////

int counter;

vector[n_observations_seen] log_density_hat;

vector[n_observations_seeing] prob_seen;

vector[n_standard] standardized_prob_seen;
//
vector[n_standard] standardized_abundance; // year effects


// observed model- hierarchichal structure of year effects
counter = 1;


for (i in 1:n_species){ // hierarchical priors on year effects by species

target += normal_lpdf(year_species_betas[counter:(counter + years_per_species[i] - 1)] | species_intercepts[i], sigma_year_species[i]);
// target += normal_lpdf(year_species_betas[counter:(counter + years_per_species[i] - 1)] | 0, sigma_year_species[i]);

counter = counter + years_per_species[i];

}

// regional hierarchical effects

counter = 1;

for (i in 1:n_species){ // hierarchical priors on region effects by species

target += normal_lpdf(region_species_betas[counter:(counter + regions_per_species[i] - 1)] | species_intercepts[i], sigma_region_species[i]);

counter = counter + regions_per_species[i];

}

log_density_hat = x_seen*betas;


// seeing model //////////////////////////////////////////////


// observed model- hierarchichal structure of year effects
counter = 1;


for (i in 1:n_species){ // hierarchical priors on year effects by species

target += normal_lpdf(seeing_year_species_betas[counter:(counter + years_per_species[i] - 1)] | seeing_species_intercepts[i], seeing_sigma_year_species[i]);
// target += normal_lpdf(year_species_betas[counter:(counter + years_per_species[i] - 1)] | 0, sigma_year_species[i]);

counter = counter + years_per_species[i];

}

// regional hierarchical effects

counter = 1;

for (i in 1:n_species){ // hierarchical priors on region effects by species

target += normal_lpdf(seeing_region_species_betas[counter:(counter + regions_per_species[i] - 1)] | seeing_species_intercepts[i], seeing_sigma_region_species[i]);

counter = counter + regions_per_species[i];

}


standardized_prob_seen = 1 ./ (1 + exp(-standard_matrix * seeing_betas));

standardized_abundance = standardized_prob_seen .* exp(year_species_betas);

prob_seen = 1 ./ (1 + exp(-x_seeing * seeing_betas));

// rest of the likelihood

//
target += cauchy_lpdf( sigma_density |0, 5);
//
target += cauchy_lpdf(sigma_year_species | 0, 5);
//
target += normal_lpdf(log_density | log_density_hat, sigma_density);

target += bernoulli_lpmf(observed | prob_seen);


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

generated quantities{

vector[n_standard] s_prob_seen;
//
vector[n_standard] standardized_abundance_index; //

s_prob_seen = 1 ./ (1 + exp(-standard_matrix * seeing_betas));
//
standardized_abundance_index = s_prob_seen .* exp(year_species_betas);
//

}


