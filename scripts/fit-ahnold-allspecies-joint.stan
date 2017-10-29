// STAN model for an integrated fit to the ahnodl model
// Goes from delta model of abundances to DiD across all species
// seen means the model fit to observed fish
// seeing means the model of the probability of observing anything
data{

  int<lower = 0> n_parameters; //number of parameters

  int<lower = 0> n_did_parameters; //number of did parameters

  int<lower = 0> n_observations_seen; //number of observations

  int<lower = 0> n_observations_seeing; //number of observations

  int<lower = 0> n_observations_did; //number of did observations

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

  matrix[n_observations_did,n_did_parameters] x_did; // covariates

  vector[n_observations_seen] log_density; // observed log densities

  int<lower=0, upper = 1> observed[n_observations_seeing]; // observed or not observed

} // close data section

parameters{

// seen fish parameters

vector[n_parameters] betas; // vector of regression coefficients for seen fish

vector<lower = 0>[n_species] sigma_year_species; // standard deviation of the year effects

vector<lower = 0>[n_species] sigma_region_species; // standard deviation of the year effects

real<lower = 0> sigma_density; // standard deviation of the observed densities

// seeing fish parameters

vector[n_parameters] seeing_betas; //vector of regression coefficients for seeing fish

vector<lower = 0>[n_species] seeing_sigma_year_species; // standard deviation of the seeing year effects

vector<lower = 0>[n_species] seeing_sigma_region_species; // standard deviation of the seeing region effects

// did parameters

vector[n_did_parameters] did_betas; //vector of regression coefficients for seeing fish

real<lower = 0> sigma_abundance; // standard deviation of the observed densities


} // close parameters

transformed parameters{
// Pull out hierarchichal parts of the betas
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

standardized_abundance = standardized_prob_seen .* exp(year_species_betas); // and finally the standardized abundance indicies

prob_seen = 1 ./ (1 + exp(-x_seeing * seeing_betas));


// did_model //////////////////////////////////////////////


target += normal_lpdf(log(standardized_abundance) | x_did * did_betas, sigma_abundance);

// rest of the likelihood

//
target += cauchy_lpdf( sigma_density |0, 5);
//
target += cauchy_lpdf( sigma_abundance |0, 5);

target += cauchy_lpdf(sigma_year_species | 0, 5);
//
target += cauchy_lpdf(sigma_region_species | 0, 5);

target += cauchy_lpdf(seeing_sigma_year_species | 0, 5);
//
target += cauchy_lpdf(seeing_sigma_region_species | 0, 5);

target += normal_lpdf(log_density | log_density_hat, sigma_density);

target += bernoulli_lpmf(observed | prob_seen);

}

generated quantities{

vector[n_standard] s_prob_seen;
//
vector[n_standard] standardized_abundance_index; //

vector[n_observations_did] log_abundance_hat; //

log_abundance_hat = x_did * did_betas;

s_prob_seen = 1 ./ (1 + exp(-standard_matrix * seeing_betas));
//
standardized_abundance_index = s_prob_seen .* exp(year_species_betas);
//

}


