// STAN model for an integrated fit to the ahnodl model
// Goes from delta model of abundances to DiD across all species
// seen means the model fit to observed fish
// seeing means the model of the probability of observing anything
data{

  int<lower = 0> n_parameters; //number of parameters

  int<lower = 0> n_did_parameters; //number of did parameters

  int<lower = 0> n_observations_seeing; //number of observations

  int<lower = 0> n_observations_did; //number of did observations

  int<lower = 0> n_species; //number of species

  int<lower = 0> n_year_species; //number of year species combinations

  int<lower = 0> n_region_species; // number of region species combinations

  int<lower = 0> n_did; // number of did effects

  int<lower = 0> n_standard; // number of standardized matrix observations

  int<lower = 0> n_non_nested; // number of non-nested parameters

  real<lower = 0> cauchy_2 ; // cauchy parameter

  int years_per_species[n_species]; // number of years for each species

  int regions_per_species[n_species]; // number of regions for each species

  int year_species_positions[n_year_species]; // the location of year-species effects

  int region_species_positions[n_region_species]; // the location of region-species effects

  // int species_intercept_positions[n_species]; // the location of each species intercept

  int did_positions[n_did]; // the location of did effects

  int non_nested_positions[n_non_nested]; // the location of did effects

  matrix[n_observations_seeing,n_parameters] x_seeing; // covariates

  matrix[n_standard,n_parameters] standard_matrix; // covariates

  matrix[n_observations_did,n_did_parameters] x_did; // covariates

  int<lower=0, upper = 1> observed[n_observations_seeing]; // observed or not observed

} // close data section

parameters{


// seeing fish parameters

vector[n_parameters] seeing_betas; //vector of regression coefficients for seeing fish

vector<lower = 0>[n_species] seeing_sigma_year_species; // standard deviation of the seeing year effects

vector<lower = 0>[n_species] seeing_sigma_region_species; // standard deviation of the seeing region effects


} // close parameters


model {
// seen model //////////////////////////////////////////////

int counter;

vector[n_observations_seeing] prob_seen;

vector[n_standard] standardized_prob_seen;

 vector[n_year_species] seeing_year_species_betas;

  vector[n_region_species] seeing_region_species_betas;

  vector[n_non_nested] non_nested_seeing_betas;

// seeing model //////////////////////////////////////////////

seeing_year_species_betas = seeing_betas[year_species_positions];

  seeing_region_species_betas = seeing_betas[region_species_positions];

  non_nested_seeing_betas = seeing_betas[non_nested_positions];

// observed model- hierarchichal structure of year effects
counter = 1;

for (i in 1:n_species){ // hierarchical priors on year effects by species

target += normal_lpdf(segment(seeing_year_species_betas, counter,years_per_species[i]) | 0, seeing_sigma_year_species[i]);

counter = counter + years_per_species[i];

}

// regional hierarchical effects

counter = 1;

for (i in 1:n_species){ // hierarchical priors on region effects by species

target += normal_lpdf(segment(seeing_region_species_betas, counter, regions_per_species[i]) | 0, seeing_sigma_region_species[i]);

counter = counter + regions_per_species[i];

}

standardized_prob_seen = 1 ./ (1 + exp(-standard_matrix * seeing_betas));

prob_seen = 1 ./ (1 + exp(-x_seeing * seeing_betas));


// did_model //////////////////////////////////////////////


target += normal_lpdf(non_nested_seeing_betas | 0, 10);

target += cauchy_lpdf(seeing_sigma_year_species | 0, cauchy_2);
//
target += cauchy_lpdf(seeing_sigma_region_species | 0, cauchy_2);

target += bernoulli_logit_lpmf(observed |(x_seeing * seeing_betas));

} // close model block



