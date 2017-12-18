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

  int<lower = 0> n_clusters; //number of geographic clusters

  int<lower = 0> n_year_species; //number of year species combinations

  int<lower = 0> n_region_clusters; // number of region species combinations

  int<lower = 0> n_did; // number of did effects

  int<lower = 0> n_non_nested; // number of non-nested parameters

  real<lower = 0> cauchy_2 ; // cauchy parameter

  int years_per_species[n_species]; // number of years for each species

  int regions_per_cluster[n_clusters]; // number of regions for each species

  int year_species_positions[n_year_species]; // the location of year-species effects

  int region_cluster_positions[n_region_clusters]; // the location of region-species effects

  // int species_intercept_positions[n_species]; // the location of each species intercept

  int did_positions[n_did]; // the location of did effects

  int non_nested_positions[n_non_nested]; // the location of did effects

  matrix[n_observations_seeing,n_parameters] x_seeing; // covariates

  matrix[n_observations_did,n_did_parameters] x_did; // covariates

  int<lower=0, upper = 1> observed[n_observations_seeing]; // observed or not observed

} // close data section

parameters{


// seeing fish parameters

vector[n_parameters] seeing_betas; //vector of regression coefficients for seeing fish

vector<lower = 0>[n_species] seeing_sigma_year_species; // standard deviation of the seeing year effects

vector<lower = 0>[n_clusters] seeing_sigma_region_cluster; // standard deviation of the seeing region effects


} // close parameters


model {
// seen model //////////////////////////////////////////////

int counter;

vector[n_observations_seeing] prob_seen;

 vector[n_year_species] seeing_year_species_betas;

vector[n_region_clusters] seeing_region_cluster_betas;

  vector[n_non_nested] non_nested_seeing_betas;

// seeing model //////////////////////////////////////////////

seeing_year_species_betas = seeing_betas[year_species_positions];

  seeing_region_cluster_betas = seeing_betas[region_cluster_positions];

  non_nested_seeing_betas = seeing_betas[non_nested_positions];

// observed model- hierarchichal structure of year effects
counter = 1;

for (i in 1:n_species){ // hierarchical priors on year effects by species

target += normal_lpdf(segment(seeing_year_species_betas, counter,years_per_species[i]) | 0, seeing_sigma_year_species[i]);

counter = counter + years_per_species[i];

}

// regional hierarchical effects

counter = 1;

for (i in 1:n_clusters){ // hierarchical priors on region effects by species

target += normal_lpdf(segment(seeing_region_cluster_betas, counter, regions_per_cluster[i]) | 0, seeing_sigma_region_cluster[i]);

counter = counter + regions_per_cluster[i];

}

// standardized_prob_seen = 1 ./ (1 + exp(-standard_matrix * seeing_betas));

// prob_seen = 1 ./ (1 + exp(-x_seeing * seeing_betas));



target += student_t_lpdf(non_nested_seeing_betas |6, 0, 2.5);

// target += normal_lpdf(non_nested_seeing_betas | 0, 1);

target += normal_lpdf(seeing_sigma_year_species | 0, 1);
//
target += normal_lpdf(seeing_sigma_region_cluster | 0, 1);


// target += cauchy_lpdf(seeing_sigma_year_species | 0, cauchy_2);
// //
// target += cauchy_lpdf(seeing_sigma_region_cluster | 0, cauchy_2);

target += bernoulli_logit_lpmf(observed |(x_seeing * seeing_betas));

} // close model block



