// STAN model for an integrated fit to the ahnodl model
// Goes from delta model of abundances to DiD across all species
// seen means the model fit to observed fish
// seeing means the model of the probability of observing anything
data{

  int<lower = 0> n_parameters; //number of parameters

  int<lower = 0> n_did_parameters; //number of did parameters

  int<lower = 0> n_observations_seen; //number of seen observations

  int<lower = 0> n_observations_did; //number of did observations

  int<lower = 0> n_species; //number of species

  int<lower = 0> n_clusters; //number of geographic clusters

  int<lower = 0> n_year_species; //number of year species combinations

  int<lower = 0> n_region_clusters; // number of region species combinations

  int<lower = 0> n_did; // number of did effects

  int<lower = 0> n_standard; // number of standardized matrix observations

  int<lower = 0> n_non_nested; // number of non-nested parameters

  real<lower = 0> cauchy_2 ; // cauchy parameter

  int years_per_species[n_species]; // number of years for each species

  int regions_per_cluster[n_clusters]; // number of regions for each species

  int year_species_positions[n_year_species]; // the location of year-species effects

  int region_cluster_positions[n_region_clusters]; // the location of region-species effects

  int did_positions[n_did]; // the location of did effects

  int non_nested_positions[n_non_nested]; // the location of did effects

  matrix[n_observations_seen,n_parameters] x_seen; // covariates

  // matrix[n_standard,n_parameters] standard_matrix; // covariates

  matrix[n_observations_did,n_did_parameters] x_did; // covariates

  vector[n_observations_seen] log_density; // observed log densities

  int log_density_species_index[n_species,2]; // matrix of start and end positions of observed densities

} // close data section

parameters{

// seen fish parameters

vector[n_parameters] betas; // vector of regression coefficients for seen fish

vector<lower = 0>[n_species] sigma_year_species; // standard deviation of the year effects

vector<lower = 0>[n_clusters] sigma_region_cluster; // standard deviation of the year effects

vector<lower = 0>[n_species] sigma_density_species; // standard deviation of the observed densities for each species

} // close parameters


model {
// seen model //////////////////////////////////////////////

int counter;

vector[n_observations_seen] log_density_hat;

vector[n_year_species] year_species_betas;

vector[n_region_clusters] region_cluster_betas;

vector[n_non_nested] non_nested_seen_betas;

// pull out specific effects

  year_species_betas = betas[year_species_positions];

  region_cluster_betas = betas[region_cluster_positions];

  non_nested_seen_betas = betas[non_nested_positions];


// observed model- hierarchichal structure of year effects
counter = 1;

for (i in 1:n_species){ // hierarchical priors on year effects by species

segment(year_species_betas,counter,years_per_species[i]) ~ normal(0, sigma_year_species[i]);

// target += normal_lpdf(segment(year_species_betas,counter,years_per_species[i]) | 0, sigma_year_species[i]);

counter = counter + years_per_species[i];

}

// regional hierarchical effects

counter = 1;

for (i in 1:n_clusters){ // hierarchical priors on region effects by species

segment(region_cluster_betas,counter, regions_per_cluster[i]) ~ normal(0, sigma_region_cluster[i]);

// target += normal_lpdf(segment(region_cluster_betas,counter, regions_per_cluster[i])| 0, sigma_region_cluster[i]);

counter = counter + regions_per_cluster[i];

}

log_density_hat = x_seen*betas;


// rest of the likelihood //////////////////////////////


non_nested_seen_betas ~ normal(0, 1);

 sigma_density_species ~ normal(0, 1);

sigma_year_species ~ normal(0, 1);

sigma_region_cluster ~ normal(0, 1);

// target += normal_lpdf(non_nested_seen_betas | 0, 1);
//
// target += normal_lpdf( sigma_density_species | 0, 1);
//
// target += normal_lpdf(sigma_year_species | 0, 1);
// //
// target += normal_lpdf(sigma_region_cluster |0, 1);



for (i in 1:n_species){ // cluster sigmas by species
// one way to do this would be ragged arrays, which are apparently complicated. For now, have to order by species abr provide indexing

log_density[log_density_species_index[i,1]:log_density_species_index[i,2]] ~ normal(log_density_hat[log_density_species_index[i,1]:log_density_species_index[i,2]], sigma_density_species[i]);

// target += normal_lpdf(log_density[log_density_species_index[i,1]:log_density_species_index[i,2]] | log_density_hat[log_density_species_index[i,1]:log_density_species_index[i,2]], sigma_density_species[i]);

}

} // close model block



