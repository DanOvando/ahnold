// STAN model for an integrated fit to the ahnodl model
// Goes from delta model of abundances to DiD across all species
// seen means the model fit to observed fish
// seeing means the model of the probability of observing anything
data{

 // common data ////////////////////

  int<lower = 0> n_parameters; //number of parameters


  int<lower = 0> n_species; //number of species

  int<lower = 0> n_clusters; //number of geographic clusters

// seen data //////////////////////

  int<lower = 0> n_observations_seen; //number of seen observations

  int<lower = 0> n_year_species; //number of year species combinations

  int<lower = 0> n_region_clusters; // number of region species combinations

  int<lower = 0> n_non_nested; // number of non-nested parameters

  int years_per_species[n_species]; // number of years for each species

  int regions_per_cluster[n_clusters]; // number of regions for each species

  int year_species_positions[n_year_species]; // the location of year-species effects

  int region_cluster_positions[n_region_clusters]; // the location of region-species effects

  int non_nested_positions[n_non_nested]; // the location of did effects

  matrix[n_observations_seen,n_parameters] x_seen; // covariates

  vector[n_observations_seen] log_density; // observed log densities

  int x_seen_species_index[n_observations_seen]; // observed log densities

  int log_density_species_index[n_species,2]; // matrix of start and end positions of observed densities

  // seeing data /////////////////////

  int<lower = 0> n_observations_seeing; //number of observations


  matrix[n_observations_seeing,n_parameters] x_seeing; // covariates

    int<lower=0, upper = 1> observed[n_observations_seeing]; // observed or not observed


} // close data section

parameters{

// seen fish parameters /////////////////

vector[n_parameters] betas; // vector of regression coefficients for seen fish

vector<lower = 0>[n_species] sigma_year_species; // standard deviation of the year effects

vector<lower = 0>[n_clusters] sigma_region_cluster; // standard deviation of the year effects

vector<lower = 0>[n_species] sigma_density_species; // standard deviation of the observed densities for each species

// seeing fish parameters /////////////////

vector[n_parameters] seeing_betas; //vector of regression coefficients for seeing fish

vector<lower = 0>[n_species] seeing_sigma_year_species; // standard deviation of the seeing year effects

vector<lower = 0>[n_clusters] seeing_sigma_region_cluster; // standard deviation of the seeing region effects



} // close parameters


model {

int counter;

vector[n_observations_seen] log_density_hat;

vector[n_year_species] year_species_betas;

vector[n_region_clusters] region_cluster_betas;

vector[n_non_nested] non_nested_seen_betas;

vector[n_observations_seeing] prob_seen;

 vector[n_year_species] seeing_year_species_betas;

vector[n_region_clusters] seeing_region_cluster_betas;

  vector[n_non_nested] non_nested_seeing_betas;


// seen model //////////////////////////////////////////////


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


//seeing model //////////////////////////////


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

prob_seen = 1 ./ (1 + exp(-x_seeing * seeing_betas));


// rest of the seen likelihood //////////////////////////////


non_nested_seen_betas ~ normal(0, 1);

sigma_density_species ~ normal(0, 1);

sigma_year_species ~ normal(0, 1);

sigma_region_cluster ~ normal(0, 1);

for (i in 1:n_species){ // cluster sigmas by species
// one way to do this would be ragged arrays, which are apparently complicated. For now, have to order by species abr provide indexing


log_density[log_density_species_index[i,1]:log_density_species_index[i,2]] ~ normal(log_density_hat[log_density_species_index[i,1]:log_density_species_index[i,2]], sigma_density_species[i]);

// target += normal_lpdf(log_density[log_density_species_index[i,1]:log_density_species_index[i,2]] | log_density_hat[log_density_species_index[i,1]:log_density_species_index[i,2]], sigma_density_species[i]);

}

// rest of seeing likelihood //////////////

 non_nested_seeing_betas ~ student_t(6, 0, 2.5);

 seeing_sigma_year_species ~normal(0, 1);

 seeing_sigma_region_cluster ~normal(0, 1);

 observed ~ bernoulli_logit((x_seeing * seeing_betas));


} // close model block

generated quantities{

vector[n_observations_seen] sigma_density_vector;

vector[n_observations_seen] predicted_log_density_hat;

vector[n_observations_seen] mean_log_density_hat;

mean_log_density_hat = x_seen*betas;

for (i in 1:n_observations_seen){

  sigma_density_vector[i] = sigma_density_species[x_seen_species_index[i]];

}


for (i in 1:n_observations_seen){ // cluster sigmas by species
// one way to do this would be ragged arrays, which are apparently complicated. For now, have to order by species abr provide indexing

predicted_log_density_hat[i] = normal_rng(mean_log_density_hat[i], sigma_density_vector[i]);

}


}


