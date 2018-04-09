
// TMB attempt at fitting ahnold
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  /////////load data/////////

  // seen data

  DATA_MATRIX(x_seen_non_nested); // non nested part of data

  DATA_MATRIX(x_seen_year_species); // year species effects

  DATA_MATRIX(x_seen_region_cluster); //region cluster data

  DATA_VECTOR(log_density); // observed log densities

  DATA_IVECTOR(seen_species_index); // index the same rows as x_seen showing what species is in that row

  DATA_IVECTOR(seen_year_species_index); // index the same length as seen_year_species_betas that shows what species that beta corresponds to

  DATA_IVECTOR(seen_region_cluster_index); // index the same length as seen_region_cluster_betas that shows what species that beta corresponds to

  // seeing data

  DATA_MATRIX(x_seeing_non_nested); // non nested part of data

  DATA_MATRIX(x_seeing_year_species); // year species effects

  DATA_MATRIX(x_seeing_region_cluster); //region cluster data

  DATA_VECTOR(any_seen); // observed log densities

  DATA_IVECTOR(seeing_species_index); // index the same rows as x_seen showing what species is in that row

  DATA_IVECTOR(seeing_year_species_index); // index the same length as seen_year_species_betas that shows what species that beta corresponds to

  DATA_IVECTOR(seeing_region_cluster_index); // index the same length as seen_region_cluster_betas that shows what species that beta corresponds to

  // did data

  DATA_MATRIX(standard_non_nested); // standardized non nested part of data

  DATA_MATRIX(standard_years); // standardized year species effects

  DATA_MATRIX(standard_regions); // standardized region cluster data

  DATA_MATRIX(x_did_non_nested); // non nested did parameters

  DATA_MATRIX(x_did_species_enviro); // species effects for did model

  DATA_MATRIX(x_did_species_intercepts); // intercepts for each species

  DATA_IVECTOR(x_did_species_effects_index); // location of each species in random effects

  DATA_INTEGER(species_enviro); // LOGICAL INDICATING WHETHER TO USE SPECIES EFFECTS

  /////////define parameters/////////

  PARAMETER_VECTOR(seen_non_nested_betas);

  PARAMETER_VECTOR(seen_year_species_betas);

  PARAMETER_VECTOR(seen_region_cluster_betas);

  PARAMETER_VECTOR(seen_year_species_sigmas);

  // PARAMETER_VECTOR(seen_region_cluster_sigmas);

  PARAMETER_VECTOR(seen_density_species_sigma);

  // seeing parameters

  PARAMETER_VECTOR(seeing_non_nested_betas);

  PARAMETER_VECTOR(seeing_year_species_betas);

  PARAMETER_VECTOR(seeing_region_cluster_betas);

  PARAMETER_VECTOR(seeing_year_species_sigmas);

  // PARAMETER_VECTOR(seeing_region_cluster_sigmas);

  // did parameters

  PARAMETER_VECTOR(did_non_nested_betas);

  PARAMETER_VECTOR(did_species_intercept_betas);

  // if (species_enviro == 1){

  PARAMETER_VECTOR(did_species_environment_betas);

  PARAMETER(species_enviro_sigma);

  // }

  PARAMETER(species_intercept_sigma);

  PARAMETER(did_sigma);


  /////////process parameters and data/////////

  Type nll; //blank storage for accumulated nll

  nll = 0;

  int i_max;
  /////////seen fish/////////

  matrix<Type> non_nested_effects = x_seen_non_nested * seen_non_nested_betas;
  //
  matrix<Type> year_species_effects = x_seen_year_species * seen_year_species_betas;
  //
  matrix<Type> region_cluster_effects = x_seen_region_cluster * seen_region_cluster_betas;
  //
  matrix<Type> log_density_hat = non_nested_effects + year_species_effects + region_cluster_effects;
  //
  i_max = seen_year_species_betas.size();

  // std::cout << i_max << "\\n";

  for (int i = 0; i < i_max; i++){

    nll -= dnorm(seen_year_species_betas(i), Type(0), exp(seen_year_species_sigmas(seen_year_species_index(i) - 1)), true);


  } // close year species effects
  //
  // i_max = seen_region_cluster_betas.size();
  // //
  // for (int i = 0; i < i_max ; i++){
  //
  //   nll -= dnorm(seen_region_cluster_betas(i), Type(0), exp(seen_region_cluster_sigmas(seen_region_cluster_index(i) - 1)), true);
  //
  //
  // } // close region cluster effects
  // //
  //
  i_max = log_density.size();

  for (int i = 0; i < i_max; i++){

    nll -= dnorm(log_density(i), log_density_hat(i), exp(seen_density_species_sigma(seen_species_index(i) - 1)), true);

  } // close log density nll


  /////////seeing fish/////////

  matrix<Type> seeing_non_nested_effects = x_seeing_non_nested * seeing_non_nested_betas;
  //
  matrix<Type> seeing_year_species_effects = x_seeing_year_species * seeing_year_species_betas;
  //
  matrix<Type> seeing_region_cluster_effects = x_seeing_region_cluster * seeing_region_cluster_betas;
  //
  matrix<Type> logit_scale_prob_seeing = seeing_non_nested_effects + seeing_year_species_effects + seeing_region_cluster_effects;
  //
  i_max = seeing_year_species_betas.size();

  for (int i = 0; i < i_max; i++){

    nll -= dnorm(seeing_year_species_betas(i), Type(0), exp(seeing_year_species_sigmas(seeing_year_species_index(i) - 1)), true);

  } // close year species effects

  // i_max = seeing_region_cluster_betas.size();
  //
  // for (int i = 0; i < i_max ; i++){
  //
  //   nll -= dnorm(seeing_region_cluster_betas(i), Type(0), exp(seeing_region_cluster_sigmas(seeing_region_cluster_index(i) - 1)), true);
  //
  // } // close region cluster effects
  // //

  vector<Type> prob_seeing =  1/ (1 + exp(-logit_scale_prob_seeing.array()));

  nll -= sum(dbinom(any_seen,Type(1),prob_seeing, true));

  vector<Type> logit_standardized_yearly_prob_seeing = standard_non_nested * seeing_non_nested_betas + standard_regions * seeing_region_cluster_betas + standard_years * seeing_year_species_betas;

  vector<Type> standardized_yearly_prob_seeing = 1 / (1 + exp(-logit_standardized_yearly_prob_seeing));

  /////////did/////////

  vector<Type> seen_abundance_index = exp(seen_year_species_betas);

  vector<Type> standardized_abundance = standardized_yearly_prob_seeing * seen_abundance_index;

  vector<Type> log_standardized_abundance = log(standardized_abundance);

  matrix<Type> did_non_nested_effects = x_did_non_nested * did_non_nested_betas;

  matrix<Type> did_intercept_effects = x_did_species_intercepts * did_species_intercept_betas;

  // if (species_enviro == 0){
  // matrix<Type> log_standardized_abundance_hat = did_non_nested_effects +did_intercept_effects;
  // }
  // if (species_enviro == 1){

    matrix<Type> did_enviro_effects = x_did_species_enviro * did_species_environment_betas;

    matrix<Type> log_standardized_abundance_hat = did_non_nested_effects +did_intercept_effects + did_enviro_effects;

  // }
  // std::cout << log_standardized_abundance_hat << "\\n";

   i_max = did_species_intercept_betas.size();
  //
   for (int i = 0; i < i_max ; i++){

     nll -= dnorm(did_species_intercept_betas(i), Type(0), exp(species_intercept_sigma), true);

   } // close species hierarchical effects

i_max = standardized_abundance.size();

  for (int i = 0; i < i_max; i++){

    nll -= dnorm(log_standardized_abundance(i),log_standardized_abundance_hat(i),exp(did_sigma), true);

  }

  // if (species_enviro == 1){

    i_max = did_species_environment_betas.size();
    //
    for (int i = 0; i < i_max ; i++){

      nll -= dnorm(did_species_environment_betas(i), Type(0), exp(species_enviro_sigma), true);

    } // close species hierarchical effects

  //}

  /////////outputs/////////

  ADREPORT(standardized_abundance);

  REPORT(log_density_hat);

  REPORT(prob_seeing);

  REPORT(standardized_yearly_prob_seeing);

  REPORT(standardized_abundance);

  REPORT(log_standardized_abundance);

  REPORT(seen_abundance_index);

  REPORT(log_standardized_abundance_hat);

  return nll;

  // return(nll);

}
