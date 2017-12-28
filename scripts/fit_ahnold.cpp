
// TMB attempt at fitting ahnold
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  /////////load data/////////

  // seen data

  DATA_MATRIX(x_seen_non_nested); // non nested part of data

  DATA_MATRIX(x_seen_did); // difference in difference terms

  DATA_IVECTOR(seen_species_index); // index the same rows as x_seen showing what species is in that row

  DATA_VECTOR(log_density); // observed log densities

  // seeing data

  DATA_MATRIX(x_seeing_non_nested); // non nested part of data

  DATA_MATRIX(x_seeing_did); // difference in difference terms

  DATA_VECTOR(any_seen); // observed log densities



  // did data

  DATA_MATRIX(standard_non_nested); // standardized non nested part of data

  DATA_MATRIX(standard_did_with_mpa); // standardized did estimators with mpa

  DATA_MATRIX(standard_did_without_mpa); // standardized did estimators without mpa

  /////////define parameters/////////

  PARAMETER_VECTOR(seen_non_nested_betas); // NON DID BETAS

  PARAMETER_VECTOR(seen_did_betas);  // DIFFERENCE IN DIFFERENCE BETAS

  PARAMETER_VECTOR(seen_density_species_sigma);

  // seeing parameters

  PARAMETER_VECTOR(seeing_non_nested_betas);

  PARAMETER_VECTOR(seeing_did_betas);

  /////////process parameters and data/////////

  Type nll; //blank storage for accumulated nll

  nll = 0;

  int i_max;
  /////////seen fish/////////

  matrix<Type> non_nested_effects = x_seen_non_nested * seen_non_nested_betas;

  matrix<Type> did_effects = x_seen_did * seen_did_betas;

  matrix<Type> log_density_hat = non_nested_effects + did_effects; // + year_species_effects + region_cluster_effects;

  i_max = log_density.size();

  for (int i = 0; i < i_max; i++){

    nll -= dnorm(log_density(i), log_density_hat(i), exp(seen_density_species_sigma(seen_species_index(i) - 1)), true);

  } // close log density nll


  /////////seeing fish/////////

  matrix<Type> seeing_non_nested_effects = x_seeing_non_nested * seeing_non_nested_betas;

  matrix<Type> seeing_did_effects = x_seeing_did * seeing_did_betas;

  matrix<Type> logit_scale_prob_seeing = seeing_non_nested_effects + seeing_did_effects; //+ seeing_year_species_effects + seeing_region_cluster_effects;

  vector<Type> prob_seeing =  1/ (1 + exp(-logit_scale_prob_seeing.array()));

  nll -= sum(dbinom(any_seen,Type(1),prob_seeing, true));

  // predict with mpa
  vector<Type> with_mpa_logit_standardized_yearly_prob_seeing = standard_non_nested * seeing_non_nested_betas + standard_did_with_mpa * seeing_did_betas; // + standard_years * seeing_year_species_betas;

  vector<Type> standardized_yearly_prob_seeing = 1 / (1 + exp(-with_mpa_logit_standardized_yearly_prob_seeing));

  vector<Type> with_mpa_standardized_abundance = standard_non_nested * seen_non_nested_betas + standard_did_with_mpa * seen_did_betas; // + standard_years * seeing_year_species_betas;

  vector<Type> log_abundance_with_mpa = standardized_yearly_prob_seeing * with_mpa_standardized_abundance;

  // predict without mpa
  vector<Type> without_mpa_logit_standardized_yearly_prob_seeing = standard_non_nested * seeing_non_nested_betas + standard_did_without_mpa * seeing_did_betas; // + standard_years * seeing_year_species_betas;

  standardized_yearly_prob_seeing = 1 / (1 + exp(-without_mpa_logit_standardized_yearly_prob_seeing));

  vector<Type> without_mpa_standardized_abundance = standard_non_nested * seen_non_nested_betas + standard_did_without_mpa * seen_did_betas; // + standard_years * seeing_year_species_betas;

  vector<Type> log_abundance_without_mpa = standardized_yearly_prob_seeing * without_mpa_standardized_abundance;

  vector<Type> net_did = log_abundance_with_mpa - log_abundance_without_mpa;
  /////////outputs/////////

  ADREPORT(net_did);

  REPORT(log_density_hat);

  REPORT(prob_seeing);

  REPORT(standardized_yearly_prob_seeing);

  REPORT(net_did);

  return nll;

  // return(nll);

}
