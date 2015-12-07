mlpa_delta_likelihood <- function(parm, Data,reg_model = 'tobit')
{
  ### Separate out parameters for observed density  ----

  beta <- parm[Data$pos_den_beta]

  parm[Data$pos_den_sigma] <- interval(parm[Data$pos_den_sigma], 1e-100, Inf)

  sigma_year <- parm[Data$pos_sigma_year]

  sigma_bi_year <- parm[Data$pos_sigma_bi_year]

  sigma_region <- parm[Data$pos_sigma_region]

  sigma_density <- parm[Data$pos_sigma_density]

  ### Log-Priors for observed ----

  year_priors <- sum(dnorm(parm[Data$pos_den_time_terms],0,sigma_year, log = T))

  bi_year_priors <- sum(dnorm(parm[Data$pos_bi_time_terms],0,sigma_bi_year, log = T))

  sigma_year_prior <- dnorm(sigma_year, .1,.2, log = T)

  sigma_bi_year_prior <- dnorm(sigma_bi_year, .1,.2, log = T)

  region_priors <- sum(dnorm(parm[Data$parm.names %in% Data$site_vars ],
                             0,sigma_region, log = T))

  sigma_region_prior <- dnorm(sigma_region, .1,.2, log = T)

  sigma_density_prior <- dnorm(sigma_density, .1,.2, log = T)

  ### Hurdle Log-Likelihood ----

  bi_beta <- parm[Data$vars_for_binom]

  bi_dat <- Data$reg_dat[,Data$beta_to_use_binom]

  #   bi_hat = tcrossprod(bi_beta, bi_dat)

  bi_hat <- bi_dat %*% bi_beta

  #   mat_bi_beta <- matrix(rep(bi_beta,Data$N),nrow  = Data$N, ncol = length(bi_beta), byrow = T)
  #
  #   bi_hat <- rowSums(bi_dat * mat_bi_beta)

  prob_hat <- exp(bi_hat)/(1 + exp(bi_hat))

  bi_loglike <-  sum(dbinom(Data$binom_dep_var,1,prob_hat, log = T))

  ### Density Log-likelihood ----

  #   mu <- tcrossprod(beta, Data$reg_dat[,Data$beta_to_use_binom == F])

  mu <- Data$reg_dat[,Data$beta_to_use_binom == F] %*% beta

  observed_density <- Data$dep_var

#     min_val <- min(observed_density)
#
#     seen <- observed_density > min_val
#
#     density_loglike <- sum(dnorm(observed_density[seen], mu[seen],sigma_density, log=TRUE)^Data$binom_dep_var)
  density_loglike <- sum(dnorm(observed_density, mu,sigma_density, log=TRUE)^(Data$binom_dep_var))

  ### Log-Posterior


#   local_vars <- ls()
#
#   post_comp <- local_vars[grepl('_loglike', local_vars) | grepl('_prior', local_vars)]
#
#   post_fmla <- paste(post_comp, sep = '+')

  LP <- density_loglike + bi_loglike  + year_priors + bi_year_priors +
    sigma_year_prior + sigma_bi_year_prior + region_priors + sigma_region_prior + sigma_density_prior

  LL <- density_loglike + bi_loglike

  Modelout <- list(LP=LP,Dev=-2*LL, Monitor=LP,yhat = 1,
                   parm=parm)

  #   Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
  #                    yhat=rnorm(length(mu), mu, sigma), parm=parm)

  #     show(proc.time() - a)

  return(Modelout)
}
