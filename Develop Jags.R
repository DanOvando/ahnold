
dat <- list( = , effort = ogdat$Effort, count = ogdat$Count,n = n,
            numpars = 10)

inits <- list(den_time_terms=rep(log(5),10), mean = rep(5,10), tau = rep(10,10),
              .RNG.name="base::Super-Duper", .RNG.seed=1, log.hypermu = 5, precision = 1)

inits <- list(lambda=rep(log(5),10), mean = rep(5,10), tau = rep(10,10),
              .RNG.name="base::Super-Duper", .RNG.seed=1, log.hypermu = 5, precision = 1)



year_priors ~ sum(dnorm(den_time_terms,0,sigma_year, log = T))

bi_year_priors ~ sum(dnorm(bi_time_terms,0,sigma_bi_year, log = T))

sigma_year_prior ~ dgamma(sigma_year, 2,0.5, log = T)

sigma_bi_year_prior ~ dgamma(sigma_bi_year, 2,0.5, log = T)

region_priors ~ sum(dnorm(den_region_terms,0,sigma_region, log = T))

sigma_region_prior ~ dgamma(sigma_region, 2,.5, log = T)

sigma_density_prior ~ dgamma(sigma_density, 2,.5, log = T)

# Hurdle Log-Likelihood

bi_beta <- parm[Data$beta_to_use_binom]

bi_dat <- Data$bi_reg_mat

bi_hat <- pmin(10,bi_dat %*% bi_beta)

prob_hat <- exp(bi_hat)/(1 + exp(bi_hat))

bi_loglike ~  sum(dbinom(binom_dep_var,1,max(1e-15,prob_hat), log = T)) # actual values of 1 return -inf if not fulfilled

# Density Log-likelihood

mu <- Data$den_reg_mat %*% beta

observed_density <- Data$dep_var

density_loglike ~ sum(dnorm(observed_density, mu,sigma_density, log=TRUE)^(binom_dep_var))
