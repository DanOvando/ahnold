library(rjags) ## For running the MCMC (Gibbs) sampler
# library(coda) ## For converting MCMC objects and diagnostics. Loads as rjags dependency.
library(dclone) ## Allows parallel updating of JAGS models
library(snow) ## Allows clusters: i.e. subsidiary R programmes running separately on your computer's different CPUs
library(jagstools) # For extracting summary statistics from MCMC chain

#----

N <- nrow(climate_df)


## Write BUGS file
model_string <- paste(
  "model{

  for(t in 1:N) {
  mu[t] <- alpha + beta*trf[t] + gamma*volc_mean[t] + delta*soi_mean[t] + eta*amo_mean[t]
  had[t]  ~ dnorm(mu[t], tau)
  y_pred[t] ~ dnorm(mu[t], tau) ## For predictions into the future
  }

  ## Noninformative Priors for all parameters
  alpha ~ dnorm(0, 0.0001)            ## intercept
  beta ~ dnorm(0, 0.0001)             ## trf coef
  gamma ~ dnorm(0, 0.0001)            ## volc coef
  delta ~ dnorm(0, 0.0001)            ## soi coef
  eta ~ dnorm(0, 0.0001)              ## amo coef
  sigma ~ dunif(0, 100)               ## Residual std dev
  tau <- pow(sigma, -2)
  had0 ~ dnorm(0.0, 1.0E-6)           ## Initialising value for prediction
  }"
    )

writeLines(model_string, con = "./BUGSfiles/noninfprior.txt")


##------------------------------------------------------------------------------
## SET UP PARALLELIZATION

load.module("lecuyer") ## JAGS module uses lecuyer random no. generator (to avoid overlap/correlation in a parallel format)

n_chains <- 3 ## i.e. 3 different chains
cl <- makeCluster(n_chains, type = "SOCK") ## SOCK is the simplest cluster formation
parLoadModule(cl, "lecuyer", quiet = T)


##------------------------------------------------------------------------------
## INTIALIZE THE CHAINS.

data_list <- list("N" = N, "had" = clim_df$had, "trf" = clim_df$trf,
                  "volc_mean" = clim_df$volc_mean,
                  "soi_mean" = clim_df$soi_mean, "amo_mean" = clim_df$amo_mean)
inits_list <- function() {
  list(alpha = 0, beta = 0, gamma = 0, delta = 0, eta = 0, sigma = 0.1)
}

##------------------------------------------------------------------------------
## RUN THE CHAINS.

parameters <- c("alpha", "beta", "gamma", "delta", "eta", "sigma", "y_pred")
par_inits <- parallel.inits(inits_list, n.chains = n_chains) ## Initialization

parJagsModel(cl, name = "jags_mod", file = bugs_file,
             data = data_list, inits = par_inits, n.chains = n_chains, n.adapt = 1000)
parUpdate(cl, "jags_mod", n.iter = 1000) ## burn-in
mod_iters <- chain_length / n_chains
mod_samples <- parCodaSamples(cl, "jags_mod", variable.names = parameters,
                              n.iter = mod_iters, n.chain = n_chains)
stopCluster(cl)


##------------------------------------------------------------------------------
## EXTRACT COEFFICIENT POSTERIORS FROM MCMC OUTPUT

## Get coefficients MCMC list into separate matrix for later. Combines all chains into one matrix.##
coefs_mat <- as.matrix(mod_samples[, c(1:6)], iters = F)

## Get summary statistics for tables ##
coefs_tab <- jagsresults(mod_samples, params = c("alpha", "beta", "gamma", "delta", "eta"))