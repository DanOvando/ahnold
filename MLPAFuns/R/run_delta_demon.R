#' Function to run MLPA demon using a
#' delta/hurdle method
#'
#' \code{run_delta_demon} uses Laplaces Demon to
#' fit a bayesian hierarchechal model to the MLPA
#' data
#' @param dat dataframe of regression variables
#' @param dep_var the name of the dependent variable
#' in the regression
#' @param pos_vars vector of variables to be potentially included
#' in the regression
#' @param delta_vars vector of variables to be included in the
#' binomial part of the regression
#' @param iterations the number of runs for the MCMC
#' @param status the fraction of iterations at which to display status
#' @param thin the thinning rate for the MCMC
#' @param burn the amount of the chain to burn off
#' @param scale_numerics T or F to center and scale numeric covariates
#' @param runpath the folder to store results
#' @param acceptance_rate the target acceptance rate
#' @param method Summon Demon for standard MCMC, Summon Reversible Demon for
#' reversible MCMC
#'

run_delta_demon <- function(dat,dep_var,pos_vars,delta_vars,iterations = 1000,status = .05,thin = 1,burn = .5,
                            scale_numerics = F,
                            runpath, acceptance_rate  = 0.234, method = 'Summon Demon') {


  # Convert data to regression format----
  observed_dat <- prep_demon(dat,pos_vars = pos_vars, scale_numerics = scale_numerics)

  binom_dat <- prep_demon(dat,pos_vars = delta_vars, scale_numerics = scale_numerics)

  vars_to_use_binom <- colnames(binom_dat)

  colnames(binom_dat) <- paste('bi',colnames(binom_dat), sep = '.')

  reg_dat <- cbind(observed_dat, binom_dat)

  J <- dim(reg_dat)[2]

  mon.names <- "LP"

  names <- colnames(reg_dat)

  beta_to_use_binom <- grepl('bi.', names, fixed = T)

  sigmas <- c('sigma_density','sigma_year','sigma_region','sigma_bi_year')

  time_vars <- names[grepl('_year.',names, fixed = T)]

  site_vars <- names[grepl('region.',names, fixed = T)]

  species_vars <- names(grepl('trophic.',names, fixed = T))

  parm.names <- c(colnames(reg_dat),sigmas)

  # Prepare reversible jump parameteres ----

  off_the_table <- (grepl('fished', parm.names, fixed = T) | grepl('mpa_applied', parm.names,fixed = T) |
                      grepl('fished_x_mpa', parm.names,fixed = T) | grepl('constant', parm.names,fixed = T) |
                      grepl('.factor', parm.names,fixed = T) | grepl('sigma', parm.names,fixed = T))

  selectable <- rep(1,length(parm.names))

  selectable[off_the_table] <- 0

  selected <- selectable

  bin.n <- J #Maximum allowable model size

  bin.p <- 0.9 #Most probable size:  bin.p x bin.n is binomial mean and median

  parm.p <- rep(1/J,J+1)

  # Prepare indices for MCMC ----

  vars_for_binom <- grepl('bi.', parm.names, fixed = T)

  pos_vars_for_binom <- which(vars_for_binom)

  parm <- rep(0,length(parm.names))

  pos_den_beta <- which(parm.names %in% colnames(reg_dat) & vars_for_binom == F)

  pos_den_sigma <- which(parm.names %in% sigmas & vars_for_binom == F)

  pos_den_time_terms <- which(parm.names %in% time_vars & vars_for_binom == F)

  pos_den_region_terms <- which(parm.names %in% site_vars & vars_for_binom == F)

  pos_bi_time_terms <- which(parm.names %in% time_vars & vars_for_binom == T)

  pos_any_betas <- which(grepl('sigma_',parm.names, fixed = T) == F)

  pos_any_sigmas <-  which(grepl('sigma_',parm.names, fixed = T))

  pos_sigma_year <- which(grepl('sigma_year',parm.names, fixed = T))

  pos_sigma_region <- which(grepl('sigma_region',parm.names, fixed = T))

  pos_sigma_density <- which(grepl('sigma_density',parm.names, fixed = T))

  pos_sigma_bi_year <- which(grepl('sigma_bi_year',parm.names, fixed = T))

  PGF <- function(Data) {
    beta <- rnorm(length(Data$pos_any_betas))
    sigma <- runif(length(Data$pos_any_sigmas))
    return(c(beta, sigma))
  }

  # Subset data to possible for regression ----

  possible <- is.na(as.matrix(reg_dat[,beta_to_use_binom == F]) %*% parm[pos_den_beta]) == F

  reg_dat <- reg_dat[possible,]

  dat <- dat[possible,]

  binom_dep_var <- as.numeric(dat[,dep_var] > min(dat[,dep_var]))

  N <- dim(dat)[1]

  # Prepare priors ----

  dep_sd <- sd(dat$log_density)

  # Prepare the Demon's snacks ----
  #

  pos_vars_for_binom <- which(vars_for_binom)

  pos_beta_to_use_binom <- which(beta_to_use_binom)

  bi_reg_mat = as.matrix(reg_dat[,beta_to_use_binom])

  den_reg_mat = as.matrix(reg_dat[,beta_to_use_binom == F])


  Data <- list(N = N,
               J=J,
               PGF=PGF,
               bi_reg_mat = bi_reg_mat,
               den_reg_mat = den_reg_mat,
               reg_dat = as.matrix(reg_dat),
               mon.names = mon.names,
               binom_dep_var = binom_dep_var,
               parm.names = parm.names,
               beta_to_use_binom = beta_to_use_binom,
               vars_for_binom = vars_for_binom,
               pos_beta_to_use_binom = pos_beta_to_use_binom,
               pos_den_beta = pos_den_beta,
               pos_den_sigma = pos_den_sigma,
               pos_den_time_terms = pos_den_time_terms,
               pos_den_region_terms = pos_den_region_terms,
               pos_any_betas = pos_any_betas,
               pos_any_sigmas = pos_any_sigmas,
               pos_sigma_year = pos_sigma_year,
               pos_sigma_density = pos_sigma_density,
               pos_sigma_region = pos_sigma_region,
               pos_sigma_bi_year = pos_sigma_bi_year,
               pos_bi_time_terms = pos_bi_time_terms,
               dep_var = as.matrix(dat[,dep_var]),
               time_vars = time_vars,
               site_vars = site_vars,
               species_vars = species_vars)
  Initial.Values <- GIV(mlpa_delta_likelihood, Data, PGF=TRUE)

  # Run Demon ----

  if (method == 'Summon Demon')
  {
    Fit <- LaplacesDemon(mlpa_delta_likelihood, Data=Data, Initial.Values = Initial.Values,
                         Covar=NULL, Iterations=iterations, Status=iterations*status, Thinning=1,
                         Algorithm = 'HARM', Specs=list(alpha.star=acceptance_rate, B = NULL),
                         parm.names = parm.names)

  }
  if (method == 'Summon Reversible Demon'){

    Fit <- LaplacesDemon(mlpa_delta_likelihood, Data=Data, Initial.Values = Initial.Values,
                         Covar=NULL, Iterations=iterations, Status=iterations*status, Thinning=1,
                         Algorithm = 'RJ', Specs=list(bin.n=bin.n, bin.p=bin.p,
                                                      parm.p=parm.p, selectable=selectable,
                                                      selected=selected), parm.names = parm.names)

  }

  return(list(demon_fit = Fit,Data = Data))
}
