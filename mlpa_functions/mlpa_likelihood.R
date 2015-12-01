#' calculate likelhood of MLPA
#'
#' \code{mlpa_likelihood} returns the negative log likelihood
#' of the MLPA data
#' @param

mlpa_likelihood <- function()
{

  #Data
  # "hurdle" model
  # First, fit a binominal process to the data to estimate the probability
  # of getting a zeros.
  #
  # The likelihood of the data will theb
  # p(binomial()) if y <= -2.322 whatever
  # p(normal()) if y > -2.322 whatever
  #
  # The priors
  # This part is fun. Let's just go with independent priors.
  # Time specific parameters will have a prior will be drawn from something
  # Region specific parameters will have have region specific hyperpriors
  # year: fixed effect. Prior will be normal/uniform with mean 0 and a big variance
  # site_side: fixed effect each drawn from a set of normal/uniform priors at the region level
  # region: fixed effects from a common distribution with
  #

}