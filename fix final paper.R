$$\p_{f,s,y} = \frac{e^{\gamma_{f,s,y}}}{1+e^{\gamma_{f,s,y}}}$$

  where we model $\gamma$ as a linear function of covariates selected to explain the probability of observing any fish at all

$$\gamma_{f,s,y} = \tau_{0} + \tau_{1}$$




  We include the DiD estimator within the hurdle portion of the model to estiamte the effect of the MLPA on the probability of zero-density observations.

For the observations with *d* > 0, we model the expected density as

$$ d_{f,s,y} = \beta_{1}fished_{f} + \beta_{2}mpa_{s,y} + \beta_{3}FxM_{f,s,y} +
  \sum\beta_{4}region_{s} + \sum\beta_{5}trophic_{f} +  \beta_{6}yearsmpa_{s,y} +... $$

  $$\beta_{7}FxYM_{f,s,y} + \beta_{7}linf_{f} + \beta_{8}vbk_{f} + \beta_{9}temp_{s,y} +
  \beta_{10}vis_{s,y} + $$

  $$ \beta_{11}templag1_{s,y} + \beta_{12}templag2_{s,y} + \beta_{13}templag3_{s,y} +
  \beta_{14}templag4_{s,y} + \beta_{15} $$

  ### Likelihoods


  For the first part we can tally the number of candidate transects *n* and the number of transects with positive densities *w* that comprised that density estimate. This translates to a binomial process where the likelihood of a given aggregation is given as

$$ Binomial likelihood $$

  where *p* is the estimated probability of the the number of positive densities.

$$ LL^{bin}_{f,s,y} = dbinom(w_{f,s,y},n_{f,s,y},p_{f,s,y}) $$

  For the observations with *d*>0, the we assume a normal likelihood.

### Priors

We construct a series of hierarchichal priors to account for clustering in the data. All parameters not listed here are assumed to have uniform priors. We model the year-term fixed effects as originating from a common normal prior

Year~y~ ~N(0,$\sigma^{y}$)

where $\sigma^{y}$ is modeled through the uninformative hyperprior

$\sigma^{y}$ ~gamma(2,0.5)

The year coefficients for the binomial portion of the model have priors constructed in the same way.

Year~y~^p^ ~ N(0,$\sigma^{y,b}$)

$\sigma^{y,b}$ ~ gamma(2,0.5)

We also assume that the prior on the regional fixed effects is distribtued normal per
Region~r~ ~ N(0,$\sigma^{r}$)

with hyperprior
$\sigma^{r}$ ~ gamma(2,0.5)

Lastly the prior on the standard deviation of the densities $\sigma$ is distributed

$\sigma$ ~ gamma(2,0.5)

## Results


You can reproduce most of the diagnostics from the Mosquito exercise. For now, let's focus on observed vs predicted to see if you're even in the same ballpark

```{r assess MCMC fit}

# load(paste(runpath,'MCMC Results.Rdata', sep = ''))
#
# dat <- bayes_reg$Data
#
# post <- bayes_reg$demon_fit$Posterior1
#
# its <- dim(post)[1]
#
# burn <- 0.5
#
# trimmed_post <- thin_mcmc(post[(burn*its):its,], thin_every = its/1000)
#
# ggmcmc(ggs(mcmc(trimmed_post)), file = paste(runpath,'mcmc diagnostics.pdf', sep = ''))
#

```


### Diagnostis

## Discussion