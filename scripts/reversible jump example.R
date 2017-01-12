rm(list = ls())
library(LaplacesDemon)
N <-  1000
J <-  100 #Number of predictors, including the intercept
X <- matrix(1,N,J)
for (j in 2:J) {X[,j] <- rnorm(N,runif(1,-3,3),runif(1,0.1,1))}
beta.orig <- runif(J,-3,3)
zero <- sample(2:J, round(J*0.9)) #Assign most parameters to be zero
beta.orig[zero] <- 0
e <- rnorm(N,0,0.1)
y <- as.vector(tcrossprod(beta.orig, X) + e)
mon.names <- "LP"
parm.names <- as.parm.names(list(beta=rep(0,J), sigma=0))
pos.beta <- grep("beta", parm.names)
pos.sigma <- grep("sigma", parm.names)
PGF <- function(Data) {
  beta <- rnorm(Data$J)
  sigma <- runif(1)
  return(c(beta, sigma))
}
MyData <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
               parm.names=parm.names, pos.beta=pos.beta, pos.sigma=pos.sigma, y=y)
### Reversible-Jump Specifications
bin.n <- J-1 #Maximum allowable model size
bin.p <- 0.4 #Most probable size:  bin.p x bin.n is binomial mean and median
parm.p <- rep(1/J,J+1)
selectable=c(0, rep(1,J-1), 0)

Model <- function(parm, Data)
{
  ### Parameters
  beta <- parm[Data$pos.beta]
  sigma <- interval(parm[Data$pos.sigma], 1e-100, Inf)
  parm[Data$pos.sigma] <- sigma
  ### Log-Prior
  beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
  sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
  ### Log-Likelihood
  mu <- tcrossprod(Data$X, t(beta))
  LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
  ### Log-Posterior
  LP <- LL + beta.prior + sigma.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=LP,
                   yhat=rnorm(length(mu), mu, sigma), parm=parm)
  return(Modelout)
}

Initial.Values <- GIV(Model, MyData, PGF=TRUE)


Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values = Initial.Values,
     Covar=NULL, Iterations=20000, Status=100, Thinning=100,
     Algorithm="RJ", Specs=list(bin.n=bin.n, bin.p=bin.p,
          parm.p=parm.p, selectable=selectable,
          selected=c(0,rep(1,J-1),0)),parm.names = parm.names)


check_post <- as.data.frame(Fit$Posterior1) %>%
  gather('var','value') %>%
  ggplot(aes(value,fill = var)) +
  scale_fill_discrete(guide = F) +
  geom_histogram() +
  facet_wrap(~var)