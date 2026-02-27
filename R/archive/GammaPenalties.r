library(RTMB)

## How penalties currenty are included.
gamma_on_tau <- function(pars){
  getAll(pars)
  sigma <- exp(logsigma)
  -dgamma(sigma^-2, shape = 0.1, scale = 1/0.1, log = TRUE)
}

obj <- MakeADFun(gamma_on_tau, parameters = list(logsigma = 0))
fit <- nlminb(0, obj$fn, obj$gr)
exp(fit$par)  ## Huge and I suspect not what you intend to do.

## Doing the same thing but just on variance.
gamma_on_var <- function(pars){
  getAll(pars)
  sigma <- exp(logsigma)
  -dgamma(sigma^2, shape = 0.1, scale = 1/0.1, log = TRUE)
}

obj2 <- MakeADFun(gamma_on_var, parameters = list(logsigma = 0))
fit2 <- nlminb(0, obj2$fn, obj2$gr)
exp(fit2$par) ## Near zero.

## What you are wanting to do based on the plot. An inverse gamma on variance.
invgamma_on_var <- function(pars){
  getAll(pars)
  sigma <- exp(logsigma)
  -dgamma(1/sigma^2, shape = 0.1, scale = 1/0.1, log = TRUE) + 2*log(sigma^2)
}

obj3 <- MakeADFun(invgamma_on_var, parameters = list(logsigma = 0))
fit3 <- nlminb(0, obj3$fn, obj3$gr)
exp(fit3$par) ## Looks much more sensible.

## Or because you can just penalize sigma directly...
## Choose parameters that relate to the mean and variance of the observed sigma.
gamma_on_sigma <- function(pars){
  getAll(pars)
  sigma <- exp(logsigma)
  -dgamma(sigma, shape = 2, scale = 1/2, log = TRUE)
}

obj4 <- MakeADFun(gamma_on_sigma, parameters = list(logsigma = 0))
fit4 <- nlminb(0, obj4$fn, obj4$gr)
exp(fit4$par) ## Looks sensible. 
