# No libraries - this is sample code


# New Sgen Solver

sGenSolver = function(loga, b){
  sMSY <- (1-gsl::lambert_W0(exp(1-loga)))/b
  fn <- function(S) { abs(log(sMSY) - log(S) - loga + b*S) }
  optimize(fn, interval = c(0, Smsy))
}

sGenOptimum <- function ( S, theta ) {
  # Function called from sGenSolver 
  loga <- theta[1]
  b <- theta[2]
  prt <- S * exp( loga - b * S)
  sMSY <- ( 1 - gsl::lambert_W0 (exp ( 1 - loga) ) ) / b
  epsilon <- log(sMSY) - log(prt)
  nlogLike <- - sum( dnorm ( epsilon, 0, 1, log = T))
  return( nlogLike )
}
sGenSolver_old <- function (loga, b) {
  # Function to estimate Sgen from loga and b Ricker parameters
  theta <- c(loga, b)
  sMSY <- (1 - gsl::lambert_W0(exp(1 - loga))) / b
  fit <- optimize(f = sGenOptimum, interval = c(0, sMSY),
                  theta = theta)
  return(fit$minimum)
}

sGenDirect <- function(loga, b){
  sMSY <- ( 1 - gsl::lambert_W0 (exp ( 1 - loga) ) ) / b
  a <- exp(loga)
  -1/b*gsl::lambert_W0(-b*sMSY/a)
}

microbenchmark::microbenchmark(
  new = sGenSolver(loga = 2, b = 0.0001),
  old = sGenSolver_old(loga = 2, b = 0.0001),
  lamw = sGenDirect(loga = 2, b = 0.0001))