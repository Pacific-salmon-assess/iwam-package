# Lambert W Internal functions

# FritschIter ####
FritschIter <- function(x, w){
  MaxEval <- 5
  CONVERGED <- FALSE
  k <- 2.0 / 3.0;
  i <- 0;
  eps <- 2.2204460492503131e-16    
  while (!CONVERGED & i < MaxEval){
    z <- log(x / w) - w
    w1 <- w + 1.0
    q <- 2.0 * w1 * (w1 + k * z)
    qmz <- q - z
    e <- z / w1 * qmz / (qmz - z)
    CONVERGED <- abs(e) <= eps
    w <- w*(1.0 + e)
    i <- i + 1
  }
  return(w)
}

# LambertW0_internal (Paul) ####
LambertW0_internal <- function(x){
  check <- 0.367879441171442334024277442949824035167694091796875 # exp(-1)
  eps <- 2.2204460492503131e-16
  if (x == Inf) {
    return(Inf);
  } else if (x < -check) {
    return(NaN);
  } else if (abs(x - check) <= eps) {
    return(-1.0);
  } else if (abs(x) <= 1e-16) {
    ## This close to 0 the W_0 branch is best estimated by its Taylor/Pade
    ## expansion whose first term is the value x and remaining terms are below
    ## machine double precision. See
    ## https://math.stackexchange.com/questions/1700919
    
    return(x);
  } else {
    w <- 0
    if (abs(x) <= 6.4e-3) {
      ## When this close to 0 the Fritsch iteration may underflow. Instead,
      ## function will use degree-6 minimax polynomial approximation of Halley
      ## iteration-based values. Should be more accurate by three orders of
      ## magnitude than Fritsch's equation (5) in this range.
      return((((((-1.0805085529250425e1 * x + 5.2100070265741278) * x -
                   2.6666665063383532) * x + 1.4999999657268301) * x -
                 1.0000000000016802) * x + 1.0000000000001752) * x +
               2.6020852139652106e-18);
      
    } else if (x <= exp(1)) {
      ## Use expansion in Corliss 4.22 to create (2, 2) Pade approximant.
      ## Equation with a few extra terms is:
      ## -1 + p - 1/3p^2 + 11/72p^3 - 43/540p^4 + 689453/8398080p^4 - O(p^5)
      ## This is just used to estimate a good starting point for the Fritsch
      ## iteration process itself.
      
      p <- sqrt(2.0 * (exp(1) * x + 1.0))
      Numer <- (0.2787037037037037 * p + 0.311111111111111) * p - 1.0;
      Denom <- (0.0768518518518518 * p + 0.688888888888889) * p + 1.0;
      w <- Numer / Denom;
    } else {
      ## Use first five terms of Corliss et al. 4.19 */
      w <- log(x)
      L_2 <- log(w)
      L_3 <- L_2 / w
      L_3_sq <- L_3 * L_3
      w <- w - L_2 + L_3 + 0.5 * L_3_sq - L_3 / w + L_3 / (w * w) - 1.5 * L_3_sq /
        w + L_3_sq * L_3 / 3.0;
    }
    return(FritschIter(x, w));
  }
} 

# dLambertW0_internal (Paul) ####
## Derivatives of LamW
dLambertW0_internal <- function(x, y, dy) {
  dy / (exp(y) * (1. + y))
}

# LambertW_internal (Original TMB - from Kasper) ####
# Code taken from https://kaskr.github.io/adcomp/lambert_8cpp_source.html
LambertW_internal <- function(x) {
  logx <- log(x)
  y <- ifelse(logx > 0, logx, 0)
  niter <- 101
  
  for (i in 1:niter){
    if (abs(logx - log(y) - y) < 1e-9) break
    y <- y - (y - exp(logx - y)) / (1 + y)
  }
  
  if (i == niter) warning("W: failed convergence")
  
  return(y)
}

# dLambertW_internal (TMB) ####
dLambertW_internal <- function(x, y, dy){
  # W <- LambertW(x)
  # DW <- 1 / (exp(W) * (1 + W))
  # list(value = W, gradient = DW)
  DW <- 1 / (exp(y) * (1 + y))
  return(DW * dy)
}

