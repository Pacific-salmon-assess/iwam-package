# Proving the relationship between logalpha and logWA

logalpha <- 1.5 # 1.5

beta <- 0.0001


Srep <- logalpha/beta

Smsy <- (1-gsl::lambert_W0(exp(1-logalpha)))/beta

a <- Smsy/Srep

logalpha.hat <- (a*(-gsl::lambert_W0(((a-1)*exp(1-1/a))/a)) + a - 1)/((a-1)*a)
  # back-simulated alpha based on predicted Srep and Smsy

beta.hat <- logalpha.hat/Srep
  # back-simulated beta

beta.hat - beta
  # back-simulated beta

logalpha.hat - logalpha
  # the difference between back-simulated alpha and alpha (Ricker from SR curve)
  # functionally the same - so we can recover our estimate of alpha from the SR curve
  # and more importantly the WA regression

  # can I simulate variable alpha's?
    # how much does the difference between recoverable alpha's change depending on the prior?


predictAlpha <- function(logWA = 4:12, b0 = c(0,0), b1 = c(0,0)){
  resalpha <- resbeta <- NULL
  for( x in logWA ){
    Srep <- b0[1] + b1[1]*x
    Smsy <- b0[2] + b1[2]*x
    
    a <- exp(Smsy-Srep)
    logalpha.hat <- (a*(-gsl::lambert_W0(((a-1)*exp(1-1/a))/a)) + a - 1)/((a-1)*a)
    beta.hat <- logalpha.hat/Srep
    resalpha <- c(resalpha, logalpha.hat)
    resbeta <- c(resbeta, beta.hat)
  }
  data.frame(logWA = logWA, logalpha=resalpha, beta=resbeta)
}

ans <- predictAlpha(logWA = 1:20, b0 = c(4.18, 3.21), b1 = c(exp(-0.409), exp(-0.411)))
plot(ans$logWA, ans$logalpha)
lm(ans$logalpha~ans$logWA)
