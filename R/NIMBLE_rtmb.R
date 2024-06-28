# NIMBLE --> RTMB Version

# Libaries ####
library(RTMB)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidybayes)
library(tmbstan) # loads rstan and StanHeaders
# library(rstan)

# Data ####

srdatwna <- read.csv(here::here("DataIn/SRinputfile.csv"))
WAbase <- read.csv(here::here("DataIn/WatershedArea.csv"))

# Target's
WAin <- read.csv(here::here("DataIn/WCVIStocks.csv"))
# WAin$WA is raw watershed areas
  # Contains CU and Inlet information - but number of rows is equal to number of stocks/rivers
mean_logtarWA <- mean(log(WAin$WA))
WAin$logtarWAshifted <- log(WAin$WA) - mean_logtarWA


# Data setup ####
srdatwna <- srdatwna %>% 
  filter(!Name %in% c("Hoko","Hoh")) 

# Remove years with NAs and renumber.
srdat <- srdatwna %>% 
  filter(!Name %in% c("Hoko","Hoh")) %>% 
  filter(Rec != "NA") %>%
  filter( !(Name == "Cowichan" & (Yr < 1985 | Yr == 1986 | Yr == 1987))) %>%
  group_by(Name, Stocknumber, Stream) %>%
  arrange(Yr) %>%
  mutate(yr_num = 0:(n()-1)) %>%
  ungroup() %>%
  arrange(Stocknumber) %>%
  mutate(lh = factor(ifelse(Stream == 0, "stream", "ocean"), levels = c("stream", "ocean")))

names <- srdat %>% 
  dplyr::select (Stocknumber, Name, lh) %>% 
  distinct()

WAbase <- WAbase %>% 
  full_join(names, by="Name") %>% 
  arrange(Stocknumber) %>%
  mutate(logWA = log(WA)) 

## Shift log WA for the mean.
mean_logWA <- mean(WAbase$logWA)
WAbase$logWAshifted <- WAbase$logWA - mean_logWA

lifehist <- srdat %>% dplyr::select(Stocknumber, Name, Stream) %>% 
  group_by(Stocknumber) %>% 
  summarize(lh=max(Stream))

## RTMB dat and par setup ####

dat <- list(srdat = srdat,
            WAbase = WAbase,
            WAin = WAin,
            # type = lifehist$lh + 1,
            logRS = log(srdat$Rec) - log(srdat$Sp))
  # what about WAbase? Do they just need to be joined?

# example for below's model
# par <- list(logLinf, logk, logsigma) # parameter values
par <- list(logAlpha_re = numeric(nrow(dat$WAbase)),
            #logAlpha = numeric(nrow(dat$WAbase)),
            logAlpha0 = 1.5,
            logAlphaSD = 10,
            logE0 = numeric(nrow(dat$WAbase)), # random
            logESD = 1,
            tauobs = 0.01 + numeric(nrow(dat$WAbase)), # Why is this 0.01 + 
            beta = c(10,0,0,0) # 4 beta's for inprod? - ocean, stream - a and b
            ) 

dat$X <- model.matrix( ~ lh*logWAshifted, data = dat$WAbase)
dat$tarX <- model.matrix( ~ lh*logtarWAshifted, data = dat$WAin) # I think these ^ are now the same in form and function

# RTMB function ####
# This is a fully MLE model parameterized for SREP (E)

f_nim <- function(par){
  getAll(dat, par)

  N_Stk = max(srdat$Stocknumber + 1)
  stk = srdat$Stocknumber + 1
  N_Obs = nrow(srdat)
  S = srdat$Sp
  N_Pred = nrow(WAin)
  
  # betaPriorMean = c(10,0,0,0)
  nbeta <- ncol(X)
  nbetatar <- ncol(tarX)
  
  # for(i in 1:nbeta) {
  #   beta[i] %~% dnorm(betaPriorMean[i], sd = sqrt(1/0.001))
  # }
  
  E <- numeric(N_Stk)
  tarE <- numeric(N_Pred)
  logAlpha <- numeric(N_Stk)
  
  #### Watershed Model
  for ( i in 1:N_Stk){
    logAlpha_re[i] %~% dnorm(0, sd = logAlphaSD) # random effect
    logAlpha[i] <- logAlpha0 + logAlpha_re[i]
    # logAlpha0 %~% dnorm(0.6, sd = 0.45) # Isn't this how it should be written?
      
    logE0[i] %~% dnorm(mean = 0, sd = logESD) # random effect
  
    # log(E[i]) <- a + b*dat[[2]]$logWAshifted + logE0[i] # **************************
    # log(E[i]) <- inprod(beta[1:nbeta], X[i,1:nbeta]) + logE0[i] # nimble function
      # sum(v1 * v2)
    logE <- sum(beta[1:nbeta] * X[i,1:nbeta]) + logE0[i]
      # Tor: should this be -sum?
    E[i] <- exp(logE)
    
    # tauobs[i] %~% dgamma(0.001, 0.001)
  }
  
  #### Ricker Model
  for (i in 1:N_Obs){
    logRS_pred <- logAlpha[stk[i]]*(1 - S[i]/E[stk[i]])
    logRS[i] %~% dnorm(logRS_pred, sd = sqrt(1/tauobs[stk[i]]))
  }
  
  #### Predictions of E for WCVI watershed areas
  for (i in 1:N_Pred){
    # create a predicted line
    # use the logE values to create a prediction
    # separate them by stream and ocean?
    # can I just evaluate a new tarlogE using a new value of X?
    logtarE <- sum(beta[1:nbeta] * tarX[i,1:nbeta]) + logE0[i]
    tarE[i] <- exp(logtarE)
    
    # How does the indexing work on the below calculations?
    # tarBeta[i] <- logAlpha[stk[i]]/tarE[i]
    # tarSMSY[i] <- (1-LambertW0(exp(1-logAlpha[stk[i]])))/tarBeta[i]
    # tarSMAX[i] <- 1/tarBeta[i]
  }
  
  # ADREPORT # gives standard errors
  ADREPORT(logRS)
  ADREPORT(E)
  ADREPORT(logAlpha)
    # Between Srep and alpha I can get beta
    # TK: Is it worth doing any calculations internally? Or leave that space for the targets?
  
  ADREPORT(tarE)
  # ADREPORT(tarBeta)
  # ADREPORT(tarSMSY)
  # ADREPORT(tarSMAX)
  
  # REPORT # just gives you the error
  
  # Create predictions
    # create a vector of: seq(min(log(WAbase$WA)), max(log(WAbase$WA)), 0.1)
    # to use to calculate a predicted line
  
  # -sum(dnorm(log(RS), logRS_pred, sigma = sqrt(1/tauobs[stk[i]]), log=TRUE))
}



# New RTMB Function - testing old method without matrix ####

# Notes from Paul
  # Starting with question about uniform distributions in TMB
# You don't need to explicitly ever include a uniform distribution unless the end points are random.
# It's just 1/(upper-lower) - the log its just -log(upper-lower)

# You aren't bayesian so you don't put a distribution on a parameter in RTMB
# But the uniform is important to know the values the parameter can take on.
# So if you had an untransformed variable then you'd want to include it to make sure it couldn't find an impossible value.

# Q: I would have to state that those 4 beta parameters in rtmb are initial parameters that are then constrained?
# No I think that's a linear model so the beta parameters don't need to be transformed.
# In a Frequentist world you just have them as par and don't add them to the likelihood.

# Well you can kind of think of all your pars being "uniformly" distributed and don't contribute to the likelihood. 
# When this is the case the posterior mode and MLE are the same.
# **When all priors are uniform then MLE and MAP are equal.

# The above notes are disjointed - but I wanted them for safe-keeping (Tor)

N_Stk <- max(srdat$Stocknumber + 1)

# parameters
par2 <- list(b0 = c(10, 10), # b0 = c(9, 9)
       bWA = c(0, 0), # bWA = c(0.83, 1)
       
       # logAlpha = numeric(N_Stk),
       logAlpha_re = numeric(nrow(dat$WAbase)),
       logE0 = numeric(N_Stk),
       tauobs = 0.01 + numeric(N_Stk), # Why can't this be zero? This doesn't run as just a string of zeros.
       logAlpha0 = 1.5, # Same
       logESD = 1, # Same
       logAlphaSD = 10 # 1 in the nimble version - no effect changing between them
       )

f_nim2 <- function(par2){
  getAll(dat, par2)
   
  N_Stk = max(srdat$Stocknumber + 1)
  stk = srdat$Stocknumber + 1
  N_Obs = nrow(srdat)
  S = srdat$Sp
  type = lifehist$lh + 1
  
  E <- numeric(N_Stk)
  log_E <- numeric(N_Stk) 
  
  logAlpha <- numeric(N_Stk) # comment on or off
  
  ## Initialize joint negative log likelihood
  nll <- 0
  
  # Prior - Now a penalty term?
  # nll <- nll - sum(dnorm(logAlpha0, 0.6, 0.45, log = TRUE))
  
  # Alternative version of alpha terms - see f_nim's translation
  
  # Slope and intercept priors
  # nll <- nll - sum(dnorm(b0[1], 0, sd = 31.6, log = TRUE))
  # nll <- nll - sum(dnorm(b0[2], 0, sd = 31.6, log = TRUE))
  # nll <- nll - sum(dnorm(bWA[1], 0, sd = 31.6, log = TRUE))
  # nll <- nll - sum(dnorm(bWA[2], 0, sd = 31.6, log = TRUE))
  
  ## Watershed Model
  for ( i in 1:N_Stk){
    # New version
    # nll <- nll - sum(dnorm(logAlpha[i], logAlpha0, sd = logAlphaSD, log = TRUE)) # random effect - is this the bayesian way?
    # Alternative version of alpha terms - see f_nim's translation
    nll <- nll - sum(dnorm(logAlpha_re[i], 0, sd = logAlphaSD, log = TRUE)) # random effect
    logAlpha[i] <- logAlpha0 + logAlpha_re[i]
    
    # logE0[i] %~% dnorm(mean = 0, sd = logESD) # random effect
    nll <- nll - sum(dnorm(logE0[i], 0, sd = logESD, log = TRUE)) # random effect
    
    # logE <- sum(beta[1:nbeta] * X[i,1:nbeta]) + logE0[i]
    log_E[i] <- b0[type[i]] + bWA[type[i]]*WAbase$logWAshifted[i] + logE0[i] ## Stock level regression
    E[i] <- exp(log_E[i])
    
    # tauobs[i] ~ dgamma(0.001, 0.001) # stock specific precision
    nll <- nll - sum(dgamma(tauobs[i], shape = 0.001, scale = 0.001))
  }
  # After this point nll is -Inf
  
  ## Ricker Model
  for (i in 1:N_Obs){
    logRS_pred <- logAlpha[stk[i]]*(1 - S[i]/E[stk[i]])
    nll <- nll - sum(dnorm(logRS[i], logRS_pred, sd = sqrt(1/tauobs[stk[i]]), log = TRUE))
    # logRS[i] %~% dnorm(logRS_pred, sd = sqrt(1/tauobs[stk[i]]))
  }

  nll
}

## MakeADFun ####
obj <- RTMB::MakeADFun(f_nim, par, random = c("logAlpha_re", "logE0")) # silent=TRUE
obj2 <- RTMB::MakeADFun(f_nim2, par2, random = c("logAlpha_re", "logE0"))

opt <- nlminb(obj$par, obj$fn, obj$gr) # THIS HAS TO BE RUN for stan?
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)

sdr <- sdreport(obj)
sdr2 <- sdreport(obj2)

sdr_est <- as.list(sdr, "Est", report=TRUE) ## ADREPORT estimates
sdr_se <- as.list(sdr, "Std", report=TRUE) ## ADREPORT standard erro

# HOW DO YOU GET OUT ALL THE NECESSARY things?
  # ADREPORT vs. REPORT
  # CI's included - which do I want them for?
  # Make a list of everything NIMBLE ouputs and replicate them

# NIMBLE_2 OUTPUTS
  # Calculate reference points - POSTERIOR PREDICTIVE
    # logalpha
    # srep
    # beta
    # smax
    # smsy
    # sgen
  # summarize reference points - ABOVE
    # ? what levels exist that need summarizing?
  # create prediction design matrix
  # predict new sites
  # make a prediction line/matrix
  # get new predictions

# Can I use the same extraction method as with the NIMBLE model?
  # are the output objects the same?
  # DIFFERENT between opt and tmbstan definitely
    # one will be posterior mean estimates (max L) and other mcmc chains



# TMB STAN ####
# ?tmbstan()
# https://github.com/kaskr/tmbstan

# below are options for parallel running WITH AN INIT FUNCTION
  # try running without a INIT
# cores <- parallel::detectCores()-1 # get number of cores
# options(mc.cores = cores) # set options for cores

# fitfirst <- tmbstan(obj, chains=1) # needs initializing values

# Tor: Weird that this is just an exact replica of par
initf1 <- function(){
  list(logAlpha_re = numeric(nrow(dat$WAbase)), # 1.5 + 
       #logAlpha = numeric(nrow(dat$WAbase)),
       logAlpha0 = 1.5,
       logAlphaSD = 10,
       logE0 = numeric(nrow(dat$WAbase)), # random
       logESD = 1,
       tauobs = 0.01 + numeric(nrow(dat$WAbase)),
       beta = c(10,0,0,0) # 4 beta's for inprod? - ocean, stream - a and b
  )
}

initf2 <- function(){
  list(b0 = c(10, 10), # b0 = c(9, 9)
       bWA = c(0, 0), # bWA = c(0.83, 1)
       logAlpha = numeric(N_Stk),
       logE0 = numeric(N_Stk),
       tauobs = 0.001 + numeric(N_Stk), # Why can't this be zero? This doesn't run as just a string of zeros.
       logAlpha0 = 1.5, # Same
       logESD = 1, # Same
       logAlphaSD = 1 # Same
  )
}

# fit <- tmbstan(obj, iter=2000, warmup=200, init=initf1,
#                chains=1, open_progress=FALSE, silent=TRUE)

fitcores <- tmbstan(obj, iter=2000, warmup=200, init=initf1,
               chains=4, open_progress=FALSE, silent=TRUE)
  # errors include: can't find getAll function
  # can't find dat

fitcores2 <- tmbstan(obj2, iter=2000, warmup=200, init=initf2,
                    chains=4, open_progress=FALSE, silent=TRUE)


# Alternative version with multiple chains
  # trying to match https://mc-stan.org/rstan/reference/stan.html#examples
# initf2 <- function(chain_id=1){
#   list(logAlpha = 1.5 + numeric(nrow(dat$WAbase)),
#        logAlpha0 = 1.5,
#        logAlphaSD = 10,
#        logE0 = numeric(nrow(dat$WAbase)),
#        logESD = 1,
#        tauobs = 0.01 + numeric(nrow(dat$WAbase)),
#        beta = c(10,0,0,0)
#   )
# }
# n_chains <- 4
# init_ll <- lapply(1:n_chains, function(id) initf2(chain_id=id))
# 
# fit2 <- tmbstan(obj, iter=2000, warmup=200, init=init_ll,
#                 chains=n_chains, open_progress=FALSE, silent=TRUE)


# Examine the pairs() plot to diagnose sampling problems
# pairs(fit, pars=names(obj$par)) # this is a massive plot
# ggpairs()

rstan::traceplot(fitcores2, pars=names(obj2$par), inc_warmup=TRUE)
# stan_trace(fit, pars=names(obj$par), inc_warmup=TRUE) # the same? just a ggplot object?
#TK: logAlpha trace is shit

library(shinystan)
launch_shinystan(fitcores)
