# Liermann Srep (E) RTMB Model with MCMC Sampling ####

# Libaries ####
library(RTMB)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(progress)
library(tmbstan)
library(TMB)
library(tidybayes)
library(bayesplot)
library(beepr) # Sounds
library(viridis)
library(latex2exp)

source(here::here("R/LambertWs.R")) # Lambert W function
source(here::here("R/helperFunctions.R")) # For bootstrapping
source(here::here("R/derived_post.R")) # For posterior extraction

# Raw data read-in ####
WAin <- c("DataIn/Parken_evalstocks.csv")
# WAin <- c("DataIn/WCVIStocks.csv")

# For re-evaluation of synoptic sets e.g. WAbase
    # Run model until setting up data section
    # Then over-write WAin <- WAbase
    # And rename logWAshifted to logWAshifted_t
    # And make sure to change type_tar to fit the expected 0:1 indexing

# New LambertW0
# See: https://github.com/pbs-assess/renewassess/blob/main/code/RTMB/PosteriorPredictiveSample.r
LambertW0 <- ADjoint(
  function(x){gsl::lambert_W0(x)},
  function(x, y, dy) {dy / (x + exp(y))}
)

# Data Manipulations ####
srdatwna <- read.csv(here::here("DataIn/SRinputfile.csv"))
WAbase <- read.csv(here::here("DataIn/WatershedArea.csv"))
WAin <- read.csv(here::here(WAin))

# Data setup ####
# Remove Hoko and Hoh stocks - consider removing SKagit OR Chehalis
# Due to sampling reasons explained in Parken 2006.
srdatwna <- srdatwna %>% 
  filter(!Name %in% c("Hoko","Hoh")) # |>
  # (!(Name == "Skagit")) # |> # Skagit is #22
  # filter( !(Name == "Chehalis")) # |> # Chehalis is #19

# Remove years with NAs and renumber.
srdat <- srdatwna %>% 
  filter(Rec != "NA") %>%
  filter( !(Name == "Cowichan" & (Yr < 1985 | Yr == 1986 | Yr == 1987))) %>%
  group_by(Name, Stocknumber, Stream) %>%
  arrange(Yr) %>%
  mutate(yr_num = 0:(n()-1)) %>%
  ungroup() %>%
  arrange(Stocknumber) %>%
  mutate(lh = factor(ifelse(Stream == 0, "stream", "ocean"), # Stream = 0, Ocean  = 1
    levels = c("stream", "ocean"))) |> 
  mutate(Stocknumber = as.integer(factor(Stocknumber)) - 1)  # Re-numbering uniquely

names <- srdat %>% 
  dplyr::select (Stocknumber, Name, lh) %>% 
  distinct()

WAbase <- WAbase %>% 
  full_join(names, by="Name") %>% 
  arrange(Stocknumber) %>%
  mutate(logWA = log(WA)) |> 
  filter(!is.na(Stocknumber))

## Shift log WA for the mean - base - makes estimation easier
mean_logWA <- mean(WAbase$logWA)
WAbase$logWAshifted <- WAbase$logWA - mean_logWA

WAin$logWA <- log(WAin$WA)
WAin$logWAshifted_t <- WAin$logWA - mean_logWA

lifehist <- srdat %>% dplyr::select(Stocknumber, Name, Stream) %>% 
  group_by(Stocknumber) %>% 
  summarize(lh=max(Stream))

## RTMB dat and par setup ####
# Dat
dat <- list(srdat = srdat,
            WAbase = WAbase,
            WAin = WAin,
            lineWA =  seq(min(WAbase$logWAshifted), max(WAbase$logWAshifted), 0.1), # Not added to likelihood
            logRS = log(srdat$Rec) - log(srdat$Sp),
            prioronly = 1) # 0 - run with data, 1 - prior prediction mode?

# External vectors
N_Stk <- max(srdat$Stocknumber + 1)
N_Obs <- nrow(srdat)

# NEW: alpha0 prior for LH specific dists. - ON/OFF
lhdiston <- T # Set true to run

# Parameters/Initial values
par <- list(b0 = c(10, 0), # Initial values for WA regression intercepts
            bWA = c(0, 0), # Inital values for WA regression slopes
            # logRS_pred = numeric(nrow(srdat)), # Zeroes - testing as a parameter
            logE_re = numeric(N_Stk), # Zeroes
            logAlpha0 = 0.6,
            logAlpha_re = numeric(nrow(dat$WAbase)), # Zeroes
            tauobs = 0.01 + numeric(N_Stk), # Constrained positive
            logESD = 1,
            logAlphaSD = 1
)
if (lhdiston) {
  par$logAlpha02 <- 0
}

f_srep <- function(par){
  getAll(dat, par)

  # logRS <- OBS(logRS) # Mark response for one-step-ahead residual calculation
    # Not currently used - with oneStepPredict()
  
  N_Stk = max(srdat$Stocknumber + 1) # number of stocks
  stk = srdat$Stocknumber + 1 # vector of stocknumbers
  N_Obs = nrow(srdat) # number of observations
  N_Pred = nrow(WAin) # number of predicted watershed areas
  
  S = srdat$Sp
  type = lifehist$lh
  type_tar = as.numeric(WAin$lh) 

  E <- numeric(N_Stk)
  logE_pred <- numeric(N_Stk)
  logE <- numeric(N_Stk)
  logAlpha <- numeric(N_Stk)
  
  # Why is logRS_pred not a parameter or vector input here?
  logRS_pred <- numeric(N_Obs) # Does this still report if not a vector?
  # logRS_sim <- numeric(N_Obs) # Posterior predictive simulated vector?
  
  E_tar <- numeric(N_Pred)
  logE_tar <- numeric(N_Pred)
  
  logAlpha_tar <- numeric(N_Pred)
  
  # Simulated line vectors
  line <- length(lineWA)
  logE_line_stream <- numeric(line)
  E_line_stream <- numeric(line)
  logE_line_ocean <- numeric(line)
  E_line_ocean <- numeric(line)
  
  nll <- 0 # Begin negative log-likelihood
  
  nll <- nll - sum(dnorm(b0[1], 10, sd = 31.6, log = TRUE)) # Prior
  nll <- nll - sum(dnorm(b0[2], 0, sd = 31.6, log = TRUE)) # Prior
  nll <- nll - sum(dnorm(bWA[1], 0, sd = 31.6, log = TRUE)) # Prior
  nll <- nll - sum(dnorm(bWA[2], 0, sd = 31.6, log = TRUE)) # Prior
  
  nll <- nll - sum(dnorm(logAlpha0, 0.6, sd = 0.45, log = TRUE)) # Prior (rM)
  if(lhdiston) nll <- nll - sum(dnorm(logAlpha02, 0, sd = 31.6, log = TRUE)) # Prior (rD)
  
  ## Second level of hierarchy - Ricker parameters:
  for (i in 1:N_Stk){
    nll <- nll - dnorm(logE_re[i], 0, sd = logESD, log = TRUE) # Unobserved
    logE[i] <- b0[1] + b0[2]*type[i] + (bWA[1] + bWA[2]*type[i]) * WAbase$logWAshifted[i] + logE_re[i] # DATA FLAG
    E[i] <- exp(logE[i])
    
    nll <- nll - dnorm(logAlpha_re[i], 0, sd  = logAlphaSD, log = TRUE) # Prior (hj)
    
    if(lhdiston) logAlpha[i] <- logAlpha0 + logAlpha02*type[i] + logAlpha_re[i] # LH specific distributions for prod.
    else logAlpha[i] <- logAlpha0 + logAlpha_re[i]

    nll <- nll - dgamma(tauobs[i], shape = 0.0001, scale = 1/0.0001, log = TRUE)
  }

  ## First level of hierarchy - Ricker model:
  for (i in 1:N_Obs){
    logRS_pred[i] <- logAlpha[stk[i]]*(1 - S[i]/E[stk[i]]) # DATA FLAG
    # nll <- nll - dnorm(logRS[i], logRS_pred[i], sd = sqrt(1/tauobs[stk[i]]), log = TRUE) # DATA FLAG

    if(!prioronly){
      nll <- nll - dnorm(logRS[i], logRS_pred[i], sd = sqrt(1/tauobs[stk[i]]), log = TRUE)
    } # If prioronly is 0, then likelihood is calculated
    # If prioronly is 1, then likelihood is not calculated
    
    # posterior predictive check
    # logRS_sim[i] <- rnorm(1, logRS_pred[i], sd = sqrt(1/tauobs[stk[i]])) # Simulated logRS for posterior predictive check
  }
  
  ## Calculate SMSY for Synoptic set - for plotting
  SMSY_r = numeric(nrow(WAbase))
  BETA_r = numeric(nrow(WAbase))
  
  for (i in 1:N_Stk){
    BETA_r[i] <- logAlpha[i] / E[i]
    SMSY_r[i] <- (1 - LambertW0(exp(1 - logAlpha[i]))) / BETA_r[i]
  }

  ## PREDICTIONS
  BETA = numeric(nrow(WAin))
  SMSY = numeric(nrow(WAin))
  SGEN = numeric(nrow(WAin))

  # Conditional posterior predictions
  for (i in 1:N_Pred){
    if(lhdiston) logAlpha_tar[i] <- logAlpha0 + logAlpha02*type_tar[i] # NEW: alpha0 prior for LH specific dists.
    else logAlpha_tar[i] <- logAlpha0

    logE_tar[i] <- b0[1] + b0[2]*type_tar[i] + (bWA[1] + bWA[2]*type_tar[i])*WAin$logWAshifted_t[i]
    E_tar[i] <- exp(logE_tar[i])
    
    # Predict BETA
    BETA[i] <- logAlpha_tar[i]/E_tar[i]
    # Predict SMSY
    SMSY[i] <- (1-LambertW0(exp(1-logAlpha_tar[i])))/BETA[i]
    # Predict SGEN
    SGEN[i] <- -1/BETA[i]*LambertW0(-BETA[i]*SMSY[i]/(exp(logAlpha_tar[i])))
  }
  
  # Create predictions on an simulated line
  for (i in 1:line){
    logE_line_ocean[i] <- b0[1] + b0[2] + (bWA[1] + bWA[2])*lineWA[i] # DATA - not in likelihood
    E_line_ocean[i] <- exp(logE_line_ocean[i])
    
    logE_line_stream[i] <- b0[1] + (bWA[1])*lineWA[i] # DATA - not in likelihood
    E_line_stream[i] <- exp(logE_line_stream[i])
  }
  
  ## ADREPORT - internal values (synoptic specific/Ricker)
  # ADREPORT(logRS) # logRS for all 501 data points
  alpha <- exp(logAlpha)
  
  REPORT(b0) # Testing simulate()
  REPORT(bWA) # Testing simulate()
  
  # ADREPORT(logRS_pred)
  REPORT(logRS_pred)
  # REPORT(logRS_sim) # Is logRS_pred as a REPORT the same as logRS_sim? (Question for Sean/Paul)
  
  ADREPORT(logE_re)
  ADREPORT(E)
  ADREPORT(logAlpha)
  ADREPORT(logAlpha_re) # random effect parameter for resampling
  ADREPORT(alpha)
  ADREPORT(SMSY_r)
  ADREPORT(BETA_r)
  ADREPORT(tauobs)
  
  # REPORT(logRS) # logRS for all 501 data points
  REPORT(logE_re)
  REPORT(E) # E (Srep) for all synoptic data set rivers (25)
  REPORT(logAlpha) # model logAlpha (25)
  REPORT(logAlpha_re) # random effect parameter for resampling
  REPORT(alpha)
  REPORT(SMSY_r)
  REPORT(BETA_r)
  REPORT(tauobs)
  
  # ADREPORT - predicted values from watershed area model
    # Mean estimate of the median (without bias correction)
  alpha_tar <- exp(logAlpha_tar)
  
  ADREPORT(E_tar) # target E (Srep) (21)
  ADREPORT(logE_tar) # exp these for the correct confidence intervals
  ADREPORT(logAlpha_tar)
  ADREPORT(alpha_tar)
  
  REPORT(E_tar)
  REPORT(logE_tar)
  REPORT(logAlpha_tar)
  REPORT(alpha_tar)

  ADREPORT(BETA)
  ADREPORT(SMSY)
  ADREPORT(SGEN)
  
  REPORT(BETA)
  REPORT(SMSY)
  REPORT(SGEN)
  
  # Simulated line values for plotting
  REPORT(E_line_stream) 
  ADREPORT(E_line_stream)
  REPORT(E_line_ocean) 
  ADREPORT(E_line_ocean)
  
  REPORT(logAlphaSD)
  REPORT(logESD)
  
  nll # output of negative log-likelihood
}

## MakeADFun ####
obj <- RTMB::MakeADFun(f_srep,
                       par,
                       random = c("logAlpha_re", "logE_re"),
                       silent=TRUE)

# Prior testing ####
    # - Pass random prior values to obj$par - then simulate new observations
    # - Requires changing par each time and running function and obj
    # - Need to write a par function that matches priors?

# priorlogRS <- obj$simulate()$logRS_pred
  # What is the difference between logRS_pred and logRS as an OBS() object?

# priorsim <- obj$simulate()

# Tor:
  # If you re-run function and MakeADFUN after creating random pars
  # e.g. par <- init() - drawn from random init starts for MCMC below
  # When you run obj$simulate() - you will have randomly generated parameters.
# par(mfrow = c(1,3))
# plot(priorsim$logRS_pred)
# plot(priorsim$logRS_pred, ylim = c(-4, 4))
# plot(dat$logRS, ylim = c(-4, 4))
# Add a line to add points() instead of re-creating the plot to see multiple sims

psimalpha <- vector("list", 100)
psimlogRS_pred <- vector("list", 100)
psimlogAlpha_re <- vector("list", 100)

for (i in 1:100){
  # Random parameter creation: Random parameter starts for prior simulation
  parp <- function() {
    
  # Can also add par <- to sequence above - see prior testing
  # listinit <- list(b0 = c(rnorm(1, 10, 1), rnorm(1, 0, 1)), # Contains negatives
  #      bWA = c(rnorm(1, 0, 1), rnorm(1, 0 , 1)), # Contains negatives
  # 
  #      # logRS_pred = rnorm(N_Obs, 0, 1),
  #      logE_re = rnorm(N_Stk, 0, 1), # Contains negatives
  #      logAlpha0 = rnorm(1, 0.6, 1), # Contains negatives
  #      logAlpha02 = rnorm(1, 0, 1) , # NEW: alpha0 prior for LH specific dists.
  #      logAlpha_re = rnorm(nrow(dat$WAbase), 0, 1), # Contains negatives
  # 
  #      tauobs = runif(N_Stk, min = 0.005, max = 0.015), # Uniform to REMAIN positive
  # 
  #      logESD = runif(1, 0.01, 3), # Positive
  #      logAlphaSD = runif(1, 0.01, 3) # Positive
  #   )
  
  logESD <- runif(1, 0, 100)
  logAlphaSD <- runif(1, 0, 100)

  listinitprior <- list(b0 = c(rnorm(1, 10, 31.6), rnorm(1, 0, 31.6)),
    bWA = c(rnorm(1, 0, 31.6), rnorm(1, 0 ,31.6)),
    logESD = logESD, # This isn't being saved internally
    logAlphaSD = logAlphaSD, # This isn't being saved internally
    
    logE_re = rnorm(N_Stk, 0, logESD),
    logAlpha0 = rnorm(1, 0.6, 0.45),
    logAlpha02 = rnorm(1, 0, 31.6),
    logAlpha_re = rnorm(nrow(dat$WAbase), 0, logAlphaSD),
    
    tauobs = rgamma(N_Stk, shape = 0.001, scale = 1/0.001)
    )
  
  # logESDalt <- runif(1, 0, 20)
  # logAlphaSDalt <- runif(1, 0, 20)
  # listinitalternate <- list(b0 = c(rnorm(1, 7, sd = 31.6), rnorm(1, 0, sd = 10)),
  #   bWA = c(rnorm(1, 0, sd = 10), rnorm(1, 0 , sd = 10)),
  #   logESD = logESDalt, # This isn't being saved internally
  #   logAlphaSD = logAlphaSDalt, # This isn't being saved internally
  #   
  #   logE_re = rnorm(N_Stk, 0, sd = logESD),
  #   logAlpha0 = rnorm(1, 0, sd = 0.71),
  #   logAlpha02 = rnorm(1, 0, sd = 31.6),
  #   logAlpha_re = rnorm(nrow(dat$WAbase), 0, sd = logAlphaSD),
  #   
  #   tauobs = rgamma(N_Stk, shape = 0.0001, scale = 1/0.0001)
  #   )
  
    return(listinitprior) # Can change between two versions - listinitprior is a direct representation of the priors
  }
  
  parpar <- parp()
  
  # enter parp into model function with MakeADFun
  objpp <- RTMB::MakeADFun(f_srep,
                       parpar,
                       random = c("logAlpha_re", "logE_re"),
                       silent=TRUE)
  
  # simulate from prior
  psimalpha[[i]] <- objpp$simulate()$alpha
  psimlogAlpha_re[[i]] <- objpp$simulate()$logAlpha_re # why logAlpha_re? 
  # ... other derived parameters?
    # psimlogAlpha0[[i]] <- objpp$simulate()$logAlpha0 - doesn't exist - not derived
    # $E ?
  
  psimlogRS_pred[[i]] <- objpp$simulate()$logRS_pred
  
  # add in observation error?
      # look at how I did it for derived_post.R?
  
}; beepr(2)

ppsimalpha <- unlist(psimalpha)
ppsimlogRS <- unlist(psimlogRS_pred)
hist(ppsimalpha, breaks = 100, freq = TRUE, xlim = c(0, 100)) # In the event of large outliers (see outputs from listinit)
hist(unlist(psimlogAlpha_re))
plot(ppsimlogRS, ylim = c(-100, 100))

# Plotting logRS by simulated logRS (INCOMPLETE)
plot(psimlogRS_pred[[1]], dat$logRS, col = rgb(0, 0, 0, 0.1), pch = 16)
for (i in 2:100) {points(psimlogRS_pred[[i]], dat$logRS, col = rgb(0, 0, 0, 0.1), pch = 16)}

# TOR: SEE ABOVE LIMITS FOR OTHER PRIOR PREDICTION SECTION *********************

  # PUSHFORWARD: simulate only the expectation from the priors
  # PREDICTIVE: simulate observations from the priors

  # Useful bayesplot info: https://mc-stan.org/bayesplot/reference/pp_check.html

# default method for stan objects - won't work oob with tmbstan unfort.
# y_rep <- example_yrep_draws() # vs. example_y_data() for y (or fitstan) ?
  # y needs to be a vector (observed)
  # y_rep needs to be a vector (generated)
# pp_check(fitstan, fun = ppc_dens_overlay) # histogram of observed vs. predicted

# LIMITS Set Upper and Lower Limits ####
upper <- numeric(length(obj$par)) + Inf
lower <- numeric(length(obj$par)) + -Inf
lower[names(obj$par) == "tauobs"] <- 0
upper[names(obj$par) == "logESD"] <- 100
lower[names(obj$par) == "logESD"] <- 0
upper[names(obj$par) == "logAlphaSD"] <- 100
lower[names(obj$par) == "logAlphaSD"] <- 0

# nlminb - MLE ####
# stan MAP: do mle via optimize - would just be using tmb obj and nlminb:

# opt <- nlminb(obj$par,
#               obj$fn,
#               obj$gr,
#               control = list(eval.max = 1e5, iter.max = 1e5, trace = 0),
#               lower = lower,
#               upper = upper
# )

# MCMC ####
# INIT FUNCTION - can also just run sampler as default random
  # Tor: Given issues with parallelization - consider that negative init
  # values may be causing issues. For example obj$fn/obj$gr can't be 
  # manually evaluted below zero.

# RANDOM INIT - MATCHING PRIORS
init <- function() {
  # Can also add par <- to sequence above - see prior testing
  listinit <- list(b0 = c(rnorm(1, 10, 1), rnorm(1, 0, 1)), # Contains negatives
       bWA = c(rnorm(1, 0, 1), rnorm(1, 0 ,1)), # Contains negatives
       
       # logRS_pred = rnorm(N_Obs, 0, 1),
       logE_re = rnorm(N_Stk, 0, 1), # Contains negatives
       logAlpha0 = rnorm(1, 0.6, 1), # Contains negatives
       # logAlpha02 = rnorm(1, 0, 1) , # NEW: alpha0 prior for LH specific dists.
       logAlpha_re = rnorm(nrow(dat$WAbase), 0, 1), # Contains negatives

       tauobs = runif(N_Stk, min = 0.005, max = 0.015), # Uniform to REMAIN positive
       
       logESD = runif(1, 0.01, 3), # Positive
       logAlphaSD = runif(1, 0.01, 3) # Positive
  )
  
  if (lhdiston) {
    listinit$logAlpha02 <- rnorm(1, 0, 1) # NEW: alpha0 prior for LH specific dists.
  }
  
  return(listinit)
}

# FIXED POINTS INIT
# initfixed <- function() {
#     listinitfixed <- list(
#        b0 = c(10, 0), # Contains negatives
#        bWA = c(0, 0), # Contains negatives
#        
#        logE_re = numeric(N_Stk), # Contains negatives
#        logAlpha0 = 0.6, # Contains negatives
#        logAlpha_re = numeric(nrow(dat$WAbase)), # Contains negatives
# 
#        tauobs = 0.01 + numeric(N_Stk), # Uniform to REMAIN positive
#        
#        logESD = 1, # Positive
#        logAlphaSD = 1 # Positive
#   )
#   
#   if (lhdiston) {
#     listinitfixed$logAlpha02 <- 0 # NEW: alpha0 prior for LH specific dists.
#   }
#   
#   return(listinitfixed)
# }

# ALTERNATIVE RANDOM INIT
# init2 <- function() {
#   list(b0 = c(rnorm(1,10,1), rnorm(1,0,1)), # Contains negatives
#        bWA = c(rnorm(1,0,1), rnorm(1,0,1)), # Contains negatives
#        
#        logE_re = rnorm(N_Stk, 0, 1), # Contains negatives
#        logAlpha0 = rnorm(1,0.6,1), # Contains negatives
#        logAlpha_re = rnorm(nrow(dat$WAbase), 0, 1), # Contains negatives
# 
#        tauobs = runif(N_Stk, min = 0.005, max = 0.015), # Uniform to REMAIN positive
#        
#        logESD = runif(1, 0.01, 3), # Positive
#        logAlphaSD = runif(1, 0.01, 3) # Positive
#   )
# }

# SAMPLE MCMC
  # Can consider in parallel - Kasper's github example doesn't currently work
# cores <- parallel::detectCores() - 2
# cores <- 1
# options(mc.cores = cores)

# Seeding
# set.seed(1) # This is now set at header of code
# rm(.Random.seed, envir=globalenv())

# The set.seed precedes this in order to have fixed runs with initfixed
set.seed(1) ; fitstan <- tmbstan(obj, iter = 2000, warmup = 1000, # default iter/2 for warmup - Typically 5000 and 1000
                   init = init, # init = init function or "random" or initfixed for fixed points
                   seed = 1, # set seed or leave out for random - now set at 1 above
                      # but not technically within sampling - unsure if each chain takes seed                  
                   # control = list(adapt_delta = 0.99, max_treedepth = 15),
                   lower = lower, upper = upper,
                   chains = 4, open_progress = FALSE, silent = TRUE); beep(2)
  # tmbstan operates by default with NUTS MCMC sampler
  # See: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0197954

# Acquire outderputs of MCMC ####
derived_obj <- derived_post(fitstan); beep(2)

# add stocknames - could add them in to each object of derived_obj?
# Stocknames <- WAin$Stock
# Srep_example <- cbind(derived_obj$deripost_summary$E_tar_adj, Stocknames) |> 
#   Srep_example[c(1, 7, 2, 3, 4, 5, 6)]

# Test and diagnostic plots ####
# names(fitstan) for complete list of parameters from stan object

  # TRACE PLOTS
mcmc_trace(as.array(fitstan)) # ALL PARAMETER TRACEPLOT
mcmc_trace(as.array(fitstan), regex_pars = "b0") # Single regex parameter traceplot e.g. b0

  # PAIRS PLOTS
# pairs_pars <- c("b0", "bWA", "logAlpha0", "logESD", "logAlphaSD")
# pairs(fitstan, pars = pairs_pars) # for specific par names from above

  # ACF
# bayesplot::mcmc_acf(fitstan, regex_pars = "b0") # ACF plot for b0

  # RHAT
# fitstan |> rhat() |> mcmc_rhat() + yaxis_text() # rhat plot for assessing rhat of each parameter

  # LOO/WAIC etc.
# Need a log-likelihood matrix to do this - would have to add a generated
  # quantity in the model to achieve this. Not generated otherwise.

# POSTERIOR PREDICTIVE CHECKS ####

    # This section has to be filled ...

# RESIDUALS? ####

# Plot generated logRS vs. logRS from data?
# plot(exp(dat$logRS), exp(derived_obj$deripost_summary$logRS_pred$Median))
plot(dat$logRS, derived_obj$deripost_summary$logRS_pred$Median)

# OR plot observed vs. generated logRS by stock? Can you see a difference?
compRS <- cbind(dat$srdat, logRS = dat$logRS)
compRS <- cbind(compRS, genMedianlogRS = derived_obj$deripost_summary$logRS_pred$Median)

par(mfrow=c(1,3))
plot(compRS$logRS, ylab = "logRS", ylim = c(-4, 4))
plot(compRS$genMedianlogRS, ylab = "genMedianlogRS", ylim = c(-4, 4))
plot(priorlogRS, ylab = "priorlogRS", ylim = c(-4, 4)) # Coming from original tmb obj$ no sampling, no nlminb

# Tor: this would make it seem like the model is not working great in terms
  # of generating observations from the Ricker model. It is more restricted
  # than the observed data.

# SAVING R OBJECTS: ####
# In script_A.R
save(my_object, file = "my_object.RData")
# In script_B.R
load("my_object.RData")