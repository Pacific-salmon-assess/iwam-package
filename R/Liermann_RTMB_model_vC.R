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
# WAin <- c("DataIn/Parken_evalstocks.csv")
WAin <- c("DataIn/WCVIStocks.csv")

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
  # ( !(Name == "Skagit")) # |> # Skagit is #22
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
            lineWA =  seq(min(WAbase$logWAshifted), max(WAbase$logWAshifted), 0.1),
            logRS = log(srdat$Rec) - log(srdat$Sp))

# External vectors
N_Stk <- max(srdat$Stocknumber + 1)
N_Obs <- nrow(srdat)

# NEW: alpha0 prior for LH specific dists. - ON/OFF
lhdiston <- T # Set true to run
# rm(logAlpha02)
# Just search for logAlpha02 and comment out per run

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

# Random parameter starts for prior simulation
# par <- function() {
#   # Can also add par <- to sequence above - see prior testing
#   listinit <- list(b0 = c(rnorm(1, 10, 1), rnorm(1, 0, 1)), # Contains negatives
#        bWA = c(rnorm(1, 0, 1), rnorm(1, 0 ,1)), # Contains negatives
#        
#        # logRS_pred = rnorm(N_Obs, 0, 1),
#        logE_re = rnorm(N_Stk, 0, 1), # Contains negatives
#        logAlpha0 = rnorm(1, 0.6, 1), # Contains negatives
#        # logAlpha02 = rnorm(1, 0, 1) , # NEW: alpha0 prior for LH specific dists.
#        logAlpha_re = rnorm(nrow(dat$WAbase), 0, 1), # Contains negatives
# 
#        tauobs = runif(N_Stk, min = 0.005, max = 0.015), # Uniform to REMAIN positive
#        
#        logESD = runif(1, 0.01, 3), # Positive
#        logAlphaSD = runif(1, 0.01, 3) # Positive
#   )
#   return(listinit)
# }
# par <- par()

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
  
  E_tar <- numeric(N_Pred)
  logE_tar <- numeric(N_Pred)
  
  logAlpha_tar <- numeric(N_Pred)
  
  # Simulated line vectors
  line <- length(lineWA)
  logE_line_stream <- numeric(line)
  E_line_stream <- numeric(line)
  logE_line_ocean <- numeric(line)
  E_line_ocean <- numeric(line)
  
  # Begin negative log-likelihood
  nll <- 0
  
  nll <- nll - sum(dnorm(b0[1], 10, sd = 31.6, log = TRUE)) # Prior
  nll <- nll - sum(dnorm(b0[2], 0, sd = 31.6, log = TRUE)) # Prior
  nll <- nll - sum(dnorm(bWA[1], 0, sd = 31.6, log = TRUE)) # Prior
  nll <- nll - sum(dnorm(bWA[2], 0, sd = 31.6, log = TRUE)) # Prior
  
  nll <- nll - sum(dnorm(logAlpha0, 0.6, sd = 0.45, log = TRUE)) # Prior (rM)
  if(lhdiston) nll <- nll - sum(dnorm(logAlpha02, 0, sd = 31.6, log = TRUE)) # Prior (rD)
  
  ## Second level of hierarchy - Ricker parameters:
  for (i in 1:N_Stk){
    nll <- nll - dnorm(logE_re[i], 0, sd = logESD, log = TRUE) # Unobserved
    logE[i] <- b0[1] + b0[2]*type[i] + (bWA[1] + bWA[2]*type[i]) * WAbase$logWAshifted[i] + logE_re[i]
    E[i] <- exp(logE[i])
    
    nll <- nll - dnorm(logAlpha_re[i], 0, sd  = logAlphaSD, log = TRUE) # Prior (hj)
    
    if(lhdiston) logAlpha[i] <- logAlpha0 + logAlpha02*type[i] + logAlpha_re[i] # LH specific distributions for prod.
    else logAlpha[i] <- logAlpha0 + logAlpha_re[i]

    nll <- nll - dgamma(tauobs[i], shape = 0.0001, scale = 1/0.0001, log = TRUE)
  }

  ## First level of hierarchy - Ricker model:
      # Where is logRS_pred coming from - and why is it not saved anywhere?
      # Added it into vector list in order to REPORT()
  for (i in 1:N_Obs){
    logRS_pred[i] <- logAlpha[stk[i]]*(1 - S[i]/E[stk[i]])
    nll <- nll - dnorm(logRS[i], logRS_pred[i], sd = sqrt(1/tauobs[stk[i]]), log = TRUE)
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
    logE_line_ocean[i] <- b0[1] + b0[2] + (bWA[1] + bWA[2])*lineWA[i]
    E_line_ocean[i] <- exp(logE_line_ocean[i])
    
    logE_line_stream[i] <- b0[1] + (bWA[1])*lineWA[i]
    E_line_stream[i] <- exp(logE_line_stream[i])
  }
  
  ## ADREPORT - internal values (synoptic specific/Ricker)
  # ADREPORT(logRS) # logRS for all 501 data points
  alpha <- exp(logAlpha)
  
  REPORT(b0) # Testing simulate()
  REPORT(bWA) # Testing simulate()
  
  # ADREPORT(logRS_pred)
  # REPORT(logRS_pred)
  REPORT(logRS_pred)
  
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

# obj$par = new thing
# obj$simulate(obj$par) a bunch of times

# LIMITS ####
# Set Upper and Lower Limits
upper <- numeric(length(obj$par)) + Inf
lower <- numeric(length(obj$par)) + -Inf
lower[names(obj$par) == "tauobs"] <- 0
upper[names(obj$par) == "logESD"] <- 100
lower[names(obj$par) == "logESD"] <- 0
upper[names(obj$par) == "logAlphaSD"] <- 100
lower[names(obj$par) == "logAlphaSD"] <- 0

# nlminb - MLE ####
# opt <- nlminb(obj$par,
#               obj$fn,
#               obj$gr,
#               control = list(eval.max = 1e5, iter.max = 1e5, trace = 0),
#               lower = lower,
#               upper = upper
# )

# PRIOR PUSHFORWARD AND PREDICTIVE CHECKS ####

# TOR: SEE ABOVE LIMITS FOR OTHER PRIOR PREDICTION SECTION *********************

  # PUSHFORWARD: simulate only the expectation from the priors
  # PREDICTIVE: simulate observations from the priors

    # This section has to be filled ...

    # Useful bayesplot info: https://mc-stan.org/bayesplot/reference/pp_check.html

  # Tor: For prior testing - would you have to work on the tmb obj$?
  # e.g. obj$simulate()
  # Does tmbstan() change the obj? Or does it just create the fitstan object?

# default method
# y_rep <- example_yrep_draws() # vs. example_y_data() for y (or fitstan) ?
  # y needs to be a vector (observed)
  # y_rep needs to be a vector (generated)
# pp_check(fitstan, fun = ppc_dens_overlay) # histogram of observed vs. predicted

# MCMC ####
# INIT FUNCTION - can also just run sampler as default random
  # Tor: Given issues with parallelization - consider that negative init
  # values may be causing issues. For example obj$fn/obj$gr can't be 
  # manually evaluted below zero.
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

initfixed <- function() {
    listinitfixed <- list(
       b0 = c(10, 0), # Contains negatives
       bWA = c(0, 0), # Contains negatives
       
       logE_re = numeric(N_Stk), # Contains negatives
       logAlpha0 = 0.6, # Contains negatives
       logAlpha_re = numeric(nrow(dat$WAbase)), # Contains negatives

       tauobs = 0.01 + numeric(N_Stk), # Uniform to REMAIN positive
       
       logESD = 1, # Positive
       logAlphaSD = 1 # Positive
  )
  
  if (lhdiston) {
    listinitfixed$logAlpha02 <- 0 # NEW: alpha0 prior for LH specific dists.
  }
  
  return(listinitfixed)
}

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

# SETTING upper and lower LIMITS for sampler - see limits above

# SAMPLE MCMC
  # Can consider in parallel - Kasper's github example doesn't currently work
# cores <- parallel::detectCores() - 2
# cores <- 1
# options(mc.cores = cores)

# Typically 5000 and 1000
set.seed(1)
fitstan <- tmbstan(obj, iter = 5000, warmup = 1000, init = initfixed, # init = init function or "random" or initfixed for fixed points
                   lower = lower, upper = upper,
                   chains = 4, open_progress = FALSE, silent = TRUE); beep(2)

  # tmbstan operates by default with NUTS MCMC sampler
  # See: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0197954

par(mfrow = c(1,1))
# Acquire outderputs of MCMC ####
derived_obj <- derived_post(fitstan); beep(2)

derived_obj$deripost_summary$SGEN_adj[15,]
derived_obj$deripost_summary$E_tar_adj[15,]
exp(derived_obj$deripost_summary$logAlpha_tar_adj$Median)
# How different are they per run?



# add stocknames
dstocks <- WAin$Stock

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

# Bootstrap Posterior to Match Original IWAM Model ####
    # The point of this is to use the Parken assumptions of productivity
    # and make bootstrapped benchmark estimates from the mcmc chains of the
    # Liermann model.

    # 1. SETUP
BS <- TRUE # default for FALSE
bsiters <- 250000 # New with Brown et al. CSAS runs
outBench <- list()
outAlpha <- list()
# set.seed <- 1
set.seed(1)
conditional <- T # Default T for conditional - FALSE for Marginal medians
prod <- c("Parken") # "LifeStageModel" or "Parken"
# bias.cor <- FALSE
  # bias.cor is also an option - but shouldn't be necessary given that these are posteriors

# library(doParallel)
# library(doRNG)
# cl <- makeCluster(detectCores() - 1)
# registerDoParallel(cl)
# registerDoRNG(seed = 1)

# Function
if (BS == TRUE) {
  # results <- foreach(k = 1:bsiters, .packages = c("dplyr", "purrr")) %dopar% {
  for (k in 1:bsiters) {
    
    # Should this use the conditional or marginal estimate of SREP?
    if (conditional == T) {
    SREP <- derived_obj$deripost_summary$E_tar_adj$Median # Conditional
    SREP_logSE <- (log(SREP) - log(derived_obj$deripost_summary$E_tar_adj$LQ_5)) / 1.96 # Conditional
    SREP_SE <- (SREP - derived_obj$deripost_summary$E_tar_adj$LQ_5) / 1.96 # Conditional
    } else {
    SREP <- derived_obj$deripost_summary$E_tar$Median # Marginal
    SREP_logSE <- (log(SREP) - log(derived_obj$deripost_summary$E_tar$LQ_5)) / 1.96 # Marginal
    SREP_SE <- (SREP - derived_obj$deripost_summary$E_tar$LQ_5) / 1.96 # Marginal
    }
  
    # Steps not included in this bootstrapping function
      # - Removing Stocks: Cypre (LIFESTAGEMODEL ONLY)
        # - Cypre is Stock # 6
        # - Could remove the 6th row per SREP and SMSY since the Stock column
        # isn't used for anything
      # - Scaling: not in Liermann at all
      # - Bias correction: not in Liermann at all
      # - Creation of the RPs dataframe
    
    # 1. LifeStageModel METHOD
      # Could USE POSTERIOR MODE INSTEAD OF MEAN - more likely to be closer to the MLE
    if (prod == "LifeStageModel") {
      Mean.Ric.A <- 1
      Ric.A <- exp(rnorm(length(SREP), Mean.Ric.A, 0))
      Sig.Ric.A <- 0.51 
     
      Ric.A <- exp(rnorm(length(SREP), Mean.Ric.A, Sig.Ric.A))
      if(min(Ric.A)<=0) Ric.A <- exp(rnorm(length(SREP), Mean.Ric.A, Sig.Ric.A))
      
      # if (bias.cor == TRUE) {
      #   sREP <- exp(rnorm(length(SREP), log(RPs$SREP) - 0.5*SREP_logSE$SE^2, SREP_logSE$SE))
      #   if(min(sREP)<=0)   sREP <- exp(rnorm(length(SREP), log(RPs$SREP) - 0.5*SREP_logSE$SE^2,
      #     SREP_logSE$SE))
      # } else {
      sREP <- exp(rnorm(length(SREP), log(SREP), SREP_logSE))
      if(min(sREP)<=0)   sREP <- exp(rnorm(length(SREP), SREP, SREP_SE))
      # }
  
      SGENcalcs <- purrr::map2_dfr (Ric.A, sREP, Sgen.fn2)
    }
    
    # 2. PARKEN METHOD
    if (prod == "Parken"){
      est_loga <- function(SMSY, SREP, shortloga=FALSE){
        loga <- nlminb(start = (0.5 - SMSY/SREP) / 0.07, 
                       objective = calc_loga, # Try to fix with a Scheurel version of LW if possible
                       SMSY= SMSY, 
                       SREP=SREP)$par
        if (shortloga) loga <- (0.5 - SMSY/SREP) / 0.07
        beta <- loga/SREP
        return( list( loga = loga , beta = beta, SMSY = SMSY, SREP = SREP) )
      }
      
      # Again: conditional or marginal estimates? 
      # Mean, Median, or Mode?
      if (conditional == T) SMSY <- derived_obj$deripost_summary$SMSY_adj$Median # Conditional
      else SMSY <- derived_obj$deripost_summary$SMSY$Median # Maringla
      
      # purrr the est_loga function for logA
      lnalpha_Parkin <- purrr::map2_dfr (SMSY, SREP, shortloga=FALSE, 
                                         est_loga)
        # Would this change if SMSY and SREP are calculated with Median vs. Mean? - Yes
        # Median should be used as it would be closer to the MLE estimate.
        # Other consideration would be NOT sampling, but just running nlminb and then 
          # bootstrapping.
      
      # Or do explicit method
      SREP_e <- SREP
      SMSY_e <- SMSY
      loga_e <- SREP_e * (SMSY_e * gsl::lambert_W0(-exp(1-SREP_e / SMSY_e) * (SREP_e-SMSY_e) / SMSY_e) +
        SREP_e - SMSY_e) / (SMSY_e * (SREP_e-SMSY_e))
      beta_e <- loga_e/SREP_e
      
      # Bias correction location 
      
      sREP <- exp(rnorm(length(SREP), log(SREP), SREP_logSE))
      if(min(sREP)<0)   sREP <- exp(rnorm(length(SREP), SREP, SREP_SE$SE))
      
      # Do SGEN calcs with new variables
      SGENcalcs <- purrr::map2_dfr (exp(median(lnalpha_Parkin$loga)), sREP, Sgen.fn2)
      SGENcalcs_e <- purrr::map2_dfr (exp(median(loga_e)), SREP_e, Sgen.fn2) # Explicit
    }
    
    # Previous bind location for RPs
    
    out <- list(bench = select(SGENcalcs, -apar, -bpar))
    outBench[[k]] <- out$bench
    
    # oute <- list(bench = select(SGENcalcs_e, -apar, -bpar))
    # outBenche[[k]] <- oute$bench
    
    if (prod == "Parken") outA <- list(alpha = exp(median(lnalpha_Parkin$loga)))
    if (prod == "LifeStageModel") outA <- list(alpha = Ric.A)
    if (prod == "Parken") outAlpha <- outA
    if (prod == "LifeStageModel") outAlpha[[k]] <- outA
    
    # result <- list(
    #   bench = select(SGENcalcs, -apar, -bpar),
    #   alpha = if (prod == "Parken") exp(median(lnalpha_Parkin$loga)) else Ric.A
    # )
    # return(result)
    
  } ; beep(2) # ; stopCluster(cl)
  
  # outBench <- lapply(results, function(x) x$bench)
  # if (prod == "Parken") {
  #   outAlpha <- lapply(results, function(x) x$alpha)
  # } else if (prod == "LifeStageModel") {
  #   outAlpha <- unique(lapply(results, function(x) x$alpha)) # Should only be one
  # }

  # 3. Outputs: Compile bootstrapped estimates of Sgen, SMSY, and SREP, and identify 5th and 
      # 95th percentiles
  SGEN.bs <- select(as.data.frame(outBench), starts_with("SGEN"))
  
  stockNames <- WAin %>% 
    # filter(Stock != "Cypre") %>% # CYPRE FLAG ## ## ## 
    pull(Stock)
  stockNames <- unique(stockNames)
  
  rownames(SGEN.bs) <- stockNames
  SGEN.boot <- data.frame(SGEN= apply(SGEN.bs, 1, quantile, 0.5), 
                          lwr=apply(SGEN.bs, 1, quantile, 0.025),
                          upr=apply(SGEN.bs, 1, quantile, 0.975) )
  
  SMSY.bs <- select(as.data.frame(outBench), starts_with("SMSY"))
  rownames(SMSY.bs) <- stockNames
  SMSY.boot <- data.frame(SMSY= apply(SMSY.bs, 1, quantile, 0.5), 
                          lwr=apply(SMSY.bs, 1, quantile, 0.025),
                          upr=apply(SMSY.bs, 1, quantile, 0.975) )
  
  SREP.bs <- select(as.data.frame(outBench), starts_with("SREP"))
  rownames(SREP.bs) <- stockNames
  SREP.boot <- data.frame(SREP= apply(SREP.bs, 1, quantile, 0.5), 
                          lwr=apply(SREP.bs, 1, quantile, 0.025),
                          upr=apply(SREP.bs, 1, quantile, 0.975) )
  
  if (prod == "LifeStageModel") {
    APAR.bs <- select(as.data.frame(outAlpha), starts_with("alpha"))
    rownames(APAR.bs) <- stockNames
    APAR.boot <- data.frame(APAR = apply(APAR.bs, 1, quantile, 0.5), 
                            lwr = apply(APAR.bs, 1, quantile, 0.025),
                            upr = apply(APAR.bs, 1, quantile, 0.975) )
  }
  
  boot <- list(SGEN.boot=SGEN.boot, SMSY.boot=SMSY.boot, 
                SREP.boot=SREP.boot) # , APAR.boot=APAR.boot)
  if (prod == "LifeStageModel") boot$APAR.boot <- APAR.boot
  
  df1 <- data.frame(boot[["SGEN.boot"]], Stock=rownames(boot[["SGEN.boot"]]), RP="SGEN") 
  df1 <- df1 %>% rename(Value=SGEN)
  df2 <- data.frame(boot[["SREP.boot"]], Stock=rownames(boot[["SREP.boot"]]), RP="SREP")
  df2 <- df2 %>% rename(Value=SREP)
  df3 <- data.frame(boot[["SMSY.boot"]], Stock=rownames(boot[["SMSY.boot"]]), RP="SMSY")
  df3 <- df3 %>% rename(Value=SMSY)
  if (prod == "LifeStageModel") {
    df4 <- data.frame(boot[["APAR.boot"]], Stock=rownames(boot[["APAR.boot"]]), RP="APAR")
    df4 <- df4 %>% rename(Value=APAR)
  }
  
  dfout <- add_row(df1, df2)
  dfout <- add_row(dfout, df3)
  if (prod == "LifeStageModel") dfout <- add_row(dfout, df4)
  rownames(dfout) <- NULL
  
  # dfout <- dfout %>% mutate(Value=signif(Value, 2)) %>% # Rounded to 2 signif digits
  #   mutate(lwr=signif(lwr,2)) %>% 
  #   mutate (upr=signif(upr,2))
  
  wasample <- WAin %>% 
        select("Stock", "WA", "lh") %>% 
        mutate(WA = round(WA, 0))
  
  BS.dfout <- merge(dfout, wasample, by="Stock", all.x=TRUE, sort=FALSE)
  if (prod == "Parken") alphaout <- outAlpha$alpha
}; beep(2)

# BS.dfout.short <- BS.dfout
# alphaout.short <- alphaout
# BS.dfout.LSM <- BS.dfout
# BS.dfout.parken <- BS.dfout

# Save one of each for analysis
# lsmout <- BS.dfout
# pout <- BS.dfout

# Prepare/load datasets for plotting ####
Parkentable1 <- read.csv(here::here("DataIn/Parken_Table1n2.csv")) # Test stocks e.g. WCVI stocks
ParkenCaseStudy <- read.csv(here::here("DataIn/Parken_evalstocks.csv")) # Case study stocks

dpars <- derived_obj$deripost_summary

targets <- WAin |>
  rename("Stock_name" = Stock) # This name can change depending on what data sets are being run

  # SMSY Estimate for TARGET STOCKS
targets1 <- cbind(targets, derived_obj$deripost_summary$SMSY) |> 
  rename("SMSY_mean" = Mean, "SMSY_median" = Median,
    "SMSY_LQ_5" = LQ_5, "SMSY_UQ_95" = UQ_95, "SMSY_Stocknum" = Stock,
    "SMSY_Mode" = PosteriorMode)
  # SGEN Estimate for TARGET STOCKS
targets2 <- cbind(targets1, derived_obj$deripost_summary$SGEN) |> 
  rename("SGEN_mean" = Mean, "SGEN_median" = Median,
    "SGEN_LQ_5" = LQ_5, "SGEN_UQ_95" = UQ_95, "SGEN_Stocknum" = Stock,
    "SGEN_Mode" = PosteriorMode)
  # Marginal Mean for SREP (E)
targets3 <- cbind(targets2, derived_obj$deripost_summary$E_tar_adj) |> 
  rename("E_tar_adj_mean" = Mean, "E_tar_adj_median" = Median,
    "E_tar_adj_LQ_5" = LQ_5, "E_tar_adj_UQ_95" = UQ_95, "E_tar_adj_Stocknum" = Stock,
    "E_tar_adj_Mode" = PosteriorMode)
  # SREP ESTIMATE FOR TARGET STOCKS
targetsAll <- cbind(targets3, derived_obj$deripost_summary$E_tar) |> 
  rename("E_tar_mean" = Mean, "E_tar_median" = Median,
    "E_tar_LQ_5" = LQ_5, "E_tar_UQ_95" = UQ_95, "E_tar_Stocknum" = Stock,
    "E_tar_Mode" = PosteriorMode)
# Ricker sigma for SYNOPTIC STOCKS 
  # Re-order Stocks to be in order of Ricker variance
# targetsAll <- cbind(targetsAll, derived_obj$deripost_summary$tauobs) |> 
#   rename("tauobs_mean" = Mean, "tauobs_median" = Median,
#     "tauobs_LQ_5" = LQ_5, "tauobs_UQ_95" = UQ_95, "tauobs_Stocknum" = Stock)
    # tauobs is based on the synoptic sets and will now have different lengths

parken <- ParkenCaseStudy |> 
  rename("SMSYp" = SMSY, "SMSYp_5" = SMSY_5, "SMSYp_95" = SMSY_95) |> 
  rename("SREPp" = SREP, "SREPp_5" = SREP_5, "SREPp_95" = SREP_95)

cols <- viridis(8, alpha=0.9, option = "mako", direction = -1)

# Point-wise Benchmark Comparison - BY LOG WA ####
ggplot() +
  
  geom_errorbar(data = parken, aes(x = fct_reorder(Stock, log(WA)), y = SREPp, ymax = SREPp_95, ymin = SREPp_5,
                                       color = "Parken",
                                       width=.1),
  ) +
  geom_point(data = parken,
             aes(x = fct_reorder(Stock, log(WA)), y = SREPp, color = "Parken")) +
  
  geom_errorbar(data = targetsAll, aes(x = fct_reorder(Stock_name, logWA),
                                     y = E_tar_mean,
                                     ymax = E_tar_UQ_95, 
                                     ymin = E_tar_LQ_5,
                                 color = "Liermann MCMC Cond.",
                                 width=.1),
                position = position_nudge(+0.1)) +
  geom_point(data = targetsAll,
             position = position_nudge(+0.1),
             aes(x = fct_reorder(Stock_name, logWA), y = E_tar_mean, color = "Liermann MCMC Cond.")) +
  
  geom_errorbar(data = targetsAll, aes(x = fct_reorder(Stock_name, logWA),
                                     y = E_tar_adj_mean,
                                     ymax = E_tar_adj_UQ_95, 
                                     ymin = E_tar_adj_LQ_5,
                                 color = "Liermann MCMC Marg.",
                                 width=.1),
                position = position_nudge(+0.2)) +
  geom_point(data = targetsAll,
             position = position_nudge(+0.2),
             aes(x = fct_reorder(Stock_name, logWA), y = E_tar_adj_mean, color = "Liermann MCMC Marg.")) +
  
  theme_classic() +
  scale_y_continuous(transform = "log", 
                     breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
  ylab(TeX("$S_{REP}$ Estimate")) +
  xlab("Stock Name") + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust = 1)) +
  scale_color_manual(name='Model',
                     breaks=c('Parken',
                              # 'RTMB MLE',
                              'Liermann MCMC Cond.',
                              'Liermann MCMC Marg.'),
                     values=c('Parken' = "black",
                              # 'RTMB MLE' = "orange",
                              'Liermann MCMC Cond.' = "skyblue",
                              'Liermann MCMC Marg.' = "darkblue"))

# Point-wise Benchmark Comparison - UNSORTED ####
# ggplot() +
#   
#   # Re-written as the Parken model is included in WAin - bound to rtmb MLE version
#   geom_errorbar(data = parken, aes(x = Stock, y = SREPp, ymax = SREPp_95, ymin = SREPp_5,
#                                        color = "Parken",
#                                        width=.1),
#                 # position = position_nudge(-0.4)
#     ) +
#   geom_point(data = parken,
#              # position = position_nudge(-0.4),
#              aes(x = Stock, y = SREPp, color = "Parken")) +
#   
#   # Add in LIERMANN from Liermann_RTMB_model.R as a global object
#   geom_errorbar(data = targetsAll, aes(x = Stock_name,
#                                      y = E_tar_median,
#                                      ymax = E_tar_UQ_95, 
#                                      ymin = E_tar_LQ_5,
#                                  color = "Liermann MCMC",
#                                  width=.1),
#                 position = position_nudge(+0.2)) +
#   geom_point(data = targetsAll,
#              position = position_nudge(+0.2),
#              aes(x = Stock_name, y = E_tar_median, color = "Liermann MCMC")) +
#   
#   theme_classic() +
#   scale_y_continuous(transform = "log", 
#                      breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
#   ylab(TeX("$S_{REP}$ Estimate")) +
#   xlab("Stock Name") + 
#   theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust = 1)) +
#   scale_color_manual(name='Model',
#                      breaks=c('Parken',
#                               # 'RTMB MLE',
#                               'Liermann MCMC'),
#                      values=c('Parken' = "black",
#                               # 'RTMB MLE' = "orange",
#                               'Liermann MCMC' = "skyblue"))

# Linear Regression: Liermann vs. Parken model ####
    # Step 1. Is the data prepared.
options(scipen = 5) # for ticks without sci. notation
col.use <- NA
for(i in 1:length(WAbase$lh)) {
  if (WAbase$lh[i]== 'stream') col.use[i] <- "forestgreen"  # stream
  else col.use[i] <- "dodgerblue3" # ocean
}

    # Step 2. PLOT base points
plot(y = dpars$E$Mean, x = WAbase$WA, pch = 20, 
     # col = ifelse(WAbase$Name == "Chehalis" | WAbase$Name == "Skagit", 'red', col.use), 
     col = col.use, cex = 1.5,
     xlab = expression("Accessible Watershed Area, km"^2), 
     ylab = expression(S[REP]), log = 'xy',
     xlim = c(50,200000) , ylim = c(200,2000000)
  )

# ADD Parken points of SYNOPTIC SET from Table 1 (Parken et al. 2006)
points(y = Parkentable1$Srep, x = Parkentable1$WA, pch = 20, col = alpha('black', 0.5))
    # Shown here because they have changed slighly in WA between model versions?

# Liermann points for TARGET ESTIMATES
    # Can change here betwen Median, Mean, and to Marginal vs. Conditional
points(y = dpars$E_tar$Median, x = WAin$WA, pch = 20, 
  col = alpha(ifelse(WAin$Stock == "Coldwater"| WAin$Stock == "Deadman", 'red', 'skyblue'), 0.5))
  # Can I write this such that for the 'red' points, the alpha is also higher?

    # Step 3. LINES
sum_pars <- summary(fitstan)
bWA1 <- sum_pars$summary[,1][3] 
bWA2 <- sum_pars$summary[,1][4] + bWA1
b01 <- sum_pars$summary[,1][1] 
b01 <- b01 - mean_logWA*bWA1 # To deal with centered/shifted watershed areas
b02 <- sum_pars$summary[,1][2] + sum_pars$summary[,1][1] - mean_logWA*bWA2

simWA <- seq(2, 13, 0.5)
Preds <- b01 + simWA*bWA1
Predso <- b02 + simWA*bWA2
lines(x = exp(simWA), y = exp(Preds), col = alpha("forestgreen", 0.5), lwd = 2, lty = 1)
lines(x = exp(simWA), y = exp(Predso), col = alpha("dodgerblue3", 0.5), lwd = 2, lty = 1)

    # Step 4. Error polygons
# Calculate quantile - uppers and lowers
  # pred_lnSMSY IS THE TARGET's
  # predlnWA should be using: WAin$logWAshifted_t
Eline_stream <- derived_obj$deripost_summary$E_line_stream |> 
  rename("line_stocks" = Stock,
    "s_mean" = Mean,
    "s_median" = Median,
    "s_LQ_5" = LQ_5,
    "s_UQ_95" = UQ_95)
Eline_ocean <- derived_obj$deripost_summary$E_line_ocean |> 
  rename("line_stocko" = Stock,
    "o_mean" = Mean,
    "o_median" = Median,
    "o_LQ_5" = LQ_5,
    "o_UQ_95" = UQ_95)
lineWA <- cbind(dat$lineWA, Eline_stream, Eline_ocean) # NEED NEW LINE VALUES

up_S <- lineWA$s_UQ_95
lo_S <- lineWA$s_LQ_5
up_O <- lineWA$o_UQ_95
lo_O <- lineWA$o_LQ_5

polygon(x=c(exp(lineWA$`dat$lineWA` + mean_logWA), exp(rev(lineWA$`dat$lineWA` + mean_logWA))), 
        y=c(up_S, rev(lo_S)), 
        col=rgb(0,0.4,0, alpha=0.2), border=NA)

polygon(x=c(exp(lineWA$`dat$lineWA` + mean_logWA), exp(rev(lineWA$`dat$lineWA` + mean_logWA))), 
        y=c(up_O, rev(lo_O)), 
        col=rgb(0,0.2,0.4, alpha=0.2), border=NA) 

    # Step 5. Grab Parken estimates for the line and add as y = mx + b
    # From Table 4. Srep Habitat Model
Preds_Parken <- 3.89 + simWA*0.693 + (0.240/2) # Stream-type
Predso_Parken <- 3.52 + simWA*0.878 + (0.133/2) # Ocean-type
lines(x = exp(simWA), y = exp(Preds_Parken), col = alpha("forestgreen", 0.5), lwd = 2, lty = 2)
lines(x = exp(simWA), y = exp(Predso_Parken), col = alpha("dodgerblue3", 0.5), lwd = 2, lty = 2)

    # Step 6. Add text to describe model equations
# q: Based on the plot coded above, how can I add text labels to each point to state the name associated with them?
# a: Use geom_text() to add text labels to the plot.
# geom_text()

# Bar plot comparison of SYNOPTIC values of SREP ####
  # Compare Parken and Liermann estimates of SREP for the SYNOPTIC STOCKS
tempEpars <- dpars$E |> 
  rename("E_stock_temp" = Stock)
bardf <- cbind(Parkentable1, tempEpars)
bardf_long <- bardf %>%
  pivot_longer(cols = c(Srep, Mean), names_to = "Category", values_to = "Value")

# Now plot using a single geom_bar()
bcompare <- ggplot(bardf_long, aes(x = Stock, y = Value, fill = Category)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  
  theme_classic() +
  scale_fill_manual(
    name = "Model",  # Custom legend title
    values = c("Srep" = "black", "Mean" = "skyblue"),  # Custom colors
    labels = c("Srep" = "Parken", "Mean" = "Liermann")  # Custom category names
  ) +
  labs(x = "Stock",
       y = expression(S[REP])) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
bcompare

# SR Curves for individual stocks - NOT FINISHED ####
Stks <- unique(srdat$Stocknumber)
NStks <- length(Stks)
par(mfrow=c(5,5), mar=c(2, 2, 1, 0.1) + 0.1, oma=c(3,3,1,1)) # To fit all the plots on one grid
# par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1, oma=c(0,0,0,0)) # To plot each individual stock

# ADDING IWAM estimates - will produce errors if it is not the standard iwam_obj name

Parken_ab <- read.csv(here::here("DataIn/Parken_Table1n2.csv")) 
  # KSR, Andrew, Lewis, Columbia - just different names - but ORDER is the same
  # Needs correction for removal of Skagit and Chehalis
if (!'Skagit' %in% WAbase$Name) Parken_ab <- Parken_ab |> filter( !(Stock == "Skagit")) 
if (!'Chehalis' %in% WAbase$Name) Parken_ab <- Parken_ab |> filter( !(Stock == "Chehalis")) 
Parken_ab <- Parken_ab |> 
  mutate(Stocknumber = as.integer(factor(Stocknumber)) - 1) # Re-order Stocknumber to be ascending

for (i in Stks){
  # Get stocknames and numbers
      # names <- pars %>% dplyr::select ("Name", "Stocknumber") %>% distinct()
      # name <- pars %>% filter (Stocknumber==i) %>% 
      #   dplyr::select ("Name") %>% distinct()
  name <- names$Name[names$Stocknumber == i]
  
  R <- srdat %>% filter (Stocknumber==i) %>% 
    dplyr::select(Rec) 
  S <- srdat %>% filter (Stocknumber==i) %>% 
    dplyr::select(Sp) 
  
  if(name != "Skagit" & name != "KSR" & name != "Humptulips" & name != "Chehalis") 
    plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)+(max(S$Sp)/3)), ylim=c(0,max(R$Rec) ) )
  if(name == "Skagit") 
    plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)*3+(max(S$Sp)*2)), ylim=c(0,max(R$Rec) ) )
  if(name == "KSR") 
    plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,1150), ylim=c(0,max(R$Rec) ) ) # xlim was 500 or 1000
  if(name == "Humptulips") # Requires extra larger xlim
    plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)+10000), ylim=c(0,max(R$Rec) ) )
  if(name == "Chehalis") # Requires extra larger xlim
    plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)+50000), ylim=c(0,max(R$Rec) ) )

  # Get alpha and beta parameters
      # NOTE: can substitute $alpha$Mean if you want to use Mean
      # Median is resistant to monotonic transformations e.g. ln, exp, ... cube
  a <- as.numeric(derived_obj$deripost_summary$logAlpha$Median[derived_obj$deripost_summary$logAlpha$Stock - 1 == i])
  b <- as.numeric(derived_obj$deripost_summary$BETA_r$Median[derived_obj$deripost_summary$BETA_r$Stock - 1 == i])
  
    # ADD IN IWAM ESTIMATES
      # Require running the IWAM model function externally - see IWAMfunctionrun.R
      # Save as a global variable
  if(exists("iwamobj")) {
    aiwam <- iwamobj$modelpars |> filter (Stocknumber==i) |>  
      filter(Param=="logA") %>%
      summarise(A=exp(Estimate)) %>%
      as.numeric()
    Sc <- iwamobj$srdat %>% filter (Stocknumber==i) %>% 
       dplyr::select(scale) %>% distinct() %>% 
       as.numeric()
    biwam <- iwamobj$modelpars %>% filter (Stocknumber==i) %>%
      filter(Param=="logB") %>%
      summarise(B=exp(Estimate)/Sc) %>%
      as.numeric()
    biwam_se <- iwamobj$modelpars %>%
      filter(Stocknumber == i, Param == "logB") %>%
      mutate(B_lower = exp(Estimate - 1.96 * Std..Error) / Sc,
        B_upper = exp(Estimate + 1.96 * Std..Error) / Sc )
  }
  
  # Parken values for skagit - from Parken et al. 2006 Table 2 (Ocean life-histories)
  skagit_alpha <- 7.74
  skagit_beta <- 0.0000657
  RR_skagit <- NA
  SS <- RR <- RR_parken <- RRiwam <- NA
  ap <- Parken_ab$Alpha # vector
  bp <- Parken_ab$Beta # vector
  
  for (j in 1:250){ # Creates a step-wise sample line by which to create a line on
    if ("Skagit" %in% WAbase$Name) if (i!=22 & i!=7) SS[j] <- j*(max(S$Sp)/100) # IF NOT SKAGIT OR KSR - 
      # When Skagit exists
    if ("Skagit" %in% WAbase$Name) if (i==22) SS[j] <- j*(max(S$Sp*3)/100) # If Skagit exists
    
    if(!"Skagit" %in% WAbase$Name) if (i!=7) SS[j] <- j*(max(S$Sp)/100) # If not KSR - When Skagit doesn't exist
    if (i==7) SS[j] <- j*(500/100) # if KSR - could also add: if ("King Salmon" %in% WAbase$Name) 
    
    RR[j] <- exp(a) * SS[j] * exp(-b * SS[j])
    
    RR_parken[j] <- ap[i+1] * SS[j] * exp(-bp[i+1] * SS[j])
    
    if ("Skagit" %in% WAbase$Name) if(i==22) {RR_skagit[j] <- skagit_alpha * SS[j] * exp(-skagit_beta * SS[j])} 
      # Skagit Line based on alpha and beta from Table 1 and 2 from Parken et al. 2006
    
    if(exists("iwamobj")) RRiwam[j] <- aiwam * SS[j] * exp(-biwam * SS[j])
  }
  
  mtext(name, side=3, cex=0.8)
  
  col.use <- "black"
  lines(x=SS, y=RR, col='black') 
  
  # For Skagit, add Parken et al. 2006 model curve
  if ("Skagit" %in% WAbase$Name) if(i==22) lines(x=SS, y=RR_skagit, lty="dashed") # }
  
  # For all stocks, added in Parken et al. 2006 model curve
  lines(x=SS, y=RR_parken, lty="dashed", col="red")

  # For all stocks, add IWAM model curve
  if(exists("iwamobj")) lines(x=SS, y=RRiwam, lty="dashed", col="forestgreen")
  
  # Calculate and Plot VERTICAL LINES FOR SMSY, SMAX, OR SREP
  # SMSY <- derived_obj$deripost_summary$SMSY_r$Median[derived_obj$deripost_summary$SMSY_r$Stock - 1 == i]
  # SMSY_ul <- derived_obj$deripost_summary$SMSY_r$UQ_95[derived_obj$deripost_summary$SMSY_r$Stock - 1 == i]
  # SMSY_ll <- derived_obj$deripost_summary$SMSY_r$LQ_5[derived_obj$deripost_summary$SMSY_r$Stock - 1 == i]
  SMAX <- 1/derived_obj$deripost_summary$BETA_r$Median[derived_obj$deripost_summary$BETA_r$Stock - 1 == i]
  SMAX_ul <- 1/derived_obj$deripost_summary$BETA_r$UQ_95[derived_obj$deripost_summary$BETA_r$Stock - 1 == i] # Rev.
  SMAX_ll <- 1/derived_obj$deripost_summary$BETA_r$LQ_5[derived_obj$deripost_summary$BETA_r$Stock - 1 == i] # Rev.
  # SREP <- derived_obj$deripost_summary$E$Median[derived_obj$deripost_summary$E$Stock - 1 == i]
  # SREP_ul <- derived_obj$deripost_summary$E$UQ_95[derived_obj$deripost_summary$E$Stock - 1 == i]
  # SREP_ll <- derived_obj$deripost_summary$E$LQ_5[derived_obj$deripost_summary$E$Stock - 1 == i]
  if(exists("iwamobj")) SMAXiwam <- 1/biwam
  if(exists("iwamobj")) SMAXiwam_ll <- 1/(biwam_se$B_lower)
  if(exists("iwamobj")) SMAXiwam_ul <- 1/(biwam_se$B_upper)
  
  # abline(v = SMSY, col=col.use, lty='dotted')
  abline(v = SMAX, col=col.use, lty='dotted')
  # abline(v = SREP, col=col.use, ly='dotted')
  
  # CI' shaded polygons - Repeat for SREP or SMSY if desired
  # IWAM SMAX CI
  if(exists("iwamobj")) polygon(x=c(SMAXiwam_ul, SMAXiwam_ll, SMAXiwam_ll, SMAXiwam_ul),
        y=c(-10000,-10000,10000+max(R$Rec),10000+max(R$Rec)),
        col=rgb(0,0.4,0, alpha=0.1), border=NA )
  
  # polygon(x=c(SMSY_ul, SMSY_ll, SMSY_ll, SMSY_ul), 
  #         y=c(-10000,-10000,10000+max(R$Rec),10000+max(R$Rec)), 
  #         col=grey(0.8, alpha=0.4), border=NA )

  polygon(x=c(SMAX_ul, SMAX_ll, SMAX_ll, SMAX_ul), 
        y=c(-10000,-10000,10000+max(R$Rec),10000+max(R$Rec)), 
        col=grey(0.8, alpha=0.4), border=NA )
  
  # Parken Smsy Estimate (vert. line) from Table 1/2 Parken et al. 2006
  # Parken_smsy <- Parken_ab$Smsy[Parken_ab$Stocknumber == Parken_ab$Stocknumber[i+1]]
  # abline(v = Parken_smsy, col="red", lty='dotted')

  Parken_smax <- 1 / Parken_ab$Beta[Parken_ab$Stocknumber == Parken_ab$Stocknumber[i+1]]
  abline(v = Parken_smax, col="red", lty='dotted')

  # Parken_srep <- Parken_ab$Srep[Parken_ab$Stocknumber == Parken_ab$Stocknumber[i+1]]
  # abline(v = Parken_srep, col="red", lty='dotted')
  
  if(exists("iwamobj")) IWAM_smax <- 1 / biwam
  if(exists("iwamobj")) abline(v = IWAM_smax, col="forestgreen", lty='dotted')
}

# Add an GLOBAL figure axis label across par()
  # x = Spawners, y = Recruitment
mtext("Spawners", side = 1, line = 1, outer = TRUE, cex = 1.3)
mtext("Recruitment", side = 2, line  = 1, outer = TRUE, cex = 1.3, las = 0)

# Compare Deviation of Marginal and Conditional Means ####
edeviation <- data.frame(re = derived_obj$deripost_summary$logE_tar_adj,
                      tar = derived_obj$deripost_summary$logE_tar,
                      dev = derived_obj$deripost_summary$logE_tar_adj$Median - derived_obj$deripost_summary$logE_tar$Median,
                      og_re = derived_obj$deripost_summary$logE_re$Median,
                      name = WAin$Stock)
  # Note: logE_re is not used - unsure what to plot this against

ggplot(edeviation, aes(x = name, y = dev)) +
  geom_hline(yintercept = 0, color = "gray") +
  geom_point(size = 4, color = "darkred") +
  labs(y = "Deviation in predicted E (More positive is\n more random effect)", x = "Stock Name") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
# Point-wise Benchmark Comparison - BY RICKER VARIANCE - ONLY for re-predicting SYNOPTIC STOCKS ####
# full <- cbind(targetsAll, parken$Stock, parken$SREPp, parken$SREPp_5, parken$SREPp_95)
# 
# ggplot() +
#   
#   # Re-written as the Parken model is included in WAin - bound to rtmb MLE version
#   geom_errorbar(data = full, aes(x = fct_reorder(parken$Stock, tauobs_mean), 
#                                  y = parken$SREPp, ymax = parken$SREPp_5, ymin = parken$SREPp_95,
#                                  color = "Parken", width=.1)) +
#   geom_point(data = full, aes(x = fct_reorder(parken$Stock, tauobs_mean), 
#                               y = parken$SREPp, color = "Parken")) +
#   
#   # Add in LIERMANN from Liermann_RTMB_model.R as a global object
#   geom_errorbar(data = full, aes(x = fct_reorder(Stock_name, tauobs_mean),
#                                  y = E_tar_median,
#                                  ymax = E_tar_UQ_95, 
#                                  ymin = E_tar_LQ_5,
#                                  color = "Liermann MCMC",
#                                  width=.1),
#                 position = position_nudge(+0.2)) +
#   geom_point(data = full,
#              position = position_nudge(+0.2),
#              aes(x = fct_reorder(Stock_name, tauobs_mean), 
#                  y = E_tar_mean, color = "Liermann MCMC")) +
#   
#   theme_classic() +
#   scale_y_continuous(transform = "log", 
#                      breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
#   ylab(TeX("$S_{REP}$ Estimate")) +
#   xlab("Stock Name") + 
#   theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust = 1)) +
#   scale_color_manual(name='Model',
#                      breaks=c('Parken',
#                               # 'RTMB MLE',
#                               'Liermann MCMC'),
#                      values=c('Parken' = "black",
#                               # 'RTMB MLE' = "orange",
#                               'Liermann MCMC' = "skyblue"))
# Testing deviations when re-predicting SYNOPTIC STOCKS ####
smsy_deviation <- derived_obj$deripost_summary$SMSY_r$Median - derived_obj$deripost_summary$SMSY$Median
smax_deviation <- (1/derived_obj$deripost_summary$BETA_r$Median) - (1/derived_obj$deripost_summary$BETA$Median)
testdf <- data.frame(name = WAin$Name, smsy_deviation = smsy_deviation, smax_deviation = smax_deviation)

# Plot of deviations in SMSY
ggplot(testdf, aes(x = name, y = smsy_deviation)) +
  geom_hline(yintercept = 0, color = "gray") +
  geom_point(size = 4, color = "darkred") +
  labs(y = "Deviation in SMSY", x = "Name") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(testdf, aes(x = name, y = smax_deviation)) +
  geom_hline(yintercept = 0, color = "gray") +
  geom_point(size = 4, color = "darkred") +
  labs(y = "Deviation in SMAX", x = "Name") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# As a regression plot
smax_r <- (1/derived_obj$deripost_summary$BETA_r$Median)
smax <- (1/derived_obj$deripost_summary$BETA$Median)
testdf2 <- data.frame(smax_r = smax_r, smax = smax)

ggplot(data = testdf2, aes(x = smax_r, y = smax)) +
  geom_point() + 
  labs(y = "Ricker Smax", x = "Predicted Smax") + 
  theme_classic() + 
  geom_abline(intercept = 0, slope = 1, color = "gray")

