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

# Shift log WA for the mean - base - makes estimation easier
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
            prioronly = 0) # 0 - run with data, 1 - prior prediction mode?

# External vectors
N_Stk <- max(srdat$Stocknumber + 1)
N_Obs <- nrow(srdat)
stk = srdat$Stocknumber + 1

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
    # half sdnormal? half t? half cauchy?
    
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
    
    # prior pushforward - turn it back into the curve that is expected based on the priors
    
    # posterior predictive check
    # logRS_sim[i] <- rnorm(1, logRS_pred[i], sd = sqrt(1/tauobs[stk[i]])) # Simulated logRS for posterior predictive check
      # do I need an rnorm or a simulate?
      # I think I could either:
      # rnorm for both parameters as the mean?
      # or simulate() and then rnorm?
    
    # outside of RTMB:
    # create an rnorm with mean logRS_pred, sd = sqrt(1/tauobs) etc.
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
  
  # QUESTION: Do you need both ADREPORT and REPORT per variable?
    # Does ADREPORT cover REPORT too?
  
  ## ADREPORT - internal values (synoptic specific/Ricker)
  alpha <- exp(logAlpha)
  
  REPORT(b0) # Testing simulate()
  REPORT(bWA) # Testing simulate()
  
  # ADREPORT(logRS_pred)
  REPORT(logRS_pred)
  # REPORT(logRS_sim) # Is logRS_pred as a REPORT the same as logRS_sim? (Question for Sean/Paul)
  
  # ADREPORT(logE_re)
  # ADREPORT(E)
  # ADREPORT(logAlpha)
  # ADREPORT(logAlpha_re) # random effect parameter for resampling
  # ADREPORT(alpha)
  # ADREPORT(SMSY_r)
  # ADREPORT(BETA_r)
  # ADREPORT(tauobs)
  
  # REPORT(logRS) # logRS for all 501 data points
  REPORT(logE_re)
  REPORT(E) # E (Srep) for all synoptic data set rivers (25)
  REPORT(logAlpha) # model logAlpha (25)
  REPORT(logAlpha_re) # random effect parameter for resampling
  REPORT(alpha)
  REPORT(SMSY_r)
  REPORT(BETA_r)
  REPORT(tauobs) # Necessary to add back in observation error?
  
  # ADREPORT - predicted values from watershed area model
    # Mean estimate of the median (without bias correction)
  alpha_tar <- exp(logAlpha_tar)
  
  # ADREPORT(E_tar) # target E (Srep) (21)
  # ADREPORT(logE_tar) # exp these for the correct confidence intervals
  # ADREPORT(logAlpha_tar)
  # ADREPORT(alpha_tar)
  
  REPORT(E_tar)
  REPORT(logE_tar)
  REPORT(logAlpha_tar)
  REPORT(alpha_tar)

  # ADREPORT(BETA)
  # ADREPORT(SMSY)
  # ADREPORT(SGEN)
  
  REPORT(BETA)
  REPORT(SMSY)
  REPORT(SGEN)
  
  # Simulated line values for plotting
  REPORT(E_line_stream) 
  # ADREPORT(E_line_stream)
  REPORT(E_line_ocean) 
  # ADREPORT(E_line_ocean)
  
  REPORT(logAlphaSD)
  REPORT(logESD)
  
  nll # output of negative log-likelihood
}



## MakeADFun ####
obj <- RTMB::MakeADFun(f_srep,
                       par,
                       random = c("logAlpha_re", "logE_re"),
                       silent=TRUE)



# Prior testing thoughts ####

# PUSHFORWARD: simulate only the expectation from the priors.
# simulations of expectations --- whatever is going into the mean part of the data likelihood
# e.g. values of expected log(R/S). You can never ‘observe’ these in reality

# PREDICTIVE: simulate observations from the priors. 
# e.g. simulating log(R/S) 
# compare these against your data if you want. The point is, 
# it’s the same type of thing as your data that you can observe.
# whatever is going into the left side of the data likelihood line (at least with Stan tilde syntax)

# Sean: It’s not in reference to what parameter is being simulated or predicted. 
# It’s about whether you’re taking draws from the expected value distribution 
# (given only the priors) or whether you’re taking draws from the expected value 
# distribution and adding on observation error too.

# Approaches:
    # 1. Exclude data from tmbstan sampling:
        # Create a logical exclusion/inclusion term in data
        # For PP - exclude nll term of observations e.g. logRS
        # run tmbstan() - will most likely fail horribly
        # Extract parameter distributions from chains (*PPF*)
        # ... (PP)
    # 2. Random generation of RTMB model:
        # Randomly sample parameters starts
        # MakeADFun
        # obj$simulate() - random generation of all likelihoods?
        # loop above niters and save them
        # This produces a distribution of expected values given only the priors (*PPF*)
        # To calculate prior predictive: add observation error to the simulated value of logRS_pred
        # by drawing an rnorm of the simulated value of logRS_pred (mu) and tauobs (sd)
        # OR rnorm drawn values (without a simulate())?
        # ** This final step is the same in both cases.

# priorlogRS <- obj$simulate()$logRS_pred

# What is the difference between logRS_pred and logRS as an OBS() object?

# Alternative approach from Sean ####
# rnorms for each of alpha prior
# a function that takes in parameters (Vectored) and creates a logRS, etc., 
# and then pipe through vector - recreating what REPORT() does
# this can then be done for all forms of prior/posterior prediction
# just avoids all the behind the scene work of the RTMB framework

# APPROACH: Random Generation ####
    # Create however many vectors of parameters you wish to test/investigate
psimalpha <- vector("list", 1000) # Vector start
psimlogRS_pred <- vector("list", 1000) # Vector start
psimlogAlpha_re <- vector("list", 1000) # Vector start

for (i in 1:100){
  # Random parameter creation: Random parameter starts for prior simulation
  parp <- function() {
  # Random parameter starts for prior simulation
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
  
  # Alternative priors from Liermann et al. 2010 - comment on/off for now
  logESDalt <- runif(1, 0, 20)
  logAlphaSDalt <- runif(1, 0, 20)
  listinitalternate <- list(b0 = c(rnorm(1, 7, sd = 31.6), rnorm(1, 0, sd = 10)),
    bWA = c(rnorm(1, 0, sd = 10), rnorm(1, 0 , sd = 10)),
    logESD = logESDalt, # This isn't being saved internally
    logAlphaSD = logAlphaSDalt, # This isn't being saved internally

    logE_re = rnorm(N_Stk, 0, sd = logESD),
    logAlpha0 = rnorm(1, 0, sd = 0.71),
    logAlpha02 = rnorm(1, 0, sd = 31.6),
    logAlpha_re = rnorm(nrow(dat$WAbase), 0, sd = logAlphaSD),

    tauobs = rgamma(N_Stk, shape = 0.0001, scale = 1/0.0001)
    )
  
    return(listinitprior) # listinitprior or listinitialternate
  }
  
  parpar <- parp()
  
  # enter parpar into model function with MakeADFun
  objpp <- RTMB::MakeADFun(f_srep,
                       parpar,
                       random = c("logAlpha_re", "logE_re"),
                       silent=TRUE)
  
  # simulate for each desired term
  psimalpha[[i]] <- objpp$simulate()$alpha
  psimlogAlpha_re[[i]] <- objpp$simulate()$logAlpha_re # why logAlpha_re? - if I want to do PP?
  # ... other derived parameters?
    # psimlogAlpha0[[i]] <- objpp$simulate()$logAlpha0 - doesn't exist - not derived
    # $E ?
  
  psimlogRS_pred[[i]] <- objpp$simulate()$logRS_pred
  
  # add in observation error?
      # look at how I did it for derived_post.R?
  
}

# Plot distributions for PRIOR PUSHFORWARD
ppsimalpha <- unlist(psimalpha) 
ppsimlogRS <- unlist(psimlogRS_pred)

hist(ppsimalpha, breaks = 100, freq = TRUE) # Graphical hiccups
hist(unlist(psimlogAlpha_re))
plot(ppsimlogRS, ylim = c(-100, 100))

# Plotting logRS by simulated logRS (INCOMPLETE)
plot(psimlogRS_pred[[1]], dat$logRS, col = rgb(0, 0, 0, 0.1), pch = 16)
for (i in 2:100) {points(psimlogRS_pred[[i]], dat$logRS, col = rgb(0, 0, 0, 0.1), pch = 16)}

# Useful bayesplot info: https://mc-stan.org/bayesplot/reference/pp_check.html

# Plot distributions for PRIOR PREDICTIVE
# First: add in observation error to logRS_pred
# ...



# LIMITS: Set Upper and Lower Limits ####
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

# PREVIOUS RANDOM INIT (DEPRECIATED)
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

# SAMPLE MCMC ####
  # Can consider in parallel - Kasper's github example doesn't currently work
# cores <- parallel::detectCores() - 2
# cores <- 1
# options(mc.cores = cores)

# Seeding
# rm(.Random.seed, envir=globalenv())

# The set.seed precedes this in order to have fixed "identical" runs with initfixed
set.seed(1) ; fitstan <- tmbstan(obj, iter = 5000, warmup = 2500, # default iter/2 for warmup - Typically 5000 and 1000
                   init = init, # init = init function or "random" or initfixed for fixed points
                   seed = 1, # set seed or leave out for random - now set at 1 above
                      # but not technically within sampling - unsure if each chain takes seed                  
                   # control = list(adapt_delta = 0.99, max_treedepth = 15),
                   lower = lower, upper = upper,
                   chains = 4, open_progress = FALSE, silent = TRUE); beep(2)
# tmbstan operates by default with NUTS MCMC sampler
# See: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0197954

# Test and diagnostic plots ####
# names(fitstan) for complete list of parameters from stan object

  # TRACE PLOTS
# mcmc_trace(as.array(fitstan)) # ALL PARAMETER TRACEPLOT
# mcmc_trace(as.array(fitstan), regex_pars = "b0") # Single regex parameter traceplot e.g. b0

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



# Acquire outputs of MCMC ####
derived_obj <- derived_post(fitstan); beep(2)
# add stocknames - see extra code from _Plots.R



# APPROACH - Prior prediction with data exclusion ####
# Make sure dat$prioronly <- 1 to not include data in the nll
# Run all the way down to derived_obj

# Pushforward tests:

# Predictive tests: *EXACTLY THE SAME AS BELOW FOR POSTERIOR*
musim <- derived_obj$deripost_full$logRS_pred
sigmasim <- derived_obj$deripost_full$tauobs
simout <- matrix(NA, nrow = 10, ncol = 501)
draws <- sample(1:10000, 10)
for (i in seq_along(draws)) {
  row <- draws[i]
  # Vectorized rnorm over 501 obs using that row
  means <- musim[row, ]                   # length 501
  sds   <- sigmasim[row, stk]             # length 501 — select right SD for each obs
  simout[i, ] <- rnorm(501, mean = means, sd = sqrt(1/sds))
}



# POSTERIOR PREDICTIVE CHECKS ####
# Run full model with dat$prioronly <- 0 for data included in nll
# Run dervied_post() to extract posterior from chains
# Randomly sample for logRS
# Compare logRS with logRS_pred (including observation error)
# e.g. rnorm(1, derived_obj$deripost_summary$logRS_pred$Median, 
           # sd  = sqrt(1/derived_obj$deripost_summary$tauobs$Median))

# each sample independent draw AND per observation
  # 501 observations x 4000 times (matrix)

# pp_logRS <- c()
# for (i in stk){
#   pp_logRS[i] <- rnorm(1, derived_obj$deripost_summary$logRS_pred$Median[i], # per sample
#                   sd = sqrt(1/derived_obj$deripost_summary$tauobs$Median[stk[i]])) # sqrt(1/tau[i])
# }
# hist(pp_logRS)

# create matrix of [10, 501] for 10 draws from 10000 iters
for (i in stk){
  ppsim[i] <- rnorm(1, derived_obj$deripost_full$logRS_pred[],
                    sd = sqrt(1/derived_obj$deripost_full$tauobs[stk[i]]))
}

# other method from derived_post
  # unsure how to work this around indexing per stock for tau
# matrices$logE_tar_adj <- apply(matrices$logE_tar, 2, FUN = function(x)rnorm(length(x), x, sd = matrices$logESD))
# matrixsim <- apply(derived_obj$deripost_full$logRS_pred, 2, FUN = function(x)rnorm(length(x), x, sd = sqrt(1/derived_obj$deripost_full$tauobs)))

# one method
musim <- derived_obj$deripost_full$logRS_pred
sigmasim <- derived_obj$deripost_full$tauobs
simout <- matrix(NA, nrow = 10, ncol = 501)
draws <- sample(1:10000, 10)
for (i in seq_along(draws)) {
  row <- draws[i]
  # Vectorized rnorm over 501 obs using that row
  means <- musim[row, ]                   # length 501
  sds   <- sigmasim[row, stk]             # length 501 — select right SD for each obs
  simout[i, ] <- rnorm(501, mean = means, sd = sqrt(1/sds))
}


# RESIDUALS? ####

# Plot generated logRS vs. logRS from data?
# plot(exp(dat$logRS), exp(derived_obj$deripost_summary$logRS_pred$Median))
plot(dat$logRS, derived_obj$deripost_summary$logRS_pred$Median)

# OR plot observed vs. generated logRS by stock? Can you see a difference?
compRS <- cbind(dat$srdat, logRS = dat$logRS)
compRS <- cbind(compRS, genMedianlogRS = derived_obj$deripost_summary$logRS_pred$Median)

par(mfrow=c(1,2))
plot(compRS$logRS, ylab = "logRS", ylim = c(-10, 10))
plot(compRS$genMedianlogRS, ylab = "genMedianlogRS", ylim = c(-10, 10))
points(
  x = rep(1:501, each = 4000),
  y = as.vector(t(derived_obj$deripost_full$logRS_pred)),
  col = rgb(0, 0, 0, alpha = 0.02),  # transparent black
  pch = 16,
  cex = 0.3
)
# plot(priorlogRS, ylab = "priorlogRS", ylim = c(-4, 4)) # Coming from original tmb obj$ no sampling, no nlminb

# Tor: this would make it seem like the model is not working great in terms
  # of generating observations from the Ricker model. It is more restricted
  # than the observed data.

# SAVING R OBJECTS: ####
# In script_A.R
save(my_object, file = "my_object.RData")
# In script_B.R
load("my_object.RData")