# Liermann Srep (E) RTMB Model with MCMC Sampling ####

# Libaries ####
library(RTMB)
library(ggplot2)
library(dplyr)
library(tidyverse) 
library(progress) # Progress bar
library(tmbstan) # MCMC tmb model sampling
# library(TMB) # Original TMB
# library(tidybayes) # bayesian visualization
library(bayesplot) # bayesian visualization
library(coda) # bayesian package
library(beepr) # Sounds
library(viridis) # Colours
library(ggridges) # Ridge plots
library(data.table) # Create data tables for pivoting
library(microbenchmark) # Timing functions

here::i_am("R/LambertWs.R") # New line for re-establishing here
							# when not working in RStudio

source(here::here("R/LambertWs.R")) # Lambert W function
source(here::here("R/helperFunctions.R")) # For bootstrapping
source(here::here("R/derived_post.R")) # For posterior extraction

# New LambertW0
# See: https://github.com/pbs-assess/renewassess/blob/main/code/RTMB/PosteriorPredictiveSample.r
LambertW0 <- ADjoint(
  function(x){gsl::lambert_W0(x)},
  function(x, y, dy) {dy / (x + exp(y))}
)

# Raw data read-in ####
WAin <- c("DataIn/Parken_evalstocks.csv") # c("DataIn/WCVIStocks.csv")
# For re-evaluation of synoptic sets e.g. WAbase
    # Run model until setting up data section
    # Then over-write WAin <- WAbase
    # And rename logWAshifted to logWAshifted_t
    # And make sure to change type_tar to fit the expected 0:1 indexing

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
  mutate(lh = factor(ifelse(Stream == 0, "stream", "ocean"),
    levels = c("stream", "ocean"))) |> # Stream = 0, Ocean  = 1
  mutate(Stocknumber = as.integer(factor(Stocknumber)) - 1) # Re-numbering uniquely

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
            lineWA =  seq(min(WAbase$logWAshifted), 
                          max(WAbase$logWAshifted), 0.1), # Not added to NLL
            logRS = log(srdat$Rec) - log(srdat$Sp),
            prioronly = 0) # 0-run with data, 1-prior prediction mode

# External vectors
N_Stk <- max(srdat$Stocknumber + 1)
N_Obs <- nrow(srdat)
stk = srdat$Stocknumber + 1

# NEW: alpha0 prior for LH specific dists.
lhdiston <- T # T = LH Specific
bias.cor <- F # T = subtract bias correction terms from expontiated mean terms

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
  
  if (bias.cor) {
	biaslogE <- -0.5*logESD^2
	biaslogAlpha <- -0.5*logAlphaSD^2
	biaslogRS <- -0.5*(sqrt(1/tauobs))^2
  } else {
	biaslogE <- 0
	biaslogAlpha <- 0
	biaslogRS <- numeric(N_Stk)
  }
  
  nll <- 0 # Begin negative log-likelihood
  
  # Can I remove the sum() from these arguments?
  nll <- nll - sum(dnorm(b0[1], 10, sd = 31.6, log = TRUE)) # Prior
  nll <- nll - sum(dnorm(b0[2], 0, sd = 31.6, log = TRUE)) # Prior
  nll <- nll - sum(dnorm(bWA[1], 0, sd = 31.6, log = TRUE)) # Prior
  nll <- nll - sum(dnorm(bWA[2], 0, sd = 31.6, log = TRUE)) # Prior
  
  nll <- nll - sum(dnorm(logAlpha0, 0.6, sd = 0.45, log = TRUE)) # Prior (rM)
  if(lhdiston) nll <- nll - sum(dnorm(logAlpha02, 0, sd = 31.6, log = TRUE)) # Prior (rD)
  
  # Half normal/half student-t/half cauchy for SD parameters
    # Half normal:
  # nll <- nll - dnorm(logAlphaSD, mean = 0, sd = 31.6, log = TRUE)
  # nll <- nll - dnorm(logESD, mean = 0, sd = 31.6, log = TRUE)
    # Half Student: Doesn't perform - particularly logESD
  # nll <- nll - dt(logAlphaSD, df = 1, log = TRUE)
  # nll <- nll - dt(logESD, df = 1, log = TRUE)
    # Half Cauchy: Still unsure how to pull - not a native TMB distribution
  
  ## Second level of hierarchy - Ricker parameters:
  for (i in 1:N_Stk){
    nll <- nll - dnorm(logE_re[i], 0, sd = 1, log = TRUE) # Median of E
	
    logE[i] <- b0[1] + b0[2]*type[i] + (bWA[1] + bWA[2]*type[i]) * WAbase$logWAshifted[i] + logE_re[i]*logESD + biaslogE
    E[i] <- exp(logE[i])
    
    # nll <- nll - dnorm(logAlpha_re[i], 0, sd  = logAlphaSD, log = TRUE) # Prior (hj)
	nll <- nll - dnorm(logAlpha_re[i], 0, sd = 1, log = TRUE) # Non-centered
    
    # Non-centered with scaling of mean by logAlphaSD
    if(lhdiston) logAlpha[i] <- logAlpha0 + logAlpha02*type[i] + logAlpha_re[i]*logAlphaSD + biaslogAlpha
    else logAlpha[i] <- logAlpha0 + logAlpha_re[i]*logAlphaSD + biaslogAlpha

    nll <- nll - dgamma(tauobs[i], shape = 0.0001, scale = 1/0.0001, log = TRUE)
  }

  ## First level of hierarchy: Ricker model:
  for (i in 1:N_Obs){
	logRS_pred[i] <- logAlpha[stk[i]]*(1 - S[i]/E[stk[i]]) + biaslogRS[stk[i]]

    if(!prioronly){
      nll <- nll - dnorm(logRS[i], logRS_pred[i], sd = sqrt(1/tauobs[stk[i]]), log = TRUE)
    } # If prioronly is 0, then likelihood is calculated
      # If prioronly is 1, then likelihood is not calculated
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
	# Conditioning on the random effect being equal to zero if bias.cor <- F
	# Conditioning on a random site, but the conditional mean
  for (i in 1:N_Pred){
    # NEW: alpha0 prior for LH specific dists.
    if(lhdiston) logAlpha_tar[i] <- logAlpha0 + logAlpha02*type_tar[i] # + biaslogAlpha + logAlphaSD^2/2
    else logAlpha_tar[i] <- logAlpha0 # + biaslogAlpha + logAlphaSD^2/2

    logE_tar[i] <- b0[1] + b0[2]*type_tar[i] + (bWA[1] + bWA[2]*type_tar[i])*WAin$logWAshifted_t[i] # + biaslogE + logESD^2/2
		# add -0.5*logESD^2 OR 
		# do the random normal when extracting the posterior
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
  alpha <- exp(logAlpha)
  
  REPORT(b0) # Testing simulate()
  REPORT(bWA) # Testing simulate()
  
  # ADREPORT(logRS_pred)
  REPORT(logRS_pred)

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
  REPORT(logAlpha0)
  REPORT(logAlpha02)
  REPORT(logAlpha_re) # random effect parameter for resampling
  REPORT(logAlphaSD)
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



# Alternative approach from Sean ####
  # rnorms for each of alpha prior
  # a function that takes in parameters (Vectored) and creates a logRS, etc., 
  # and then pipe through vector - recreating what REPORT() does
  # this can then be done for all forms of prior/posterior prediction
  # just avoids all the behind the scene work of the RTMB framework

# PRIOR PUSHFORWARD APPROACH: Random Generation ####
randomgenmethod <- F # Turn on or off for now
if(randomgenmethod == T){
  # Create however many vectors of parameters you wish to test/investigate
  geniters <- 1000
  psimalpha <- vector("list", geniters) # Vector start
  psimlogRS_pred <- vector("list", geniters) # Vector start
  psimlogAlpha_re <- vector("list", geniters) # Vector start
  
  for (i in 1:geniters){
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
  # Second: plot the distribution of logRS_pred with observation error
}



# LIMITS: Set Upper and Lower Limits ####
upper <- numeric(length(obj$par)) + Inf
lower <- numeric(length(obj$par)) + -Inf
lower[names(obj$par) == "tauobs"] <- 0
upper[names(obj$par) == "logESD"] <- 100 # Turn off for half dists.
lower[names(obj$par) == "logESD"] <- 0
upper[names(obj$par) == "logAlphaSD"] <- 100 # Turn off for half dists.
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
       
       # Should these be 0.01 to 100s? Would that be more accurate?
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



# SAMPLE MCMC ####
  # Can consider in parallel - Kasper's github example doesn't currently work
# cores <- parallel::detectCores() - 2
# cores <- 1
# options(mc.cores = cores)

# The set.seed precedes this in order to have fixed "identical" runs with initfixed
set.seed(1) ; fitstan <- tmbstan(obj, iter = 5000, warmup = 2500, # default iter/2 for warmup 
                                                                  # - Typically 5000 and 1000
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

# Matching Liermann et al. statistics
# mcmc_chains <- As.mcmc.list(fitstan)
# autocorr.plot(mcmc_chains)
# geweke_results <- geweke.diag(mcmc_chains)
# print(geweke_results)
# par(mfrow = c(2, 2))  # adjust depending on number of chains
# for (i in 1:length(mcmc_chains)) {
#   geweke.plot(mcmc_chains[[i]], main = paste("Chain", i))
# }
# heidel_results <- heidel.diag(mcmc_chains)
# print(heidel_results)
# effectiveSize(mcmc_chains)
# gelman.diag(mcmc_chains)


# Acquire outputs of MCMC ####
derived_obj <- derived_post(fitstan); beep(2)
# add stocknames - see extra code from _Plots.R

# SAVING R OBJECTS: ####
# In script_A.R
save(derived_obj, file = "derived_obj.RData")

if(dat$prioronly == 1){print("Prior Prediction Mode")} else {print("Posterior Prediction Mode")}


# PRIOR PUSHFORWARD/PREDICTIVE CHECKS ####
# APPROACH - Prior prediction with data exclusion
# Make sure dat$prioronly <- 1 to not include data in the nll
# Run all the way down to derived_obj

# Pushforward tests:
  # Just take the Prior Prediction Mode (dat$prioronly == 1)
  # Then plot desired priors/parameters
  # E.g. logAlpha
  # geom_hist() or geom_density() are potential options

# Predictive tests: 
  # Exactly the same as for Posterior - just depends
  # on whether data has been included in the likelihood or not.
  # Check:
  
  # Note: prior predictive checks (if taking the method of excluding data)
  # Can be repeated using the below posterior predictive checks.
  # Specifically: logRS, logAlpha



# POSTERIOR PREDICTIVE CHECKS ####
  # Run full model with dat$prioronly <- 0 for data included in nll
  # Run dervied_post() to extract posterior from chains
  # Randomly sample for logRS
  # Compare logRS with logRS_pred (including observation error)
  # e.g. rnorm(1, derived_obj$deripost_summary$logRS_pred$Median, 
             # sd  = sqrt(1/derived_obj$deripost_summary$tauobs$Median))

slogRS_pred <- derived_obj$deripost_full$logRS_pred
stauobs <- derived_obj$deripost_full$tauobs
simlogRS <- matrix(NA, nrow = dim(slogRS_pred)[1], ncol = dim(slogRS_pred)[2])
for (i in 1:dim(slogRS_pred)[1]){
  simlogRS[i, ] <- rnorm(dim(slogRS_pred)[2], # 501
                        mean = slogRS_pred[i, ], 
                        sd = sqrt(1/stauobs[i, stk]))
} # dim of simlogRS is 10

# Retrieve 9 samples out of the 10000 iterations
nsim <- 9
draws <- sample(1:dim(slogRS_pred)[1], nsim)
savedsims <- simlogRS[draws,] # [1:9,]

# Plot
par(mfrow = c(5, 2)) # 5 x 2 grid
plot(dat$logRS, xlab = "", ylab = "") # , col = stk - if you want stocks colored
    # Do I want to set xlim, ylim for visualization?
for (i in 1:nsim) {
  plot(simlogRS[i, ], xlab = "", ylab = "")
} 
mtext("Observation Index", side = 1, outer = TRUE, line = -2)
mtext("log(R/S)", side = 2, outer = TRUE, line = -1.5)
par(mfrow = c(1, 1))

# POSTERIOR P-VALUES
    # Take the mean and sd of each of the iterations e.g. do more than 9 in this case
    # Then plot the histogram of these means against the mean of
    # the observed data
    # Then the p-value is the proportion of simulated means that are
    # ABOVE the mean of the observations.
    # Report this value.

# Mean value of simulations: simlogRS
# Mean value of observed data: dat$logRS
mean_simlogRS <- apply(simlogRS, 1, mean) # Mean of all simulated logRS
  # this should produce 10,000 means
hist(mean_simlogRS, xlab = "Mean of Log(R/S)", 
     main = "Visualization of Posterior Predictive P-Value")
abline(v = mean(dat$logRS), col = "red", lwd = 3) # Mean of observed logRS

pvalpp <- sum(mean_simlogRS > mean(dat$logRS))/length(mean_simlogRS) # Proportion of simulated means above observed mean
print(paste0("Posterior Predictive P-Value = ", pvalpp))
  # Therefore a slightly higher bias than observed predictions



# Posterior Predictive Distribution: logAlpha
    # Manually write out and add rnorms for hyper parameters.
# I think I am interested in the following equation parts:
      # logAlpha[i] <- logAlpha0 + logAlpha02*type[i] + logAlpha_re[i]
      # and dnorm(logAlpha_re[i], 0, sd  = logAlphaSD, log = TRUE)
    # For posterior predictive it would be:
      # rnorms for each of the 3? assuming median value?
      # Liermann says "for populations with no SR data"
# 1. pull draws for logAlpha0
# 2. pull draws for logAlphaSD
# 3. iter for sims
slogAlpha0 <- derived_obj$deripost_full$logAlpha0 # dim(10000, 1)
slogAlpha02 <- derived_obj$deripost_full$logAlpha02 # dim(10000, 1)
slogAlphaSD <- derived_obj$deripost_full$logAlphaSD # dim(10000, 1)
simlogAlpha_s <- matrix(NA, nrow = dim(slogAlpha0)[1], ncol = dim(slogAlpha0)[2])
simlogAlpha_o <- matrix(NA, nrow = dim(slogAlpha0)[1], ncol = dim(slogAlpha0)[2])
for (i in 1:dim(slogAlpha0)[1]){
  simlogAlpha_s[i] <- slogAlpha0[i] + rnorm(1, mean = 0, sd = slogAlphaSD[i])
  simlogAlpha_o[i] <- simlogAlpha_s[i] + slogAlpha02[i]
}

ggplot() +
  geom_density(aes(x = simlogAlpha_s), # For a new STREAM observation
               fill = "skyblue", alpha = 0.4, color = "forestgreen", linewidth = 1.2) +
  geom_density(aes(x = simlogAlpha_o), # For a new OCEAN observation
              fill = "skyblue", alpha = 0.4, color = "skyblue", linewidth = 1.2) +
  theme_classic() +
  # labs(x = "Mean of uncentered logAlpha Posterior Predictive Distribution (Stream and Ocean)", 
         # y = "Density")
  labs(x = "Mean of uncentered logAlpha Prior Predictive Distribution (Stream and Ocean)", 
       y = "Density")

# Pushforward (under prior predictive): logAlpha ####
ggplot() +
  geom_density(aes(x = slogAlpha0), # For a pushforward alpha?
               fill = "skyblue", alpha = 0.4, color = "forestgreen", linewidth = 1.2) +
  theme_classic() +
  labs(x = "logAlpha Prior Pushforward Distribution", y = "Density")

# Posterior OR Prior Distribution ridge plot: logAlpha ####
  # IF Prior - then its pushforward - its a per stock value - NOT a new observ.
        # This is trash code - needs revision
dfalpharidge <- derived_obj$deripost_full$logAlpha[, 1:25]
Stocknames <- WAin$Stock
colnames(dfalpharidge) <- Stocknames
alpharidgetable <- as.data.table(dfalpharidge)
alpharidgetable_long <- melt(alpharidgetable, measure.vars = Stocknames,
                variable.name = "Stock", value.name = "Value")
TypeLabels <- ifelse(lifehist$lh == 0, "S", "O")
alpharidgetable_long[, Type := TypeLabels[match(Stock, Stocknames)]]
n_S <- sum(TypeLabels == "S")
n_O <- sum(TypeLabels == "O")

ggplot(alpharidgetable_long, aes(x = Value, y = Stock, fill = interaction(Type, Stock))) +
  geom_density_ridges(color = "gray20", alpha = 0.8, scale = 1.2) +
  theme_classic() +
  labs(x = "logAlpha", y = "Stock") +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      setNames(colorRampPalette(c("#a8e6a3", "#0b6e0b"))(n_S),
               paste("S.", Stocknames[TypeLabels == "S"], sep = "")),
      setNames(colorRampPalette(c("#a3d5ff", "#084a8b"))(n_O),
               paste("O.", Stocknames[TypeLabels == "O"], sep = ""))
    )
  )

# ggplot(alpharidge_long, aes(x = Value, y = factor(Index), fill = factor(Index))) +
#   geom_density_ridges(alpha = 0.5, scale = 1.2, color = "forestgreen") +
#   theme_classic() +
#   labs(x = "logAlpha Posterior", y = "Stock ID") +
#   theme(legend.position = "none")
# Change ordering of stocks e.g. by WA size as with point-wise comp. plots
# Can I add in ocean/stream-type identifiers and colour via?

# beta?

# Posterior Predictive Calculation for: SREP (E)
	# logE[i] <- b0[1] + b0[2]*type[i] + (bWA[1] + bWA[2]*type[i]) * WAbase$logWAshifted[i] + logE_re[i]*logESD
	# Need b0, bWA, WA, and random effects
	# Will get one for stream and one for ocean
	

# Posterior Distributions of b0 and bWA
ggplot() +
  geom_density(aes(x = derived_obj$deripost_full$b0[,1]), 
               fill = "skyblue", alpha = 0.4, color = "forestgreen", linewidth = 1.2) +
  geom_density(aes(x = (derived_obj$deripost_full$b0[,2] + derived_obj$deripost_full$b0[,1])),
              fill = "skyblue", alpha = 0.4, color = "skyblue", linewidth = 1.2) +
  theme_classic() +
  labs(x = "b0", y = "Density")

ggplot() +
  geom_density(aes(x = derived_obj$deripost_full$bWA[,1]), 
               fill = "skyblue", alpha = 0.4, color = "forestgreen", size = 1.2) +
  geom_density(aes(x = (derived_obj$deripost_full$bWA[,2] + derived_obj$deripost_full$bWA[,1])),
              fill = "skyblue", alpha = 0.4, color = "skyblue", size = 1.2) +
  theme_classic() +
  labs(x = "bWA", y = "Density")



# Plot curves?
  # E.g. instead of logRS, plot R vs. S
  # or Residuals (see below work) to compare observed logRS to 
  # In order to get back to R - you can use srdat$Sp and the base Ricker
  # equation - assuming some kind of randomness
# 1. Get raw data for R and S
	# srdat$Rec and srdat$Sp
# 2. Get out posterior predictive parameters for Ricker model PER STOCK
	# simlogRS is 501 random observations 
	# Alpha: simlogAlpha_s and simlogAlpha_o (each a draw for a single value for the different alpha's)
	# SREP (E): logE[i] <- b0[1] + b0[2]*type[i] + (bWA[1] + bWA[2]*type[i]) * WAbase$logWAshifted[i] + logE_re[i]*logESD
	# May want to sample 1000 from the 10,000 iterations for ease of viewing?
# X. USE POSTERIORS for STOCKS - WITHOUT RANDOM ERROR - JUST AS MODEL SEES IT?
	# derived_obj$deripost_summary$E
	# derived_obj$deripost_summary$logAlpha
# lineSREPmedian <- derived_obj$deripost_summary$E$Median
# lineAlphamedian <- exp(derived_obj$deripost_summary$logAlpha$Median)
# 3. Create an imaginary line of spawners
# 4. Use imaginary line to solve for recruits PER STOCK
	# e.g. IWAM: RR[j] <- a * SS[j] * exp(-b * SS[j])
	# logRS_pred[i] <- logAlpha[stk[i]]*(1 - S[i]/E[stk[i]])
	# logR/S <- loga * (1 - S/E)
	# R <-  S e ^ (a (1 - S/E) + w) # I think w is incorporated already?
	# RR[j] <- SS[j] * exp(alpha * (1 - S / SREP))
# create an object for the alphas - if it is random draw, it has to match the same number draw for SREP
# create an object for the SREPs
# SS <- c()
# for (j in 1:100) {SS[j] <- j*(max(srdat$Sp[srdat$Name == 'Harrison'])/100)} # where j is 1:100 (so a 100 point line)
# RR <- c()
# for (i in 1:25){
	# for (j in 1:100) {
		# SS[j] <- *(max(srdat$Sp[srdat$Name == unique(srdat$Name)[i]])/100)
		# RR[j] <- SS[j] * lineAlpha[1] ^ (1 - SS[j] / lineSREP[1])
	# }
# }
# for (j in 1:100) {
	# RR[j] <- SS[j] * exp(lineAlpha[1] * (1 - SS[j] / lineSREP[1]))
	# RR[j] <- SS[j] * lineAlpha[1] ^ (1 - SS[j] / lineSREP[1])
# }
# 5. Plot spawners by recruits PER STOCK
	# Do it for one stock first - and then iterate over
	# and will need a grid by parmfrow or with ggplot
# Base plot e.g. plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)), ylim=c(0,max(R$Rec) ) )
# plot(x = srdat$Sp[srdat$Name == 'Harrison'], y = srdat$Rec[srdat$Name == 'Harrison']) # Stock 1 of synoptic set
# lines(x = SS, y = RR)

lineSREPdraws <- derived_obj$deripost_full$E # 10000, 25
lineAlphadraws <- exp(derived_obj$deripost_full$logAlpha) # 10000, 25
SSdraws <- matrix(NA, nrow = 10000, ncol = 100) # this needs to be 10,000 by 100?
RRdraws <- matrix(NA, nrow = 10000, ncol = 100)
RRmedian <- lineAlphamedian <- lineSREPmedian <- SSmedian <- NA

rowsample <- sample(1:10000, 1000)

par(mfrow=c(5,5), mar=c(2, 2, 1, 0.1) + 0.1, oma=c(3,3,1,1))

for (i in 1:25){
	for (k in 1:10000){ # length(rowsample)
		for (j in 1:100){
			SSdraws[k, j] <- j*(max(srdat$Sp[srdat$Name == unique(srdat$Name)[i]])/100) # srdat$Name == 'Harrison'
			RRdraws[k, j] <- SSdraws[k, j] * lineAlphadraws[k, i] ^ (1 - SSdraws[k, j] / lineSREPdraws[k, i])
		}
	}

	lineAlphamedian <- median(lineAlphadraws[,i]) # Per stock median (single value)
	lineSREPmedian <- median(lineSREPdraws[,i])
	
	for (j in 1:100){
		SSmedian[j] <- median(SSdraws[,j]) # single over-written value
		RRmedian[j] <- SSmedian[j] * lineAlphamedian ^ (1 - SSmedian[j] / lineSREPmedian) # I need 100 of these for a line
	}

	plot(x = srdat$Sp[srdat$Name == unique(srdat$Name)[i]], 
		y = srdat$Rec[srdat$Name == unique(srdat$Name)[i]],
		xlim=c(0,max(srdat$Sp[srdat$Name == unique(srdat$Name)[i]])), 
		ylim=c(0,max(srdat$Rec[srdat$Name == unique(srdat$Name)[i]]) ) ) # Stock 1 of synoptic set
	
	mtext(unique(srdat$Name[srdat$Name == unique(srdat$Name)[i]]), side=3, cex=0.8)
	
	# subset matrix of draws
	SSsubdraws <- SSdraws[rowsample,]
	RRsubdraws <- RRdraws[rowsample,]
	
	for (k in 1:1000){ # subset for 1000 out of 10,000
		lines(SSsubdraws[k,], RRsubdraws[k,], col=rgb(0, 0, 0, alpha=0.1))
		lines(x = SSmedian, y = RRmedian, col = c("red")) # Median line
	}
	
	mtext("Spawners", side = 1, line = 1, outer = TRUE, cex = 1.3)
	mtext("Recruitment", side = 2, line  = 1, outer = TRUE, cex = 1.3, las = 0)
}

# Optimized version of the above
lineSREPdraws <- derived_obj$deripost_full$E # 10000, 25
lineAlphadraws <- exp(derived_obj$deripost_full$logAlpha) # 10000, 25
rowsample <- sample(1:10000, 1000)
par(mfrow = c(5, 5), mar = c(2, 2, 1, 0.1) + 0.1, oma = c(3, 3, 1, 1))

for (i in 1:25) {
  stock_name <- unique(srdat$Name)[i]
  spawners   <- srdat$Sp[srdat$Name == stock_name]
  recruits   <- srdat$Rec[srdat$Name == stock_name]
  Smax       <- max(spawners)

  SSseq <- seq(Smax/100, Smax, length.out = 100)
  alpha_draws <- lineAlphadraws[, i]
  srep_draws  <- lineSREPdraws[, i]

  SSmat <- matrix(SSseq, nrow = 10000, ncol = 100, byrow = TRUE)
  RRmat <- SSmat * alpha_draws^(1 - SSmat / srep_draws)

  RRmed <- SSseq * median(alpha_draws)^(1 - SSseq / median(srep_draws))

  plot(spawners, recruits, xlim = c(0, Smax), ylim = c(0, max(recruits)))
  mtext(stock_name, side = 3, cex = 0.8)
  matlines(t(SSmat[rowsample, ]), t(RRmat[rowsample, ]),
           col = rgb(0, 0, 0, 0.1), lty = 1)
  lines(SSseq, RRmed, col = "red", lwd = 2)
}

mtext("Spawners", side = 1, line = 1, outer = TRUE, cex = 1.3)
mtext("Recruitment", side = 2, line = 1, outer = TRUE, cex = 1.3)



# RESIDUALS? ####

# Plot generated logRS vs. logRS from data?
# plot(exp(dat$logRS), exp(derived_obj$deripost_summary$logRS_pred$Median))
plot(dat$logRS, derived_obj$deripost_summary$logRS_pred$Median)
points(dat$logRS, simlogRS[1,])

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
# plot(priorlogRS, ylab = "priorlogRS", ylim = c(-4, 4)) 
# Coming from original tmb obj$ no sampling, no nlminb

# Tor: this would make it seem like the model is not working great in terms
  # of generating observations from the Ricker model. It is more restricted
  # than the observed data.
