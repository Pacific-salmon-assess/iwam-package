# Liermann Srep (E) RTMB Model with MCMC Sampling ####

# Libaries ####
library(RTMB)
library(ggplot2)
library(gridExtra)
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

here::i_am("R/LambertWs.R") # For non-RStudio functionality
source(here::here("R/LambertWs.R")) # Lambert W function
source(here::here("R/helperFunctions.R")) # For bootstrapping
source(here::here("R/derived_post.R")) # For posterior extraction

# New LambertW0 see: https://github.com/pbs-assess/renewassess/blob/main/code/RTMB/PosteriorPredictiveSample.r
LambertW0 <- ADjoint(
  function(x){gsl::lambert_W0(x)},
  function(x, y, dy) {dy / (x + exp(y))}
)

# Raw data read-in ####
WAin <- c("DataIn/Parken_evalstocks.csv") # c("DataIn/WCVIStocks.csv")

# For predicting/re-evaluation of synoptic sets: WAbase
    # 1. Run model until setting up data section
    # 2. Then over-write WAin <- WAbase
    # 3. And rename logWAshifted to logWAshifted_t
    # 4. And make sure to change type_tar to fit the expected 0:1 indexing

# Data Manipulations ####
srdatwna <- read.csv(here::here("DataIn/SRinputfile.csv")) # Consider _TK ************************************
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
                          max(WAbase$logWAshifted), 0.1),
			mean_logWA = mean_logWA,
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
            logSMAX_re = numeric(N_Stk), # Zeroes # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            logAlpha0 = 0.6,
            logAlpha_re = numeric(nrow(dat$WAbase)), # Zeroes
            tauobs = 0.01 + numeric(N_Stk), # Constrained positive
            logSMAX_sd = 1, # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            logAlpha_sd = 1
)
if (lhdiston) {
  par$logAlpha02 <- 0
}

f_smax <- function(par){
  getAll(dat, par)
  
  N_Stk = max(srdat$Stocknumber + 1) # number of stocks
  stk = srdat$Stocknumber + 1 # vector of stocknumbers
  N_Obs = nrow(srdat) # number of observations
  N_Pred = nrow(WAin) # number of predicted watershed areas
  
  S = srdat$Sp
  type = lifehist$lh
  type_tar = as.numeric(WAin$lh) 

  SMAX <- numeric(N_Stk) # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  logSMAX <- numeric(N_Stk) # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  logAlpha <- numeric(N_Stk)
  
  # Why is logRS_pred not a parameter or vector input here?
  logRS_pred <- numeric(N_Obs) # Does this still report if not a vector?

  SMAX_tar <- numeric(N_Pred) # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  logSMAX_tar <- numeric(N_Pred) # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  logAlpha_tar <- numeric(N_Pred)
  
  # Simulated line vectors
  line <- length(lineWA)
  logSMAX_line_stream <- numeric(line) # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  SMAX_line_stream <- numeric(line) # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  logSMAX_line_ocean <- numeric(line) # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  SMAX_line_ocean <- numeric(line) # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  if (bias.cor) {
	biaslogSMAX <- -0.5*logSMAX_sd^2 # Global # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	biaslogAlpha <- -0.5*logAlpha_sd^2 # Global
	biaslogRS <- -0.5*(sqrt(1/tauobs))^2 # Stock-specific
  } else {
	biaslogSMAX <- 0 # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	biaslogAlpha <- 0
	biaslogRS <- numeric(N_Stk)
  }
  
  nll <- 0
  
  # Can I remove the sum() from these arguments?
  nll <- nll - sum(dnorm(b0[1], 10, sd = 31.6, log = TRUE)) # Prior b0
  nll <- nll - sum(dnorm(b0[2], 0, sd = 31.6, log = TRUE)) # Prior b0
  nll <- nll - sum(dnorm(bWA[1], 0, sd = 31.6, log = TRUE)) # Prior bWA
  nll <- nll - sum(dnorm(bWA[2], 0, sd = 31.6, log = TRUE)) # Prior bWA
  
  nll <- nll - sum(dnorm(logAlpha0, 0.6, sd = 0.45, log = TRUE)) # Prior (rM)
  if(lhdiston) nll <- nll - sum(dnorm(logAlpha02, 0, sd = 31.6, log = TRUE)) # Prior (rD)
  
  ## Second level of hierarchy - Ricker parameters:
  for (i in 1:N_Stk){
    nll <- nll - dnorm(logSMAX_re[i], 0, sd = 1, log = TRUE) 
	
    logSMAX[i] <- b0[1] + b0[2]*type[i] + (bWA[1] + bWA[2]*type[i]) * WAbase$logWAshifted[i] + logSMAX_re[i]*logSMAX_sd + biaslogSMAX # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SMAX[i] <- exp(logSMAX[i]) # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
	nll <- nll - dnorm(logAlpha_re[i], 0, sd = 1, log = TRUE)
    if(lhdiston) logAlpha[i] <- logAlpha0 + logAlpha02*type[i] + logAlpha_re[i]*logAlpha_sd + biaslogAlpha
    else logAlpha[i] <- logAlpha0 + logAlpha_re[i]*logAlpha_sd + biaslogAlpha

    nll <- nll - dgamma(tauobs[i], shape = 0.0001, scale = 1/0.0001, log = TRUE)
  }

  ## First level of hierarchy: Ricker model:
  for (i in 1:N_Obs){
	# logRS_pred[i] <- logAlpha[stk[i]]*(1 - S[i]/SREP[stk[i]]) + biaslogRS[stk[i]] # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	# logRS_pred[i] <- alpha_pred[stk[i]]*(1 - S[i]/SREP[stk[i]]) + biaslogRS[stk[i]] # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	alpha_pred <- exp(logAlpha) 
	# logRS_pred[i] <- exp(logAlpha[stk[i]]) - S[i]/SMAX[stk[i]]
	logRS_pred[i] <- alpha_pred[stk[i]] - S[i]/SMAX[stk[i]] + biaslogRS[stk[i]]
	# logRS_pred[i] <- logAlpha[stk[i]] - S[i]/SMAX[stk[i]] + biaslogRS[stk[i]]


    if(!prioronly){ # If prioronly is 1, then likelihood is not calculated, if 0 then it is
      nll <- nll - dnorm(logRS[i], logRS_pred[i], sd = sqrt(1/tauobs[stk[i]]), log = TRUE) # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    } 
  }
  
  ## Calculate SMSY for Synoptic set - for plotting
  SMSY_r = numeric(nrow(WAbase))
  BETA_r = numeric(nrow(WAbase))
  SREP_r = numeric(nrow(WAbase))
  
  for (i in 1:N_Stk){
    # BETA_r[i] <- logAlpha[i] / SREP[i] # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	BETA_r[i] <- 1/SMAX[i] # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SMSY_r[i] <- (1 - LambertW0(exp(1 - alpha_pred[i]))) / BETA_r[i]
	SREP_r[i] <- alpha_pred[i]/BETA_r[i] # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  }

  ## PREDICTIONS
  BETA = numeric(nrow(WAin))
  SMSY = numeric(nrow(WAin))
  SGEN = numeric(nrow(WAin))
  SREP = numeric(nrow(WAin))

  for (i in 1:N_Pred){
    # NEW: alpha0 prior for LH specific dists.
    if(lhdiston) logAlpha_tar[i] <- logAlpha0 + logAlpha02*type_tar[i] + biaslogAlpha # + biaslogAlpha + logAlpha_sd^2/2
    else logAlpha_tar[i] <- logAlpha0 + biaslogAlpha # + biaslogAlpha + logAlpha_sd^2/2

    logSMAX_tar[i] <- b0[1] + b0[2]*type_tar[i] + (bWA[1] + bWA[2]*type_tar[i])*WAin$logWAshifted_t[i] + biaslogSMAX # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    SMAX_tar[i] <- exp(logSMAX_tar[i]) # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
	alpha_tar <- exp(logAlpha_tar)
    # Predict BETA
    # BETA[i] <- logAlpha_tar[i]/SREP_tar[i] # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	BETA[i] <- 1/SMAX_tar[i] # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    # Predict SMSY
    SMSY[i] <- (1-LambertW0(exp(1-alpha_tar[i])))/BETA[i]
    # Predict SGEN
    SGEN[i] <- -1/BETA[i]*LambertW0(-BETA[i]*SMSY[i]/(exp(alpha_tar[i])))
	# Predict SREP
	SREP[i] <- alpha_tar[i]/BETA[i]
  }
  
  # Create predictions on an simulated line
  for (i in 1:line){
    logSMAX_line_ocean[i] <- b0[1] + b0[2] + (bWA[1] + bWA[2])*lineWA[i] + biaslogSMAX
    SMAX_line_ocean[i] <- exp(logSMAX_line_ocean[i])
    
    logSMAX_line_stream[i] <- b0[1] + (bWA[1])*lineWA[i] + biaslogSMAX 
    SMAX_line_stream[i] <- exp(logSMAX_line_stream[i])
  }
  
  ## ADREPORT - internal values (synoptic specific/Ricker) 
  
  REPORT(b0) # Testing simulate()
  REPORT(bWA) # Testing simulate()
  
  REPORT(logRS_pred)
  
  # alpha <- exp(logAlpha)
  # REPORT(logRS) # logRS for all 501 data points
  REPORT(logSMAX_re)
  REPORT(logSMAX_sd)
  REPORT(SMAX) # E (Srep) for all synoptic data set rivers (25)
  REPORT(logSMAX)
  REPORT(logAlpha) # model logAlpha (25)
  REPORT(logAlpha0)
  REPORT(logAlpha02)
  REPORT(logAlpha_re) # random effect parameter for resampling
  REPORT(logAlpha_sd)
  REPORT(alpha_pred) # Also called alpha_pred
  REPORT(SMSY_r)
  REPORT(BETA_r)
  REPORT(SREP_r) # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  REPORT(tauobs) # Necessary to add back in observation error?
  
  # ADREPORT - predicted values from watershed area model
    # Mean estimate of the median (without bias correction)
  # alpha_tar <- exp(logAlpha_tar)
  REPORT(SMAX_tar)
  REPORT(logSMAX_tar)
  REPORT(logAlpha_tar)
  REPORT(alpha_tar)
  
  REPORT(BETA)
  REPORT(SMSY)
  REPORT(SGEN)
  REPORT(SREP) # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  # Simulated line values for plotting
  REPORT(SMAX_line_stream) 
  REPORT(logSMAX_line_stream) 
  REPORT(SMAX_line_ocean) 
  REPORT(logSMAX_line_ocean)
  
  nll
}

## MakeADFun ####
obj <- RTMB::MakeADFun(f_smax, par,
    random = c("logAlpha_re", "logSMAX_re"), # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    silent=TRUE)

# LIMITS: Set Upper and Lower Limits ####
upper <- numeric(length(obj$par)) + Inf
lower <- numeric(length(obj$par)) + -Inf
lower[names(obj$par) == "tauobs"] <- 0
upper[names(obj$par) == "logSMAX_sd"] <- 100 # Turn off for half dists. # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
lower[names(obj$par) == "logSMAX_sd"] <- 0 # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
upper[names(obj$par) == "logAlpha_sd"] <- 100 # Turn off for half dists.
lower[names(obj$par) == "logAlpha_sd"] <- 0

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
# RANDOM INIT - MATCHING PRIORS
init <- function() {
  # Can also add par <- to sequence above - see prior testing
  listinit <- list(b0 = c(rnorm(1, 10, 1), rnorm(1, 0, 1)), # Contains negatives
       bWA = c(rnorm(1, 0, 1), rnorm(1, 0 ,1)), # Contains negatives
       
       # logRS_pred = rnorm(N_Obs, 0, 1),
       logSMAX_re = rnorm(N_Stk, 0, 1), # Contains negatives # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       logAlpha0 = rnorm(1, 0.6, 1), # Contains negatives
       # logAlpha02 = rnorm(1, 0, 1) , # NEW: alpha0 prior for LH specific dists.
       logAlpha_re = rnorm(nrow(dat$WAbase), 0, 1), # Contains negatives

       tauobs = runif(N_Stk, min = 0.005, max = 0.015), # Uniform to REMAIN positive
       
       # Should these be 0.01 to 100s? Would that be more accurate?
       logSMAX_sd = runif(1, 0.01, 3), # Positive # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
       logAlpha_sd = runif(1, 0.01, 3) # Positive
  )
  
  if (lhdiston) {
    listinit$logAlpha02 <- rnorm(1, 0, 1) # alpha0 prior for LH specific dists.
  }
  
  return(listinit)
}

# SAMPLE MCMC ####
# The set.seed precedes this in order to have fixed "identical" runs with initfixed
set.seed(1) ; fitstan <- tmbstan(obj, iter = 5000, warmup = 2500, # default iter/2 for warmup 
                                                                  # - Typically 5000 and 1000
                   init = init, # init = init function or "random" or initfixed for fixed points
                   seed = 1, # set seed or leave out for random - now set at 1 above
                      # but not technically within sampling - unsure if each chain takes seed                  
                   # control = list(adapt_delta = 0.99, max_treedepth = 15),
                   lower = lower, upper = upper,
                   chains = 4, open_progress = FALSE, silent = TRUE); beep(2)
# tmbstan operates by default with NUTS MCMC sampler see: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0197954

# Save stan fit object
# save(fitstan, file = "fitstan.RData")
# if(dat$prioronly == 1) {save(fitstan, file = "fitstan_prioronly.RData")} else {save(fitstan, file = "fitstan.RData")}

# Test and diagnostic plots ####
# names(fitstan) for complete list of parameters from stan object

  # TRACE PLOTS
# mcmc_trace(as.array(fitstan)) # ALL PARAMETER TRACEPLOT
# mcmc_trace(as.array(fitstan), regex_pars = "b0") # Single regex parameter traceplot e.g. b0

  # PAIRS PLOTS
# pairs_pars <- c("b0", "bWA", "logAlpha0") # "logAlpha_sd" "logSREP_sd"
# pairs(fitstan, pars = pairs_pars) # for specific par names from above
# bayesplot::mcmc_pairs(fitstan, regex_pars = pairs_pars)
# bayesplot::mcmc_pairs(fitstan, regex_pars = c("tauobs"))
# bayesplot::mcmc_pairs(fitstan, regex_pars = c("logSREP_re"))
# bayesplot::mcmc_pairs(fitstan, regex_pars = c("logAlpha_re"))

  # ACF
# bayesplot::mcmc_acf(fitstan, regex_pars = "b0") # ACF plot for b0
# bayesplot::mcmc_acf(fitstan) # For every and all - should be a faster way to save them

  # RHAT
# fitstan |> rhat() |> mcmc_rhat() + yaxis_text() # rhat plot for assessing rhat of each parameter

# Acquire outputs of MCMC ####
derived_obj <- derived_post(fitstan, model = 'SMAX'); beep(2)
	# add stocknames - see extra code from _Plots.R
	dsmax <- derived_obj
	fitsmax <- fitstan

# SAVING R OBJECTS: ####
# save(derived_obj, file = "derived_obj.RData")
# if(dat$prioronly == 1) {save(derived_obj, file = "derived_obj_prioronly.RData")} else {save(derived_obj, file = "derived_obj.RData")}
if(dat$prioronly == 1){print("Prior Prediction Mode")} else {print("Posterior Prediction Mode")}


