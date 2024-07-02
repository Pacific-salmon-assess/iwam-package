# IWAMsrep RTMB Model ####

# The following code includes the wrapper function and processing for the 
# IWAM srep model based on: NIMBLE_rtmb.R
# There are two possible model variants depending on:
  # Using a design matrix, or,
  # using explicit betas (intercepts and slopes)
# The following wrapper function will initially focus only on the explicit model.

# Libaries ####
library(RTMB)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidybayes)
library(tmbstan) 

# Wrapper function start ####
compiler::enableJIT(0) # Run first without just to see if bug is fixed yet

# Internal running saves
WAin <- c("DataIn/WCVIStocks.csv")

IWAMsrep_rtmb <- function(WAin = c("DataIn/WCVIStocks.csv") # ,
                      # remove.EnhStocks = FALSE,
                      # run.bootstraps = TRUE, # to turn on or off the bootstrap function added at the end
                      # bs_seed = 1, # seed for bootstrapping
                      # bs_nBS = 10, # trials for bootstrapping
                      # plot = FALSE # whether or not to create plots stored in DataOut/
)
{

# Data ####
srdatwna <- read.csv(here::here("DataIn/SRinputfile.csv"))
WAbase <- read.csv(here::here("DataIn/WatershedArea.csv"))
# WAin <- read.csv(here::here("DataIn/WCVIStocks.csv"))
WAin <- read.csv(here::here(WAin))

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

# Dat
dat <- list(srdat = srdat,
            WAbase = WAbase,
            WAin = WAin,
            # type = lifehist$lh + 1,
            logRS = log(srdat$Rec) - log(srdat$Sp))

# External vectors
N_Stk <- max(srdat$Stocknumber + 1)

# Parameters
par <- list(b0 = c(10, 10), # b0 = c(9, 9)
             bWA = c(0, 0), # bWA = c(0.83, 1)
             # logAlpha = numeric(N_Stk),
             logAlpha_re = numeric(nrow(dat$WAbase)),
             logE0 = numeric(N_Stk),
             tauobs = 0.01 + numeric(N_Stk), # Why can't this be zero? This doesn't run as just a string of zeros.
             logAlpha0 = 1.5,
             logESD = 1,
             logAlphaSD = 10
)

f_srep <- function(par){
  getAll(dat, par)
  
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
  # if (penalized = TRUE){
  #   nll <- nll - sum(dnorm(logAlpha0, 0.6, 0.45, log = TRUE))
  # }
  
  # Alternative version of alpha terms - see f_nim's translation
  # Slope and intercept priors
  # nll <- nll - sum(dnorm(b0[1], 0, sd = 31.6, log = TRUE))
  # nll <- nll - sum(dnorm(b0[2], 0, sd = 31.6, log = TRUE))
  # nll <- nll - sum(dnorm(bWA[1], 0, sd = 31.6, log = TRUE))
  # nll <- nll - sum(dnorm(bWA[2], 0, sd = 31.6, log = TRUE))
  
  ## Watershed Model
  for ( i in 1:N_Stk){
    # if (penalized = TRUE){
    # Run logAlpha0 as a penalty term/prior
    # New version
    # nll <- nll - sum(dnorm(logAlpha[i], logAlpha0, sd = logAlphaSD, log = TRUE)) # random effect - is this the bayesian way?
    # }
    
    nll <- nll - sum(dnorm(logAlpha_re[i], 0, sd = logAlphaSD, log = TRUE)) # random effect
    logAlpha[i] <- logAlpha0 + logAlpha_re[i]
    
    nll <- nll - sum(dnorm(logE0[i], 0, sd = logESD, log = TRUE)) # random effect
    log_E <- b0[type[i]] + bWA[type[i]]*WAbase$logWAshifted[i] + logE0[i] ## Stock level regression
    E[i] <- exp(log_E)
    
    nll <- nll - sum(dgamma(tauobs[i], shape = 0.001, scale = 0.001))
  }

  ## Ricker Model
  for (i in 1:N_Obs){
    logRS_pred <- logAlpha[stk[i]]*(1 - S[i]/E[stk[i]])
    nll <- nll - sum(dnorm(logRS[i], logRS_pred, sd = sqrt(1/tauobs[stk[i]]), log = TRUE))
  }
  
  # ADREPORT - estimates with standard error
  ADREPORT(logRS)
  ADREPORT(E)
  ADREPORT(logAlpha)

  nll
  
}

## MakeADFun ####
obj <- RTMB::MakeADFun(f_srep, par, random = c("logAlpha_re", "logE0"))
opt <- nlminb(obj$par, obj$fn, obj$gr)

sdr <- sdreport(obj)
sdr_full <- summary(RTMB::sdreport(obj))

sdr_est <- as.list(sdr, "Est", report=TRUE) ## ADREPORT estimates
sdr_se <- as.list(sdr, "Std", report=TRUE) ## ADREPORT standard error

return(list(opt = opt,
            obj = obj,
            sdr = sdr,
            sdr_full = sdr_full,
            sdr_est = sdr_est,
            sdr_se = sdr_se
))

} # End of IWAMsrep_rtmb

test <- IWAMsrep_rtmb() # default test run for outputs
