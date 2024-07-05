# IWAMsrep RTMB Model ####

# The following code includes the wrapper function and processing for the 
# IWAM srep model based on: NIMBLE_rtmb.R
# There are two possible model variants depending on:
  # Using a design matrix, or,
  # using explicit betas (intercepts and slopes)
# The following wrapper function will initially focus only on the explicit model.

# MODEL DESCRIPTION ####

# This is the model from Liermann et al. (CITATION YEAR).
# ....



# Libaries ####
library(RTMB)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidybayes)
library(tmbstan) 

# Imported LambertW0 ####
# FritschIter <- function(x, w){
#   MaxEval <- 5
#   CONVERGED <- FALSE
#   k <- 2.0 / 3.0;
#   i <- 0;
#   eps <- 2.2204460492503131e-16    
#   while (!CONVERGED & i < MaxEval){
#     z <- log(x / w) - w
#     w1 <- w + 1.0
#     q <- 2.0 * w1 * (w1 + k * z)
#     qmz <- q - z
#     e <- z / w1 * qmz / (qmz - z)
#     CONVERGED <- abs(e) <= eps
#     w <- w*(1.0 + e)
#     i <- i + 1
#   }
#   return(w)
# }
# 
# LambertW0_internal <- function(x){
#   check <- 0.367879441171442334024277442949824035167694091796875 # exp(-1)
#   eps <- 2.2204460492503131e-16
#   if (x == Inf) {
#     return(Inf);
#   } else if (x < -check) {
#     return(NaN);
#   } else if (abs(x - check) <= eps) {
#     return(-1.0);
#   } else if (abs(x) <= 1e-16) {
#     ## This close to 0 the W_0 branch is best estimated by its Taylor/Pade
#     ## expansion whose first term is the value x and remaining terms are below
#     ## machine double precision. See
#     ## https://math.stackexchange.com/questions/1700919
#     
#     return(x);
#   } else {
#     w <- 0
#     if (abs(x) <= 6.4e-3) {
#       ## When this close to 0 the Fritsch iteration may underflow. Instead,
#       ## function will use degree-6 minimax polynomial approximation of Halley
#       ## iteration-based values. Should be more accurate by three orders of
#       ## magnitude than Fritsch's equation (5) in this range.
#       return((((((-1.0805085529250425e1 * x + 5.2100070265741278) * x -
#                    2.6666665063383532) * x + 1.4999999657268301) * x -
#                  1.0000000000016802) * x + 1.0000000000001752) * x +
#                2.6020852139652106e-18);
#       
#     } else if (x <= exp(1)) {
#       ## Use expansion in Corliss 4.22 to create (2, 2) Pade approximant.
#       ## Equation with a few extra terms is:
#       ## -1 + p - 1/3p^2 + 11/72p^3 - 43/540p^4 + 689453/8398080p^4 - O(p^5)
#       ## This is just used to estimate a good starting point for the Fritsch
#       ## iteration process itself.
#       
#       p <- sqrt(2.0 * (exp(1) * x + 1.0))
#       Numer <- (0.2787037037037037 * p + 0.311111111111111) * p - 1.0;
#       Denom <- (0.0768518518518518 * p + 0.688888888888889) * p + 1.0;
#       w <- Numer / Denom;
#     } else {
#       ## Use first five terms of Corliss et al. 4.19 */
#       w <- log(x)
#       L_2 <- log(w)
#       L_3 <- L_2 / w
#       L_3_sq <- L_3 * L_3
#       w <- w - L_2 + L_3 + 0.5 * L_3_sq - L_3 / w + L_3 / (w * w) - 1.5 * L_3_sq /
#         w + L_3_sq * L_3 / 3.0;
#     }
#     return(FritschIter(x, w));
#   }
# }
# 
# ## Derivatives of LamW
# dLambertW0_internal <- function(x, y, dy) {
#   dy / (exp(y) * (1. + y))
# }
# 
# LambertW0 <- RTMB:::ADjoint(LambertW0_internal, dLambertW0_internal)
 
# Wrapper Function ####
# compiler::enableJIT(0) # Run first without just to see if bug is fixed yet

# Internal running
WAin <- c("DataIn/WCVIStocks.csv")

IWAMsrep_rtmb <- function(WAin = c("DataIn/WCVIStocks.csv") # ,
                      # remove.EnhStocks = FALSE,
                      # run.bootstraps = TRUE, # to turn on or off the bootstrap function added at the end
                      # bs_seed = 1, # seed for bootstrapping
                      # bs_nBS = 10, # trials for bootstrapping
                      # plot = FALSE # whether or not to create plots stored in DataOut/
)
{

  # Just in case atm
  compiler::enableJIT(0) # Run first without just to see if bug is fixed yet

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
  
  ## Derivatives of LamW
  dLambertW0_internal <- function(x, y, dy) {
    dy / (exp(y) * (1. + y))
  }
  
  LambertW0 <- RTMB:::ADjoint(LambertW0_internal, dLambertW0_internal)
  
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

## Shift log WA for the mean - base
mean_logWA <- mean(WAbase$logWA)
WAbase$logWAshifted <- WAbase$logWA - mean_logWA

## Shift log WA for the mean - targets
mean_logtarWA <- mean(log(WAin$WA))
WAin$logtarWAshifted <- log(WAin$WA) - mean_logtarWA

lifehist <- srdat %>% dplyr::select(Stocknumber, Name, Stream) %>% 
  group_by(Stocknumber) %>% 
  summarize(lh=max(Stream))


## RTMB dat and par setup ####

# Dat
dat <- list(srdat = srdat,
            WAbase = WAbase,
            WAin = WAin,
            logRS = log(srdat$Rec) - log(srdat$Sp))

# External vectors
N_Stk <- max(srdat$Stocknumber + 1)

# Parameters
  # TK: Add in comments to say what the pars are - similar to ModelBook comparisons so it is easier to read first time
par <- list(b0 = c(10, 10), # WA regression intercept initial value
              # b0 = c(9, 9)
             bWA = c(0, 0), 
              # bWA = c(0.83, 1)
             # logAlpha = numeric(N_Stk), # comment on or off depending if using "nll" or not
             logAlpha_re = numeric(nrow(dat$WAbase)),
             logAlpha_re_pred = numeric(nrow(dat$WAin)),
             logE0 = numeric(N_Stk),
             tauobs = 0.01 + numeric(N_Stk), # Why can't this be zero? This doesn't run as just a string of zeros.
             logAlpha0 = 1.5,
             logESD = 1,
             logAlphaSD = 10
)

# TK/CH: Add a version with a log-normal back-transformation adjustment. Unless there is a citable reference
  # to explain why.

f_srep <- function(par){
  getAll(dat, par)
  
  N_Stk = max(srdat$Stocknumber + 1)
  stk = srdat$Stocknumber + 1
  N_Obs = nrow(srdat)
  N_Pred = nrow(WAin) # number of predicted watershed areas
  S = srdat$Sp
  type = lifehist$lh + 1
  
  E <- numeric(N_Stk)
  # log_E <- numeric(N_Stk) 
  E_tar <- numeric(N_Pred)
  
  logAlpha <- numeric(N_Stk) # comment on or off
  logAlpha_tar <- numeric(N_Pred)
  
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
  
  # PREDICTIONS
    # Predict target watershed areas
  BETA = numeric(nrow(WAin))
  SMSY = numeric(nrow(WAin))
  SGEN = numeric(nrow(WAin))
  
  for (i in 1:N_Pred){
    # Do I need to create a new alpha parameter to estimate here
      # And would it need to be added to the nll
    nll <- nll - sum(dnorm(logAlpha_re_pred[i], 0, sd = logAlphaSD, log = TRUE)) # new random effect?
    logAlpha_tar[i] <- logAlpha0 + logAlpha_re_pred[i]
    # Carrie: Use the MLE of the hyper-logalpha
    
    # Predict E for target watershed areas
    log_E_tar <- b0[type[i]] + bWA[type[i]]*WAin$logtarWAshifted[i] + logE0[i] ## Stock level regression
    E_tar[i] <- exp(log_E_tar)
    
    # Predict BETA
    BETA[i] <- logAlpha_tar[i]/E_tar[i]
    # Predict SMSY
    SMSY[i] <- (1-LambertW0(exp(1-logAlpha_tar[i])))/BETA[i]
    # Predict SGEN
    SGEN[i] <- -1/BETA[i]*LambertW0(-BETA[i]*SMSY[i]/(exp(logAlpha_tar[i])))
  }
  
  # Predict along a preset line for plotting
  
  # ADREPORT - internals
  # ADREPORT(logRS) # logRS for all 501 data points
    # I think it might be worth it not to report for model speed?
    # Question: Does ADREPORT slow down RTMB? It is more to calculate.
  ADREPORT(E) # E (Srep) for all synotopic data set rivers (25)
  ADREPORT(logAlpha) # model logAlpha (25)
  
  # ADREPORT - predicted
  ADREPORT(E_tar) # target E (Srep) (21)
  ADREPORT(logAlpha_tar)
  ADREPORT(BETA)
  ADREPORT(SMSY)
  ADREPORT(SGEN)
  # you need to do everything else internally so as to make sure of ADREPORT and se's
  
  nll # nll must be the final line of the function
}

## MakeADFun ####
obj <- RTMB::MakeADFun(f_srep, par, random = c("logAlpha_re", "logAlpha_re_pred", "logE0"))
opt <- nlminb(obj$par, obj$fn, obj$gr)
  # 463.0525 - reference target

sdr <- sdreport(obj)
sdr_full <- summary(RTMB::sdreport(obj))

sdr_est <- as.list(sdr, "Est", report=TRUE) ## ADREPORT estimates
sdr_se <- as.list(sdr, "Std", report=TRUE) ## ADREPORT standard error

## Outputs and Plotting Preparations (interior)
  # IWAMsmax produces the following:
  # - SR curve
  # - WA regression plots
  # - csv files of target WA and associated target estimates of SMSY and SREP
  # * parameter lists for targets supplied
    # * E
    # * logAlpha
  # * parameters to be calculated for targets
    # * BETA = logalpha/E
    # * SMAX = 1/beta
    # * SMSY = (1-LambertW0(exp(1-logalpha)))/beta # EXPLICIT
    # * SGEN = -1/beta*LambertW0(-beta*Smsy[i,j]/(exp(logalpha))) # EXPLICIT

# Step 1: Create a dataframe for the target populations by joining WAin 
  # with $E_tar
WArefpoints <- cbind(WAin, 
                     E = sdr_est$E_tar,
                     E_se = sdr_se$E_tar,
                     logalpha = sdr_est$logAlpha_tar,
                     logalpha_se = sdr_se$logAlpha_tar,
                     BETA = sdr_est$BETA,
                     BETA_se = sdr_se$BETA,
                     SMSY = sdr_est$SMSY,
                     SMSY_se = sdr_se$SMSY,
                     SGEN = sdr_est$SGEN,
                     SGEN_se = sdr_se$SGEN)

# Step 2: Calculate the new columns of BETA, SMAX, SMSY, and SGEN to populate
  # df with
  # * Are we using a mean/median global alpha?
  # ** Or should that be a new parameter to be assessed in the model itself?

# Step 3: Print out

return(list(opt = opt,
            obj = obj,
            sdr = sdr,
            sdr_full = sdr_full,
            sdr_est = sdr_est,
            sdr_se = sdr_se,
            refpoints = WArefpoints
))

} # End of IWAMsrep_rtmb

test <- IWAMsrep_rtmb() # default test run for outputs
