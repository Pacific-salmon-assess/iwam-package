# TMB --> RTMB

# Libaries ####

library(RTMB)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidybayes)
library(tmbstan)

# Imported LambertW0 ####

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

WAin <- read.csv(here::here("DataIn/WCVIStocks.csv"))

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

# Bring back scaling
source (here::here("R/helperFunctions.R"))
srdat <- digit_scaling(srdat)
srdat_scale <- srdat %>% dplyr::select(Stocknumber, scale) %>% distinct()
srdat_scale <- srdat_scale$scale
scale_TMB <- srdat$scale

# Bring back life history reference
lifehist <- srdat %>% dplyr::select(Stocknumber, Name, Stream) %>% 
  group_by(Stocknumber) %>% 
  summarize(lh=max(Stream)) %>% 
  arrange (Stocknumber)

WAbase <- WAbase %>% 
  full_join(names, by="Name") %>% 
  arrange(Stocknumber) %>%
  mutate(logWA = log(WA)) 

## Shift log WA for the mean.
mean_logWA <- mean(WAbase$logWA)
WAbase$logWAshifted <- WAbase$logWA - mean_logWA

## RTMB dat and par setup ####
# Create Dat list
dat <- list(srdat = srdat,
            WAbase = WAbase,
            
            logRS = log((srdat$Rec/srdat$scale)/(srdat$Sp/srdat$scale)),
            lifehist = lifehist$lh,
            
            scale = srdat_scale,
            scale_TMB = srdat$scale,
            S = srdat$Sp/scale_TMB,
            
            Tau_dist = 0.1,
            Tau_D_dist = 1,
            
            logMuA_stream_mean = 1.5,
            logMuA_stream_sig = 2,
            logMuA_ocean_mean = 0,
            logMuA_ocean_sig = 2
            )

# Create a data list for the targets whether they are stream or ocean
dat$pred_lnWA <- seq(min(log(WAbase$WA)), max(log(WAbase$WA)), 0.1)

dat$target_lnWA_ocean <- WAin %>% 
  mutate (lnWA=log(WA)) %>%
  filter(lh==1) %>% 
  pull(lnWA)

dat$target_lnWA_stream <- WAin %>%
  mutate (lnWA=log(WA)) %>%
  filter(lh==0) %>%
  pull(lnWA)

# Create Par list
B <- srdat %>% group_by(Stocknumber) %>% 
  summarise(m = -lm(log(Rec/Sp) ~ Sp)$coef[2])

par <- list(logA = (srdat %>% group_by (Stocknumber) %>% 
                          summarise(yi = lm(log(Rec / Sp) ~ Sp)$coef[1]))$yi, # random effect
            logB = log ( 1/ ( (1/B$m)/dat$scale )), # fixed effect
            logSigma = rep(-2, length(unique(srdat$Name))),
            logSigmaA = -2,
            
            logMuA_stream = 1.5,
            logMuA_ocean = 0,
            
            logDelta1 = 3,
            logDelta1_ocean = 0,
            logDelta2 = log(0.72),
            Delta2_ocean = 0,
            
            logNu1 = 3,
            logNu1_ocean = 0,
            logNu2 = log(0.72),
            Nu2_ocean = 0,
            
            logNuSigma = -0.412,
            logDeltaSigma = -0.412
)

# RTMB function ####

f_tmb <- function(par){
  getAll(dat, par)
  
  N_stk = max(srdat$Stocknumber + 1) # TMB starts at zero, R doesn't
  N_obs = nrow(srdat)
  N_pred = length(pred_lnWA) # (srdat$Stocknumber) # How long is this supposed to be?
  stk = srdat$Stocknumber + 1 # TMB starts at zero, R doesn't
  # S = srdat$Sp/scale_TMB
  # yr = srdat$yr_num # Not used in TMB code

  sigma_delta = exp(logDeltaSigma)
  sigma_nu = exp(logNuSigma)
  
  sigmaA = exp(logSigmaA)
  sigma = exp(logSigma)
  
  SMSY = rep(0, length(unique(srdat$Name)))
  pred_lnSMSY = rep(0, length(unique(srdat$Name)))
  pred_lnSREP = rep(0, length(unique(srdat$Name)))
  logRS_pred = rep(0, N_obs)
  
  ## Initialize joint negative log likelihood
  nll <- 0
  
  # PENALTY TERMS???
  ## Add priors for hyperpars:
    ## MuA prior for stream type
  # logMuA_stream %~% dnorm(logMuA_stream_mean, sd = logMuA_stream_sig)
  nll <- nll - sum(dnorm(logMuA_stream, mean = logMuA_stream_mean, sd = logMuA_stream_sig, log=TRUE)) # nll V.
    ## MuA prior for ocean type
  # logMuA_ocean %~% dnorm(logMuA_ocean_mean, sd = logMuA_ocean_sig)
  nll <- nll - sum(dnorm(logMuA_ocean, mean = logMuA_ocean_mean, sd = logMuA_ocean_sig, log=TRUE)) # nll V.
  
  ## sigmaA prior
  # Tau_sigmaA <- 1/sigmaA^2
  # Tau_sigmaA %~% dgamma(shape = Tau_dist, scale = 1/Tau_dist) 
  nll <- nll - sum(dgamma(1/sigmaA^2, shape = Tau_dist, scale = 1/Tau_dist, log=TRUE)) # nll V.
  # 1/sigmaA^2 %~% dgamma(shape = Tau_dist, scale = 1/Tau_dist)

  ## Add hierarchical structure to A:
  # Tau_sigma <- 1/sigma^2
  for (i in 1:N_stk){
    # logA[i] %~% dnorm(logMuA_stream + logMuA_ocean * lifehist[i], Tau_sigmaA)
    nll <- nll - sum(dnorm(logA[i], logMuA_stream + logMuA_ocean * lifehist[i], sd=sigmaA, log=TRUE)) # nll V.
    # logA[i] %~% dnorm(logMuA_stream + logMuA_ocean * lifehist[i], sigmaA)
    # Tau_sigma[i] %~% dgamma(shape = Tau_dist, scale = 1/Tau_dist)
    nll <- nll - sum(dgamma(1/sigma[i]^2, shape = Tau_dist, scale = 1/Tau_dist, log=TRUE)) # nll V.
    # 1/sigma[i]^2 %~% dgamma(shape = Tau_dist, scale = 1/Tau_dist)
  }
  
  ## Standard Ricker model: 
  for (i in 1:N_obs){
    logRS_pred[i] <- logA[stk[i]] - exp(logB[stk[i]]) * S[i] - sigma[stk[i]]^2/2 # BIAS CORRECTION
    # logRS_pred[i] <- logA[stk[i]] - exp(logB[stk[i]]) * S[i] - sigma[stk[i]]^2/2 # BIAS CORRECTION
      # TK: Should above term be Tau_sigma or sigma?
      # TK: Tau_sigma seems to preform worse
    # logRS_pred <- logA[stk[i]] - exp(logB[stk[i]]) * S[i] # NO BIAS CORRECTION
    # logRS[i] %~% dnorm(logRS_pred[i], Tau_sigma[stk[i]])
    
    # logRS_pred[i] %~% dnorm(logRS[i], Tau_sigma[stk[i]]) # attempting to reverse the order for RTMB
    nll <- nll - sum(dnorm(logRS_pred[i], logRS[i], sd=sigma[stk[i]], log=TRUE)) # nll V.
    # logRS_pred[i] %~% dnorm(logRS[i], sigma[stk[i]]) # attempting to reverse the order for RTMB
  }
  
  ## Calculate SMSY and SREP
  for(i in 1:N_stk){
    SMSY[i] =  (1 - LambertW0(exp(1-logA[i]))) / exp(logB[i]) # Using Paul's new function
  }
  SREP = logA / exp(logB)
  
  ## Inverse gamma prior on sigma_delta and sigma_nu
  # Tau_delta <- 1/sigma_delta^2
  # Tau_delta %~% dgamma(shape = Tau_D_dist, scale = 1/Tau_D_dist)
  nll <- nll - sum(dgamma(1/sigma_delta^2, shape = Tau_D_dist, scale = 1/Tau_D_dist, log=TRUE)) # nll V.
  # 1/sigma_delta^2 %~% dgamma(shape = Tau_D_dist, scale = 1/Tau_D_dist)# test gamma - what do you need to use rate or scale --> only one, one is shape
  
  # Tau_nu <- 1/sigma_nu^2
  # Tau_nu %~% dgamma(shape = Tau_D_dist, scale = 1/Tau_D_dist)
  nll <- nll - sum(dgamma(1/sigma_nu^2, shape = Tau_D_dist, scale = 1/Tau_D_dist, log=TRUE))# nll V.
  # 1/sigma_nu^2 %~% dgamma(shape = Tau_D_dist, scale = 1/Tau_D_dist)
  
  ## Watershed Model
  for (i in 1:N_stk){
    pred_lnSMSY[i] <- logDelta1 + logDelta1_ocean * lifehist[i] + (exp(logDelta2) + Delta2_ocean * lifehist[i]) * log(WAbase$WA[i])
    pred_lnSREP[i] <- logNu1 + logNu1_ocean * lifehist[i] + (exp(logNu2) + Nu2_ocean * lifehist[i]) * log(WAbase$WA[i])
    
    # lnSMSY <- log(SMSY[i] * scale[i])
    # lnSREP <- log(SREP[i] * scale[i])
    # lnSMSY %~% dnorm(pred_lnSMSY, sd = Tau_delta)
    # lnSREP %~% dnorm(pred_lnSREP, sd = Tau_nu)
    # pred_lnSMSY[i] %~% dnorm(log(SMSY[i]*scale[i]), Tau_delta) # attempting to reverse the order for RTMB
    # pred_lnSREP[i] %~% dnorm(log(SREP[i]*scale[i]), Tau_nu) # attempting to reverse the order for RTMB
    
    nll <- nll - sum(dnorm(pred_lnSMSY[i], log(SMSY[i]*scale[i]), sd = sigma_delta, log=TRUE))
    nll <- nll - sum(dnorm(pred_lnSREP[i], log(SREP[i]*scale[i]), sd = sigma_nu, log=TRUE))
    
    # pred_lnSMSY[i] %~% dnorm(log(SMSY[i]*scale[i]), sigma_delta) # attempting to reverse the order for RTMB
    # pred_lnSREP[i] %~% dnorm(log(SREP[i]*scale[i]), sigma_nu) # attempting to reverse the order for RTMB
  }
  
  
  # TARGET PREDICTIONS
  ## Get predicted values for plotting WA regresssion with CIs
    # Create their vectors first
  pred_lnSMSY_stream_CI = rep(0, N_pred)
  pred_lnSMSY_ocean_CI = rep(0, N_pred)
  pred_lnSREP_stream_CI = rep(0, N_pred)
  pred_lnSREP_ocean_CI = rep(0, N_pred)
  
  for (i in 1:N_pred){
    pred_lnSMSY_stream_CI[i] = logDelta1 + exp(logDelta2) * pred_lnWA[i]
    pred_lnSMSY_ocean_CI[i] = logDelta1 + logDelta1_ocean + (exp(logDelta2) + Delta2_ocean) * pred_lnWA[i]
    pred_lnSREP_stream_CI[i] = logNu1 + exp(logNu2) * pred_lnWA[i]
    pred_lnSREP_ocean_CI[i] = logNu1 + logNu1_ocean + (exp(logNu2) + Nu2_ocean) * pred_lnWA[i]
  }
  
  ## Get predicted values for stream-type target stocks with CIs
  N_target_stream = length(target_lnWA_stream) # length of 0 is still 1
  target_lnSMSY_stream = rep(0, N_target_stream)
  target_lnSREP_stream = rep(0, N_target_stream)
  if (length(dat$target_lnWA_stream) > 1){
    for (i in 0:N_target_stream){ # What is N_target_stream
      target_lnSMSY_stream[i] = logDelta1 + exp(logDelta2) * target_lnWA_stream[i]
      target_lnSREP_stream[i] = logNu1 + exp(logNu2) * target_lnWA_stream[i]
    }
    }
  
  ## Get predicted values for ocean-type target stocks with CIs
  N_target_ocean = length(target_lnWA_ocean) 
  target_lnSMSY_ocean = rep(0, N_target_ocean)
  target_lnSREP_ocean = rep(0, N_target_ocean)
  if (length(dat$target_lnWA_ocean) > 1){
    for (i in 0:N_target_ocean){
      target_lnSMSY_ocean[i] = logDelta1 + logDelta1_ocean + (exp(logDelta2) + Delta2_ocean) * target_lnWA_ocean[i]
      target_lnSREP_ocean[i] = logNu1 + logNu1_ocean + (exp(logNu2) + Nu2_ocean) * target_lnWA_ocean[i]
    }
    }

  # OUPUTS
  lnSMSY = log(SMSY*scale)
  lnSREP = log(SREP*scale)

  ADREPORT(SMSY)
  ADREPORT(SREP)
  ADREPORT(logRS_pred)
  ADREPORT(pred_lnSMSY)
  ADREPORT(pred_lnSREP)
  ADREPORT(lnSMSY)
  ADREPORT(lnSREP)
  
  ADREPORT(pred_lnSMSY_stream_CI)
  ADREPORT(pred_lnSMSY_ocean_CI)
  ADREPORT(pred_lnSREP_stream_CI)
  ADREPORT(pred_lnSREP_ocean_CI)
  
  ADREPORT(target_lnSMSY_ocean)
  ADREPORT(target_lnSREP_ocean)
  ADREPORT(target_lnSMSY_stream)
  ADREPORT(target_lnSREP_stream)
  
  # REPORT(nLL)
  
  ## Return
  nll
  
}


## MakeADFun ####

# obj_test <- RTMB::MakeADFun(f_tmb, par, silent = TRUE) # runs
obj <- RTMB::MakeADFun(f_tmb, par, random = c("logA"), silent = TRUE) # create the rtmb object

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 5)) # optimization

# OLD TMB model optimization code
# upper <- unlist(obj$par)
# upper[1:length(upper)]<- Inf
# 
# lower <- unlist(obj$par)
# lower[1:length(lower)]<- -Inf
# 
# opt <- nlminb(obj$par,
#               obj$fn,
#               obj$gr,
#               control = list(eval.max = 1e5, iter.max = 1e5),
#               lower=lower,
#               upper=upper)

# Summarize
sdr <- sdreport(obj)
# sdr_full <- summary(sdreport(obj)) # Includes all ADREPORTS

# Example for how to get ADREPORTS out
sdr_est <- as.list(sdr, "Est", report=TRUE) ## ADREPORT estimates
sdr_se <- as.list(sdr, "Std", report=TRUE) ## ADREPORT standard error

# Best to input a list into the the above example to get out all parameters quickly


## MCMC run through tmbstan ####
# Create a par function
initf1 <- function(){
  list(logA = (srdat %>% group_by (Stocknumber) %>% 
                 summarise(yi = lm(log(Rec / Sp) ~ Sp)$coef[1]))$yi, # random effect
       logB = log ( 1/ ( (1/B$m)/dat$scale )), # fixed effect
       logSigma = rep(-2, length(unique(srdat$Name))),
       logSigmaA = -2,
       
       logMuA_stream = 1.5,
       logMuA_ocean = 0,
       
       logDelta1 = 3,
       logDelta1_ocean = 0,
       logDelta2 = log(0.72),
       Delta2_ocean = 0,
       
       logNu1 = 3,
       logNu1_ocean = 0,
       logNu2 = log(0.72),
       Nu2_ocean = 0,
       
       logNuSigma = -0.412,
       logDeltaSigma = -0.412
  )
}

fitcores <- tmbstan(obj, iter=2000, warmup=200, init=initf1,
                    chains=4, open_progress=FALSE, silent=TRUE)

## Traceplots ####



# Everything else ####
