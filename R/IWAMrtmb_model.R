# IWAM Model

# Made in order to create one unified function with a matching set of
# SREP and SMAX models written in RTMB.
# This will include the following features to be developed:
  # 1. Switch between SREP and SMAX
  # 2. Swtich between Frequentist and Bayesian estimation methods
    # a. Switch for priors on/off
      # e.g. an add prior function instead of an on/off switch within the rtmb(f)

# Libaries ####
library(RTMB)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(progress) # progress bar for iterative loops
library(tidybayes)
library(tmbstan)


# Wrapper Function ####
compiler::enableJIT(0) # Run first without just to see if bug is fixed yet

# Test pieces
WAin = c("DataIn/WCVIStocks.csv")
mod = c("srep")
bias.cor = FALSE
priors = FALSE

# IWAM FUNCTION ####

IWAM_r <- function(WAin = c("DataIn/WCVIStocks.csv"),
                   mod = c("srep"), # or smax - changes what model type is run
                   priors = FALSE, # turn on/off bayesian priors via fpriors 
                   bias.cor = FALSE, # turn on/off usage of a "bias correction" term
                   run.sim = TRUE, # run simulations/bootstraps
                    # TK: Does simulate() provide the same results as IWAM's old bootstrap function?
                    # TK: How do they differ?
                    # TK: Maybe at first just use one of them
                   nsim = 10,
                   sim.seed = 1 # seed for bootstrapping
                   # plot = FALSE
                   # remove.EnhStocks = FALSE # Old data manipulation section for SMAX version
                    # TK: To remove for now
                   ){

  # Progress bar ####
  pb <- progress_bar$new(total = nsim) # nsim can be subbed for however many simulations/bootstraps are going to be run

  # LAMBERT W INTERNAL FUNCTION ####
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
  
  ## Derivatives of LamW ####
  dLambertW0_internal <- function(x, y, dy) {
    dy / (exp(y) * (1. + y))
  }
  
  LambertW0 <- RTMB:::ADjoint(LambertW0_internal, dLambertW0_internal)
  
  
  

  # DATA ENTRY ####
  srdatwna <- read.csv(here::here("DataIn/SRinputfile.csv"))
  WAbase <- read.csv(here::here("DataIn/WatershedArea.csv"))
  # WAin <- read.csv(here::here("DataIn/WCVIStocks.csv")) # For testing
  # WAin <- read.csv(here::here(c("DataIn/Ordered_backcalculated_noagg.csv"))) # For testing
  WAin <- read.csv(here::here(WAin))
  
  
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
  
  # TK: IN SMAX VERSION - WOULD PREVIOUSLY USE SCALING
    # See IWAMsmax_RTMB_model.R
  
  WAbase <- WAbase %>% 
    full_join(names, by="Name") %>% 
    arrange(Stocknumber) %>%
    mutate(logWA = log(WA)) 
  
  # Shift log WA for the mean - base - makes estimation easier
  mean_logWA <- mean(WAbase$logWA)
  WAbase$logWAshifted <- WAbase$logWA - mean_logWA
  
  WAin$logWA <- log(WAin$WA)
  WAin$logWAshifted_t <- WAin$logWA - mean_logWA
  
  lifehist <- srdat %>% dplyr::select(Stocknumber, Name, Stream) %>% 
    group_by(Stocknumber) %>% 
    summarize(lh=max(Stream))
  
  # DAT AND PAR ####
  
  ## * SREP ####
  # Fully covered by dat in SMAX section
  # dat <- list(srdat = srdat,
  #             WAbase = WAbase,
  #             WAin = WAin,
  #             logRS = log(srdat$Rec) - log(srdat$Sp))
  
  # External vectors
  # N_Stk <- max(srdat$Stocknumber + 1)
  
  # Parameters
  # TK: Add in comments to say what the pars are - similar to ModelBook comparisons so it is easier to read first time
  # par <- list(b0 = c(10, 10), # Initial values for WA regression intercepts
  #             # b0 = c(9, 9)
  #             bWA = c(0, 0), # Inital values for WA regression slopes
  #             # bWA = c(0.83, 1)
  #             # logAlpha = numeric(N_Stk), # comment on or off depending if using "nll" or not
  #             logAlpha_re = numeric(nrow(dat$WAbase)), 
  #             logAlpha_re_pred = numeric(nrow(dat$WAin)), 
  #             logE0 = numeric(max(srdat$Stocknumber + 1)),
  #             logE0_ = numeric(nrow(dat$WAin)),
  #             tauobs = 0.01 + numeric(max(srdat$Stocknumber + 1)), # Why can't this be zero? This doesn't run as just a string of zeros.
  #             logAlpha0 = 1.5,
  #             logESD = 1,
  #             logAlphaSD = 10
  # )
  
  ## * SMAX ####
  dat <- list(srdat = srdat,
              WAbase = WAbase,
              WAin = WAin,
              logRS = log(srdat$Rec) - log(srdat$Sp), 
              
              # Aim is to eliminate scale usage
              logRS_scaled = log((srdat$Rec/srdat$scale)/(srdat$Sp/srdat$scale)),
              lifehist = lifehist$lh, # On/off switch in regression
              
              # Scales - may be able to remove all
              scale = srdat_scale, # DNE
              scale_TMB = srdat$scale,
              S_scaled = srdat$Sp/scale_TMB,
              # S = srdat$Sp,
              
              # Constants
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
  
  # Create pars
  if (mod == "srep"){
    par <- list( #### SREP
      b0 = c(10, 10), # Initial values for WA regression intercepts
      # b0 = c(9, 9)
      bWA = c(0, 0), # Inital values for WA regression slopes
      # bWA = c(0.83, 1)
      # logAlpha = numeric(N_Stk), # comment on or off depending if using "nll" or not
      logAlpha_re = numeric(nrow(dat$WAbase)), 
      logAlpha_re_pred = numeric(nrow(dat$WAin)), 
      logE0 = numeric(max(srdat$Stocknumber + 1)),
      logE0_ = numeric(nrow(dat$WAin)),
      tauobs = 0.01 + numeric(max(srdat$Stocknumber + 1)), # Why can't this be zero? This doesn't run as just a string of zeros.
      logAlpha0 = 1.5,
      logESD = 1,
      logAlphaSD = 10)
  } else {
    B <- srdat %>% group_by(Stocknumber) %>% 
      summarise(m = -lm(log(Rec/Sp) ~ Sp)$coef[2])
    par <- list( #### SMAX
      logA = (srdat %>% group_by (Stocknumber) %>% 
                summarise(yi = lm(log(Rec / Sp) ~ Sp)$coef[1]))$yi, # random effect
      logB = log ( 1/ ( (1/B$m)/dat$scale )), # fixed effect
      # logB = log(1/(1/B$m)),
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
  
  # PRIOR FUNCTION ####
  
  fpriors <- function(mod = mod # srep (Liermann priors), and smax (Holt priors)
                        # TK: mod should hopefully supplied by the over-riding IWAM_r function
                      ){
    # Create a nll function for priors to be used in fiwam depending on call
    getAll(dat, par) # Is this required? If it is - should par be going into the function call?
    
    # Constants
    
    nll <- 0
    
    # Priors
    if (mod == "srep"){ # SREP: Liermann
      
    } else { # SMAX: Holt
      
    }
    
    return(nll)
  }
  
  # RTMB FUNCTION ####
  
  fiwam <- function(par){ # THIS IS WRITTEN FIRST AS MLE
    getAll(dat, par)
    
    # Constants: SREP
    N_Stk = max(srdat$Stocknumber + 1)
    stk = srdat$Stocknumber + 1
    N_Obs = nrow(srdat)
    S = srdat$Sp
    type = lifehist + 1 # Used to be lifehist$lh
    E <- numeric(N_Stk)
    logAlpha <- numeric(N_Stk)
    
    # N_Pred = nrow(WAin) # PREDICTION
    # type_ = WAin$lh + 1 # PREDICTION
    # E_tar <- numeric(N_Pred) # PREDICTION
    # log_E_tar <- numeric(N_Pred) # PREDICTION
    # logAlpha_tar <- numeric(N_Pred) # PREDICTION
    
    # Constants: SMAX
    if (mod == "smax"){
      sigma_delta = exp(logDeltaSigma)
      sigma_nu = exp(logNuSigma)
      sigmaA = exp(logSigmaA)
      sigma = exp(logSigma)
      SMSY = rep(0, length(unique(srdat$Name)))
      pred_lnSMSY = rep(0, length(unique(srdat$Name)))
      pred_lnSREP = rep(0, length(unique(srdat$Name)))
    }

    logRS_pred = rep(0, N_Obs)
    
    
    # nll init
    nll <- 0
    
    ### * SREP ####
    if (mod == "srep"){
      ## HYPER-PRIORS ####
        # None for SREP model
      
      ## WATERSHED MODEL ####
      for ( i in 1:N_Stk){
        # nll <- nll - sum(dnorm(logAlpha0, 0.6, sd = 0.45, log = TRUE)) # Test
        # nll <- nll - sum(dunif(logAlphaSD, 0, 100, log = TRUE)) # Test
        # nll <- nll - sum(dunif(logESD, 0, 100, log = TRUE)) # Test
        # nll <- nll - sum(dnorm(logAlpha[i], logAlpha0, sd = logAlphaSD)) # Test
        # nll <- nll - sum(dnorm(logE0[i], 0, sd = logESD)) # Test
        
        nll <- nll - sum(dnorm(logAlpha_re[i], 0, sd = logAlphaSD, log = TRUE)) # Random effect on logAlpha
        logAlpha[i] <- logAlpha0 + logAlpha_re[i] # Moving the scale
        
        nll <- nll - sum(dnorm(logE0[i], 0, sd = logESD, log = TRUE)) # Random effect
        log_E <- b0[type[i]] + bWA[type[i]]*WAbase$logWAshifted[i] + logE0[i] # Stock level regression
        E[i] <- exp(log_E)
        
        nll <- nll - sum(dgamma(tauobs[i], shape = 0.001, scale = 0.001)) 
        # Should be removed as a prior - given an initial value instead?
      }
      
      ## RICKER MODEL ####
      for (i in 1:N_Obs){
        if (bias.cor == FALSE){ # Will this work if it comes from outside this function?
          logRS_pred[i] <- logAlpha[stk[i]]*(1 - S[i]/E[stk[i]])
        } else {
          logRS_pred[i] <- logAlpha[stk[i]]*(1 - S[i]/E[stk[i]]) - 1/tauobs[stk[i]]/2 # BIAS CORRECTION
        }

        nll <- nll - sum(dnorm(logRS[i], logRS_pred[i], sd = sqrt(1/tauobs[stk[i]]), log = TRUE))
      }
      
      ## PREDICTION ####
        # TK: Add later
      
      ## REPORTING ####
      ADREPORT(E)
      ADREPORT(logAlpha)
      ADREPORT(logRS_pred) # Why the hell not
      
      nll
    }
    
    ### * SMAX ####
    if (mod == "smax"){
      ## HYPER-PRIORS ####
        # TK: Now supplied by the prior function
      if (priors == TRUE){
        nll <- nll - fpriors(mod = mod) # input the function sum of nll from priors function
          # TK: This should have a really simple input/default setting
      } else {
        ## Add MLE priors/penalities for hyperpars:
          ## MuA prior for stream type
        nll <- nll - sum(dnorm(logMuA_stream, mean = logMuA_stream_mean, sd = logMuA_stream_sig, log=TRUE))
          ## MuA prior for ocean type
        nll <- nll - sum(dnorm(logMuA_ocean, mean = logMuA_ocean_mean, sd = logMuA_ocean_sig, log=TRUE))
        
          ## sigmaA prior
        nll <- nll - sum(dgamma(1/sigmaA^2, shape = Tau_dist, scale = 1/Tau_dist, log=TRUE))

        ## Add hierarchical structure to A:
        for (i in 1:N_Stk){
          nll <- nll - sum(dnorm(logA[i], logMuA_stream + logMuA_ocean * lifehist[i], sd=sigmaA, log=TRUE))
          nll <- nll - sum(dgamma(1/sigma[i]^2, shape = Tau_dist, scale = 1/Tau_dist, log=TRUE))
        }
        
        ## Inverse gamma prior on sigma_delta and sigma_nu
        nll <- nll - sum(dgamma(1/sigma_delta^2, shape = Tau_D_dist, scale = 1/Tau_D_dist, log=TRUE))
        
        nll <- nll - sum(dgamma(1/sigma_nu^2, shape = Tau_D_dist, scale = 1/Tau_D_dist, log=TRUE))
      }

      ## WATERSHED MODEL ####
      for (i in 1:N_Stk){
        pred_lnSMSY[i] <- logDelta1 + logDelta1_ocean * lifehist[i] + (exp(logDelta2) + Delta2_ocean * lifehist[i]) * log(WAbase$WA[i])
        pred_lnSREP[i] <- logNu1 + logNu1_ocean * lifehist[i] + (exp(logNu2) + Nu2_ocean * lifehist[i]) * log(WAbase$WA[i])
        
          # Errors - lacking scaling
        nll <- nll - sum(dnorm(pred_lnSMSY[i], log(SMSY[i]*scale[i]), sd = sigma_delta, log=TRUE))
        nll <- nll - sum(dnorm(pred_lnSREP[i], log(SREP[i]*scale[i]), sd = sigma_nu, log=TRUE))
      }
      
      ## RICKER MODEL ####
        # TK: Is this actually faster - or should the TRUE FALSE IF statement be outside the actual loop?
      for (i in 1:N_Obs){
        if (bias.cor == FALSE){
          logRS_pred[i] <- logA[stk[i]] - exp(logB[stk[i]]) * S_scaled[i]
        } else {
          logRS_pred[i] <- logA[stk[i]] - exp(logB[stk[i]]) * S_scaled[i] - sigma[stk[i]]^2/2 # BIAS CORRECTION
        }
        nll <- nll - sum(dnorm(logRS_scale[i], logRS_pred[i], sd=sigma[stk[i]], log=TRUE))
      }
      
      ## Calculate SMSY and SREP
      for(i in 1:N_Stk){
        SMSY[i] =  (1 - LambertW0(exp(1-logA[i]))) / exp(logB[i]) # Using Paul's new function
      }
      SREP = logA / exp(logB)
      
      ## PREDICTION ####
        # TK: Add later
      
      ADREPORT(SMSY)
      ADREPORT(SREP)
      ADREPORT(logRS_pred)
      ADREPORT(pred_lnSMSY)
      ADREPORT(pred_lnSREP)
      ADREPORT(lnSMSY)
      ADREPORT(lnSREP)
    
      nll
    }
    
    return(nll)
  }
  
  # RTMB OBJ and OPT ####
  
  if (mod == "srep") {# determine what the random variables are depending on the model being run
    # random_vars <- c("logAlpha_re", "logAlpha_re_pred", "logE0", "logE0_")
    random_vars <- c("logAlpha_re", "logE0")
  } else {
      random_vars <- c("logA")
    }
  
  obj <- RTMB::MakeADFun(fiwam, par, random = random_vars, silent = TRUE) # create the rtmb object
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 0)) # optimization
  
  # SUMMARIES ####
  sdr <- sdreport(obj)
  sdr_full <- summary(RTMB::sdreport(obj))
  sdr_est <- as.list(sdr, "Est", report=TRUE) ## ADREPORT estimates
  sdr_se <- as.list(sdr, "Std", report=TRUE) ## ADREPORT standard error
  
  # OUTPUT MANIPULATION ####
  
  # SIMULATE/bootstrap ####
  
  # RETURN ####
  return(list(opt = opt,
              obj = obj,
              sdr = sdr,
              sdr_full = sdr_full,
              sdr_est = sdr_est,
              sdr_se = sdr_se)
         )
  
}

# Test function
# IWAM_r <- function(WAin = c("DataIn/WCVIStocks.csv"),
#                    mod = c("srep"), # or smax - changes what model type is run
#                    priors = FALSE, # turn on/off bayesian priors via fpriors
#                    bias.cor = FALSE, # turn on/off usage of a "bias correction" term
#                    run.sim = TRUE, # run simulations/bootstraps
#                    # TK: Does simulate() provide the same results as IWAM's old bootstrap function?
#                    # TK: How do they differ?
#                    # TK: Maybe at first just use one of them
#                    nsim = 10,
#                    sim.seed = 1 # seed for bootstrapping
#                    # plot = FALSE
#                    # remove.EnhStocks = FALSE # Old data manipulation section for SMAX version
#                    # TK: To remove for now
# )

test <- IWAM_r()
test <- IWAM_r(mod = c("smax"))
