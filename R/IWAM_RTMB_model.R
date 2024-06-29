# IWAM 2 ####

# This is the rtmb version of the IWAM model.
# The primary reason of this re-compilation is to check if the rtmb model can be processed
# using the same code as is in "IWAM_model.R".

# Libaries ####

library(RTMB)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidybayes)
library(tmbstan)

source (here::here("R/helperFunctions.R"))
source(here::here("R/PlotFunctions.R"))
source(here::here("R/Get_LRP_bs.R")) 
  # Also attaches TMB's package 
  # Something about this is bricking the RTMB function code

## Turn off byte compiler ####
compiler::enableJIT(0)

## Wrapper function Beginning ####

# remove.EnhStocks = FALSE # Just to make the bottom copied code run
# WAin <- read.csv(here::here("DataIn/WCVIStocks.csv"))
# WAin <- c("DataIn/WCVIStocks.csv")

IWAM_rtmb <- function(WAin = c("DataIn/WCVIStocks.csv"),
                      remove.EnhStocks = FALSE,
                      run.bootstraps = TRUE, # to turn on or off the bootstrap function added at the end
                      bs_seed = 1, # seed for bootstrapping
                      bs_nBS = 10, # trials for bootstrapping
                      plot = FALSE, # whether or not to create plots stored in DataOut/
)
{
  
  ## Imported LambertW0 ####
  
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
  
  
  
  ## Data ####
  srdatwna <- read.csv(here::here("DataIn/SRinputfile.csv"))
  WAbase <- read.csv(here::here("DataIn/WatershedArea.csv"))
  WAin <- read.csv(here::here(WAin))
  
  ## Data setup ####
  srdatwna <- srdatwna %>% 
    filter(!Name %in% c("Hoko","Hoh")) 
  
  # Remove years with NAs and re-numerate.
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
  
  ## RTMB function ####

  f_tmb <- function(par){
    getAll(dat, par)
  
    N_stk = max(srdat$Stocknumber + 1)
    N_obs = nrow(srdat)
    N_pred = length(pred_lnWA)
    stk = srdat$Stocknumber + 1
  
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
  
    ## Add priors for hyperpars:
    ## MuA prior for stream type
    nll <- nll - sum(dnorm(logMuA_stream, mean = logMuA_stream_mean, sd = logMuA_stream_sig, log=TRUE)) # nll V.
    ## MuA prior for ocean type
    nll <- nll - sum(dnorm(logMuA_ocean, mean = logMuA_ocean_mean, sd = logMuA_ocean_sig, log=TRUE)) # nll V.
  
    ## sigmaA prior
    nll <- nll - sum(dgamma(1/sigmaA^2, shape = Tau_dist, scale = 1/Tau_dist, log=TRUE)) # nll V.
      # Should this be in log space?
  
    ## Add hierarchical structure to A:
    for (i in 1:N_stk){
      nll <- nll - sum(dnorm(logA[i], logMuA_stream + logMuA_ocean * lifehist[i], sd=sigmaA, log=TRUE)) # nll V.
      nll <- nll - sum(dgamma(1/sigma[i]^2, shape = Tau_dist, scale = 1/Tau_dist, log=TRUE)) # nll V.
        # Should this be in log space?
    }
  
    ## Standard Ricker model:
    for (i in 1:N_obs){
      logRS_pred[i] <- logA[stk[i]] - exp(logB[stk[i]]) * S[i] - sigma[stk[i]]^2/2 # BIAS CORRECTION
  
      nll <- nll - sum(dnorm(logRS_pred[i], logRS[i], sd=sigma[stk[i]], log=TRUE)) # nll V.
    }
  
    ## Calculate SMSY and SREP
    for(i in 1:N_stk){
      SMSY[i] =  (1 - LambertW0(exp(1-logA[i]))) / exp(logB[i]) # Using Paul's new function
    }
    SREP = logA / exp(logB)
  
    ## Inverse gamma prior on sigma_delta and sigma_nu
    nll <- nll - sum(dgamma(1/sigma_delta^2, shape = Tau_D_dist, scale = 1/Tau_D_dist, log=TRUE)) # nll V.
  
    nll <- nll - sum(dgamma(1/sigma_nu^2, shape = Tau_D_dist, scale = 1/Tau_D_dist, log=TRUE))# nll V.
  
    ## Watershed Model
    for (i in 1:N_stk){
      pred_lnSMSY[i] <- logDelta1 + logDelta1_ocean * lifehist[i] + (exp(logDelta2) + Delta2_ocean * lifehist[i]) * log(WAbase$WA[i])
      pred_lnSREP[i] <- logNu1 + logNu1_ocean * lifehist[i] + (exp(logNu2) + Nu2_ocean * lifehist[i]) * log(WAbase$WA[i])
  
      nll <- nll - sum(dnorm(pred_lnSMSY[i], log(SMSY[i]*scale[i]), sd = sigma_delta, log=TRUE))
      nll <- nll - sum(dnorm(pred_lnSREP[i], log(SREP[i]*scale[i]), sd = sigma_nu, log=TRUE))
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
  
    ## Return
    return(nll)
  }
  
  ## MakeADFun ####
    # NOTE: If you run f_tmb - before you run obj - there will be math errors
    # TK: I have no idea why
  obj <- RTMB::MakeADFun(f_tmb, par, random = c("logA"), silent = TRUE) # create the rtmb object
  
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(trace = 5)) # optimization
  
  # sdr <- summary(sdreport(obj)) # This is missing something compared to the TMB version
  
  # Must be renamed before being passed into function - otherwise summary doesn't work
  # obj_rtmb <- obj
  
  # Return function
  # return(list(opt = opt,
  #             obj_rtmb = obj_rtmb
  #             )) # This can be added to whenever
  
  ## Compile model outputs ####
  # Create Table of outputs
  all_pars <- data.frame(summary(RTMB::sdreport(obj)))
  all_pars$Param <- row.names(all_pars)
  all_pars$Param <- sapply(all_pars$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
  
  pars <- data.frame()
  pars <- all_pars %>% filter (Param %in% c("logA", 
                                            "logB", 
                                            "logSigma",  
                                            "SMSY", 
                                            "SREP",
                                            "lnSMSY",
                                            "lnSREP"))
  
  stnum <- unique(srdat$Stocknumber) 
  pars$Stocknumber <- rep(stnum) 
  pars <- left_join(pars, unique(srdat[, c("Stocknumber", "Name")]))
  
  logDeltaSigma <- all_pars %>% filter (Param %in% c("logDeltaSigma")) 
  DeltaSigmaUCL <- exp(logDeltaSigma$Estimate + logDeltaSigma$Std..Error*1.96)
  DeltaSigmaLCL <- exp(logDeltaSigma$Estimate - logDeltaSigma$Std..Error*1.96) 
  DeltaSigma <- exp(logDeltaSigma$Estimate)
  
  # Combine again and rename
  # pars <- All_Est <- All_Ests_std
  pars$Param <- sapply(pars$Param, function(x) (unlist(strsplit(x, "[_]"))[[1]]))
  pars <- pars %>% left_join(lifehist, by="Stocknumber")
  
  
  all_Deltas <- data.frame()
  all_Deltas <- all_pars %>% filter (Param %in% c("logDelta1", 
                                                  "logDelta2",
                                                  "sigma_delta", 
                                                  "Delta2_bounded", 
                                                  "logDelta1_ocean", 
                                                  #"logDelta2ocean", # does not exist
                                                  "Delta2_ocean", 
                                                  "logNu1", 
                                                  "logNu2", 
                                                  "sigma_nu", 
                                                  "logNu1_ocean", 
                                                  "Nu2_ocean"))
  
  
  ## Calculate diagnostics and plot SR curves, etc. ####
  # Get predicted values and calculate r2
    # 
  
  pred_RS <- data.frame()
  pred_RS <- all_pars %>% filter (Param %in% c("logRS_pred"))
  
  # * all_pred previously Preds
  all_pred <- srdat %>% dplyr::select("Stocknumber",
                                      "yr_num", 
                                      "Sp", 
                                      "Rec", 
                                      "scale", 
                                      "Name") %>% 
              add_column(Pred=pred_RS$Estimate)
  
  # Check that length of vector of pred_RS is same as the length of spawner abundances in TMB data element
  # for spawners
  if (length(pred_RS$Estimate) == length(dat$S)) {
    print("Lengths checked passed.")
  } else {
    print("WARNING: The output and inputs are not the same length.")
  }
  
  # mutate the predicted values with scale 
  # RE-SCALED VALUES
  # These Preds_stds are not used for plotting
  all_pred <- all_pred %>% mutate(ObslogRS = log ( (Rec / scale) / (Sp/scale) ) )
  r2 <- all_pred %>% group_by(Stocknumber) %>% summarize(r2=cor(ObslogRS,Pred)^2)
  
  # Get predicted values and their SEs to plot CIs
  # *These are not re-scaled*
  # They are used in the plotting functions and scaled within
  # pred_lnSMSY_S and pred_lnSMSY_O don't occur anywhere ************************************************************
  pred_lnSMSY <- data.frame() 
  pred_lnSMSY <- all_pars %>% filter (Param %in% c("pred_lnSMSY_stream", # not included
                                                   "pred_lnSMSY_ocean", # not included
                                                   "pred_lnSMSY_CI", # not included
                                                   "pred_lnSMSY_stream_CI", 
                                                   "pred_lnSMSY_ocean_CI"))
  # pred_lnSREP_S and pred_lnSREP_O don't occur elsewhere ***********************************************************
  pred_lnSREP <- data.frame() 
  pred_lnSREP <- all_pars %>% filter (Param %in% c("pred_lnSREP_stream", # not included
                                                   "pred_lnSREP_ocean", # not included
                                                   "pred_lnSREP_CI", # not included
                                                   "pred_lnSREP_stream_CI", 
                                                   "pred_lnSREP_ocean_CI"))
  
  # Calculate standardized residuals
  # These are RE-SCALED values
  SRes <- all_pred %>% mutate ( Res = ObslogRS- Pred) #%>% mutate (StdRes = Res/??)
  sigma <- pars %>% filter(Param=="logSigma") %>% dplyr::select(Stocknumber, Estimate, Name)
  SRes <- SRes %>% left_join(sigma) %>% rename(logSig = Estimate)
  # Slight problem naming with Estimate - can produce clang errors due to overlap
  # with function name.
  SRes <- SRes %>% mutate (StdRes = Res/exp(logSig))
  
  ## PLOTTING ####
  
  mod <- "IWAM_Liermann" 
  
  #### * Plot SR Curves ----------------------------------------------------------
  if (plot==TRUE){
    png(paste("DataOut/SR_rtmb_", mod, ".png", sep=""), width=7, height=7, units="in", res=500)
    PlotSRCurve(srdat=srdat, pars=pars, r2=r2, removeSkagit = FALSE, mod=mod)
    dev.off()
    png(paste("DataOut/SRLin_rtmb_", mod, ".png", sep=""), width=7, height=7, units="in", res=1000)
    PlotSRLinear(srdat=srdat, pars=pars, r2=r2, removeSkagit = FALSE)
    dev.off()
    png(paste("DataOut/StdResid_rtmb_", mod, ".png", sep=""), width=7, height=7, units="in", res=1000)
    PlotStdResid(SRes)
    dev.off()
    png(paste("DataOut/ACF_rtmb_", mod, ".png", sep=""), width=7, height=7, units="in", res=1000)
    Plotacf(SRes)
    dev.off()
  }
  
  #### * Plot WA Regression ------------------------------------------------------
  if(plot==TRUE){
    png(paste("DataOut/WAregSMSY_rtmb_", mod, "_wBC.png", sep=""), width=7, height=7, units="in", res=500)
    par(mfrow=c(1,1), mar=c(4, 4, 4, 2) + 0.1)
    title_plot <- "Prior Ricker sigma and prior WA regression sigma"
    plotWAregressionSMSY (pars, all_Deltas, srdat, lifehist, WAbase, pred_lnSMSY, 
                          pred_lnWA = dat$pred_lnWA, title1=title_plot, mod)
    dev.off()
    
    png(paste("DataOut/WAregSREP_rtmb_", mod, "_wBC.png", sep=""), width=7, height=7, units="in", res=500)
    par(mfrow=c(1,1), mar=c(4, 4, 4, 2) + 0.1)
    title_plot <- "Prior Ricker sigmas and prior on WA regression sigma"
    plotWAregressionSREP (pars, all_Deltas, srdat, lifehist, WAbase, pred_lnSREP, 
                          pred_lnWA = dat$pred_lnWA, title1=title_plot, mod)
    dev.off()
  }
  
  ## Calculate prediction intervals for SMSY and SREP for additional stocks ####
  
  # Get predicted values to estimate prediction intervals
  # These values are RE-SCALED to raw estimates during outputting in the TMB code
  pred_lnSMSY_pi <- data.frame()
  pred_lnSMSY_pi <- all_pars %>% filter (Param %in% c("pred_lnSMSY", "lnSMSY"))
  pred_lnSREP_pi <- data.frame()
  pred_lnSREP_pi <- all_pars %>% filter (Param %in% c("pred_lnSREP", "lnSREP"))
  
  pred_lnSMSY_pi$Stocknumber <- rep(stnum)
  pred_lnSREP_pi$Stocknumber <- rep(stnum)
  
  # To calculate prediction intervals, first get predicted and observed logSMSY 
  # and logSREP values for synoptic data set
  #   (actually only need observed logSMSY and logSREP values)
  
  #  First need to get the scale for each stock
  scale_pi <- srdat %>% dplyr::select(Stocknumber, scale) %>% distinct()
  pred_lnSMSY_pi <- pred_lnSMSY_pi %>% left_join(unique(srdat[, c("Stocknumber", "Name")])) %>% 
    left_join(scale_pi)
  pred_lnSREP_pi <- pred_lnSREP_pi %>% left_join(unique(srdat[, c("Stocknumber", "Name")])) %>% 
    left_join(scale_pi)
  
  # Then need to separate observed stream vs ocean type
  
  # pred_lSMSY_stream = predicted log SMSY for stream type
  # pred_lSMSY_ocean = predicted log SMSY for ocean type
  # obs_lSMSY_stream = observed log SMSY for stream type
  # obs_lSMSY_ocean = observed log SMSY for ocean type
  # and same for SREP
  
  pred_lSMSY_stream <- pred_lnSMSY_pi %>% filter(Param=="pred_lnSMSY") %>% left_join(lifehist) %>% 
    filter(lh == 0) %>% pull(Estimate) #Predicted lnSMSY from WA regression- stream - PlSMSYs
  pred_lSMSY_ocean <- pred_lnSMSY_pi %>% filter(Param=="pred_lnSMSY") %>% left_join(lifehist) %>% 
    filter(lh == 1) %>% pull(Estimate) #Predicted lnSMSY from WA regression- ocean
  obs_lSMSY_stream <- pred_lnSMSY_pi %>% filter(Param=="lnSMSY") %>% left_join(lifehist) %>% 
    filter( lh== 0) %>% pull(Estimate) # "observed" lnSMSY data output from SR models- stream
  obs_lSMSY_ocean <- pred_lnSMSY_pi %>% filter(Param=="lnSMSY") %>% left_join(lifehist) %>% 
    filter(lh == 1) %>% pull(Estimate) # "observed" lnSMSY data output from SR models- ocean
  
  pred_SREP_stream <- pred_lnSREP_pi %>% filter(Param=="pred_lnSREP") %>% left_join(lifehist) %>% 
    filter(lh == 0) %>% pull(Estimate) #Predicted lnSMSY from WA regression- stream
  pred_SREP_ocean <- pred_lnSREP_pi %>% filter(Param=="pred_lnSREP") %>% left_join(lifehist) %>% 
    filter(lh == 1) %>% pull(Estimate) #Predicted lnSMSY from WA regression- ocean
  obs_SREP_stream <- pred_lnSREP_pi %>% filter(Param=="lnSREP") %>% left_join(lifehist) %>% 
    filter(lh == 0) %>% pull(Estimate) # "observed" lnSMSY data output from SR models- stream
  obs_SREP_ocean <- pred_lnSREP_pi %>% filter(Param=="lnSREP") %>% left_join(lifehist) %>% 
    filter(lh == 1) %>% pull(Estimate) # "observed" lnSMSY data output from SR models- ocean
  
  
  # Get watershed areas for synoptic data set to calculate PIs for stream and ocean
    # Doesn't work because WAbase$lh is now a factor ****
  wa_stream <- WAbase %>% left_join(lifehist, by = c("Stocknumber")) %>% filter(lh.y == 0) %>% pull(WA)
  wa_ocean <- WAbase %>% left_join(lifehist, by = c("Stocknumber")) %>% filter(lh.y == 1) %>% pull(WA)
  
  
  # TK: Ignoring CU and inlet's for now!
  
  
  stocknames_stream <- c(as.vector(WAin$Stock[WAin$lh == 0]))
  stocknames_ocean <- c(as.vector(WAin$Stock[WAin$lh == 1]))
  
  # Get Predicted SMSY and SREP values for the target stocks and their Prediction Intervals
  # For single life-history events (stream OR ocean targets)
  # targetSMSY <- data.frame() 
  # targetSREP <- data.frame()
  
  # For instances of both life histories (stream AND ocean targets)
  target_SMSY_ocean <- data.frame()
  target_SREP_ocean <- data.frame()
  
  target_SMSY_stream <- data.frame()
  target_SREP_stream <- data.frame()
  
  if (any(WAin$lh == 1)) {
    # Step 1
    target_SMSY_ocean <- all_pars %>% filter (Param %in% c("target_lnSMSY_ocean")) %>% add_column(Stock = stocknames_ocean)
    # with latest dataset -> backcalced:
    # stocknames has 116 values - of these 116 values, 85 are ocean, and 32 are stream
    # need to make sure that this difference is understood
    # consider making a stocknames_ocean and _stream?
    target_SREP_ocean <- all_pars %>% filter (Param %in% c("target_lnSREP_ocean")) %>% add_column(Stock = stocknames_ocean)
    target_SMSY_pull_ocean <- target_SMSY_ocean %>% pull(Estimate)
    target_SREP_pull_ocean <- target_SREP_ocean %>% pull(Estimate)
    # Step 2
    target_SMSY_pi_ocean <- PredInt(x = log(wa_ocean), y = obs_lSMSY_ocean, Predy = target_SMSY_pull_ocean, Newx = dat$target_lnWA_ocean)
    target_SREP_pi_ocean <- PredInt(x = log(wa_ocean), y = obs_SREP_ocean, Predy = target_SREP_pull_ocean, Newx = dat$target_lnWA_ocean)
    # Step 3
    target_SMSY_ocean <- target_SMSY_ocean %>% add_column(LL = exp(target_SMSY_pi_ocean$lwr), UL = exp(target_SMSY_pi_ocean$upr))
    target_SREP_ocean <- target_SREP_ocean %>% add_column(LL = exp(target_SREP_pi_ocean$lwr), UL = exp(target_SREP_pi_ocean$upr))
    # Step 4
    target_SMSY_ocean <- target_SMSY_ocean %>% mutate (Estimate = exp(Estimate)) %>% dplyr::select(-Std..Error, - Param) %>% 
      add_column(Param = "SMSY")
    target_SREP_ocean <- target_SREP_ocean %>% mutate (Estimate = exp(Estimate)) %>% dplyr::select(-Std..Error, - Param) %>% 
      add_column(Param = "SREP")
    # Step 5
    target_estimates_SMSY_ocean <- target_SMSY_ocean %>% mutate(Estimate = round(Estimate, 0), LL = round(LL,0), UL = round(UL,0))
    target_estimates_SREP_ocean <- target_SREP_ocean %>% mutate(Estimate = round(Estimate, 0), LL = round(LL,0), UL = round(UL,0))
    
    data_out_ocean <- target_estimates_SMSY_ocean %>% bind_rows(target_estimates_SREP_ocean)
    
  }
  
  # WARNING
  length_check_params_ocean <- all_pars %>% filter (Param %in% c("target_lnSMSY_ocean"))
  if (length(stocknames_ocean) == length(length_check_params_ocean$Estimate)) { # originally this checked the full stocknames list
    print("Lengths checked passed - Ocean life histories.")
  } else {
    print("WARNING: The output and inputs are not the same length.")
  }
  
  # Does not currently work * no data for stream
  if (any(WAin$lh == 0)){
    # Step 1
    target_SMSY_stream <- all_pars %>% filter (Param %in% c("target_lnSMSY_stream")) %>% add_column(Stock = stocknames_stream)
    target_SREP_stream <- all_pars %>% filter (Param %in% c("target_lnSREP_stream")) %>% add_column(Stock = stocknames_stream)
    target_SMSY_pull_stream <- target_SMSY_stream %>% pull(Estimate)
    target_SREP_pull_stream <- target_SREP_stream %>% pull(Estimate)
    # Step 2
    target_SMSY_pi_stream <- PredInt(x = log(wa_stream), y = obs_lSMSY_stream, Predy = target_SMSY_pull_stream, Newx = dat$target_lnWA_stream)
    target_SREP_pi_stream <- PredInt(x = log(wa_stream), y = obs_SREP_stream, Predy = target_SREP_pull_stream, Newx = dat$target_lnWA_stream)
    # Step 3
    target_SMSY_stream <- target_SMSY_stream %>% add_column(LL = exp(target_SMSY_pi_stream$lwr), UL = exp(target_SMSY_pi_stream$upr))
    target_SREP_stream <- target_SREP_stream %>% add_column(LL = exp(target_SREP_pi_stream$lwr), UL = exp(target_SREP_pi_stream$upr))
    # Step 4
    target_SMSY_stream <- target_SMSY_stream %>% mutate (Estimate = exp(Estimate)) %>% dplyr::select(-Std..Error, - Param) %>%
      add_column(Param = "SMSY")
    target_SREP_stream <- target_SREP_stream %>% mutate (Estimate = exp(Estimate)) %>% dplyr::select(-Std..Error, - Param) %>%
      add_column(Param = "SREP")
    # Step 5
    target_estimates_SMSY_stream <- target_SMSY_stream %>% mutate(Estimate = round(Estimate, 0), LL = round(LL,0), UL = round(UL,0))
    target_estimates_SREP_stream <- target_SREP_stream %>% mutate(Estimate = round(Estimate, 0), LL = round(LL,0), UL = round(UL,0))
    
    data_out_stream <- target_estimates_SMSY_stream %>% bind_rows(target_estimates_SREP_stream)
    
  }
  
  # WARNING
  length_check_params_stream <- all_pars %>% filter (Param %in% c("target_lnSMSY_stream"))
  if (length(stocknames_stream) == length(length_check_params_stream$Estimate)) { # originally this checked the full stocknames list
    print("Lengths checked passed - Stream life histories.")
  } else {
    print("WARNING: The output and inputs are not the same length.")
  }
  
  # Final combination of SMSY and SREP estimates into final df
  # if statement for stream or ocean presence
  if (all(c(0, 1) %in% WAin$lh)) {
    print("Both life histories are present and will be combined.")
    data_out_combined <- data_out_ocean %>% bind_rows(data_out_stream)
  } else {
    print("Only one life history detected.")
  }
  
  ## Data Printing ####
  
  # For NO ENHANCEMENT datasets
  if (remove.EnhStocks) {
    if (all(WAin$lh == 0)) { # stream only
      if(remove.EnhStocks) write.csv(data_out_stream, here::here("DataOut/rtmb_dataout_target_stream_noEnh.csv"))
      if(remove.EnhStocks) datain <- c("DataOut/rtmb_dataout_target_stream_noEnh.csv")
      print("Stream life histories detected and moved forward.")
      print("DataOut: rtmb_dataout_target_stream_noEnh")
      
    }
    if (all(WAin$lh == 1)) { # ocean only
      if(remove.EnhStocks) write.csv(data_out_ocean, here::here("DataOut/rtmb_dataout_target_ocean_noEnh.csv"))
      if(remove.EnhStocks) datain <- c("DataOut/rtmb_dataout_target_ocean_noEnh.csv")
      print("Ocean life histories detected and moved forward.")
      print("DataOut: rtmb_dataout_target_ocean_noEnh")
      
    }
    if (all(c(0, 1) %in% WAin$lh)) { # stream and ocean 
      if(remove.EnhStocks) write.csv(data_out_combined, here::here("DataOut/rtmb_dataout_target_noEnh.csv"))
      if(remove.EnhStocks) datain <- c("DataOut/rtmb_dataout_target_noEnh.csv")
      print("Stream and ocean life histories detected, combined, and moved forward.")
      print("DataOut: rtmb_dataout_target_noEnh")
      
    }
  }
  
  # For ENHANCEMENT datasets
  if (!remove.EnhStocks) {
    if (all(WAin$lh == 0)) { # stream
      if(!remove.EnhStocks) write.csv(data_out_stream, here::here("DataOut/rtmb_dataout_target_stream_wEnh.csv"))
      if(!remove.EnhStocks) datain <- c("DataOut/rtmb_dataout_target_stream_wEnh.csv")
      print("Stream life histories detected and moved forward.")
      print("DataOut: rtmb_dataout_target_stream_wEnh")
      
    } 
    if (all(WAin$lh == 1)) { # ocean
      if(!remove.EnhStocks) write.csv(data_out_ocean, here::here("DataOut/rtmb_dataout_target_ocean_wEnh.csv"))
      if(!remove.EnhStocks) datain <- c("DataOut/rtmb_dataout_target_ocean_wEnh.csv")
      print("Ocean life histories detected and moved forward.")
      print("DataOut: rtmb_dataout_target_ocean_wEnh")
      
    }
    if (all(c(0, 1) %in% WAin$lh)) { # combined
      if(!remove.EnhStocks) write.csv(data_out_combined, here::here("DataOut/rtmb_dataout_target_wEnh.csv"))
      if(!remove.EnhStocks) datain <- c("DataOut/rtmb_dataout_target_wEnh.csv")
      print("Stream and ocean life histories detected, combined, and moved forward.")
      print("DataOut: rtmb_dataout_target_wEnh")
      
    } 
  }
  
  # Data is now printed in the same form as IWAM_model.R
  
  ## Bootstrap function ####
    # Run the below for internal runs ********************************************
  # run.bootstraps <- TRUE
  # bs_seed <- 1
  # bs_nBS <- 10
  # datain <- c("DataOut/rtmb_dataout_target_ocean_wEnh.csv")
  
  if (run.bootstraps == TRUE){
    # set.seed(1) #10#12#13 (work for 1000), for 100, 200, 300, (for 5000trials), 1, 2, 3 (for 20000trials)
    set.seed(bs_seed)
    # nBS <- 10 # number trials for bootstrapping (original 20000), for testing use 10
    nBS <- bs_nBS
    outBench <- list()
    
    for (k in 1:nBS) {
      # datain must match the above Step 6 for writing output target estimates
      out <- Get.LRP.bs(datain = datain, # "DataOut/FUNCTIONTEST_dataout_target_ocean_noEnh.csv"
                        Bern_logistic=FALSE, 
                        prod="LifeStageModel",
                        LOO = NA, 
                        run_logReg=FALSE) 
      outBench[[k]] <- out$bench
    }
    # Get.LRP.bs has WCVIstocks cooked in
    # need to set it so there is an input dataset
    # datain is not enough - or maybe certain loops just need re-objectification
    
    # # Is 200 enough trials? Yes
    # running.mean <- cumsum(LRP.bs$fit) / seq_along(LRP.bs$fit) 
    # plot(running.mean)
    
    # Compile bootstrapped estimates of Sgen, SMSY, and SREP, and identify 5th and 
    # 95th percentiles
    SGEN.bs <- select(as.data.frame(outBench), starts_with("SGEN"))
    
    stockNames <- read.csv(here::here(datain)) %>% 
      #  filter(Stock != "Cypre") %>% 
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
    
    boot <- list(SGEN.boot=SGEN.boot, SMSY.boot=SMSY.boot, 
                 SREP.boot=SREP.boot)
    
    df1 <- data.frame(boot[["SGEN.boot"]], Stock=rownames(boot[["SGEN.boot"]]), RP="SGEN") 
    df1 <- df1 %>% rename(Value=SGEN)
    df2 <- data.frame(boot[["SREP.boot"]], Stock=rownames(boot[["SREP.boot"]]), RP="SREP")
    df2 <- df2 %>% rename(Value=SREP)
    df3 <- data.frame(boot[["SMSY.boot"]], Stock=rownames(boot[["SMSY.boot"]]), RP="SMSY")
    df3 <- df3 %>% rename(Value=SMSY)
    
    dfout <- add_row(df1, df2)
    dfout <- add_row(dfout, df3)
    rownames(dfout) <- NULL
    # now round to 2 signif digits
    dfout <- dfout %>% mutate(Value=signif(Value, 2)) %>% 
      mutate(lwr=signif(lwr,2)) %>% 
      mutate (upr=signif(upr,2))
    
    # Add a function for IWAM_func to rename outputs?
    write.csv(dfout, here::here("DataOut/rtmb_getLRP-BootstrappedRPs.csv"))
  }
  
  # This is where the old table function used to be. RIP.
  
  #### Return function ####
  return(list(opt = opt,
              dataname = datain,
              dfout = dfout,
              iwampars = pars,
              all_Deltas = all_Deltas,
              srdat = srdat,
              lh = lifehist,
              WAbase = WAbase,
              pred_lnSREP = pred_lnSREP,
              pred_lnSMSY = pred_lnSMSY,
              pred_lnSREP_pi = pred_lnSREP_pi,
              pred_lnSMSY_pi = pred_lnSMSY_pi,
              pred_lnWA = dat$pred_lnWA,
              SRes = SRes,
              r2 = r2)) # This can be added to whenever
  
} # End IWAM_rtmb function

# Test run IWAM_rtmb func
test <- IWAM_rtmb() # default run
  # confirmed same objective value - code matches - we are good






