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
library(progress) # progress bar for iterative loops
library(tmbstan)
library(tidybayes) # Plotting usage

source(here::here("R/LambertWs.R"))

# Wrapper Function ####
compiler::enableJIT(0) # Run first without just to see if bug is fixed yet
# seems not to run on the first attempt - answers [3] instead of [0]

# Saves for internal runs without the wrapper function
WAin <- c("DataIn/Parken_evalstocks.csv")
nsim <- 10

lier_rtmb <- function(WAin = c("DataIn/WCVIStocks.csv"),
                          nsim = 10 # default nsim for bootstrapping
                          # remove.EnhStocks = FALSE,
                          # run.bootstraps = TRUE, # to turn on or off the bootstrap
                              # function added at the end
                          # bs_seed = 1, # seed for bootstrapping
                          # bs_nBS = 10, # trials for bootstrapping
                          # plot = FALSE # whether or not to create plots stored in DataOut/
)
{
  # compiler::enableJIT(0) # Run first without just to see if bug is fixed yet
  
  pb <- progress_bar$new(total = nsim)
  
  # Original LambertW0
  LambertW0 <- RTMB:::ADjoint(LambertW0_internal, dLambertW0_internal)
  # New LambertW0
  # See: https://github.com/pbs-assess/renewassess/blob/main/code/RTMB/PosteriorPredictiveSample.r
  # LambertW0 <- ADjoint(
  #   function(x){gsl::lambert_W0(x)},
  #   function(x, y, dy) {dy / (x + exp(y))}
  # )
  
  # Data ####
  srdatwna <- read.csv(here::here("DataIn/SRinputfile.csv"))
  WAbase <- read.csv(here::here("DataIn/WatershedArea.csv"))
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
  
  WAbase <- WAbase %>% 
    full_join(names, by="Name") %>% 
    arrange(Stocknumber) %>%
    mutate(logWA = log(WA)) 
  
  ## Shift log WA for the mean - base - makes estimation easier
  mean_logWA <- mean(WAbase$logWA)
  WAbase$logWAshifted <- WAbase$logWA - mean_logWA
  
  WAin$logWA <- log(WAin$WA)
  WAin$logWAshifted_t <- WAin$logWA - mean_logWA
  
  lifehist <- srdat %>% dplyr::select(Stocknumber, Name, Stream) %>% 
    group_by(Stocknumber) %>% 
    summarize(lh=max(Stream))
  
  ## RTMB data and par setup ####
  # Data
  dat <- list(srdat = srdat,
              WAbase = WAbase,
              WAin = WAin,
              logRS = log(srdat$Rec) - log(srdat$Sp))
  
  # External vectors
  N_Stk <- max(srdat$Stocknumber + 1) # 25
  
  # Parameters
  par <- list(b0 = c(10, 10), # Initial values for WA regression intercepts
              bWA = c(0, 0), # Inital values for WA regression slopes
              
              # logE = numeric(N_Stk), # Was logE0
              logE0 = numeric(N_Stk),
              # logE_tar = numeric(nrow(dat$WAin)), # Targets, was logE0_tar
              logE0_tar = numeric(nrow(dat$WAin)),
              rm = 1,
              # hj = numeric(nrow(dat$WAbase)), # numeric(nrow(dat$WAbase)) or numeric(N_Stk)
              # hj_pred = numeric(nrow(dat$WAin)),

              tauobs = 0.01 + numeric(N_Stk),

              logESD = 1,
              logAlphaSD = 1
  )

  f_srep <- function(par){
    getAll(dat, par)
    
    N_Stk = max(srdat$Stocknumber + 1)
    stk = srdat$Stocknumber + 1
    N_Obs = nrow(srdat)
    N_Pred = nrow(WAin) # number of predicted watershed areas
    
    S = srdat$Sp
    type = lifehist$lh + 1
    
    E <- numeric(N_Stk)
    logE <- numeric(N_Stk)
    logAlpha <- numeric(N_Stk) # comment on or off
    
    E_tar <- numeric(N_Pred)
    logE_tar <- numeric(N_Pred)
    logAlpha_tar <- numeric(N_Pred)
    
    # logE_pred <- numeric(N_Stk)
    # logE_pred_tar <- numeric(N_Pred)
    
    nll <- 0
    
    nll <- nll - sum(dnorm(b0[1], 10, sd = 31.6, log = TRUE))
    nll <- nll - sum(dnorm(b0[2], 0, sd = 31.6, log = TRUE))
    nll <- nll - sum(dnorm(bWA[1], 0, sd = 31.6, log = TRUE)) # Was 10, 31.6
    nll <- nll - sum(dnorm(bWA[2], 0, sd = 31.6, log = TRUE))
    
    nll <- nll - sum(dnorm(rm, 0.6, sd = 0.45, log = TRUE)) # New position.
    
    ## Second level of hierarchy - Ricker parameters:
    for ( i in 1:N_Stk){
      # logE_pred[i] <- b0[type[i]] + bWA[type[i]]*WAbase$logWAshifted[i] # + logE0[i]
      # nll <- nll - sum(dnorm(logE[i], logE_pred[i], sd = logESD, log = TRUE))
      # E[i] <- exp(logE[i])
      
      nll <- nll - sum(dnorm(logE0[i], 0, sd = logESD, log = TRUE))
      logE <- b0[type[i]] + bWA[type[i]]*WAbase$logWAshifted[i] + logE0[i]
      E[i] <- exp(logE)

            
      nll <- nll - sum(dnorm(hj[i], 0, sd  = logAlphaSD, log = TRUE)) # Liermann prior table
      logAlpha[i] <- rm + hj[i]
      
      # nll <- nll - sum(dnorm(logAlpha[i], rm, sd = logAlphaSD, log = TRUE)) # Was hj[i] --> logAlpha[i]
      
      
      nll <- nll - sum(dgamma(tauobs[i], shape = 0.0001, scale = 1/0.0001, log = TRUE))
    }

    ## First level of hierarchy - Ricker model:
    for (i in 1:N_Obs){
      logRS_pred <- logAlpha[stk[i]]*(1 - S[i]/E[stk[i]])
      nll <- nll - sum(dnorm(logRS[i], logRS_pred, sd = sqrt(1/tauobs[stk[i]]), log = TRUE))
    }
    
    
    
    ## PREDICTIONS
    BETA = numeric(nrow(WAin))
    SMSY = numeric(nrow(WAin))
    SGEN = numeric(nrow(WAin))
    
    for (i in 1:N_Pred){
      # nll <- nll - sum(dnorm(hj_tar[i], 0, sd = logAlphaSD, log  = TRUE))
      # logAlpha_tar[i] <- rm + hj_tar[i]
      # nll <- nll - sum(dnorm(hj_tar[i], rm, sd = logAlphaSD, log = TRUE))
              # hj_pred[i] goes onwards instead of logAlpha_tar[i]
      
      logAlpha_tar[i] <- dnorm(logAlpha_tar[i], rm, sd = logAlphaSD, log = TRUE)
              # What is this alpha?

      # Predict E for target watershed areas
      # logE_pred_tar[i] <- b0[type[i]] + bWA[type[i]]*WAin$logWAshifted_t[i] # + logE0_tar[i]
      # logE_tar[i] <- dnorm(logE_tar[i], logE_pred_tar[i], sd = logESD, log = TRUE)
      # E_tar[i] <- exp(logE_tar[i])

      logE0_tar[i] <- dnorm(logE0_tar[i], 0, sd = logESD, log = TRUE)
      logE_tar[i] <- b0[type[i]] + bWA[type[i]]*WAin$logWAshifted_t[i] + logE0_tar[i]
      E_tar[i] <- exp(logE_tar[i])

      # Predict BETA
      BETA[i] <- logAlpha_tar[i]/E_tar[i]
      # Predict SMSY
      SMSY[i] <- (1-LambertW0(exp(1-logAlpha_tar[i])))/BETA[i]
      # Predict SGEN
      SGEN[i] <- -1/BETA[i]*LambertW0(-BETA[i]*SMSY[i]/(exp(logAlpha_tar[i])))
    }
    
    ## ADREPORT - internals
      # ALWAYS THE NUMBER OF SYNOPTIC SETS - RN 25
    # ADREPORT(logRS) # logRS for all 501 data points
    ADREPORT(E)
    ADREPORT(logAlpha)
  
    REPORT(E) # E (Srep) for all synoptic data set rivers (25)
    REPORT(logAlpha) # model logAlpha (25)
    
    # ADREPORT - predicted
      # NUMBER OF PREDICTED STOCKS WILL CHANGE
    # Mean estimate of the median (without bias correction)
    ADREPORT(E_tar) # target E (Srep) (21)
    ADREPORT(logE_tar) # exp these for the correct confidence intervals
    ADREPORT(logAlpha_tar)
    
    REPORT(E_tar)
    REPORT(logAlpha_tar)
    
    ADREPORT(BETA)
    ADREPORT(SMSY)
    ADREPORT(SGEN)
    # Symmetrical - would need to get the confidence interval on the log-scale and then exp()
    
    REPORT(BETA)
    REPORT(SMSY)
    REPORT(SGEN)
    
    nll
  }
  
  ## MakeADFun ####
  obj <- RTMB::MakeADFun(f_srep,
                         par,
                         silent=TRUE)
  # obj <- RTMB::MakeADFun(f_srep,
  #                        par,
  #                        # random = c("hj", "hj_pred", "logE0", "logE0_pred"),
  #                        # random = c("hj", "logE0"),
  #                        random = c("logE0", "rm"), # ?
  #                        silent=TRUE)
  
  # New limits
  upper <- numeric(length(obj$par)) + Inf
  lower <- numeric(length(obj$par)) + -Inf
  lower[names(obj$par) == "tauobs"] <- 0
  upper[names(obj$par) == "logESD"] <- 100
  lower[names(obj$par) == "logESD"] <- 0
  upper[names(obj$par) == "logAlphaSD"] <- 100
  lower[names(obj$par) == "logAlphaSD"] <- 0
  # Tor: Limits (upper/lower) for nlminb do not seem to effect the convergence
    # or objective function
  
  opt <- nlminb(obj$par, 
                obj$fn, 
                obj$gr, 
                # control = list(trace = 0) # ,
                control = list(eval.max = 1e5, iter.max = 1e5, trace = 0),
                lower = lower,
                upper = upper
  )
                
  # osdr <- sdreport(obj)
  
  # SAVE MLE ####
  sdr <- RTMB::sdreport(obj)
  sdr_full <- summary(RTMB::sdreport(obj))

  sdr_est <- as.list(sdr, "Est", report=TRUE) ## ADREPORT estimates
  sdr_se <- as.list(sdr, "Std", report=TRUE) ## ADREPORT standard error

  # Create the correct quantiles for E
  # log_E_tar for sdr_est and sdr_se
  # Must first calculate their quantile in log-space AND THEN exponentiate
  E_quant <- cbind(WAin,
    E_tar = exp(sdr_est$log_E_tar), # Mean
    E_LQ = exp(sdr_est$log_E_tar - (sdr_se$log_E_tar*1.96)), # E LQ
    E_UQ = exp(sdr_est$log_E_tar + (sdr_se$log_E_tar*1.96)), # E UQ
    SGEN.Mean = sdr_est$SGEN, # SGEN
    SGEN.LQ = sdr_est$SGEN - (sdr_se$SGEN*1.96),
    SGEN.UQ = sdr_est$SGEN + (sdr_se$SGEN*1.96),
    SMSY.Mean = sdr_est$SMSY, # SMSY
    SMSY.LQ = sdr_est$SMSY - (sdr_se$SMSY*1.96),
    SMSY.UQ = sdr_est$SMSY + (sdr_se$SMSY*1.96)
  )
  
  rtmb_full <- E_quant
  # rtmb <- E_quant # Saved as version with no nll for predictions
  
  # Simulate #### THIS SECTION IS COMMENTED OUT ####
  # sgen = smsy = beta = NULL
  #     # obj$simulate() # 1000 times in a for loop and then track it - as a bootstrap
  # nsim <- nsim # number of sims
  # 
  # for (i in 1:nsim){
  #       # pb$tick()
  #   temp <- obj$simulate()
  #       # error in logAlpha[i] <- rm + hj[i]
  #       # Simulate one or more responses from the distribution corresponding to a fitted model object.
  #   sgen <- rbind(sgen, temp$SGEN)
  #   beta <- rbind(beta, temp$BETA)
  #   smsy <- rbind(smsy, temp$SMSY)
  # }
  # # quantiles by row - apply
  # # dim(x) must have a positive length
  # SGEN_boot <- apply(sgen, 2, FUN = function(x){c(mean(x), quantile(x, c(0.0275,0.5,0.975)))})
  # SMSY_boot <- apply(smsy, 2, FUN = function(x){c(mean(x), quantile(x, c(0.0275,0.5,0.975)))})
  # BETA_boot <- apply(beta, 2, FUN = function(x){c(mean(x), quantile(x, c(0.0275,0.5,0.975)))})
  # 
  # # bind together bootstrapped values
  # # transpose so that it is 115 rows instead of columns
  # # and then add in an identifier
  # # and the rbind them together - since the columns will all be the same
  # SGEN_boot_ <- t(as.data.frame(SGEN_boot)) # transpose SGEN_boot (once as a df)
  # SMSY_boot_ <- t(as.data.frame(SMSY_boot))
  # BETA_boot_ <- t(as.data.frame(BETA_boot))
  # # rename columns
  # colnames(SGEN_boot_) <- c("Median","LQ","Mean","UQ")
  # colnames(SMSY_boot_) <- c("Median","LQ","Mean","UQ")
  # colnames(BETA_boot_) <- c("Median","LQ","Mean","UQ")
  # # now rows are just numbered - can just cbind into main df since they should be in same order
  # # and are of length (nrow(t(test))) for e.g.
  # # nrow(SGEN_boot_)
  # 
  # # CURRENTLY NaN Std. Error ***************************************************
  # sdr <- RTMB::sdreport(obj)
  # sdr_full <- summary(RTMB::sdreport(obj))
  # 
  # sdr_est <- as.list(sdr, "Est", report=TRUE) ## ADREPORT estimates
  # sdr_se <- as.list(sdr, "Std", report=TRUE) ## ADREPORT standard error
  # 
  # # Create the correct quantiles for E
  # # log_E_tar for sdr_est and sdr_se
  # # Must first calculate their quantile in log-space AND THEN exponentiate
  # E_quant <- cbind(WAin,
  #                  E_tar = exp(sdr_est$log_E_tar), # Mean
  #                  E_LQ = exp(sdr_est$log_E_tar - (sdr_se$log_E_tar*1.96)), # E LQ
  #                  E_UQ = exp(sdr_est$log_E_tar + (sdr_se$log_E_tar*1.96)) # E UQ
  # )

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

  # WArefpoints <- cbind(E_quant, # Now contains WAin's information
  #                      # E = sdr_est$E_tar,
  #                      # E_se = sdr_se$E_tar,
  #                      logalpha = sdr_est$logAlpha_tar,
  #                      logalpha_se = sdr_se$logAlpha_tar,
  #                      SGEN = SGEN_boot_,
  #                      SMSY = SMSY_boot_,
  #                      BETA = BETA_boot_
  # )
  
  
  # MCMC ####
    # INIT FUNCTION
  init <- function(){
    list(b0 = c(10, 10), # Initial values for WA regression intercepts
         bWA = c(0, 0), # Inital values for WA regression slopes
         
         logE0 = numeric(N_Stk),
         logE0_tar = numeric(nrow(dat$WAin)),
         hj = numeric(nrow(dat$WAbase)), # numeric(nrow(dat$WAbase)) or numeric(N_Stk)
         hj_pred = numeric(nrow(dat$WAin)),
         # WAin or WAbase
         rm = 1,
         
         tauobs = 0.01 + numeric(N_Stk),
         
         logESD = 1,
         logAlphaSD = 1
    )
  }
    # SETTING LIMITS
  upper <- numeric(length(obj$par)) + Inf
  lower <- numeric(length(obj$par)) + -Inf
  
  lower[names(obj$par) == "tauobs"] <- 0
  upper[names(obj$par) == "logESD"] <- 100
  lower[names(obj$par) == "logESD"] <- 0
  upper[names(obj$par) == "logAlphaSD"] <- 100
  lower[names(obj$par) == "logAlphaSD"] <- 0
  
    # SAMPLE
  fitstan <- tmbstan(obj, iter = 5000, warmup = 500, init = init,
                     lower = lower, upper = upper,
                     chains = 4, open_progress = FALSE, silent = TRUE)
  
  return(list(opt = opt,
              obj = obj,
              fitstan = fitstan # MCMC object for processing
              # sdr = sdr,
              # sdr_full = sdr_full,
              # sdr_est = sdr_est,
              # sdr_se = sdr_se,
              # refpoints = WArefpoints
  ))

}

# Run the function test ####
# test_srep <- lier_rtmb(WAin = c("DataIn/Parken_evalstocks.csv"),
#                            # Parken_evalstocks.csv
#                            # WCVIStocks.csv
#                            nsim = 10) # default test run for outputs

source(here::here("R/derived_post.R")) # ~ 4 minutes total run time
derived_obj <- derived_post(fitstan)

# Test plots
traceplot(fitstan, pars=names(obj$par), inc_warmup=TRUE)
pairs(fitstan, pars=names(obj$par))


# Saving for plotting
# rtmb_full <- rtmb_full |> 
#    rename("Stock_name" = Stock_num)
# testfull <- cbind(rtmb_full, derived_obj$deripost_summary$SMSY)

# ******************************************************************************
## MCMC run through tmbstan - OUTSIDE OF FUNCTION ####

# Create a initial value function 
# EXAMPLE from TMB_rtmb.R using the SMAX model version
# initf1_EXAMPLE <- function(){
#   list(logA = (srdat %>% group_by (Stocknumber) %>%
#                  summarise(yi = lm(log(Rec / Sp) ~ Sp)$coef[1]))$yi, # random effect
#        logB = log ( 1/ ( (1/B$m)/dat$scale )), # fixed effect
#        logSigma = rep(-2, length(unique(srdat$Name))),
#        logSigmaA = -2,
#        logMuA_stream = 1.5,
#        logMuA_ocean = 0,
#        logDelta1 = 3,
#        logDelta1_ocean = 0,
#        logDelta2 = log(0.72),
#        Delta2_ocean = 0,
#        logNu1 = 3,
#        logNu1_ocean = 0,
#        logNu2 = log(0.72),
#        Nu2_ocean = 0,
#        logNuSigma = -0.412,
#        logDeltaSigma = -0.412
#   )
# }

# initf1 <- function(){
#   list( # list out pars
#     b0 = c(10, 10),
#     bWA = c(0, 0),
#     logAlpha_re = numeric(nrow(dat$WAbase)), 
#     logAlpha_re_pred = numeric(nrow(dat$WAin)),
#     logE0 = numeric(N_Stk),
#     logE0_ = numeric(nrow(dat$WAin)),
#     tauobs = 0.01 + numeric(N_Stk),
#     # logAlpha0 = 1.5,
#     logESD = 1,
#     logAlphaSD = 10
#   )
# }

# initf2 <- function(){
#   list(b0 = c(10, 10), # Initial values for WA regression intercepts
#        bWA = c(0, 0), # Inital values for WA regression slopes
#        
#        logE0 = numeric(N_Stk),
#        logE0_pred = numeric(nrow(dat$WAin)),
#        hj = numeric(nrow(dat$WAbase)), # numeric(nrow(dat$WAbase)) or numeric(N_Stk)
#        hj_pred = numeric(nrow(dat$WAin)),
#           # WAin or WAbase
#        rm = 1,
#        
#        tauobs = 0.01 + numeric(N_Stk),
#        
#        logESD = 1,
#        logAlphaSD = 1
#        )
# }
# 
# # Run the cores using:
#   # obj <- RTMB object
#   # init <- initial value function defined above
#   # rest of usual chain and iter parameters
# 
# # Create bounds - moved upwards to creation of rtmb obj
# upper <- numeric(length(obj$par)) + Inf
# lower <- numeric(length(obj$par)) + -Inf
# 
# lower[names(obj$par) == "tauobs"] <- 0
# upper[names(obj$par) == "logESD"] <- 100
# lower[names(obj$par) == "logESD"] <- 0
# upper[names(obj$par) == "logAlphaSD"] <- 100
# lower[names(obj$par) == "logAlphaSD"] <- 0
# 
# fitstan <- tmbstan(obj, iter = 5000, warmup = 500, init = initf2,
#                    lower = lower, upper = upper,
#                    chains = 4, open_progress = FALSE, silent = TRUE)
#   # Can I add in a beepr?
# # https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# # Examine the pairs() plot to diagnose sampling problems
# 
# source(here::here("R/derived_post.R")) # ~ 4 minutes total run time
# derived_obj <- derived_post(fitstan)
#   # Will contain df and matrices of derived posteriors
# 
# # Test plots
# traceplot(fitstan, pars=names(obj$par), inc_warmup=TRUE)
# pairs(fitstan, pars=names(obj$par))






# ******************************************************************************
# ********************** TESTING ***********************************************
# ******************************************************************************
# Getting things out 
# Probably use tidybayes or something
## Can extract marginal posteriors easily
post <- as.matrix(fitstan)
  # post is 133 by 18000 as a matrix
  # which is 133: the number of parameters (for each iteration)
  # 18000 is 4500 x 4
hist(post[,'b0[1]'])
hist(post[,'bWA[1]'])

## What if you want a posterior for derived quantities in the report? Just
  ## loop through each posterior sample (row) and call the report function
  ## which returns a list. The last column is the log-posterior density (lp__)
  ## and needs to be dropped

# This only works for single values e.g. sd0 is only element
# obj$report(post[1,-ncol(post)])
  # Can be subset by obj$report(post)$NAME e.g. SGEN or something
  # obj$report(post[18000,133]) - but 133 is a TRUE/FALSE issue
# outs <- rep(NA, length.out=nrow(post))
  # This is a list of 18000 long

# I need a matrix that is:
  # [18000, 132] ? because of the issue of lp__ as above?
  # [rows, columns]

# What is the contents of obj$report(post[,X])
  # Are these the standard error's? Or the false reports of the parameters
  # hence why [,1] and [,2] are massive

# Therefore I don't want to deal with [,X] and only loop through the rows of 
  # samples
  # But I want to loop into a matrix that is 
  # 18000 long and length(names(obj$report(post))) <- 7
  # This will be a long dataframe of all derived posteriors?
  # AND another dimension for all of the 25 stocks
    # e.g. length(obj$report(post)$E)

# deripost <- matrix(0, nrow = nrow(post), ncol = length(obj$report(post)))
# deripost <- data.frame()
  # create a df that has 18000 rows, 7 normal columns, and an 
  # 8th column for stock that is 1:25 repeated
  # This will mean that the total number of rows is 18000 X 25
  # The other way to do this would be 18000 X 7 - and then have 25 + 1 columns

  # obj$report(post[1,])[1][[1]][1] refers to the:
    # 1st value - ie. stock #1
    # in the list of values in the 1st named derived parameter out of 7
    # therefore one could generalize this statement
    # by: obj$report(post[1:18000,])[1:7][[1]][[1:25]]
    # and that would result in a single value out of 18000 X 7 X 25 total values

  # Does it make sense to create a series of lists for each stock or par?

# Kasper's method?
# for(i in 1:nrow(post)){ # 18000 iterations
#   r <- obj$report(post[i,]) # loop through all the obj$report names
#   for (j in 1:length(obj$report(post[1,])$E)) { # 25 stocks
#     
#   }
#     for (k in 1:length(obj$report(post))) { # 7 parameters
#       deripost[i, k] <- r[k]
#         # creating an error b/c r[k] has 25 objects within it
#     }
# }

# Alternate method of 7 matrices of 25 by 18000
  # Create matrices
    # Define dimensions
# n_rows <- 18000
# # n_cols <- 25 # This will be subject to change depending on the variable
# n_cols <- sapply(obj$report(post), function(x) length(x))
# n_matrices <- 7

# This takes A LONG time - requires speeding up
# for (k in 1:length(obj$report(post))) { # 7 Total
#   for (j in 1:length(obj$report(post[1,])[1][[1]])) { # 25 Total
#     for (i in 1:nrow(post)) { # 18000 Total
#       matrices[k][[k]][j] <- obj$report(post[i,])[k][[1]][[j]]
#     }
#   }
# }

# Vectorized form of above
# matrices <- lapply(1:n_matrices, function(x) matrix(NA, nrow = n_rows, ncol = n_cols))
# matrices <- lapply(seq_len(n_matrices), function(k) {
#   matrix(NA, nrow = n_rows, ncol = n_cols[k])
# })
#   # This needs to be changed such that the matrix made is matched to the length
#   # of the parameter it is describing
# names(matrices) <- names(obj$report(post[1,]))
#     # Perform the nested loop
# 
# pb <- txtProgressBar(min = 0, max = n_rows, style = 3)
# report_i <- NULL
# for (i in 1:n_rows) {  # Outer loop over rows of 'post'
#   report_i <- obj$report(post[i,])  # Extract report once per row
#   for (k in 1:n_matrices) {  # Loop over matrices
#     for (j in 1:n_cols) {  # Loop over columns
#       matrices[[k]][i, j] <- report_i[[k]][[j]]  # Fill pre-allocated matrices
#     }
#   }
#   setTxtProgressBar(pb, i)
# }
# close(pb)
#   # This produces 7 different matrices - 1 for each derived parameter
#   # Each matrix is 25 columns by 18000 rows
# 
# pb <- txtProgressBar(min = 0, max = n_rows, style = 3)
# for (i in seq_len(n_rows)) {  
#   report_i <- obj$report(post[i,])  
#   for (k in seq_len(n_matrices)) {  # This is pretty instant
#     matrices[[k]][i, seq_len(n_cols[k])] <- report_i[[k]]
#   }
#   setTxtProgressBar(pb, i)
# }
# close(pb)

# The next step is to make these more usable for final data visualization etc.
  # Therefore we need: 
    # - Mean
    # - Median
    # - Quantiles (5% and 95%)
    # - Std. Error?
  # In the form of columns: Stockname, PAR X mean, PAR X median, PAR X LQ, PAR X UQ, etc...
    # Therefore for each matrix I need the following transformations:
      # - Extract mean, median, and quantile
      # - Create a new dataframe PER PARAMETER w/ Stock, Mean, Median, LQ, UQ columns
      # - Merge these dataframe's together if needed
      # - Or make a generalized statement for them to be added to

# summary(coda::mcmc(matrices$E))
# quantile(matrices$E, probs = c(0.05, 0.95))
# apply(matrices$E, 2 , median) or for mean - 2 refers to BY COLUMN

# The following code needs to be generalized to be looped over
# testdf <- data.frame(
#   Stock = character(25),
#   Mean = numeric(25),
#   Median = numeric(25),
#   LQ = numeric(25),
#   UQ = numeric(25)
# )
# 
# testdf$Stock <- c(1:ncol(matrices$E))
# testdf$Mean <- apply(matrices$E, 2 , mean)
# testdf$Median <- apply(matrices$E, 2 , median)
# testdf$LQ <- apply(matrices$E, 2, quantile, probs = c(0.05, 0.95))[1,] # 0.05
# testdf$UQ <- apply(matrices$E, 2, quantile, probs = c(0.05, 0.95))[2,] # 0.95

# Loop above
  # 1. Create a series of blank df's like the matrices
  # 2. loop through the above creations
# Use replicate to create 7 identical dataframes
# dataframes <- lapply(seq_len(n_matrices), function(k) {
#   data.frame(
#     Stock = character(n_cols[k]),
#     Mean = numeric(n_cols[k]),
#     Median = numeric(n_cols[k]),
#     LQ_5 = numeric(n_cols[k]),
#     UQ_95 = numeric(n_cols[k])
#   )
# })
# names(dataframes) <- names(matrices)
# 
# for (i in 1:n_matrices) {
#   dataframes[[i]]$Stock <- c(1:ncol(matrices[[i]]))
#   dataframes[[i]]$Mean <- apply(matrices[[i]], 2 , mean)
#   dataframes[[i]]$Median <- apply(matrices[[i]], 2 , median)
#   dataframes[[i]]$LQ_5 <- apply(matrices[[i]], 2, quantile, probs = c(0.05, 0.95))[1,] # 0.05
#   dataframes[[i]]$UQ_95 <- apply(matrices[[i]], 2, quantile, probs = c(0.05, 0.95))[2,] # 0.95
# }

# Now have dataframes that can be joined by stock id for E_tar (SREP), SMSY, and SGEN
# ******************************************************************************

# Next step is to plot them against the LoopedIWAMruns.R point comparison plot
  # to see how the MCMC and Liermann model stand up against others
  # Most importantly - how the quantiles appear.






  
# ******************************************************************************
## Let's get posterior samples for all our REPORTS.
# https://github.com/pbs-assess/renewassess/blob/main/material/RTMB%20Recap%20and%20Nimble%20Demo/Jacobians.Rmd
# post <- as.matrix(fit.stan)
# report <- obj$report()
# post.report <- matrix(0,  nrow = nrow(post), ncol = length(report))
# colnames(post.report) <- names(report)
# 
# for (i in 1:nrow(post)) { # 18000 iterations 
#   # i is length of post
#   # need a second one that is k for length of report
#   for (k in 1:length(report)) {
#     post.report[i, k] <- obj$report(post[i,k])
#   }
# }
#   post.report[i,] <- do.call('c', obj$report(post[i,seq_along(obj$par)]))
  
summary(coda::mcmc(post.report))


# ******************************************************************************
# Posterior package usage
library(posterior)
summarise_draws(fitstan)


# Rstan and tidybayes
library(rstan)

str(rstan::extract(fitstan)) # this is a list format of post
wstanout <- rstan::extract(fitstan) # a subsettable version of the above

# Drawing Point summaries and intervals
  # This is only for parameter estimates themselves
fitstan |> 
  spread_draws(logESD) |> 
  median_qi()
