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
library(progress) # progress bar for iterative loops
library(tmbstan)
library(tidybayes) # Plotting usage

source(here::here("R/LambertWs.R"))

# Wrapper Function ####
compiler::enableJIT(0) # Run first without just to see if bug is fixed yet
# seems not to run on the first attempt - answers [3] instead of [0]

# Internal running
# WAin <- c("DataIn/WCVIStocks.csv")
# WAin <- c("DataIn/Ordered_backcalculated_noagg.csv")
WAin <- c("DataIn/Parken_evalstocks.csv")
nsim <- 10

IWAMsrep_rtmb <- function(WAin = c("DataIn/WCVIStocks.csv"),
                          nsim = 10 # default nsim for bootstrapping
                          # remove.EnhStocks = FALSE,
                          # run.bootstraps = TRUE, # to turn on or off the bootstrap function added at the end
                          # bs_seed = 1, # seed for bootstrapping
                          # bs_nBS = 10, # trials for bootstrapping
                          # plot = FALSE # whether or not to create plots stored in DataOut/
)
{
  
  # Just in case atm
  # compiler::enableJIT(0) # Run first without just to see if bug is fixed yet
  
  pb <- progress_bar$new(total = nsim)
  
  LambertW0 <- RTMB:::ADjoint(LambertW0_internal, dLambertW0_internal)
  
  # Data ####
  srdatwna <- read.csv(here::here("DataIn/SRinputfile.csv"))
  WAbase <- read.csv(here::here("DataIn/WatershedArea.csv"))
  # WAin <- read.csv(here::here("DataIn/WCVIStocks.csv"))
  WAin <- read.csv(here::here(WAin))
  # WAin2 = read.csv(here::here(c("DataIn/Ordered_backcalculated_noagg.csv")))
  
  
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
  
  ## RTMB dat and par setup ####
  
  # Dat
  dat <- list(srdat = srdat,
              WAbase = WAbase,
              WAin = WAin,
              logRS = log(srdat$Rec) - log(srdat$Sp))
  
  # External vectors
  N_Stk <- max(srdat$Stocknumber + 1) # 25
  
  # Parameters/Initial values
  par <- list(b0 = c(10, 10), # Initial values for WA regression intercepts
              bWA = c(0, 0), # Inital values for WA regression slopes
              
              logE0 = numeric(N_Stk),
              logE0_pred = numeric(nrow(dat$WAin)), # Targets
              hj = numeric(nrow(dat$WAbase)), # numeric(nrow(dat$WAbase)) or numeric(N_Stk)
              rm = 1,
              hj_pred = numeric(nrow(dat$WAbase)),

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
    log_E <- numeric(N_Stk)
    logAlpha <- numeric(N_Stk) # comment on or off
    
    E_tar <- numeric(N_Pred)
    log_E_tar <- numeric(N_Pred)
    logAlpha_tar <- numeric(N_Pred)
    
    
    nll <- 0
    
    nll <- nll - sum(dnorm(b0[1], 10, sd = 31.6, log = TRUE))
    nll <- nll - sum(dnorm(b0[2], 0, sd = 31.6, log = TRUE))
    nll <- nll - sum(dnorm(bWA[1], 10, sd = 31.6, log = TRUE)) 
    nll <- nll - sum(dnorm(bWA[2], 0, sd = 31.6, log = TRUE))
    
    ## Second level of hierarchy - Ricker parameters:
    for ( i in 1:N_Stk){
      nll <- nll - sum(dnorm(logE0[i], 0, sd = logESD, log = TRUE))

      log_E <- b0[type[i]] + bWA[type[i]]*WAbase$logWAshifted[i] + logE0[i]
      E[i] <- exp(log_E)
      
      logAlpha[i] <- rm + hj[i]
      nll <- nll - sum(dnorm(rm, 0.6, sd = 0.45, log = TRUE)) # Liermann prior table
      nll <- nll - sum(dnorm(hj[i], 0, sd  = logAlphaSD, log = TRUE)) # Liermann prior table
      
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
    
    # MAKE SURE THIS MIRRORS ABOVE - with terms - match terms to param list too *****************************************
    for (i in 1:N_Pred){
      nll <- nll - sum(dnorm(hj_pred[i], 0, sd = logAlphaSD, log  = TRUE)) # **********************************************
      # nll <- nll - sum(dnorm(rm, 0.6, sd = 0.45, log = TRUE)) # IS THIS NECESSARY HERE?
      logAlpha_tar[i] <- rm + hj_pred[i] # ******************************************************************************

      # Predict E for target watershed areas
      nll <- nll - sum(dnorm(logE0_pred[i], 0, sd = logESD, log = TRUE))
      log_E_tar[i] <- b0[type[i]] + bWA[type[i]]*WAin$logWAshifted_t[i] + logE0_pred[i]
      E_tar[i] <- exp(log_E_tar[i])

      # Predict BETA
      BETA[i] <- logAlpha_tar[i]/E_tar[i]
      # Predict SMSY
      SMSY[i] <- (1-LambertW0(exp(1-logAlpha_tar[i])))/BETA[i]
      # Predict SGEN
      SGEN[i] <- -1/BETA[i]*LambertW0(-BETA[i]*SMSY[i]/(exp(logAlpha_tar[i])))
    }
    
    
    ## ADREPORT - internals
    # ADREPORT(logRS) # logRS for all 501 data points
    # I think it might be worth it not to report for model speed?
    # Question: Does ADREPORT slow down RTMB? It is more to calculate.
    REPORT(E) # E (Srep) for all synoptic data set rivers (25)
    REPORT(logAlpha) # model logAlpha (25)
    
    # ADREPORT - predicted
    # Mean estimate of the median (without bias correction)
    # ADREPORT(E_tar) # target E (Srep) (21)
    # ADREPORT(log_E_tar) # exp these for the correct confidence intervals
    # ADREPORT(logAlpha_tar)
    REPORT(E_tar)
    REPORT(logAlpha_tar)
    
    # ADREPORT(BETA)
    # ADREPORT(SMSY)
    # ADREPORT(SGEN)
    # Symmetrical - would need to get the confidence interval on the log-scale and then exp()
    REPORT(BETA)
    REPORT(SMSY)
    REPORT(SGEN)
    
    nll
  }
  
  ## MakeADFun ####
  # obj <- RTMB::MakeADFun(f_srep,
  #                        par,
  #                        random = c("logAlpha_re", "logAlpha_re_pred", "logE0", "logE0_"),
  #                        silent=TRUE)
  obj <- MakeADFun(f_srep,
                         par,
                         random = c("hj", "hj_pred", "logE0", "logE0_pred"),
                         silent=TRUE)
  opt <- nlminb(obj$par, 
                obj$fn, 
                obj$gr, 
                control = list(trace = 0)) 
  
  # osdr <- sdreport(obj)
  
  # Simulate ####
  sgen = smsy = beta = NULL
  # obj$simulate() # 1000 times in a for loop and then track it - as a bootstrap
  nsim <- nsim # number of sims
  
  for (i in 1:nsim){
    # pb$tick()
    temp <- obj$simulate()
    # Simulate one or more responses from the distribution corresponding to a fitted model object.
    sgen <- rbind(sgen, temp$SGEN)
    beta <- rbind(beta, temp$BETA)
    smsy <- rbind(smsy, temp$SMSY)
  }
  # quantiles by row - apply
  # dim(x) must have a positive length
  SGEN_boot <- apply(sgen, 2, FUN = function(x){c(mean(x), quantile(x, c(0.0275,0.5,0.975)))})
  SMSY_boot <- apply(smsy, 2, FUN = function(x){c(mean(x), quantile(x, c(0.0275,0.5,0.975)))})
  BETA_boot <- apply(beta, 2, FUN = function(x){c(mean(x), quantile(x, c(0.0275,0.5,0.975)))})
  
  # bind together bootstrapped values
  # transpose so that it is 115 rows instead of columns 
  # and then add in an identifier
  # and the rbind them together - since the columns will all be the same
  SGEN_boot_ <- t(as.data.frame(SGEN_boot)) # transpose SGEN_boot (once as a df)
  SMSY_boot_ <- t(as.data.frame(SMSY_boot))
  BETA_boot_ <- t(as.data.frame(BETA_boot))
  # rename columns
  colnames(SGEN_boot_) <- c("Median","LQ","Mean","UQ")
  colnames(SMSY_boot_) <- c("Median","LQ","Mean","UQ")
  colnames(BETA_boot_) <- c("Median","LQ","Mean","UQ")
  # now rows are just numbered - can just cbind into main df since they should be in same order
  # and are of length (nrow(t(test))) for e.g.
  # nrow(SGEN_boot_)
  
  sdr <- sdreport(obj)
  sdr_full <- summary(RTMB::sdreport(obj))
  
  sdr_est <- as.list(sdr, "Est", report=TRUE) ## ADREPORT estimates
  sdr_se <- as.list(sdr, "Std", report=TRUE) ## ADREPORT standard error
  
  # Create the correct quantiles for E
  # log_E_tar for sdr_est and sdr_se
  # Must first calculate their quantile in log-space AND THEN exponentiate
  E_quant <- cbind(WAin,
                   E_tar = exp(sdr_est$log_E_tar), # Mean
                   E_LQ = exp(sdr_est$log_E_tar - (sdr_se$log_E_tar*1.96)), # E LQ
                   E_UQ = exp(sdr_est$log_E_tar + (sdr_se$log_E_tar*1.96)) # E UQ
  )
  
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
  
  WArefpoints <- cbind(E_quant, # Now contains WAin's information
                       # E = sdr_est$E_tar,
                       # E_se = sdr_se$E_tar,
                       logalpha = sdr_est$logAlpha_tar,
                       logalpha_se = sdr_se$logAlpha_tar,
                       SGEN = SGEN_boot_,
                       SMSY = SMSY_boot_,
                       BETA = BETA_boot_
  )
  
  return(list(opt = opt,
              obj = obj,
              sdr = sdr,
              sdr_full = sdr_full,
              sdr_est = sdr_est,
              sdr_se = sdr_se,
              refpoints = WArefpoints
  ))
  
}

# Testing ####
test_srep <- IWAMsrep_rtmb(WAin = c("DataIn/Parken_evalstocks.csv"),
                           # Parken_evalstocks.csv
                           # WCVIStocks.csv
                           nsim = 1000) # default test run for outputs


## MCMC run through tmbstan ####

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

initf2 <- function(){
  list(b0 = c(10, 10), # Initial values for WA regression intercepts
       bWA = c(0, 0), # Inital values for WA regression slopes
       
       logE0 = numeric(N_Stk),
       logE0_pred = numeric(nrow(dat$WAin)),
       hj = numeric(nrow(dat$WAbase)), # numeric(nrow(dat$WAbase)) or numeric(N_Stk)
       hj_pred = numeric(nrow(dat$WAbase)),
       rm = 1,
       
       tauobs = 0.01 + numeric(N_Stk),
       
       logESD = 1,
       logAlphaSD = 1
       )
}

# Run the cores using:
# obj <- RTMB object
# init <- initial value function defined above
# rest of usual chain and iter parameters
upper <- numeric(length(obj$par)) + Inf
lower <- numeric(length(obj$par)) + -Inf

lower[names(obj$par) == "tauobs"] <- 0
upper[names(obj$par) == "logESD"] <- 100
lower[names(obj$par) == "logESD"] <- 0
upper[names(obj$par) == "logAlphaSD"] <- 100
lower[names(obj$par) == "logAlphaSD"] <- 0

fitstan <- tmbstan(obj, iter=2000, warmup=200, init=initf2,
                   lower = lower, upper = upper,
                   chains=4, open_progress=FALSE, silent=TRUE)
# https://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup
# Examine the pairs() plot to diagnose sampling problems

# Getting things out 
# Probably use tidybayes or something
traceplot(fitstan, pars=names(obj$par), inc_warmup=TRUE)
pairs(fitstan, pars=names(obj$par))


## Can extract marginal posteriors easily
post <- as.matrix(fitstan)
# hist(post[,'b0[1]'])
# hist(post[,'logAlpha0'])
## What if you want a posterior for derived quantities in the report? Just
## loop through each posterior sample (row) and call the report function
## which returns a list. The last column is the log-posterior density (lp__)
## and needs to be dropped
# This only works for single values
obj$report(post[1,-ncol(post)])         # sd0 is only element
SMSY_mc <- rep(NA, length.out=nrow(post))
for(i in 1:nrow(post)){
  r <- obj$report(post[i,-ncol(post)])
  SMSY_mc[i] <- r$SMSY_mc
}
hist(sd0)



## Let's get posterior samples for all our REPORTS.post <- as.matrix(fit.stan)
post <- as.matrix(fitstan)
report <- obj$report() # REPORTS ONLY IF REPORT IS IN TMB FUNCTION
post.report <- matrix(0,  nrow = nrow(post), ncol = length(report))
colnames(post.report) <- names(report)
for( i in 1:nrow(post) )
  post.report[i,] <- do.call('c', 
                             obj$report(post[i,seq_along(obj$par)]))


## Let's get posterior samples for all our REPORTS.
post <- as.matrix(fitstan)
post.report <- NULL
for( i in 1:nrow(post) )
  post.report <- rbind(post.report, 
                       unlist(obj$report(post[i,seq_along(obj$par)])))


# shinystan example:
# library(shinystan)
# launch_shinystan(fitstan)





