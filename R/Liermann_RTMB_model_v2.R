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
library(bayesplot)

library(viridis)
library(latex2exp)

source(here::here("R/LambertWs.R"))

# Wrapper Function ####
compiler::enableJIT(0) # Run first without just to see if bug is fixed yet
# seems not to run on the first attempt - answers [3] instead of [0]

# Saves for internal runs without the wrapper function
WAin <- c("DataIn/Parken_evalstocks.csv")
nsim <- 10

# lier_rtmb <- function(WAin = c("DataIn/WCVIStocks.csv"),
#                           nsim = 10 # default nsim for bootstrapping
#                           # remove.EnhStocks = FALSE,
#                           # run.bootstraps = TRUE, # to turn on or off the bootstrap
#                               # function added at the end
#                           # bs_seed = 1, # seed for bootstrapping
#                           # bs_nBS = 10, # trials for bootstrapping
#                           # plot = FALSE # whether or not to create plots stored in DataOut/
# )
# {

  # Just in case atm
  # compiler::enableJIT(0) # Run first without just to see if bug is fixed yet
  
  pb <- progress_bar$new(total = nsim)
  
  # Original LambertW0
  # LambertW0 <- RTMB:::ADjoint(LambertW0_internal, dLambertW0_internal)
  
  # New LambertW0
  # See: https://github.com/pbs-assess/renewassess/blob/main/code/RTMB/PosteriorPredictiveSample.r
  LambertW0 <- ADjoint(
    function(x){gsl::lambert_W0(x)},
    function(x, y, dy) {dy / (x + exp(y))}
  )
  
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
              
              logE_re = numeric(N_Stk),
              # logE = numeric(N_Stk),
              logAlpha0 = 0.6,
              logAlpha_re = numeric(nrow(dat$WAbase)), # numeric(nrow(dat$WAbase)) or numeric(N_Stk)
              # logAlpha = numeric(nrow(dat$WAbase)), # Turn on or off if using a non-zero prior
              # logE0_pred = numeric(nrow(dat$WAin)), # Targets
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
    logE_pred <- numeric(N_Stk)
    logE <- numeric(N_Stk)
    logAlpha <- numeric(N_Stk)
    
    E_tar <- numeric(N_Pred)
    logE_tar <- numeric(N_Pred)
    
    logAlpha_tar <- numeric(N_Pred) # ???? Is this the right length ????
    
    # logE0_pred <- numeric(nrow(dat$WAin))
    # hj_pred <- numeric(nrow(dat$WAin)) # ???? Is this the right length ????
    
    nll <- 0
    
    nll <- nll - sum(dnorm(b0[1], 10, sd = 31.6, log = TRUE)) # Prior
    nll <- nll - sum(dnorm(b0[2], 0, sd = 31.6, log = TRUE)) # Prior
    nll <- nll - sum(dnorm(bWA[1], 10, sd = 31.6, log = TRUE)) # Prior
    nll <- nll - sum(dnorm(bWA[2], 0, sd = 31.6, log = TRUE)) # Prior
    
    nll <- nll - sum(dnorm(logAlpha0, 0.6, sd = 0.45, log = TRUE)) # Prior
    
    ## Second level of hierarchy - Ricker parameters:
    for (i in 1:N_Stk){
      nll <- nll - dnorm(logE_re[i], 0, sd = logESD, log = TRUE) # Unobserved
      logE[i] <- b0[type[i]] + bWA[type[i]]*WAbase$logWAshifted[i] + logE_re[i]
      # nll <- nll - dnorm(logE[i], logE_pred[i], sd = logESD, log = TRUE)
      E[i] <- exp(logE[i])
      
      nll <- nll - dnorm(logAlpha_re[i], 0, sd  = logAlphaSD, log = TRUE) # CHANGED **************************
      logAlpha[i] <- logAlpha0 + logAlpha_re[i] # CHANGED ****************************************************
      # nll <- nll - dnorm(logAlpha[i], logAlpha0, sd = logAlphaSD, log = TRUE) # FIX ************************
      
      nll <- nll - dgamma(tauobs[i], shape = 0.0001, scale = 1/0.0001, log = TRUE)
    }

    ## First level of hierarchy - Ricker model:
    for (i in 1:N_Obs){
      logRS_pred <- logAlpha[stk[i]]*(1 - S[i]/E[stk[i]])
      nll <- nll - dnorm(logRS[i], logRS_pred, sd = sqrt(1/tauobs[stk[i]]), log = TRUE)
    }
    
    
    ## PREDICTIONS
    BETA = numeric(nrow(WAin))
    SMSY = numeric(nrow(WAin))
    SGEN = numeric(nrow(WAin))

    for (i in 1:N_Pred){
      # Posterior of the mean NOT the posterior predicted mean
          # Are you trying to predict onto a new site vs.
          # Predict the mean sitesq
      
        # Are you included the uncertainty of the random effect or not?
      
      # nll <- nll - sum(dnorm(hj_pred[i], 0, sd = logAlphaSD, log  = TRUE))
      # dnorm(hj_pred[i], 0, sd = logAlphaSD, log  = TRUE) 
      # logAlpha_tar[i] <- logAlpha0 + hj_pred[i]
      
      logAlpha_tar[i] <- logAlpha0
      
      # dnorm(logAlpha_tar[i], logAlpha0, sd = logAlphaSD, log = TRUE)
      
      # nll <- nll - sum(dnorm(hj_pred[i], logAlpha0, sd = logAlphaSD, log = TRUE))
      # dnorm(hj_pred[i], logAlpha0, sd = logAlphaSD, log = TRUE)

      # Predict E for target watershed areas
      # nll <- nll - sum(dnorm(logE0_pred[i], 0, sd = logESD, log = TRUE))
      # dnorm(logE0_pred[i], 0, sd = logESD, log = TRUE) 
      # logE_tar[i] <- b0[type[i]] + bWA[type[i]]*WAin$logWAshifted_t[i] + logE0_pred[i]
      logE_tar[i] <- b0[type[i]] + bWA[type[i]]*WAin$logWAshifted_t[i]
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
    REPORT(logE_tar)
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
  # obj <- RTMB::MakeADFun(f_srep,
  #                        par,
  #                        random = c("logAlpha_re", "logAlpha_re_pred", "logE_re", "logE0_"),
  #                        silent=TRUE)
  obj <- RTMB::MakeADFun(f_srep,
                         par,
                         # random = c("logAlpha_re", "hj_pred", "logE_re", "logE0_pred"), # ORIGINAL
                         random = c("logAlpha_re", "logE_re"),
                         # random = c("logAlpha", "logE_re"),
                         silent=TRUE)
  
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
  
# } # Whatever

# MCMC ####
  # INIT FUNCTION
init <- function(){
  list(b0 = c(10, 10),
       bWA = c(0, 0),
       
       logE_re = numeric(N_Stk),
       logAlpha0 = 0.6,
       logAlpha_re = numeric(nrow(dat$WAbase)),
       # logAlpha = numeric(nrow(dat$WAbase)), # Turn off for non-zeroed parameterization
       
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
                   # consider adapt_delta or max_treedepth
                   chains = 4, open_progress = FALSE, silent = TRUE)

# return(list(opt = opt,
#             obj = obj,
#             fitstan = fitstan # MCMC object for processing
#             # sdr = sdr,
#             # sdr_full = sdr_full,
#             # sdr_est = sdr_est,
#             # sdr_se = sdr_se,
#             # refpoints = WArefpoints
# ))

# Run the function test ####
# test_srep <- lier_rtmb(WAin = c("DataIn/Parken_evalstocks.csv"),
#                            # Parken_evalstocks.csv
#                            # WCVIStocks.csv
#                            nsim = 10) # default test run for outputs

source(here::here("R/derived_post.R")) # ~ 4 minutes total run time
derived_obj <- derived_post(fitstan)

# Test plots ####
traceplot(fitstan, pars=names(obj$par), inc_warmup=TRUE)
# names(obj$par)
pairs_pars <- c("b0", "bWA", "logAlpha0", "logESD", "logAlphaSD")
# pairs(fitstan, pars=names(obj$par))
pairs(fitstan, pars = pairs_pars) # for specific par names from above
fitstan |> rhat() |> mcmc_rhat() + yaxis_text() # rhat plot for assessing rhat of each parameter

# Saving for plotting ####

# bind in WAin for stocknames

#  E_quant <- cbind(WAin, #
targets <- WAin |> 
  rename("Stock_name" = Stock)

targets1 <- cbind(targets, derived_obj$deripost_summary$SMSY) |> 
  rename("SMSY_mean" = Mean, "SMSY_median" = Median,
    "SMSY_LQ_5" = LQ_5, "SMSY_UQ_95" = UQ_95, "SMSY_Stocknum" = Stock)
  # Do I have to do this every time?
targets2 <- cbind(targets1, derived_obj$deripost_summary$SGEN) |> 
  rename("SGEN_mean" = Mean, "SGEN_median" = Median,
    "SGEN_LQ_5" = LQ_5, "SGEN_UQ_95" = UQ_95, "SGEN_Stocknum" = Stock)
targetsAll <- cbind(targets2, derived_obj$deripost_summary$E_tar) |> 
  rename("E_tar_mean" = Mean, "E_tar_median" = Median,
    "E_tar_LQ_5" = LQ_5, "E_tar_UQ_95" = UQ_95, "E_tar_Stocknum" = Stock)
  # This method seems slow and inefficient

# Pointwise comparison plot ####
parken <- read.csv(here::here("DataIn/Parken_evalstocks.csv"))

parken <- parken |> 
  rename("SMSY_i" = SMSY, "lwr" = SMSY_5, "upr" = SMSY_95)

cols <- viridis(8, alpha=0.9, option = "mako", direction = -1)

ggplot() +
  
  # Re-written as the Parken model is included in WAin - bound to rtmb MLE version
  geom_errorbar(data = parken, aes(x = Stock, y = SMSY_i, ymax = lwr, ymin = upr,
                                       color = "Parken",
                                       width=.1),
                # position = position_nudge(-0.4)
    ) +
  geom_point(data = parken,
             # position = position_nudge(-0.4),
             aes(x = Stock, y = SMSY_i, color = "Parken")) +
  
  # # Add in RTMB from IWAMsrep_RTMB_model.R as a global object (run internal - function is broken)
  # geom_errorbar(data = testfull, aes(x = Stock_name, y = SMSY.Mean, ymax = SMSY.UQ, ymin = SMSY.LQ,
  #                                  color = "RTMB MLE",
  #                                  width=.1),
  #               position = position_nudge(-0.2)) +
  # geom_point(data = testfull,
  #            position = position_nudge(-0.2),
  #            aes(x = Stock_name, y = SMSY.Mean, color = "RTMB MLE")) +
  
  # Add in LIERMANN from Liermann_RTMB_model.R as a global object
  geom_errorbar(data = targetsAll, aes(x = Stock_name,
                                     y = SMSY_mean,
                                     ymax = SMSY_UQ_95, 
                                     ymin = SMSY_LQ_5,
                                 color = "Liermann MCMC",
                                 width=.1),
                position = position_nudge(+0.2)) +
  geom_point(data = targetsAll,
             position = position_nudge(+0.2),
             aes(x = Stock_name, y = SMSY_mean, color = "Liermann MCMC")) +
  
  theme_classic() +
  scale_y_continuous(transform = "log", 
                     breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
  ylab(TeX("$S_{MSY}$ Estimate")) +
  xlab("Stock Name") + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust = 1)) +
  scale_color_manual(name='Model',
                     breaks=c('Parken',
                              # 'RTMB MLE',
                              'Liermann MCMC'),
                     values=c('Parken' = "black",
                              # 'RTMB MLE' = "orange",
                              'Liermann MCMC' = "skyblue"))
