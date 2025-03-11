# IWAMsrep RTMB Model with MCMC Sampling ####

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
source(here::here("R/derived_post.R")) # ~ 4 minutes total run time

# Wrapper Function ####
compiler::enableJIT(0) # Run first without just to see if bug is fixed yet
# seems not to run on the first attempt - answers [3] instead of [0]

# Saves for internal runs without the wrapper function
# WAin <- c("DataIn/Parken_evalstocks.csv")
WAin <- c("DataIn/WCVIStocks.csv") # For use in comparing StockRecruit relationships
  # Unless you don't want Parken estimates alongside
# WAin <- c("DataIn/WCVIStocks_NanPunt.csv") # For testing

# nsim <- 10
# pb <- progress_bar$new(total = nsim)

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
# Remove Hoko and Hoh stocks - consider removing SKagit OR Chehalis
    # Due to sampling reasons explained in Parken 2006.
srdatwna <- srdatwna %>% 
  filter(!Name %in% c("Hoko","Hoh")) # |> 
  # filter( !(Name == "Skagit")) # |> # Skagit is #22
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
    levels = c("stream", "ocean"))) |> 
  mutate(Stocknumber = as.integer(factor(Stocknumber)) - 1)  # Re-numbering uniquely

names <- srdat %>% 
  dplyr::select (Stocknumber, Name, lh) %>% 
  distinct()

# Remove Skagit or Chehalis from WAbase
# WAbase <- WAbase |> 
#   filter(!(Name == "Skagit"))
  # filter(!(Name == "Chehalis"))

WAbase <- WAbase %>% 
  full_join(names, by="Name") %>% 
  arrange(Stocknumber) %>%
  mutate(logWA = log(WA)) |> 
  filter(!is.na(Stocknumber))

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
            lineWA =  seq(min(WAbase$logWAshifted), max(WAbase$logWAshifted), 0.1),
            logRS = log(srdat$Rec) - log(srdat$Sp))

# External vectors
N_Stk <- max(srdat$Stocknumber + 1) # 25
# N_Stk <- length(unique(srdat$Stocknumber)) #

# Parameters/Initial values
par <- list(b0 = c(10, 0), # Initial values for WA regression intercepts
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
  type = lifehist$lh # + 1
  type_tar = WAin$lh # used to be lh_new --> now changed "lh to lh_old"
  
  E <- numeric(N_Stk)
  logE_pred <- numeric(N_Stk)
  logE <- numeric(N_Stk)
  logAlpha <- numeric(N_Stk)
  
  E_tar <- numeric(N_Pred)
  logE_tar <- numeric(N_Pred)
  
  logAlpha_tar <- numeric(N_Pred)
  
  # Imaginary line vectors
  line <- length(lineWA)
  logE_line_stream <- numeric(line)
  E_line_stream <- numeric(line)
  logE_line_ocean <- numeric(line)
  E_line_ocean <- numeric(line)
  
  # Begin negative log-likelihood
  nll <- 0
  
  nll <- nll - sum(dnorm(b0[1], 10, sd = 31.6, log = TRUE)) # Prior
  nll <- nll - sum(dnorm(b0[2], 0, sd = 31.6, log = TRUE)) # Prior - OG 0
  nll <- nll - sum(dnorm(bWA[1], 0, sd = 31.6, log = TRUE)) # Prior - OG 10
  nll <- nll - sum(dnorm(bWA[2], 0, sd = 31.6, log = TRUE)) # Prior
  
  nll <- nll - sum(dnorm(logAlpha0, 0.6, sd = 0.45, log = TRUE)) # Prior
  
  ## Second level of hierarchy - Ricker parameters:
  for (i in 1:N_Stk){
    nll <- nll - dnorm(logE_re[i], 0, sd = logESD, log = TRUE) # Unobserved
    # logE[i] <- b0[type[i]] + bWA[type[i]]*WAbase$logWAshifted[i] + logE_re[i]
    logE[i] <- b0[1] + b0[2]*type[i] + (bWA[1] + bWA[2]*type[i]) * WAbase$logWAshifted[i] + logE_re[i]
    # nll <- nll - dnorm(logE[i], logE_pred[i], sd = logESD, log = TRUE)
    E[i] <- exp(logE[i])
    
    nll <- nll - dnorm(logAlpha_re[i], 0, sd  = logAlphaSD, log = TRUE)
    logAlpha[i] <- logAlpha0 + logAlpha_re[i] 
    # nll <- nll - dnorm(logAlpha[i], logAlpha0, sd = logAlphaSD, log = TRUE) 
    
    nll <- nll - dgamma(tauobs[i], shape = 0.0001, scale = 1/0.0001, log = TRUE)
  }

  ## First level of hierarchy - Ricker model:
  for (i in 1:N_Obs){
    logRS_pred <- logAlpha[stk[i]]*(1 - S[i]/E[stk[i]])
    nll <- nll - dnorm(logRS[i], logRS_pred, sd = sqrt(1/tauobs[stk[i]]), log = TRUE)
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
    
    ## logE_tar[i] <- b0[type[i]] + bWA[type[i]]*WAin$logWAshifted_t[i]
    # logE_tar[i] <- b0[1] + b0[2]*type[i] + (bWA[1] + bWA[2]*type[i])*WAin$logWAshifted_t[i]
    logE_tar[i] <- b0[1] + b0[2]*type_tar[i] + (bWA[1] + bWA[2]*type_tar[i])*WAin$logWAshifted_t[i]
    E_tar[i] <- exp(logE_tar[i])

    # Predict BETA
    BETA[i] <- logAlpha_tar[i]/E_tar[i]
    # Predict SMSY
    SMSY[i] <- (1-LambertW0(exp(1-logAlpha_tar[i])))/BETA[i]
    # Predict SGEN
    SGEN[i] <- -1/BETA[i]*LambertW0(-BETA[i]*SMSY[i]/(exp(logAlpha_tar[i])))
  }
  
  for (i in 1:line){
    # Create predictions on an imaginary line
    logE_line_ocean[i] <- b0[1] + b0[2] + (bWA[1] + bWA[2])*lineWA[i]
    E_line_ocean[i] <- exp(logE_line_ocean[i])
    
    logE_line_stream[i] <- b0[1] + (bWA[1])*lineWA[i]
    E_line_stream[i] <- exp(logE_line_stream[i])
  }
  
  ## ADREPORT - internals
    # ALWAYS THE NUMBER OF SYNOPTIC SETS - RN 25
  # ADREPORT(logRS) # logRS for all 501 data points
  ADREPORT(logE_re)
  REPORT(logE_re)
  
  ADREPORT(E)
  ADREPORT(logAlpha)
  ADREPORT(SMSY_r)
  ADREPORT(BETA_r)
  ADREPORT(tauobs)
  
  REPORT(E) # E (Srep) for all synoptic data set rivers (25)
  REPORT(logAlpha) # model logAlpha (25)
  REPORT(SMSY_r)
  REPORT(BETA_r)
  REPORT(tauobs)
  
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
  
  REPORT(E_line_stream) # Imaginary line values for plotting
  ADREPORT(E_line_stream)
  REPORT(E_line_ocean) # Imaginary line values for plotting
  ADREPORT(E_line_ocean)
  
  nll
}

## MakeADFun ####
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

opt <- nlminb(obj$par, 
              obj$fn, 
              obj$gr, 
              control = list(eval.max = 1e5, iter.max = 1e5, trace = 0),
              lower = lower,
              upper = upper
)

# MCMC ####
# INIT FUNCTION
init <- function(){
  list(b0 = c(10, 0),
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

# Acquire outputs of MCMC ####
derived_obj <- derived_post(fitstan)

# Test and diagnositic plots ####
# traceplot(fitstan, pars=names(obj$par), inc_warmup=TRUE)
# pairs_pars <- c("b0", "bWA", "logAlpha0", "logESD", "logAlphaSD")
# pairs(fitstan, pars = pairs_pars) # for specific par names from above
# fitstan |> rhat() |> mcmc_rhat() + yaxis_text() # rhat plot for assessing rhat of each parameter

# Bootstrap Posterior to Match Original IWAM Model ####
    # The point of this is to use the Parken assumptions of productivity
    # and make bootstrapped benchmark estimates from the mcmc chains of the
    # Liermann model.

    # Follow the code of Get_LRP_bs.R prod == "Parken"

    # 1. SETUP ** (replicated below - DEPRECIATED)
# Read in data
# Filter out SMSY into and object per stock and SREP
# Join above
# Calculate scale **
# Select out SE for SREP 
source(here::here("R/helperFunctions.R"))
bsiters <- 20000
BS <- TRUE # avoid if you don't want to run bootstraps
outBench <- list()
set.seed <- 1
prod <- c("LifeStageModel") # "LifeStageModel"
#bias.cor is also an option - but I don't needed given that these
  # are done on posteriors

    # 1. LifeStageModel METHOD
# USE POSTERIOR MODE INSTEAD OF MEAN - more likely to be closer to the MLE

    # 2. PARKEN METHOD
# Function
if (BS == TRUE) {
  for (k in 1:bsiters) {
        # E_tar, SMSY
    SREP <- derived_obj$deripost_summary$E_tar$Median # Mean or Median?
        # Calculate SE
    # SREP_logSE <- RPs %>% mutate(SE = (log(RPs$SREP) - log(RPs$SREPLL)) / 1.96) # OLD VERSION
    SREP_logSE <- (log(SREP) - log(derived_obj$deripost_summary$E_tar$LQ_5)) / 1.96
    SREP_SE <- (SREP - derived_obj$deripost_summary$E_tar$LQ_5) / 1.96
    
    if (prod == "LifeStageModel") {
      Mean.Ric.A <- 1
      Ric.A <- exp(rnorm(length(SREP), Mean.Ric.A, 0))
      Sig.Ric.A <- 0.51 
     
      Ric.A <- exp(rnorm(length(SREP), Mean.Ric.A, Sig.Ric.A))
      if(min(Ric.A)<=0) Ric.A <- exp(rnorm(length(SREP), Mean.Ric.A, Sig.Ric.A))
      
      # if (bias.cor == TRUE) {
      #   sREP <- exp(rnorm(length(SREP), log(RPs$SREP) - 0.5*SREP_logSE$SE^2, SREP_logSE$SE))
      #   if(min(sREP)<=0)   sREP <- exp(rnorm(length(SREP), log(RPs$SREP) - 0.5*SREP_logSE$SE^2,
      #     SREP_logSE$SE))
      # } else {
        sREP <- exp(rnorm(length(SREP), log(SREP), SREP_logSE))
        if(min(sREP)<=0)   sREP <- exp(rnorm(length(SREP), SREP, SREP_SE))
      # }
  
      SGENcalcs <- purrr::map2_dfr (Ric.A, sREP, Sgen.fn2)
    }
    
    if (prod == "Parken"){
      est_loga <- function(SMSY, SREP, shortloga=FALSE){
        loga <- nlminb(start = (0.5 - SMSY/SREP) / 0.07, 
                       objective = calc_loga, # Try to fix with a Scheurel version of LW if possible
                       SMSY= SMSY, 
                       SREP=SREP)$par
        if(shortloga) loga <- (0.5 - SMSY/SREP) / 0.07
        beta <- loga/SREP
        return( list( loga = loga , beta = beta, SMSY = SMSY, SREP = SREP) )
      }
      
      SMSY <- derived_obj$deripost_summary$SMSY$Median # Mean or Median?
      # RPs_long <- read.csv(here::here(datain))
      # SMSY <- RPs_long %>% filter(Param == "SMSY") %>% select(Estimate) 
      # SREP <- RPs_long %>% filter(Param == "SREP") %>% select(Estimate)
      
      # purrr the est_loga function for logA
      lnalpha_Parkin <- purrr::map2_dfr (SMSY, SREP, shortloga=FALSE, 
                                         est_loga)
      
      # Or do explicit method
      # SREP_e <- SREP
      # SMSY_e <- SMSY
      # loga_e <- SREP_e*(SMSY_e*gsl::lambert_W0(-exp(1-SREP_e/SMSY_e)*(SREP_e-SMSY_e)/SMSY_e) + SREP_e - SMSY_e)/(SMSY_e*(SREP_e-SMSY_e))
      # beta_e <- loga_e/SREP_e
      
      # Bias correction could go here - see earlier versions  
      sREP <- exp(rnorm(length(SREP), log(SREP), SREP_logSE))
      if(min(sREP)<0)   sREP <- exp(rnorm(length(SREP), SREP, SREP_SE$SE))
      
      # Do SGEN calcs with new variables
      SGENcalcs <- purrr::map2_dfr (exp(median(lnalpha_Parkin$loga)), sREP, Sgen.fn2)
      # SGENcalcs_e <- purrr::map2_dfr (exp(median(loga_e)), sREP_e, Sgen.fn2)
    }
    
    # Bind everything together
    # outlist_RPs <- list(Stock = WAin$Stock,
    #   SGEN = SGENcalcs$SGEN,
    #   SMSY = SGENcalcs$SMSY,
    #   SMSYLL = derived_obj$deripost_summary$SMSY$LQ_5,
    #   SMSYUL = derived_obj$deripost_summary$SMSY$UQ_95,
    #   SREP = derived_obj$deripost_summary$E_tar$Mean,
    #   SREPLL = derived_obj$deripost_summary$E_tar$LQ_5,
    #   SREPUL = derived_obj$deripost_summary$E_tar$UQ_95,
    #   a.par = SGENcalcs$apar
    #   )
    
    out <- list(bench = select(SGENcalcs, -apar, -bpar))
    
    # outBench[[k]] <- out # This doesn't translate into the SGEN.bs statement below
    outBench[k] <- out
  }
  # beep(sound = 2)
      # 3. Outputs: Compile bootstrapped estimates of Sgen, SMSY, and SREP, and identify 5th and 
      # 95th percentiles
  SGEN.bs <- select(as.data.frame(outBench), starts_with("SGEN"))
  
  stockNames <- WAin %>% 
    # filter(Stock != "Cypre") %>% # **** CYPRE FLAG
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
  
  wasample <- WAin %>% 
        select("Stock", "WA", "lh") %>% 
        mutate(WA = round(WA, 0))
  
  # This should be identical to the dfout standard of the original IWAM_model.R
  BS.dfout <- merge(dfout, wasample, by="Stock", all.x=TRUE, sort=FALSE)
  
  beep(sound = 2)
}

# Saving for plotting ####
  # Re-order Stocks to be in ascending order by logWA
  # Re-order Stocks to be in order of Ricker variance - which is what term? Is it extracted ???

# Parken Table 1 and 2
Parkentable1 <- read.csv(here::here("DataIn/Parken_Table1n2.csv")) # Test stocks e.g. WCVI stocks
ParkenCaseStudy <- read.csv(here::here("DataIn/Parken_evalstocks.csv")) # Case study stocks

dpars <- derived_obj$deripost_summary

targets <- WAin |> 
  rename("Stock_name" = Stock)

  # Do I have to do this every time? This method seems slow and inefficient
# SMSY Estimate for TARGET STOCKS
targets1 <- cbind(targets, derived_obj$deripost_summary$SMSY) |> 
  rename("SMSY_mean" = Mean, "SMSY_median" = Median,
    "SMSY_LQ_5" = LQ_5, "SMSY_UQ_95" = UQ_95, "SMSY_Stocknum" = Stock,
    "SMSY_Mode" = PosteriorMode)
# SGEN Estimate for TARGET STOCKS
targets2 <- cbind(targets1, derived_obj$deripost_summary$SGEN) |> 
  rename("SGEN_mean" = Mean, "SGEN_median" = Median,
    "SGEN_LQ_5" = LQ_5, "SGEN_UQ_95" = UQ_95, "SGEN_Stocknum" = Stock,
    "SGEN_Mode" = PosteriorMode)
# SREP ESTIMATE FOR TARGET STOCKS
targetsAll <- cbind(targets2, derived_obj$deripost_summary$E_tar) |> 
  rename("E_tar_mean" = Mean, "E_tar_median" = Median,
    "E_tar_LQ_5" = LQ_5, "E_tar_UQ_95" = UQ_95, "E_tar_Stocknum" = Stock,
    "E_tar_Mode" = PosteriorMode)

# Ricker sigma for SYNOPTIC STOCKS
# targets3 <- cbind(targets, derived_obj$deripost_summary$tauobs) |> 
#   rename("tauobs_mean" = Mean, "tauobs_median" = Median,
#     "tauobs_LQ_5" = LQ_5, "tauobs_UQ_95" = UQ_95, "tauobs_Stocknum" = Stock)
    # tauobs is based on the synoptic sets and will now have different lengths

# Pointwise comparison plot data prep ####
parken <- read.csv(here::here("DataIn/Parken_evalstocks.csv"))

parken <- parken |> 
  rename("SMSYp" = SMSY, "SMSYp_5" = SMSY_5, "SMSYp_95" = SMSY_95) |> 
  rename("SREPp" = SREP, "SREPp_5" = SREP_5, "SREPp_95" = SREP_95)

cols <- viridis(8, alpha=0.9, option = "mako", direction = -1)

# ORIGINAL ####
ggplot() +
  
  # Re-written as the Parken model is included in WAin - bound to rtmb MLE version
  geom_errorbar(data = parken, aes(x = Stock, y = SREPp, ymax = SREPp_95, ymin = SREPp_5,
                                       color = "Parken",
                                       width=.1),
                # position = position_nudge(-0.4)
    ) +
  geom_point(data = parken,
             # position = position_nudge(-0.4),
             aes(x = Stock, y = SREPp, color = "Parken")) +
  
  # Add in LIERMANN from Liermann_RTMB_model.R as a global object
  geom_errorbar(data = targetsAll, aes(x = Stock_name,
                                     y = E_tar_mean,
                                     ymax = E_tar_UQ_95, 
                                     ymin = E_tar_LQ_5,
                                 color = "Liermann MCMC",
                                 width=.1),
                position = position_nudge(+0.2)) +
  geom_point(data = targetsAll,
             position = position_nudge(+0.2),
             aes(x = Stock_name, y = E_tar_mean, color = "Liermann MCMC")) +
  
  theme_classic() +
  scale_y_continuous(transform = "log", 
                     breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
  ylab(TeX("$S_{REP}$ Estimate")) +
  xlab("Stock Name") + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust = 1)) +
  scale_color_manual(name='Model',
                     breaks=c('Parken',
                              # 'RTMB MLE',
                              'Liermann MCMC'),
                     values=c('Parken' = "black",
                              # 'RTMB MLE' = "orange",
                              'Liermann MCMC' = "skyblue"))

# ORDERED BY LOG WA ####
ggplot() +
  
  # Re-written as the Parken model is included in WAin - bound to rtmb MLE version
  geom_errorbar(data = parken, aes(x = fct_reorder(Stock, log(WA)), y = SREPp, ymax = SREPp_95, ymin = SREPp_5,
                                       color = "Parken",
                                       width=.1),
                # position = position_nudge(-0.4)
    ) +
  geom_point(data = parken,
             # position = position_nudge(-0.4),
             aes(x = fct_reorder(Stock, log(WA)), y = SREPp, color = "Parken")) +
  
  # Add in LIERMANN from Liermann_RTMB_model.R as a global object
  geom_errorbar(data = targetsAll, aes(x = fct_reorder(Stock_name, logWA),
                                     y = E_tar_mean,
                                     ymax = E_tar_UQ_95, 
                                     ymin = E_tar_LQ_5,
                                 color = "Liermann MCMC",
                                 width=.1),
                position = position_nudge(+0.2)) +
  geom_point(data = targetsAll,
             position = position_nudge(+0.2),
             aes(x = fct_reorder(Stock_name, logWA), y = E_tar_mean, color = "Liermann MCMC")) +
               # shape = factor(lh_new))) + 
  # scale_shape_manual(values = c(1, 16)) + 
  
  theme_classic() +
  scale_y_continuous(transform = "log", 
                     breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
  ylab(TeX("$S_{REP}$ Estimate")) +
  xlab("Stock Name") + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust = 1)) +
  scale_color_manual(name='Model',
                     breaks=c('Parken',
                              # 'RTMB MLE',
                              'Liermann MCMC'),
                     values=c('Parken' = "black",
                              # 'RTMB MLE' = "orange",
                              'Liermann MCMC' = "skyblue"))

# ORDERED BY RICKER VARIANCE - NOT WORKING (length diffs) ####
# First re-order targetsAll - and then match the order of Parken stock?
  # Or just cbind parken to targetsAll and then re-factor them as you plot?
full <- cbind(targetsAll, parken$Stock, parken$SREPp, parken$SREPp_5, parken$SREPp_95)

ggplot() +
  
  # Re-written as the Parken model is included in WAin - bound to rtmb MLE version
  geom_errorbar(data = full, aes(x = fct_reorder(parken$Stock, tauobs_mean), 
                                 y = parken$SREPp, ymax = parken$SREPp_5, ymin = parken$SREPp_95,
                                 color = "Parken", width=.1)) +
  geom_point(data = full, aes(x = fct_reorder(parken$Stock, tauobs_mean), 
                              y = parken$SREPp, color = "Parken")) +
  
  # Add in LIERMANN from Liermann_RTMB_model.R as a global object
  geom_errorbar(data = full, aes(x = fct_reorder(Stock_name, tauobs_mean),
                                 y = E_tar_mean,
                                 ymax = E_tar_UQ_95, 
                                 ymin = E_tar_LQ_5,
                                 color = "Liermann MCMC",
                                 width=.1),
                position = position_nudge(+0.2)) +
  geom_point(data = full,
             position = position_nudge(+0.2),
             aes(x = fct_reorder(Stock_name, tauobs_mean), 
                 y = E_tar_mean, color = "Liermann MCMC")) +
  
  theme_classic() +
  scale_y_continuous(transform = "log", 
                     breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
  ylab(TeX("$S_{REP}$ Estimate")) +
  xlab("Stock Name") + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust = 1)) +
  scale_color_manual(name='Model',
                     breaks=c('Parken',
                              # 'RTMB MLE',
                              'Liermann MCMC'),
                     values=c('Parken' = "black",
                              # 'RTMB MLE' = "orange",
                              'Liermann MCMC' = "skyblue"))

# Linear Regression: Liermann vs. Parken model ####

# TOR's THOUGHT: maybe Liermann MCMC version is just actually lower and higher - and the points 
  # furthest from the lines are the ones with issues?
  # Why else are the two lines on the outside of the datasets?

    # Step 1. Is the data prepared.
options(scipen = 5) # for ticks without sci. notation
col.use <- NA
for(i in 1:length(WAbase$lh)) {
  if (WAbase$lh[i]== 'stream') col.use[i] <- "forestgreen"  # stream
  else col.use[i] <- "dodgerblue3" # ocean
}

    # Step 2. PLOT
plot(y = dpars$E$Mean, x = WAbase$WA, pch = 20, 
     col = ifelse(WAbase$Name == "Chehalis", 'red', col.use), 
     xlab = expression("Accessible Watershed Area, km"^2), 
     ylab = expression(S[REP]), log = 'xy',
     xlim = c(50,200000) , ylim = c(200,2000000)
  )

points(y =  dpars$E$Mean, x = WAbase$WA, pch = 20, 
       col = ifelse(WAbase$Name == "Chehalis", 'red', col.use), cex = 1.5)

# Target points
# points(y = dpars$E_tar$Mean, x = WAin$WA, pch = 1, col = 'black')

# Parken points of SYNOPTIC SET
# Parkentable1
points(y = Parkentable1$Srep, x = Parkentable1$WA, pch = 1, col = 'black')


    # Step 3. LINES
sum_pars <- summary(fitstan)
bWA1 <- sum_pars$summary[,1][3] 
bWA2 <- sum_pars$summary[,1][4] + bWA1
b01 <- sum_pars$summary[,1][1] 
b01 <- b01 - mean_logWA*bWA1 # To deal with centered/shifted watershed areas
b02 <- sum_pars$summary[,1][2] + sum_pars$summary[,1][1] - mean_logWA*bWA2

simWA <- seq(2, 13, 0.5)
Preds <- b01 + simWA*bWA1
Predso <- b02 + simWA*bWA2
lines(x = exp(simWA), y = exp(Preds), col="forestgreen", lwd = 2, lty = 1)
lines(x = exp(simWA), y = exp(Predso), col="dodgerblue3", lwd = 2, lty = 1)


    # Step 4. Error polygons
# Calculate quantile - uppers and lowers
  # derived_obj$deripost_summary$E has LQ_5 and UQ_95

# pred_lnSMSY IS THE TARGET's
  # predlnWA should be using: WAin$logWAshifted_t
Eline_stream <- derived_obj$deripost_summary$E_line_stream |> 
  rename("line_stocks" = Stock,
    "s_mean" = Mean,
    "s_median" = Median,
    "s_LQ_5" = LQ_5,
    "s_UQ_95" = UQ_95)
Eline_ocean <- derived_obj$deripost_summary$E_line_ocean |> 
  rename("line_stocko" = Stock,
    "o_mean" = Mean,
    "o_median" = Median,
    "o_LQ_5" = LQ_5,
    "o_UQ_95" = UQ_95)
lineWA <- cbind(dat$lineWA, Eline_stream, Eline_ocean) # NEED NEW LINE VALUES

up_S <- lineWA$s_UQ_95
lo_S <- lineWA$s_LQ_5
up_O <- lineWA$o_UQ_95
lo_O <- lineWA$o_LQ_5

polygon(x=c(exp(lineWA$`dat$lineWA` + mean_logWA), exp(rev(lineWA$`dat$lineWA` + mean_logWA))), 
        y=c(up_S, rev(lo_S)), 
        col=rgb(0,0.4,0, alpha=0.2), border=NA)

polygon(x=c(exp(lineWA$`dat$lineWA` + mean_logWA), exp(rev(lineWA$`dat$lineWA` + mean_logWA))), 
        y=c(up_O, rev(lo_O)), 
        col=rgb(0,0.2,0.4, alpha=0.2), border=NA) 

    # Step 5. Grab Parken estimates for the line and add as y = mx + b
# From Table 4. Srep Habitat Model
# stream-type:
    # y = 3.89 + 0.693 * WA + (0.240/2) # variance added
# abline(a = 3.89, b = 0.693, col="forestgreen", lwd=2)
# ocean-type:
    # y = 3.52 + 0.878 * WA + (0.133/2) # variance added
# abline(a = 3.52, b = 0.878, col="dodgerblue3", lwd=2)

Preds_Parken <- 3.89 + simWA*0.693 + (0.240/2)
Predso_Parken <- 3.52 + simWA*0.878 + (0.133/2)
lines(x = exp(simWA), y = exp(Preds_Parken), col="forestgreen", lwd = 2, lty = 2)
lines(x = exp(simWA), y = exp(Predso_Parken), col="dodgerblue3", lwd = 2, lty = 2)

    # Step 6. Add text to describe model equations


# Bar plot comparison of SYNOPTIC values of SREP - NOT WORKING (length diffs) ####
  # Compare Parken and Liermann estimates of SREP for the SYNOPTIC STOCKS
tempEpars <- dpars$E |> 
  rename("E_stock_temp" = Stock)
bardf <- cbind(Parkentable1, tempEpars)

# bcompare <- ggplot(bardf, aes(x = Stock, fill = factor(Stock))) + 
#   geom_bar(aes(y = Srep, fill = "Parken"), stat = 'identity', position = 'dodge', alpha = 0.5) + # Parken
#   geom_bar(aes(y = Mean, fill = "Liermann"), stat = 'identity', position = 'dodge', alpha = 0.5) + # Liermann
#   theme_classic() +
#   scale_fill_manual(name = "Model", values = c("Parken" = "black", "Liermann" = "skyblue")) +
#   # coord_cartesian(ylim = c(0, 100)) +  # Zooms in without removing data
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# bcompare

# Alternative
# Sample Data (Make sure bardf is in long format)
bardf_long <- bardf %>%
  pivot_longer(cols = c(Srep, Mean), names_to = "Category", values_to = "Value")

# Now plot using a single geom_bar()
bcompare <- ggplot(bardf_long, aes(x = Stock, y = Value, fill = Category)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  
  theme_classic() +
  scale_fill_manual(
    name = "Model",  # Custom legend title
    values = c("Srep" = "black", "Mean" = "skyblue"),  # Custom colors
    labels = c("Srep" = "Parken", "Mean" = "Liermann")  # Custom category names
  ) +
  labs(x = "Stock",
       y = expression(S[REP])) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
bcompare

# SR Curves for individual stocks - NOT FINISHED ####
    # Taken from PlotFunctions.R - PlotSRCurve

# Function inputs: srdat, pars, SMSY_std = NULL, StksNum_ar = NULL
  # stdkNum_surv = NULL, stks_surv, r2, removeSkagit (T/F), mod

Stks <- unique(srdat$Stocknumber)
NStks <- length(Stks)
# par(mfrow=c(5,5), mar=c(2, 2, 1, 0.1) + 0.1)
par(mfrow=c(5,5), mar=c(2, 2, 1, 0.1) + 0.1, oma=c(3,3,1,1))

Parken_ab <- read.csv(here::here("DataIn/Parken_Table1n2.csv")) 

for (i in Stks){
  # Get stocknames and numbers
      # names <- pars %>% dplyr::select ("Name", "Stocknumber") %>% distinct()
      # name <- pars %>% filter (Stocknumber==i) %>% 
      #   dplyr::select ("Name") %>% distinct()
  name <- names$Name[names$Stocknumber == i]
  
  R <- srdat %>% filter (Stocknumber==i) %>% 
    dplyr::select(Rec) 
  S <- srdat %>% filter (Stocknumber==i) %>% 
    dplyr::select(Sp) 
      # Sc <- srdat %>% filter (Stocknumber==i) %>% 
      #   dplyr::select(scale) %>% distinct() %>% 
      #   as.numeric()
  
  if(name != "Skagit" & name != "KSR") 
    plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)), ylim=c(0,max(R$Rec) ) )
  if(name == "Skagit") 
    plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)*3), ylim=c(0,max(R$Rec) ) )
  if(name == "KSR") 
    plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,500), ylim=c(0,max(R$Rec) ) )
  
  # Get alpha and beta parameters
  a <- as.numeric(derived_obj$deripost_summary$logAlpha$Median[derived_obj$deripost_summary$logAlpha$Stock - 1 == i])
  b <- as.numeric(derived_obj$deripost_summary$BETA_r$Median[derived_obj$deripost_summary$BETA_r$Stock-1 == i])
  
  # a <- pars %>% filter (Stocknumber==i) %>% 
  #   filter(Param=="logA") %>% 
  #   summarise(A=exp(Estimate)) %>% 
  #   as.numeric()
  # # Divide b by scale
  # b <- pars %>% filter (Stocknumber==i) %>% 
  #   filter(Param=="logB") %>% 
  #   summarise(B=exp(Estimate)/Sc) %>% 
  #   as.numeric()
  
  # Parken values for skagit
    # These are from Parken et al. 2006 Table 2 (Ocean life-histories)
    # The complete table can now be found in DataIn/Parken_Table1n2.csv
  skagit_alpha <- 7.74
  skagit_beta <- 0.0000657
  RR_skagit <- NA
  SS <- RR <- RR_parken <- NA
  #RR_std <- NA
  ap <- Parken_ab$Alpha
  bp <- Parken_ab$Beta
  # ap <- ParkenCaseStudy$ # Do the case study stocks have alpha and beta estimates?
  # bp <- ParkenCaseStudy$ # Do the case study stocks have alpha and beta estimates?
  
  for (j in 1:100){ # Creates a step-wise sample line by which to create a line on
    if (i!=22 & i!=7) SS[j] <- j*(max(S$Sp)/100) # IF NOT SKAGIT OR KSR
    if (i==22) SS[j] <- j*(max(S$Sp*3)/100) # Skagit
    if (i==7) SS[j] <- j*(500/100) # KSR
    
    RR[j] <- a * SS[j] * exp(-b * SS[j])
    RR_parken[j] <- ap[i+1] * SS[j] * exp(-bp[i+1] * SS[j])
    
    if(i==22) {RR_skagit[j] <- skagit_alpha * SS[j] * exp(-skagit_beta * SS[j])} # Skagit Line based on
      # alpha and beta from Table 1 and 2 from Parken et al. 2006
  }
  
  col.use <- "black"
  lines(x=SS, y=RR, col=col.use) 
  
  # For Skagit, add Parken et al. 2006 model curve
  #if(removeSkagit==FALSE) {
  if(i==22) lines(x=SS, y=RR_skagit, lty="dashed")# }
  
  # For all stocks, added in Parken et al. 2006 model curve
  lines(x=SS, y=RR_parken, lty="dashed", col="red")
  
  mtext(name, side=3, cex=0.8)
  
  # Plot SMSY_stream (black for std, red for AR(1),  # *************************
    # and dashed for Parken et al. 2006)
  # SMSY <- pars %>% filter (Stocknumber==i) %>% 
  #   filter(Param=="SMSY") %>% 
  #   summarise(SMSY = Estimate * Sc) %>% 
  #   as.numeric()
  # SMSY_ul <- pars %>% filter (Stocknumber==i) %>% 
  #   filter(Param=="SMSY") %>% 
  #   summarise(SMSY_ul = Estimate * Sc + 1.96 * Std..Error * Sc ) %>% 
  #   as.numeric()
  # SMSY_ll <- pars %>% filter (Stocknumber==i) %>% 
  #   filter(Param=="SMSY") %>% 
  #   summarise(SMSY_ul = Estimate * Sc - 1.96 * Std..Error * Sc ) %>% 
  #   as.numeric()
  
  SMSY <- derived_obj$deripost_summary$SMSY_r$Median[derived_obj$deripost_summary$SMSY_r$Stock - 1 == i]
  SMSY_ul <- derived_obj$deripost_summary$SMSY_r$UQ_95[derived_obj$deripost_summary$SMSY_r$Stock - 1 == i]
  SMSY_ll <- derived_obj$deripost_summary$SMSY_r$LQ_5[derived_obj$deripost_summary$SMSY_r$Stock - 1 == i]
  
  abline(v = SMSY, col=col.use, lty='dotted')
  
  polygon(x=c(SMSY_ul, SMSY_ll, SMSY_ll, SMSY_ul), 
          y=c(-10000,-10000,10000+max(R$Rec),10000+max(R$Rec)), 
          col=grey(0.8, alpha=0.4), border=NA )

  # if(!is.null(SMSY_std)) {
  #   SMSY_std <- SMSY_std %>% right_join(names) %>% filter(Name==name$Name) # filter(Stocknumber != 22) 
  #   if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedCombined"|mod=="Ricker_AllMod"){ 
  #     if(i %in% stksNum_ar) abline(v=SMSY_std$Estimate[which(SMSY_std$Stocknumber==i)]*scale.stock[i+1] , col="black")
  #     if(i %in% stksNum_surv) abline(v=SMSY_std$Estimate[which(SMSY_std$Stocknumber==i)]*scale.stock[i+1] , col="black")
  #   }
  # }
  
  # Parken Smsy Estimate (vert. line) from Table 1/2 Parken et al. 2006
    # Stocks not ordered the same way as other files - alphabetical instead
  Parken_smsy <- Parken_ab$Smsy[Parken_ab$Stocknumber == Parken_ab$Stocknumber[i+1]] 
    # ordered by stocknumber - but starting at 1 instead of 0 (instead of 0:24 as per Stks - 1:25 --> i+1)
  abline(v = Parken_smsy, col="red", lty='dotted')
  
  # Add r2 coefficient values
  # if(is.data.frame(r2)==TRUE) {
  #   lab <-  r2 %>% 
  #     filter(Stocknumber==i) %>% 
  #     dplyr::select(r2) %>% 
  #     as.numeric() %>% 
  #     round(2)
  #   legend("topright", legend = "", title= paste0("r2=",lab), bty="n")
  # }
}

# Add an GLOBAL figure axis label across par()
  # x = Spawners
  # y = Recruitment
mtext("Spawners", side = 1, line = 1, outer = TRUE, cex = 1.3)
mtext("Recruitment", side = 2, line  = 1, outer = TRUE, cex = 1.3, las = 0)

# END ###