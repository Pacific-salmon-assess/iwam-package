# BOOTSTRAPPING: Liermann Srep (E) RTMB Model with MCMC Sampling ####

# Libaries ####
library(RTMB)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(progress)
library(tmbstan)
library(TMB)
library(tidybayes)
library(bayesplot)
library(beepr) # Sounds
library(viridis)
library(latex2exp)

source(here::here("R/LambertWs.R")) # Lambert W function
source(here::here("R/helperFunctions.R")) # For bootstrapping

# SAVING R OBJECTS: ####
# In script_A.R
# save(derived_obj, file = "derived_obj.RData")
# In script_B.R
load(here::here("derived_obj.RData"))
WAin <- c("DataIn/Parken_evalstocks.csv")
WAin <- read.csv(here::here(WAin))

# Bootstrap Posterior to Match Original IWAM Model ####
    # The point of this is to use the Parken assumptions of productivity
    # and make bootstrapped benchmark estimates from the mcmc chains of the
    # Liermann model.

    # 1. SETUP
BS <- TRUE # default for FALSE
bsiters <- 100 # New with Brown et al. CSAS runs
outBench <- list()
outAlpha <- list()
# set.seed <- 1
# set.seed(1) # Now set at header of code
conditional <- TRUE # Default T for conditional - FALSE for Marginal medians
MCMC <- TRUE # Default F to use original method, T takes logA from MCMC samples to estimate SGEN
prod <- c("Parken") # "LifeStageModel" or "Parken"
# bias.cor <- FALSE
  # bias.cor is also an option - but shouldn't be necessary given that these are posteriors

# library(doParallel)
# library(doRNG)
# cl <- makeCluster(detectCores() - 1)
# registerDoParallel(cl)
# registerDoRNG(seed = 1)

# Function
set.seed(1) ; if (BS == TRUE) {
  
  pb <- txtProgressBar(min = 1, max = bsiters, style = 3, title = "Bootstrap")
  start_time <- Sys.time()
  
  # results <- foreach(k = 1:bsiters, .packages = c("dplyr", "purrr")) %dopar% {
  for (k in 1:bsiters) {
    
    # Should this use the conditional or marginal estimate of SREP?
    if (conditional == T) {
    SREP <- derived_obj$deripost_summary$E_tar_adj$Median # Conditional
    SREP_logSE <- (log(SREP) - log(derived_obj$deripost_summary$E_tar_adj$LQ_5)) / 1.96 # Conditional
    SREP_SE <- (SREP - derived_obj$deripost_summary$E_tar_adj$LQ_5) / 1.96 # Conditional
    } else {
    SREP <- derived_obj$deripost_summary$E_tar$Median # Marginal
    SREP_logSE <- (log(SREP) - log(derived_obj$deripost_summary$E_tar$LQ_5)) / 1.96 # Marginal
    SREP_SE <- (SREP - derived_obj$deripost_summary$E_tar$LQ_5) / 1.96 # Marginal
    }
  
    # Steps not included in this bootstrapping function
      # - Removing Stocks: Cypre (LIFESTAGEMODEL ONLY)
        # - Cypre is Stock # 6
        # - Could remove the 6th row per SREP and SMSY since the Stock column
        # isn't used for anything
      # - Scaling: not in Liermann at all
      # - Bias correction: not in Liermann at all
      # - Creation of the RPs dataframe
    
    # 1. LifeStageModel METHOD
      # Could USE POSTERIOR MODE INSTEAD OF MEAN - more likely to be closer to the MLE
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
    
    # 2. PARKEN METHOD
    if (prod == "Parken"){
      est_loga <- function(SMSY, SREP, shortloga=FALSE){
        loga <- nlminb(start = (0.5 - SMSY/SREP) / 0.07, 
                       objective = calc_loga, # Try to fix with a Scheurel version of LW if possible
                       SMSY = SMSY, 
                       SREP = SREP)$par
        if (shortloga) loga <- (0.5 - SMSY/SREP) / 0.07
        beta <- loga/SREP
        return( list( loga = loga , beta = beta, SMSY = SMSY, SREP = SREP) )
      }
      
      # Skip above 
		# - take loga, SREP, SMSY, and BETA from MCMC
	  if (conditional == T) lnalphamc <- derived_obj$deripost_summary$logAlpha_tar_adj$Median else # Conditional
		lnalphamc <- derived_obj$deripost_summary$logAlpha_tar$Median # Marginal
      
      # Again: conditional or marginal estimates? 
      # Mean, Median, or Mode?
      if (conditional == T) SMSY <- derived_obj$deripost_summary$SMSY_adj$Median else # Conditional
		SMSY <- derived_obj$deripost_summary$SMSY$Median # Maringal
      
      # purrr the est_loga function for logA
      if (MCMC == F) lnalpha_Parkin <- purrr::map2_dfr (SMSY, SREP, shortloga=FALSE, est_loga)
        # Would this change if SMSY and SREP are calculated with Median vs. Mean? - Yes
        # Median should be used as it would be closer to the MLE estimate.
        # Other consideration would be NOT sampling, but just running nlminb and then 
          # bootstrapping.
      
      # Or do explicit method
			# Concerned that this method slightly undercalculates alpha
      # SREP_e <- SREP
      # SMSY_e <- SMSY
      # loga_e <- SREP_e * (SMSY_e * gsl::lambert_W0(-exp(1-SREP_e / SMSY_e) * (SREP_e-SMSY_e) / SMSY_e) +
        # SREP_e - SMSY_e) / (SMSY_e * (SREP_e-SMSY_e))
      # beta_e <- loga_e/SREP_e
      
      # Bias correction location #
      
      sREP <- exp(rnorm(length(SREP), log(SREP), SREP_logSE))
      if(min(sREP)<0)   sREP <- exp(rnorm(length(SREP), SREP, SREP_SE$SE))
      
      # Do SGEN calcs with new variables
      if (MCMC == F) SGENcalcs <- purrr::map2_dfr (exp(median(lnalpha_Parkin$loga)), sREP, Sgen.fn2) else
		SGENcalcs <- purrr::map2_dfr (exp(median(lnalphamc)), sREP, Sgen.fn2)
      # SGENcalcs_e <- purrr::map2_dfr (exp(median(loga_e)), SREP_e, Sgen.fn2) # Explicit
    }
    
    setTxtProgressBar(pb, k)
    # Calculate and display elapsed time dynamically
    current_time <- Sys.time()
    elapsed_time <- difftime(current_time, start_time, units = "secs")
    # Print running timer
    cat(sprintf("\rElapsed time: %s seconds", round(as.numeric(elapsed_time), 2)))
    
    # Previous bind location for RPs
    
    out <- list(bench = select(SGENcalcs, -apar, -bpar))
    outBench[[k]] <- out$bench
			# This produces: SGEN (new), SMSY (new), and feeds out used SREP
    
    # oute <- list(bench = select(SGENcalcs_e, -apar, -bpar))
    # outBenche[[k]] <- oute$bench
    
    if (prod == "Parken" & MCMC == F) outA <- list(alpha = exp(median(lnalpha_Parkin$loga))) else
		outA <- list(alpha = exp(median(lnalphamc)))
    if (prod == "LifeStageModel") outA <- list(alpha = Ric.A)
    if (prod == "Parken") outAlpha <- outA
    if (prod == "LifeStageModel") outAlpha[[k]] <- outA
    
  } # ; stopCluster(cl)
  
  close(pb)
  # End timing
  end_time <- Sys.time()
  # Final elapsed time
  total_elapsed_time <- end_time - start_time
  cat("\nTotal time taken:", total_elapsed_time, "\n")
  
  # outBench <- lapply(results, function(x) x$bench)
  # if (prod == "Parken") {
  #   outAlpha <- lapply(results, function(x) x$alpha)
  # } else if (prod == "LifeStageModel") {
  #   outAlpha <- unique(lapply(results, function(x) x$alpha)) # Should only be one
  # }

  # 3. Outputs: Compile bootstrapped estimates of Sgen, SMSY, and SREP, and identify 5th and 
      # 95th percentiles
  SGEN.bs <- select(as.data.frame(outBench), starts_with("SGEN"))
  
  stockNames <- WAin %>% 
    # filter(Stock != "Cypre") %>% # CYPRE FLAG ## ## ## 
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
  
  if (prod == "LifeStageModel") {
    APAR.bs <- select(as.data.frame(outAlpha), starts_with("alpha"))
    rownames(APAR.bs) <- stockNames
    APAR.boot <- data.frame(APAR = apply(APAR.bs, 1, quantile, 0.5), 
                            lwr = apply(APAR.bs, 1, quantile, 0.025),
                            upr = apply(APAR.bs, 1, quantile, 0.975) )
  }
  
  boot <- list(SGEN.boot=SGEN.boot, SMSY.boot=SMSY.boot, 
                SREP.boot=SREP.boot) # , APAR.boot=APAR.boot)
  if (prod == "LifeStageModel") boot$APAR.boot <- APAR.boot
  
  df1 <- data.frame(boot[["SGEN.boot"]], Stock=rownames(boot[["SGEN.boot"]]), RP="SGEN") 
  df1 <- df1 %>% rename(Value=SGEN)
  df2 <- data.frame(boot[["SREP.boot"]], Stock=rownames(boot[["SREP.boot"]]), RP="SREP")
  df2 <- df2 %>% rename(Value=SREP)
  df3 <- data.frame(boot[["SMSY.boot"]], Stock=rownames(boot[["SMSY.boot"]]), RP="SMSY")
  df3 <- df3 %>% rename(Value=SMSY)
  if (prod == "LifeStageModel") {
    df4 <- data.frame(boot[["APAR.boot"]], Stock=rownames(boot[["APAR.boot"]]), RP="APAR")
    df4 <- df4 %>% rename(Value=APAR)
  }
  
  dfout <- add_row(df1, df2)
  dfout <- add_row(dfout, df3)
  if (prod == "LifeStageModel") dfout <- add_row(dfout, df4)
  rownames(dfout) <- NULL
  
  # dfout <- dfout %>% mutate(Value=signif(Value, 2)) %>% # Rounded to 2 signif digits
  #   mutate(lwr=signif(lwr,2)) %>% 
  #   mutate (upr=signif(upr,2))
  
  wasample <- WAin %>% 
        select("Stock", "WA", "lh") %>% 
        mutate(WA = round(WA, 0))
  
  BS.dfout <- merge(dfout, wasample, by="Stock", all.x=TRUE, sort=FALSE)
  if (prod == "Parken") alphaout <- outAlpha$alpha
}; beep(2)

# BS.dfout.LSM <- BS.dfout
# BS.dfout.parken <- BS.dfout

