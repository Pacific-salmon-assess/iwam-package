# BOOTSTRAPPING: Simulating IWAM Model Posteriors under new Prior Assumptions ####

	# Original code from first iteration of IWAM Model by Holt et al. 2023 (CSAS)
	# Second edit for: Liermann Srep (E) RTMB Model with MCMC Sampling
	# Current model usage is: Liermann Smax (1/Beta) RTMB Model

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
library(mvtnorm)

source(here::here("R/LambertWs.R")) # Lambert W function
source(here::here("R/helperFunctions.R")) # For bootstrapping

options(scipen = 999)

    # 1. Data 
# WAin <- c("DataIn/Parken_evalstocks.csv") # c("DataIn/Ordered_backcalculated_noagg.csv")
# WAin <- read.csv(here::here(WAin))

	# 2. Original Model Inputs
# BS <- TRUE # default for FALSE
# bsiters <- 1000 # New with Brown et al. CSAS runs is 20,000 or 250,000
# outBench <- list()
# outAlpha <- list()

# adj <- F # Default T for conditional - FALSE for Marginal medians - FLIPPED ****
	# This references whether or not predictions have included
	# stock-level error (through rnorm)
	# Naming is under investigation
# MCMC <- TRUE # Default F to use original method, T takes logA from MCMC samples to estimate SGEN FOR PARKEN METHOD ****
# prod <- c("LifeStageModel") # "LifeStageModel" or "Parken" or "RunReconstruction"
# bias.cor <- F # True to add a bias correction term for sREP
  # bias.cor is also an option - but shouldn't be necessary given that these are posteriors
# model <- c("SMAX") # Options are SMAX or SREP
# Ricprior <- c(1, 0.51) # mu, sd



# HOW TO SOURCE:
# source(here::here("R/Liermann_RTMB_model_Bootstrap.R")) 
# In your desired script - this will load the required packages, sources, and run the function
# Which you can then adjust to your desires

dobootstrap <- function(BS = TRUE, # Can be removed
						bsiters = 1000, # Simulation iterations
						adj = FALSE, # F = Conditional, T = Marginal FOR posteriors
						bias.cor = FALSE, # Apply bias-correction
						prod = c("LifeStageModel"), # Assumption
						MCMC = TRUE, # RENAME - TRUE use posterior alpha in parken assumption
						model = c("SMAX"), # IWAM Input model version
						Ricprior = c(1, 0.3), # LifeStageModel assumed alpha prior
						round = FALSE, # Round final reported values
						WAinname = c("DataIn/Parken_evalstocks.csv") # Comparison stocks
						) {
	# Begin function
	WAin <- read.csv(here::here(WAinname))	
	outBench <- list()
	outAlpha <- list()
	
	set.seed(1) ; if (BS == TRUE) {
  
    pb <- txtProgressBar(min = 1, max = bsiters, style = 3, title = "Bootstrap")
    start_time <- Sys.time()
  
    for (k in 1:bsiters) {    
		# Use mean SREP on realscale
		# Conditional if bias corrected, Marginal if rnorm (?)
		if (model == 'SREP'){
		if (adj == T) { # These should be means - not medians!!!
			SREP <- dsrep$deripost_summary$SREP_tar_adj$Median 
			SREP_logSE <- (log(SREP) - log(dsrep$deripost_summary$SREP_tar_adj$LQ_5)) / 1.96 
			SREP_SE <- (SREP - dsrep$deripost_summary$SREP_tar_adj$LQ_5) / 1.96
		} else {
			SREP <- dsrep$deripost_summary$SREP_tar$Median
			SREP_logSE <- (log(SREP) - log(dsrep$deripost_summary$SREP_tar$LQ_5)) / 1.96 
			SREP_SE <- (SREP - dsrep$deripost_summary$SREP_tar$LQ_5) / 1.96 
		}}
	
		# SREP summaries FROM the SMAX model)
			# Useful in showing the affect of the global mean alpha on beta
		# if (model == 'SREP'){
		# if (adj == T) { # These should be means - not medians!!!
			# SREP <- dsmax$deripost_summary$SREP_adj$Median 
			# SREP_logSE <- (log(SREP) - log(dsmax$deripost_summary$SREP_adj$LQ_5)) / 1.96 
			# SREP_SE <- (SREP - dsmax$deripost_summary$SREP_adj$LQ_5) / 1.96
		# } else {
			# SREP <- dsmax$deripost_summary$SREP$Median
			# SREP_logSE <- (log(SREP) - log(dsmax$deripost_summary$SREP$LQ_5)) / 1.96 
			# SREP_SE <- (SREP - dsmax$deripost_summary$SREP$LQ_5) / 1.96 
		# }}
	
		if (model == 'SMAX') {
			if (adj == T) { # Looks very skewed in distribution
				SMAX <- dsmax$deripost_summary$SMAX_tar_adj$Median # 25 val
				SMAX_logSE <- (log(SMAX) - log(dsmax$deripost_summary$SMAX_tar_adj$LQ_5)) / 1.96 # 25 val
				SMAX_SE <- (SMAX - dsmax$deripost_summary$SMAX_tar_adj$LQ_5) / 1.96 # 25 val
				SMAX_SD <- sd(SMAX) # 1 val
			} else {
				SMAX <- dsmax$deripost_summary$SMAX_tar$Median # 25 val per pop
				logSMAX <- log(SMAX)
				SMAX_logSE <- (log(SMAX) - log(dsmax$deripost_summary$SMAX_tar$LQ_5)) / 1.96 # 25 val
				SMAX_SE <- (SMAX - dsmax$deripost_summary$SMAX_tar$LQ_5) / 1.96 # 25 val
				logSMAX_SD <- apply(log(dsmax$deripost_full$SMAX_tar), 2, sd) # 25 val per pop
				
				# Pull alpha posteriors in order to get alpha_SD
				logALPHA <- dsmax$deripost_summary$logalpha_tar$Median 
				logALPHA_SD <- apply(dsmax$deripost_full$logalpha_tar, 2, sd)
				
				# Calculate rho cor() for alpha and beta (1/SMAX) as a list PER POPULATION
				rho <- sapply(1:ncol(dsmax$deripost_full$logalpha_tar), 
					function(i) cor(dsmax$deripost_full$logalpha_tar[,i], log(dsmax$deripost_full$SMAX_tar[,i])))
				# rho_re <- sapply(1:nrow(dsmax$deripost_full$Alpha_re), 
					# function(i) cor(dsmax$deripost_full$Alpha_re[i,], log(dsmax$deripost_full$SMAX_tar[i,])))
					# random effects aren't correlated and neither are their posteriors
		}}
	
		# Establish the correlation/covariance matrix
		sigma_matrix <- list()
		for (i in 1:length(rho)) {
			sigma_matrix[[i]] <- matrix(c(logSMAX_SD[i]^2, 
							 rho[i] * logSMAX_SD[i] * logALPHA_SD[i],
							 rho[i] * logSMAX_SD[i] * logALPHA_SD[i], 
							 logALPHA_SD[i]^2),
							 nrow = 2)
		}
	  
		if (model == 'SREP' & prod == 'RunReconstruction'){
		inSREP <- derived_obj$deripost_summary$SREP_tar_adj |>
			mutate(CU = WAin$CU) |>
			mutate(Inlet = WAin$Inlet) |>
			mutate(Stock = WAin$Stock) |>
			# filter (Stock != "Cypre") |>
			rename(inlets=Inlet)
		}
    
		# 1. LifeStageModel METHOD
		  # Could USE POSTERIOR MODE INSTEAD OF MEAN - more likely to be closer to the MLE
		if (prod == "LifeStageModel") {
			Mean.Ric.A <- Ricprior[1] # Default = 1
			Sig.Ric.A <- Ricprior[2] # Default = 0.51
		  
			# if (model == 'SREP') Ric.A <- exp(rnorm(length(SREP), Mean.Ric.A, 0))
			# if (model == 'SMAX') Ric.A <- exp(rnorm(length(SMAX), Mean.Ric.A, 0))
			if (model == 'SREP') Ric.A <- exp(rnorm(length(SREP), Mean.Ric.A, Sig.Ric.A)) 
			if (model == 'SMAX') Ric.A <- exp(rnorm(length(SMAX), Mean.Ric.A, Sig.Ric.A))
			# if (model == 'SREP') if(min(Ric.A)<=0) Ric.A <- exp(rnorm(length(SREP), Mean.Ric.A, Sig.Ric.A))
			# if (model == 'SMAX') if(min(Ric.A)<=0) Ric.A <- exp(rnorm(length(SMAX), Mean.Ric.A, Sig.Ric.A))
		  
			# Bias correction - SREP
			if (model == 'SREP'){
			if (bias.cor == TRUE) { # Only if mean
				sREP <- exp(rnorm(length(SREP), log(SREP) - 0.5*SREP_logSE^2, SREP_logSE))
				if(min(sREP)<=0) sREP <- exp(rnorm(length(SREP), log(SREP) - 0.5*SREP_logSE^2, SREP_logSE))
			} else {
				sREP <- exp(rnorm(length(SREP), log(SREP), SREP_logSE))
				if(min(sREP)<=0)   sREP <- exp(rnorm(length(SREP), SREP, SREP_SE))
			}}
		  
			# Bias correction - SMAX
			if (model == 'SMAX'){
			if (bias.cor == TRUE) { # Only if mean
				sMAX <- exp(rnorm(length(SMAX), log(SMAX) - 0.5*SMAX_logSE^2, SMAX_logSE))
				if(min(sMAX)<=0) sMAX <- exp(rnorm(length(SMAX), log(SMAX) - 0.5*SMAX_logSE^2, SMAX_logSE))
			} else {
				sMAX <- exp(rnorm(length(SMAX), log(SMAX), SMAX_logSE))
				if(min(sMAX)<=0)   sMAX <- exp(rnorm(length(SMAX), SMAX, SMAX_SE))
			}}
		  
			delta <- c()
			adjustment_mean <- c()
			adjustment_sd <- c()
			  
			new_samples <- matrix(NA, nrow = 2, ncol = 25)
			new_beta <- c()
			for (i in 1:length(rho)){
				delta[i] <- log(Ric.A[i]) - logALPHA[i]
				
				adjustment_mean[i] <- rho[i] * (logSMAX_SD[i] / logALPHA_SD[i]) * delta[i]
				adjustment_sd[i] <- sqrt(logSMAX_SD[i]^2 - (rho[i] * logSMAX_SD[i])^2)
				
				new_beta[i] <- rnorm(1, logSMAX[i] + adjustment_mean[i], adjustment_sd[i]) # logSMAX_SD[i]
				# conditional on alpha draw a new beta
				
				# Old mvnorm function:
				# new_samples[,i] <- rmvnorm(n = c(1,1), 
					# mean = c(logSMAX[i] + adjustment[i], log(Ric.A[i])), # ifelse(log.a<0,0,log.a) #  + adjustment[i]
					# sigma = sigma_matrix[[i]])
			}
		  
			# METHODS of Forcing alpha above log(0) AFTER BEING DRAWN
			# new_samples[1,] <- ifelse(new_samples[1,] < 0, 0.01, new_samples[1,]) # Not letting it go negative either 0 or 0.01
			# Ric.A <- ifelse(Ric.A < 1, 1.01, Ric.A)
		  
			if (model == 'SREP') SGENcalcs <- purrr::map2_dfr (Ric.A, sREP, Sgen.fn2)
			# if (model == 'SMAX') SGENcalcs <- purrr::map2_dfr (Ric.A, sMAX, Sgen.fn4)
		  
			# NEW SGENcalcs for mvnorm version:
			# if (model == 'SMAX') SGENcalcs <- purrr::map2_dfr (exp(new_samples[2,]), exp(new_samples[1,]), Sgen.fn4)
			if (model == 'SMAX') SGENcalcs <- purrr::map2_dfr (Ric.A, exp(new_beta), Sgen.fn4)
			# The main difference from .fn2 is that the beta par here is unaffected by the NEW assumed alpha
			# In .fn2 - the b.par is affected by the residual information in SREP AND the NEW assumed alpha - which warps it
			# When you use the SMAX model - but calculate b.par THROUGH SREP - you can get to similar answers
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
		  
			# Skip above - take loga, SREP, SMSY, and BETA from MCMC
			if (adj == T) lnalphamc <- derived_obj$deripost_summary$logAlpha_tar_adj$Median else 
				lnalphamc <- derived_obj$deripost_summary$logAlpha_tar$Median 
		  
			# Again: conditional or marginal estimates? 
				# Mean, Median, or Mode?
			if (adj == T) SMSY <- derived_obj$deripost_summary$SMSY_adj$Median else 
				SMSY <- derived_obj$deripost_summary$SMSY$Median 
		  
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
				# Note the < instead of the <= in the above if statement for LSM assumption
			if (bias.cor == TRUE) {
				sREP <- exp(rnorm(length(SREP), log(SREP) - 0.5*SREP_logSE^2, SREP_logSE))
				if(min(sREP) < 0) sREP <- exp(rnorm(length(SREP), log(SREP) - 0.5*SREP_logSE^2, SREP_logSE))
			} else {
				sREP <- exp(rnorm(length(SREP), log(SREP), SREP_logSE))
				if(min(sREP) < 0)   sREP <- exp(rnorm(length(SREP), SREP, SREP_SE))
			}
			# sREP <- exp(rnorm(length(SREP), log(SREP), SREP_logSE))
			# if(min(sREP)<0)   sREP <- exp(rnorm(length(SREP), SREP, SREP_SE$SE))
		  
			# Do SGEN calcs with new variables
			if (MCMC == F) SGENcalcs <- purrr::map2_dfr (exp(median(lnalpha_Parkin$loga)), sREP, Sgen.fn2) else
				SGENcalcs <- purrr::map2_dfr (exp(median(lnalphamc)), sREP, Sgen.fn2)
			# SGENcalcs_e <- purrr::map2_dfr (exp(median(loga_e)), SREP_e, Sgen.fn2) # Explicit
		}
    
		# 3. RUN RECONSTRUCTION
		if(prod == "RunReconstruction"){ # Requires Inlet devisions per stocks
			lnalpha_inlet <- read.csv(here::here("DataIn/CUPars_wBC.csv")) %>% 
				select(alpha,stkName) %>% 
				rename(inlets=stkName, lnalpha=alpha)
			lnalpha_nBC_inlet <- read.csv(here::here("DataIn/CUPars_nBC.csv")) %>% 
				select(alpha,stkName) %>% 
				rename(inlets=stkName, lnalpha_nBC=alpha)
			targetstocks <- read.csv(here::here("DataIn/WCVIStocks.csv")) %>% # Previously WCVIStocks - Changed to: Parken_evalstocks
				# filter (Stock != "Cypre") %>%
				rename(inlets=Inlet)
			Ric.A <- lnalpha_inlet %>% 
				left_join(targetstocks, by="inlets") %>% 
				select(c(lnalpha,inlets,CU,Stock))
			
			inSREP <- inSREP %>% 
				left_join(Ric.A, by = join_by(Stock, CU, inlets)) %>% 
				mutate(a.RR = exp(lnalpha))
			inSREP[inSREP$Stock=="Nitinat",]$a.RR <- exp(1)
			inSREP[inSREP$Stock=="San Juan",]$a.RR <- exp(1)
			inSREP[inSREP$Stock=="Nitinat",]$a.RR <- exp(1)
			
			# Below needs to create a new row for each inlet AS A STOCK and insert a a.RR value
			
			inSREP[inSREP$Stock=="Barkley",]$a.RR <- inSREP[inSREP$inlets=="Barkley",]$a.RR[1]
			inSREP <- inSREP |> 
				add_row(Stock = "Barkley", 
					inlets = "Barkley", 
					a.RR = inSREP[inSREP$inlets=="Barkley",]$a.RR[1])
			
			inSREP[inSREP$Stock=="Clayoquot",]$a.RR <- inSREP[inSREP$inlets=="Clayoquot",]$a.RR[1]
			inSREP <- inSREP |> 
				add_row(Stock = "Clayoquot", 
					inlets = "Clayoquot", 
					a.RR = inSREP[inSREP$inlets=="Clayoquot",]$a.RR[1])

			inSREP[inSREP$Stock=="Kyuquot",]$a.RR <- inSREP[inSREP$inlets=="Kyuquot",]$a.RR[1]
			inSREP <- inSREP |> 
				add_row(Stock = "Kyuquot", 
					inlets = "Kyuquot", 
					a.RR = inSREP[inSREP$inlets=="Kyuquot",]$a.RR[1])

			inSREP[inSREP$Stock=="Nootka/Esperanza",]$a.RR <- inSREP[inSREP$inlets=="Nootka/Esperanza",]$a.RR[1]
			inSREP <- inSREP |> 
				add_row(Stock = "Nootka/Esperanza", 
					inlets = "Nootka/Esperanza", 
					a.RR = inSREP[inSREP$inlets=="Nootka/Esperanza",]$a.RR[1])

			inSREP[inSREP$Stock=="Quatsino",]$a.RR <- inSREP[inSREP$inlets=="Quatsino",]$a.RR[1]
			inSREP <- inSREP |> 
				add_row(Stock = "Quatsino", 
					inlets = "Quatsino", 
					a.RR = inSREP[inSREP$inlets=="Quatsino",]$a.RR[1])

			inSREP[inSREP$Stock=="WCVI South",]$a.RR <- inSREP[inSREP$inlets=="Barkley",]$a.RR[1]
			inSREP <- inSREP |> 
				add_row(Stock = "WCVI South", 
					inlets = "Barkley", 
					a.RR = inSREP[inSREP$inlets=="Barkley",]$a.RR[1])

			inSREP[inSREP$Stock=="WCVI Nootka & Kyuquot",]$a.RR <- inSREP[inSREP$inlets=="Nootka/Esperanza",]$a.RR[1]
			inSREP <- inSREP |> 
				add_row(Stock = "WCVI Nootka & Kyuquot", 
					inlets = "Nootka/Esperanza", 
					a.RR = inSREP[inSREP$inlets=="Nootka/Esperanza",]$a.RR[1])

			inSREP[inSREP$Stock=="WCVI North",]$a.RR <- inSREP[inSREP$inlets=="Quatsino",]$a.RR[1]
			inSREP <- inSREP |> 
				add_row(Stock = "WCVI North", 
					inlets = "Quatsino", 
					a.RR = inSREP[inSREP$inlets=="Quatsino",]$a.RR[1])
			
			inSREP <- inSREP %>% select(-c(inlets, CU, lnalpha)) %>% rename(a.par=a.RR)
			
			Sig.Ric.A <- 0.51 # 0.255 for a narrower plausible bound
			
			Ric.A.hi <- exp(rnorm(length(inSREP$Median), log(inSREP$a.par), Sig.Ric.A))
			if(min(Ric.A.hi)<0) Ric.A <- exp(rnorm(length(inSREP$Median), inSREP$a.RR, Sig.Ric.A))
			
			sREP <- exp(rnorm(length(inSREP$Median), log(inSREP$Median), SREP_logSE))
			if(min(sREP)<0)   sREP <- exp(rnorm(length(inSREP$Median), inSREP$Median, 
												SREP_SE))
			
			SGENcalcs <- purrr::map2_dfr (Ric.A.hi, sREP, Sgen.fn2)

		} 
	
		# Timer tick
		setTxtProgressBar(pb, k)
		current_time <- Sys.time()
		elapsed_time <- difftime(current_time, start_time, units = "secs")
		cat(sprintf("\rElapsed time: %s seconds", round(as.numeric(elapsed_time), 2)))
		
		out <- list(bench = select(SGENcalcs, -apar, -bpar))
		outBench[[k]] <- out$bench
		
		# oute <- list(bench = select(SGENcalcs_e, -apar, -bpar))
		# outBenche[[k]] <- oute$bench
		
		if (prod == "Parken" & MCMC == F) outA <- list(alpha = exp(median(lnalpha_Parkin$loga))) # else
		if (prod == "Parken" & MCMC ==T) outA <- list(alpha = exp(median(lnalphamc)))
		if (prod == "LifeStageModel") outA <- list(alpha = Ric.A)
		if (prod == "RunReconstruction") outA <- list(alpha = Ric.A.hi)
		if (prod == "Parken") outAlpha <- outA
		if (prod == "LifeStageModel") outAlpha[[k]] <- outA
		if (prod == "RunReconstruction") outAlpha[[k]] <- outA
	
	}
  
	# Timer end
	close(pb)
	end_time <- Sys.time()
	total_elapsed_time <- end_time - start_time
	cat("\nTotal time taken:", total_elapsed_time, "\n")

	# 3. Outputs: Compile bootstrapped estimates of Sgen, SMSY, and SREP, 5th and 95th percentiles  
	stockNames <- WAin %>% 
		# filter(Stock != "Cypre") %>% # CYPRE FLAG ## ## ## 
		pull(Stock)
	stockNames <- unique(stockNames)
  
	SGEN.bs <- select(as.data.frame(outBench), starts_with("SGEN"))
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

	SMAX.bs <- select(as.data.frame(outBench), starts_with("SMAX"))
	rownames(SMAX.bs) <- stockNames
	SMAX.boot <- data.frame(SMAX= apply(SMAX.bs, 1, quantile, 0.5), 
                          lwr=apply(SMAX.bs, 1, quantile, 0.025),
                          upr=apply(SMAX.bs, 1, quantile, 0.975) )
  
	if (prod == "LifeStageModel" | prod == "RunReconstruction") {
		APAR.bs <- select(as.data.frame(outAlpha), starts_with("alpha"))
		rownames(APAR.bs) <- stockNames
		APAR.boot <- data.frame(APAR = apply(APAR.bs, 1, quantile, 0.5), 
                            lwr = apply(APAR.bs, 1, quantile, 0.025),
                            upr = apply(APAR.bs, 1, quantile, 0.975) )
	}
  
	boot <- list(SGEN.boot=SGEN.boot, SMSY.boot=SMSY.boot, 
                SREP.boot=SREP.boot, SMAX.boot=SMAX.boot) # , APAR.boot=APAR.boot)
	if (prod == "LifeStageModel" | prod == "RunReconstruction") boot$APAR.boot <- APAR.boot
  
	df1 <- data.frame(boot[["SGEN.boot"]], Stock=rownames(boot[["SGEN.boot"]]), RP="SGEN") 
	df1 <- df1 %>% rename(Value=SGEN)
	df2 <- data.frame(boot[["SREP.boot"]], Stock=rownames(boot[["SREP.boot"]]), RP="SREP")
	df2 <- df2 %>% rename(Value=SREP)
	df3 <- data.frame(boot[["SMSY.boot"]], Stock=rownames(boot[["SMSY.boot"]]), RP="SMSY")
	df3 <- df3 %>% rename(Value=SMSY)
	df4 <- data.frame(boot[["SMAX.boot"]], Stock=rownames(boot[["SMAX.boot"]]), RP="SMAX") 
	df4 <- df4 %>% rename(Value=SMAX)
	if (prod == "LifeStageModel" | prod == "RunReconstruction") {
		df5 <- data.frame(boot[["APAR.boot"]], Stock=rownames(boot[["APAR.boot"]]), RP="APAR")
		df5 <- df5 %>% rename(Value=APAR)
	}
  
	dfout <- add_row(df1, df2)
	dfout <- add_row(dfout, df3)
	dfout <- add_row(dfout, df4)
	if (prod == "LifeStageModel" | prod == "RunReconstruction") dfout <- add_row(dfout, df5)
	rownames(dfout) <- NULL
  
	if (round == TRUE) {
		dfout <- dfout %>% mutate(Value=signif(Value, 2)) %>% # Rounded to 2 signif digits
			mutate(lwr=signif(lwr,2)) %>% 
		    mutate (upr=signif(upr,2))
	}
	
	wasample <- WAin %>% 
        select("Stock", "WA", "lh") %>% 
        mutate(WA = round(WA, 0))
  
	BS.dfout <- merge(dfout, wasample, by="Stock", all.x=TRUE, sort=FALSE)
	if (prod == "Parken") alphaout <- outAlpha$alpha
	}; beep(2)

bsout <- list(BS.dfout = BS.dfout)

if (prod == "Parken") {
  bsout$alphaout <- alphaout
}

return(bsout)
}