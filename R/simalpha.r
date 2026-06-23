# BOOTSTRAPPING: Simulating IWAM Model Posteriors under new Prior Assumptions ####

#### Libaries ####
# library(RTMB)
# library(ggplot2)
library(dplyr)
library(tidyverse)
library(progress)
# library(tmbstan)
library(TMB)
# library(tidybayes)
# library(bayesplot)
library(beepr) # Sounds
# library(viridis)
# library(latex2exp)
# library(mvtnorm)
library(MASS)

source(here::here("R/LambertWs.R")) # Lambert W function
source(here::here("R/helperFunctions.R")) # For bootstrapping

options(scipen = 999)

simalpha <- function(bsiters = 1000, # Simulation iterations
					newalpha = c(1, 0.3), # LifeStageModel assumed alpha prior
					prior_rho = c(-0.4), # default rho value for covariance matrix
					roundfinal = FALSE, # Round final reported values
					WAinname = c("DataIn/Parken_evalstocks.csv") # Comparison stocks
					) {
						
	WAin <- read.csv(here::here(WAinname))	
	outBench <- list()
	outAlpha <- list()
	outBpar <- list()
  
    pb <- txtProgressBar(min = 1, max = bsiters, style = 3, title = "Bootstrapping")
    start_time <- Sys.time()
  
	poplen <- ncol(dsmax$deripost_full$Alpha_tar_adj) # 25 stock predictions to be made
	
	logSMAX_sd <- dsmax$deripost_summary$logSMAX_sd$Mean
	
	b0 <- dsmax$deripost_full$b0 # stream, +ocean
	bWA <- dsmax$deripost_full$bWA
	Alpha0 <- exp(dsmax$deripost_full$Alpha0)
	
	post <- cbind(Alpha0, b0[,1], bWA[,1]) # Full posteriors for all
	post2 <- cbind(exp(dsmax$deripost_full$Alpha0 + dsmax$deripost_full$Alpha02), b0[,1] + b0[,2], bWA[,1] + bWA[,2])
	# Then build a covariance matrix based on the above matrix 
	covmatrix <- cov(post)
	covmatrix2 <- cov(post2)
	mu <- apply(post, 2, mean)
	mu2 <- apply(post2, 2, mean)

	# cov is sqrt(cor)*variable
	# correlation is: covmatrix[1,2]/(sqrt(covmatrix[1,1]*covmatrix[2,2]))
	# if you have correlation just solve for covaraince
	
	# e.g. cov2cor() to find previous correlation matrix
		# NOTE: for SLOPES (bWA) - stream is very small negative and ocean is very small positive
			# therefore the slopes are going in different directions
			# as alpha goes down - stream slope goes up (Steeper)
			# as alpha goes down - ocean slope goes down (Flatter)
			# which effectively draws the lines together
	# covmatrix[1,2] <- covmatrix[2,1] <- prior_rho * sqrt(covmatrix[1,1])*sqrt(covmatrix[2,2]) 
	# covmatrix2[1,2] <- covmatrix2[2,1] <- prior_rho * sqrt(covmatrix2[1,1])*sqrt(covmatrix2[2,2])

    for (k in 1:bsiters){
	## From https://statproofbook.github.io/P/mvn-cond.html
		logalphak <- rnorm(1, newalpha[1], newalpha[2])
	  	
		# stream type (base case)
		mu_new <- mu[-1] + covmatrix[-1,1]/covmatrix[1,1]*(logalphak - mu[1])
		var_new <- covmatrix[-1,-1] - (covmatrix[-1,1, drop = FALSE] / covmatrix[1,1]) %*% covmatrix[1,-1, drop = FALSE]
		bnew <- MASS::mvrnorm(1, mu = mu_new, var_new)
		
		# ocean type (+ additive case)
		mu2_new <- mu2[-1] + covmatrix2[-1,1]/covmatrix2[1,1]*(logalphak - mu2[1])
		var2_new <- covmatrix2[-1,-1] - (covmatrix2[-1,1, drop = FALSE] / covmatrix2[1,1]) %*% covmatrix2[1,-1, drop = FALSE]
		b2new <- MASS::mvrnorm(1, mu = mu2_new, var2_new)
		
		newlogSMAX <- c()
		type_tar <- as.numeric(WAin$lh)
		logWAshifted <- log(WAin$WA) - mean(WAbase$logWA)
		# calculate newlogSMAX with 2 different lines depending on lifehistory type:
		newlogSMAX <- ifelse(WAin$lh == 0, bnew[1] + bnew[2] * logWAshifted + rnorm(poplen, 0, sd = logSMAX_sd),
											b2new[1] + b2new[2] * logWAshifted + rnorm(poplen, 0, sd = logSMAX_sd))

		# newlogSMAX <- bnew[1] + bnew[2] * logWAshifted + rnorm(poplen, 0, sd = logSMAX_sd)
		# newlogSMAX2 <- b2new[1] + b2new[2] * logWAshifted + rnorm(poplen, 0, sd = logSMAX_sd)
		
		SGENcalcs <- purrr::map2_dfr (exp(rep(logalphak, poplen)), exp(newlogSMAX), Sgen.fn4)
		# SGENcalco <- purrr::map2_dfr (exp(rep(logalphak, poplen)), exp(newlogSMAX2), Sgen.fn4)
		
		out <- list(bench = dplyr::select(SGENcalcs, -apar, -bpar),
					bpar = dplyr::select(SGENcalcs, bpar))
		outBench[[k]] <- out$bench
		outBpar[[k]] <- out$bpar
		outA <- list(alpha = exp(rep(logalphak, poplen)))
		outAlpha[[k]] <- outA
		
		# Timer tick
		setTxtProgressBar(pb, k)
		current_time <- Sys.time()
		elapsed_time <- difftime(current_time, start_time, units = "secs")
		cat(sprintf("\rElapsed time: %s seconds", round(as.numeric(elapsed_time), 2)))
	}
  
	# Timer end
	close(pb)
	end_time <- Sys.time()
	total_elapsed_time <- end_time - start_time
	cat("\nTotal time taken:", total_elapsed_time, "\n")

	# 3. Outputs: Compile bootstrapped estimates of Sgen, SMSY, and SREP, 5th and 95th percentiles  
	stockNames <- WAin %>% 
		pull(Stock)
	stockNames <- unique(stockNames)
  
	SGEN.bs <- dplyr::select(as.data.frame(outBench), starts_with("SGEN"))
	rownames(SGEN.bs) <- stockNames
	SGEN.boot <- data.frame(SGEN = apply(SGEN.bs, 1, quantile, 0.5), 
                          lwr = apply(SGEN.bs, 1, quantile, 0.025),
                          upr = apply(SGEN.bs, 1, quantile, 0.975) )
  
	SMSY.bs <- dplyr::select(as.data.frame(outBench), starts_with("SMSY"))
	rownames(SMSY.bs) <- stockNames
	SMSY.boot <- data.frame(SMSY = apply(SMSY.bs, 1, quantile, 0.5), 
                          lwr = apply(SMSY.bs, 1, quantile, 0.025),
                          upr = apply(SMSY.bs, 1, quantile, 0.975) )
  
	SREP.bs <- dplyr::select(as.data.frame(outBench), starts_with("SREP"))
	rownames(SREP.bs) <- stockNames
	SREP.boot <- data.frame(SREP = apply(SREP.bs, 1, quantile, 0.5), 
                          lwr = apply(SREP.bs, 1, quantile, 0.025),
                          upr = apply(SREP.bs, 1, quantile, 0.975) )

	SMAX.bs <- dplyr::select(as.data.frame(outBench), starts_with("SMAX"))
	rownames(SMAX.bs) <- stockNames
	SMAX.boot <- data.frame(SMAX = apply(SMAX.bs, 1, quantile, 0.5), 
                          lwr = apply(SMAX.bs, 1, quantile, 0.025),
                          upr = apply(SMAX.bs, 1, quantile, 0.975) )
  

	APAR.bs <- dplyr::select(as.data.frame(outAlpha), starts_with("alpha"))
	rownames(APAR.bs) <- stockNames
	APAR.boot <- data.frame(APAR = apply(APAR.bs, 1, quantile, 0.5), 
						lwr = apply(APAR.bs, 1, quantile, 0.025),
						upr = apply(APAR.bs, 1, quantile, 0.975) )


	boot <- list(SGEN.boot = SGEN.boot, SMSY.boot = SMSY.boot, 
                SREP.boot = SREP.boot, SMAX.boot = SMAX.boot)
	boot$APAR.boot <- APAR.boot
  
	df1 <- data.frame(boot[["SGEN.boot"]], Stock = rownames(boot[["SGEN.boot"]]), RP = "SGEN") 
	df1 <- df1 %>% rename(Value = SGEN)
	df2 <- data.frame(boot[["SREP.boot"]], Stock = rownames(boot[["SREP.boot"]]), RP = "SREP")
	df2 <- df2 %>% rename(Value = SREP)
	df3 <- data.frame(boot[["SMSY.boot"]], Stock = rownames(boot[["SMSY.boot"]]), RP = "SMSY")
	df3 <- df3 %>% rename(Value = SMSY)
	df4 <- data.frame(boot[["SMAX.boot"]], Stock = rownames(boot[["SMAX.boot"]]), RP = "SMAX") 
	df4 <- df4 %>% rename(Value = SMAX)

	df5 <- data.frame(boot[["APAR.boot"]], Stock = rownames(boot[["APAR.boot"]]), RP = "APAR")
	df5 <- df5 %>% rename(Value = APAR)

  
	dfout <- add_row(df1, df2)
	dfout <- add_row(dfout, df3)
	dfout <- add_row(dfout, df4)
	dfout <- add_row(dfout, df5)
	rownames(dfout) <- NULL
  
	if (roundfinal == TRUE) {
		dfout <- dfout %>% dplyr::mutate(Value = signif(Value, 2)) %>% # Rounded to 2 signif digits
			mutate(lwr = signif(lwr,2)) %>% 
		    mutate (upr = signif(upr,2))
	}
	
	wasample <- WAin %>% 
        dplyr::select("Stock", "WA", "lh") %>% 
        dplyr::mutate(WA = round(WA, 0))
  
	BS.dfout <- merge(dfout, wasample, by="Stock", all.x = TRUE, sort = FALSE)

	bsout <- list(BS.dfout = BS.dfout, bpar = outBpar)
	return(bsout)
}