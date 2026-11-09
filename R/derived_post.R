# The following is for getting out derived posterior samples from RTMB ####

# Library calls
library(RTMB)
library(MCMCglmm) # For posterior.mode()
library(beepr)
library(HDInterval) # High density interval for posteriors with skew/multimodality

options(scipen = 999)

derived_post <- function(x, model) {
	post <- as.matrix(x)
	
	# MATRIX EXTRACTION of full samples
	n_rows <- dim(x)[1] * dim(x)[2] # dim(fitstan) - 4500 post-warmup draws * 4 chains = 18000
	n_cols <- sapply(obj$report(post), function(x) length(x)) # Number of columns/stocks PER parameter
	n_matrices <- length(n_cols)
	  
	matrices <- lapply(seq_len(n_matrices), function(k) {
		matrix(NA, nrow = n_rows, ncol = n_cols[k])
	})
	names(matrices) <- names(obj$report(post[1,]))
  
  	pb <- txtProgressBar(min = 0, max = n_rows, style = 3, title = "Matrix Extraction")
	start_time <- Sys.time()
  
	for (i in seq_len(n_rows)) {  
		report_i <- obj$report(post[i,])
		for (k in seq_len(n_matrices)) { 
			matrices[[k]][i, seq_len(n_cols[k])] <- report_i[[k]]
		}
		
		setTxtProgressBar(pb, i)
		current_time <- Sys.time()
		elapsed_time <- difftime(current_time, start_time, units = "secs")
   		cat(sprintf("\rElapsed time: %s seconds", round(as.numeric(elapsed_time), 2)))
	}
  
	close(pb)
	end_time <- Sys.time()
	total_elapsed_time <- end_time - start_time
	cat("\nTotal time taken:", total_elapsed_time, "\n")
  
	# Random effects - per row of the matrix for logE_tar and logAlpha_tar
		# you would add in their respective rnorms
	#   > names(matrices)
	#  [1] "logE_re"       "E"             "logAlpha"      "alpha"         "SMSY_r"       
	#  [6] "BETA_r"        "tauobs"        "E_tar"         "logE_tar"      "logAlpha_tar" 
	# [11] "alpha_tar"     "BETA"          "SMSY"          "SGEN"          "E_line_stream"
	# [16] "E_line_ocean"  "logAlphaSD"    "logESD"  
		# so logAlpha_tar exists, and logE_tar exists HOWEVER _re object's don't exist
		# they come from the following distributions:
			# logAlpha_tar <- dnorm(logAlpha_re, mean = 0, sd = logAlphaSD, log = TRUE)
			# logE_tar <- dnorm(logE_re, mean = 0, sd = logESD, log = TRUE)
		# matrices$logE_re <- dnorm(...) # ? is this necessary?
  
	if (dat$prioronly == 0) { # Just to be able to run prioronly without these additional items
		# LambertW0 doesn't work in prioronly settings
		if (model == 'SREP'){
		# monte carlo integration - WRONG
		# This creates a prediction interval - including new site variability
			# If you take the mean of logE_tar_pp --> marginal mean - WRONG
		matrices$logSREP_tar_pp <- apply(matrices$logSREP_tar, 2, 
			FUN = function(x)rnorm(length(x), x, sd = matrices$logSREP_sd))

		matrices$logSREP_line_ocean_pp <- apply(matrices$logSREP_line_ocean, 2, 
			FUN = function(x)rnorm(length(x), x, sd = matrices$logSREP_sd))
		
		matrices$logSREP_line_stream_pp <- apply(matrices$logSREP_line_stream, 2, 
			FUN = function(x)rnorm(length(x), x, sd = matrices$logSREP_sd))

		# mean(exp(matrices$logE_tar_pp[,1])) # Should be the same as below - larger
		# exp(mean(matrices$logE_tar_pp[,1] + matrices$logESD^2/2)) # Should be the same as above - smaller
		# NOT CURRENTLY THE SAME - SOMETHING IS WRONG

		matrices$logAlpha_tar_pp <- apply(matrices$logAlpha_tar, 2, 
			FUN = function(x)rnorm(length(x), x, sd = matrices$Alpha_sd))

		# Instead of doing this false conditional estimate
		# This is where we should calculate the marginal mean estimate with bias correction subtracted
		# E.g.:
		# matrices$logE_tarbias <- matrices$logE_tar - (0.5 * (sqrt(1/matrices$tauobs))^2)
		# matrices$logAlpha_tarbias <- matrices$logAlpha_tar - (0.5 * ()^2)

		# Now we have logE and logAlpha for our targets
			# These need to be transformed - and ALSO need to re-calculate 
			# BETA, SMSY, and SGEN with these new values
			# Would this be done on EVERY iteration for EACH chain?
		# e.g.
		matrices$Alpha_tar_pp <- exp(matrices$logAlpha_tar_pp) # matrices$logAlpha_tar_pp is actually loglog - needs to be exp
	  
		matrices$SREP_tar_pp <- exp(matrices$logSREP_tar_pp) # for transformed E
		matrices$BETA_pp <- matrices$Alpha_tar_pp / matrices$SREP_tar_pp
		matrices$SMSY_pp <- (1 - LambertW0(exp(1 - matrices$Alpha_tar_pp))) / matrices$BETA_pp
		matrices$SGEN_pp <- -1/ matrices$BETA_pp * 
			LambertW0(- matrices$BETA_pp * matrices$SMSY_pp / (exp(matrices$Alpha_tar_pp)))

		matrices$SREP_line_ocean_pp <- exp(matrices$logSREP_line_ocean_pp) 
		matrices$SREP_line_stream_pp <- exp(matrices$logSREP_line_stream_pp)
		}
	  
		if (model == 'SMAX'){
			matrices$logSMAX_tar_pp <- apply(matrices$logSMAX_tar, 2, 
				FUN = function(x)rnorm(length(x), x, sd = matrices$logSMAX_sd))

			matrices$logSMAX_line_ocean_pp <- apply(matrices$logSMAX_line_ocean, 2, 
				FUN = function(x)rnorm(length(x), x, sd = matrices$logSMAX_sd))
			
			matrices$logSMAX_line_stream_pp <- apply(matrices$logSMAX_line_stream, 2, 
				FUN = function(x)rnorm(length(x), x, sd = matrices$logSMAX_sd))

			# Does this calculate the marginal posterior OR the posterior predictive distribution?
			matrices$logAlpha_tar_pp <- apply(matrices$logAlpha_tar, 2, 
				FUN = function(x)rnorm(length(x), x, sd = matrices$Alpha_sd))

			matrices$Alpha_tar_pp <- exp(matrices$logAlpha_tar_pp)

			matrices$SMAX_tar_pp <- exp(matrices$logSMAX_tar_pp) 
			matrices$BETA_pp <- 1 / matrices$SMAX_tar_pp
			matrices$SREP_pp <- matrices$Alpha_tar_pp/matrices$BETA_pp
			matrices$SMSY_pp <- (1 - LambertW0(exp(1 - matrices$Alpha_tar_pp))) / matrices$BETA_pp
			matrices$SGEN_pp <- -1/ matrices$BETA_pp * 
				LambertW0(- matrices$BETA_pp * matrices$SMSY_pp / (exp(matrices$Alpha_tar_pp)))

			matrices$SMAX_line_ocean_pp <- exp(matrices$logSMAX_line_ocean_pp) 
			matrices$SMAX_line_stream_pp <- exp(matrices$logSMAX_line_stream_pp)
		}
	}

	# # Prior predictions and posterior predictions - adding in observation error
	# matrices$logRS_pred_obs <- apply(matrices$logRS_pred, 2,
	#   FUN = function(x) rnorm(length(x), x, sd = sqrt(1/matrices$tauobs)))
	# # sqrt(1/tauobs[stk[i]])
  
	# Re-calculate lengths for added parameters
	n_matrices <- length(matrices)
	n_cols <- sapply(matrices, function(x) dim(x)) # Now with 2 dimensions
  
	# DATAFRAME CREATION of summarized derived posteriors
	dataframes <- lapply(seq_len(n_matrices), function(k) {
		data.frame(
			Stock = character(n_cols[2,k]),
			Mean = numeric(n_cols[2,k]),
			Median = numeric(n_cols[2,k]),
			PosteriorMode = numeric(n_cols[2,k]),
			LQ_025 = numeric(n_cols[2,k]),
			UQ_975 = numeric(n_cols[2,k])
		)
	})
  
	names(dataframes) <- names(matrices)
  
	for (i in 1:n_matrices) {
		hdi_list <- apply(matrices[[i]], 2, function(x) hdi(x, credMass = 0.9)) # 0.9 = 90% HDI
		dataframes[[i]]$Stock <- c(1:ncol(matrices[[i]]))
		dataframes[[i]]$Mean <- apply(matrices[[i]], 2 , mean)
		dataframes[[i]]$Median <- apply(matrices[[i]], 2 , median)
		dataframes[[i]]$PosteriorMode <- apply(matrices[[i]], 2, 
			function(col) {posterior.mode(as.mcmc(col))}) 
		dataframes[[i]]$LQ_025 <- apply(matrices[[i]], 2, quantile, probs = c(0.025, 0.975))[1,] # 0.05
		dataframes[[i]]$UQ_975 <- apply(matrices[[i]], 2, quantile, probs = c(0.025, 0.975))[2,] # 0.95
		# Consider adding HDPI intervals here
			# Using and apply function with hdi(x, credMass = 0.89)
		# dataframes[[i]]$HDIlwr <- sapply(hdi_list, `[`, 1) # lower bound
		# dataframes[[i]]#HDIupr <- sapply(hdi_list, `[`, 2) # upper bound
	}

	return(list(deripost_summary = dataframes, # consider derivpostsum and derivpostfull or just dpostsum and dpostfull
            deripost_full = matrices)
	)
  
}
