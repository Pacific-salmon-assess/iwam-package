# The following is for getting out derived posterior samples from RTMB

# Library calls
library(RTMB)
library(MCMCglmm) # For posterior.mode()
library(beepr)

# Function start
derived_post <- function(x) {
  post <- as.matrix(x) # Where x is the stan fit object
  
  # mcmc_dens(post, regex_pars = "b0")
  
  # MATRIX EXTRACTION of full samples
  n_rows <- dim(x)[1] * dim(x)[2] # dim(fitstan) - 4500 post-warmup draws * 4 chains = 18000
  n_cols <- sapply(obj$report(post), function(x) length(x)) # Number of columns/stocks PER parameter
  n_matrices <- length(n_cols)
  
  matrices <- lapply(seq_len(n_matrices), function(k) {
    matrix(NA, nrow = n_rows, ncol = n_cols[k])
  })
  names(matrices) <- names(obj$report(post[1,]))
  
  
  # Current version
  pb <- txtProgressBar(min = 0, max = n_rows, style = 3, title = "Matrix Extraction")
  start_time <- Sys.time()
  
  # res = NULL
  for (i in seq_len(n_rows)) {  
    report_i <- obj$report(post[i,])
    for (k in seq_len(n_matrices)) { 
      matrices[[k]][i, seq_len(n_cols[k])] <- report_i[[k]]
      # This should be correctly adapted for a different number of columns per parameters
    }
    setTxtProgressBar(pb, i)
    
    # Calculate and display elapsed time dynamically
    current_time <- Sys.time()
    elapsed_time <- difftime(current_time, start_time, units = "secs")
    
    # Print running timer
    cat(sprintf("\rElapsed time: %s seconds", round(as.numeric(elapsed_time), 2)))
  }
  
  close(pb)
  # End timing
  end_time <- Sys.time()
  # Final elapsed time
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
  
  matrices$logE_tar_adj <- apply(matrices$logE_tar, 2, FUN = function(x)rnorm(length(x), x, sd = matrices$logESD))
  matrices$logAlpha_tar_adj <- apply(matrices$logAlpha_tar, 2, FUN = function(x)rnorm(length(x), x, sd = matrices$logAlphaSD))

  # Now we have logE and logAlpha for our targets
    # These need to be transformed - and ALSO need to re-calculate 
    # BETA, SMSY, and SGEN with these new values
    # Would this be done on EVERY iteration for EACH chain?
  # e.g. 
  matrices$E_tar_adj <- exp(matrices$logE_tar_adj) # for transformed E
  matrices$BETA_adj <- matrices$logAlpha_tar_adj / matrices$E_tar_adj
  matrices$SMSY_adj <- (1 - LambertW0(exp(1 - matrices$logAlpha_tar_adj))) / matrices$BETA_adj
  matrices$SGEN_adj <- -1/ matrices$BETA_adj * 
    LambertW0(- matrices$BETA_adj * matrices$SMSY_adj / (exp(matrices$logAlpha_tar_adj)))

  # Re-calculate lengths for added parameters
  n_matrices <- length(matrices)
  n_cols <- sapply(matrices, function(x) dim(x)) # Now with 2 dimensions
  
  # DATAFRAME CREATION of summarized derived posteriors
  dataframes <- lapply(seq_len(n_matrices), function(k) {
    data.frame(
      Stock = character(n_cols[2,k]),
      Mean = numeric(n_cols[2,k]),
      Median = numeric(n_cols[2,k]),
      PosteriorMode = numeric(n_cols[2,k]), # New addition of posterior Mode
      LQ_5 = numeric(n_cols[2,k]),
      UQ_95 = numeric(n_cols[2,k])
    )
  })
  
  names(dataframes) <- names(matrices)
  
  for (i in 1:n_matrices) {
    dataframes[[i]]$Stock <- c(1:ncol(matrices[[i]]))
    dataframes[[i]]$Mean <- apply(matrices[[i]], 2 , mean)
    dataframes[[i]]$Median <- apply(matrices[[i]], 2 , median)
    dataframes[[i]]$PosteriorMode <- apply(matrices[[i]], 2, function(col) {posterior.mode(as.mcmc(col))}) 
      # New addition of posterior Mode
    dataframes[[i]]$LQ_5 <- apply(matrices[[i]], 2, quantile, probs = c(0.05, 0.95))[1,] # 0.05
    dataframes[[i]]$UQ_95 <- apply(matrices[[i]], 2, quantile, probs = c(0.05, 0.95))[2,] # 0.95
  }

  return(list(deripost_summary = dataframes,
              deripost_full = matrices)
  )
  
} # ; beep(2)
