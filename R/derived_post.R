# The following is for getting out derived posterior samples from RTMB

# Library calls
library(RTMB)
library(MCMCglmm) # For posterior.mode()
library(beepr)

# Function start
derived_post <- function(x) {
  post <- as.matrix(x) # Where x is the stan fit object
  
  # MATRIX EXTRACTION of full samples
  n_rows <- 18000
  n_cols <- sapply(obj$report(post), function(x) length(x))
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
  
  # DATAFRAME CREATION of summarized derived posteriors
  dataframes <- lapply(seq_len(n_matrices), function(k) {
    data.frame(
      Stock = character(n_cols[k]),
      Mean = numeric(n_cols[k]),
      Median = numeric(n_cols[k]),
      PosteriorMode = numeric(n_cols[k]), # New addition of posterior Mode
      LQ_5 = numeric(n_cols[k]),
      UQ_95 = numeric(n_cols[k])
    )
  })
  
  names(dataframes) <- names(matrices)
  
  for (i in 1:n_matrices) {
    dataframes[[i]]$Stock <- c(1:ncol(matrices[[i]]))
    dataframes[[i]]$Mean <- apply(matrices[[i]], 2 , mean)
    dataframes[[i]]$Median <- apply(matrices[[i]], 2 , median)
    dataframes[[i]]$PosteriorMode <- apply(matrices[[i]], 2, function(col) {posterior.mode(as.mcmc(col))}) # New addition of posterior Mode
    dataframes[[i]]$LQ_5 <- apply(matrices[[i]], 2, quantile, probs = c(0.05, 0.95))[1,] # 0.05
    dataframes[[i]]$UQ_95 <- apply(matrices[[i]], 2, quantile, probs = c(0.05, 0.95))[2,] # 0.95
  }

  return(list(deripost_summary = dataframes,
              deripost_full = matrices)
  )
  beep(2)
}
