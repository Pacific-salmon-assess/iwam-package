# The following is for getting out derived posterior samples from RTMB

# Library calls
library(RTMB)

# Function start
derived_post <- function(x) {
  post <- as.matrix(x) # Where x is the stan fit object
  
  
  # Old version
  # pb <- txtProgressBar(min = 0, max = n_rows, style = 3, title = "Matrix Extraction")
  # for (i in 1:n_rows) {  # Outer loop over rows of 'post'
  #   report_i <- obj$report(post[i,])  # Extract report once per row (18000)
  #   for (k in 1:n_matrices) {  # Loop over matrices (7)
  #     for (j in 1:n_cols) {  # Loop over columns (25)
  #       matrices[[k]][i, j] <- report_i[[k]][[j]]  # Fill pre-allocated matrices
  #     }
  #   }
  #   setTxtProgressBar(pb, i)
  # }
  # close(pb)
  
  # Current version
  start_time <- Sys.time()
  
  fn1  = function(){
    res = NULL
    
  for (i in 1:100) {  # seq_len(n_rows)
    report_i <- obj$report(post[i,])  
    res <- rbind(res, unlist(report_i))
  }
}

  library(microbenchmark)  
  microbenchmark(tor = fn2(), 
                 paul = fn1(), times = 20)
  
  
  # pb <- txtProgressBar(min = 0, max = n_rows, style = 3, title = "Matrix Extraction")
  start_time <- Sys.time()
  

    fn2 = function(){
      # MATRIX EXTRACTION of full samples
      n_rows <- 100
      n_cols <- sapply(obj$report(post), function(x) length(x))
      n_matrices <- 7
      
      matrices <- lapply(seq_len(n_matrices), function(k) {
        matrix(NA, nrow = n_rows, ncol = n_cols[k])
      })
      names(matrices) <- names(obj$report(post[1,]))
      
      for (i in 1:100) { # seq_len(n_rows)
        
      report_i <- obj$report(post[i,])
    # res <- rbind(res, unlist(report_i))
    # res[,grepl("SGEN", colnames(res))]
    for (k in seq_len(n_matrices)) {
      matrices[[k]][i, seq_len(n_cols[k])] <- report_i[[k]]
      # This should be correctly adapted for a different number of columns per parameters
    }
      }
    }

    
    
    
    
      # close(pb)
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
      LQ_5 = numeric(n_cols[k]),
      UQ_95 = numeric(n_cols[k])
    )
  })
  names(dataframes) <- names(matrices)
  
  ## PAUL ADDED
  namec <- c("SGEN", "logAlpha", "SMSY", "E", "E_tar", "BETA", "logAlpha_tar")
  namesres <- gsub("*[0-9]", "", colnames(res))
  res[, namesres == namec[2]]
  ans <- NULL
  for(namei in namec){
    ans <- rbind(ans, data.frame(Parameter = namei, stock = 1:sum(namesres == namei), Mean = apply(res[,namesres == namei], 2, mean)))
    # or could cbind
    # ans <- cbind()
  }
  ## END PAUL ADDED
  
  for (i in 1:n_matrices) {
    dataframes[[i]]$Stock <- c(1:ncol(matrices[[i]]))
    dataframes[[i]]$Mean <- apply(matrices[[i]], 2 , mean)
    dataframes[[i]]$Median <- apply(matrices[[i]], 2 , median)
    dataframes[[i]]$LQ_5 <- apply(matrices[[i]], 2, quantile, probs = c(0.05, 0.95))[1,] # 0.05
    dataframes[[i]]$UQ_95 <- apply(matrices[[i]], 2, quantile, probs = c(0.05, 0.95))[2,] # 0.95
  }
  
  return(list(deripost_summary = dataframes,
              deripost_full = matrices)
  )
  
}
