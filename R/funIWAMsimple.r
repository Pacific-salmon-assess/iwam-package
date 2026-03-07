# Liermann SMAX (B) RTMB Model with MCMC Sampling ####
# Begin IWAM model function: ####
	# The purpose of this function is a simple - no input - function for running IWAM
	# under the default conditions/priors/settings.
	# This is for speed testing.
	# For a complete function see: funIWAM
	# For a function designed for single WA see: funIWAMsingle
	

# Libaries ####
library(RTMB)
library(dplyr)
library(tidyverse) 
library(progress) 
library(tmbstan) 
library(beepr) 

here::i_am("R/LambertWs.R") # For non-RStudio functionality
source(here::here("R/LambertWs.R")) # Lambert W function
source(here::here("R/helperFunctions.R")) # For bootstrapping

options(scipen = 999)

## Turn off byte compiler ####
	# This the only way to get RTMB functions to run within wrapper functions
compiler::enableJIT(0)

IWAMsimple <- function() {

	LambertW0 <- ADjoint(
		function(x){gsl::lambert_W0(x)},
		function(x, y, dy) {dy / (x + exp(y))}
	)

	# For predicting/re-evaluation of synoptic sets: WAbase
		# 1. Run model until setting up data section
		# 2. Then over-write WAin <- WAbase
		# 3. And rename logWAshifted to logWAshifted_t
		# 4. And make sure to change type_tar to fit the expected 0:1 indexing

	# Raw data read-in ####
	WAin <- c("DataIn/Parken_evalstocks.csv")
	# c("DataIn/WCVIStocks.csv") or c("DataIn/Parken_evalstocks.csv") 
	# WAin <- c("DataIn/Ordered_backcalculated_noagg.csv")

	srdatwna <- read.csv(here::here("DataIn/SRinputfile.csv")) # Consider _TK **** as alternative with longer dataset
	WAbase <- read.csv(here::here("DataIn/WatershedArea.csv"))
	WAin <- read.csv(here::here(WAin))

	srdatwna <- srdatwna %>% 
		filter(!Name %in% c("Hoko","Hoh")) # |>
		# (!(Name == "Skagit")) # |> # Skagit is #22
		# filter( !(Name == "Chehalis")) # |> # Chehalis is #19

	srdat <- srdatwna %>% 
		filter(Rec != "NA") %>%
		filter( !(Name == "Cowichan" & (Yr < 1985 | Yr == 1986 | Yr == 1987))) %>%
		group_by(Name, Stocknumber, Stream) %>%
		arrange(Yr) %>%
		mutate(yr_num = 0:(n()-1)) %>%
		ungroup() %>%
		arrange(Stocknumber) %>%
		mutate(lh = factor(ifelse(Stream == 0, "stream", "ocean"),
			levels = c("stream", "ocean"))) |> # Stream = 0, Ocean  = 1
		mutate(Stocknumber = as.integer(factor(Stocknumber)) - 1) # Re-numbering uniquely

	names <- srdat %>% 
		dplyr::select (Stocknumber, Name, lh) %>% 
		distinct()

	WAbase <- WAbase %>% 
		full_join(names, by="Name") %>% 
		arrange(Stocknumber) %>%
		mutate(logWA = log(WA)) |> 
		filter(!is.na(Stocknumber))

	mean_logWA <- mean(WAbase$logWA)
	WAbase$logWAshifted <- WAbase$logWA - mean_logWA

	WAin$logWA <- log(WAin$WA)
	WAin$logWAshifted_t <- WAin$logWA - mean_logWA

	lifehist <- srdat %>% dplyr::select(Stocknumber, Name, Stream) %>% 
		group_by(Stocknumber) %>% 
		summarize(lh=max(Stream))

	dat <- list(srdat = srdat,
				WAbase = WAbase,
				WAin = WAin,
				lifehist = lifehist,
				lineWA =  seq(min(WAbase$logWAshifted), 
							  max(WAbase$logWAshifted), 0.1),
				mean_logWA = mean_logWA,
				logRS = log(srdat$Rec) - log(srdat$Sp),
				prioronly = 0) # 0-run with data, 1-prior prediction mode

	if(dat$prioronly == 1){print("Prior Prediction Mode")} else 
		{print("Posterior Prediction Mode")}

	N_Stk <- max(srdat$Stocknumber + 1)
	lhdiston <- T # T = LH Specific
	bias.cor <- F # T = subtract bias correction terms from expontiated mean terms

	par <- list(b0 = c(10, 0), # Initial values for WA regression intercepts
				bWA = c(0, 0), # Inital values for WA regression slopes
				# logRS_pred = numeric(nrow(srdat)), # Zeroes - testing as a parameter
				logSMAX_re = numeric(N_Stk), # Zeroes 
				Alpha0 = 0.6,
				Alpha_re = numeric(nrow(dat$WAbase)), # Zeroes
				tauobs = 0.01 + numeric(N_Stk), # Constrained positive
				logSMAX_sd = 1,
				Alpha_sd = 1
	)
	if (lhdiston) {
	  par$Alpha02 <- 0
	}
	
	f_smax <- function(par){
	  getAll(dat, par)
	  
	  N_Stk = max(srdat$Stocknumber + 1) # number of stocks
	  stk = srdat$Stocknumber + 1 # vector of stocknumbers
	  N_Obs = nrow(srdat) # number of observations
	  N_Pred = nrow(WAin) # number of predicted watershed areas
	  
	  S = srdat$Sp
	  type = lifehist$lh
	  type_tar = as.numeric(WAin$lh) 

	  SMAX <- numeric(N_Stk)
	  logSMAX <- numeric(N_Stk)
	  logAlpha <- numeric(N_Stk) 
	  
	  logRS_pred <- numeric(N_Obs) # Does this still report if not a vector?

	  SMAX_tar <- numeric(N_Pred)
	  logSMAX_tar <- numeric(N_Pred) 
	  logAlpha_tar <- numeric(N_Pred)
	  
	  line <- length(lineWA)
	  logSMAX_line_stream <- numeric(line) 
	  SMAX_line_stream <- numeric(line) 
	  logSMAX_line_ocean <- numeric(line) 
	  SMAX_line_ocean <- numeric(line) 
	  
	  if (bias.cor) {
		biaslogSMAX <- -0.5*logSMAX_sd^2 # Global 
		biaslogAlpha <- -0.5*Alpha_sd^2 # Global 
		biaslogRS <- -0.5*(sqrt(1/tauobs))^2 # Stock-specific
	  } else {
		biaslogSMAX <- 0 
		biaslogAlpha <- 0 
		biaslogRS <- numeric(N_Stk)
	  }
	  
	  nll <- 0
	  
	  nll <- nll - sum(dnorm(b0[1], 10, sd = 31.6, log = TRUE)) # Prior b0
	  nll <- nll - sum(dnorm(b0[2], 0, sd = 31.6, log = TRUE)) # Prior b0
	  nll <- nll - sum(dnorm(bWA[1], 0, sd = 31.6, log = TRUE)) # Prior bWA
	  nll <- nll - sum(dnorm(bWA[2], 0, sd = 31.6, log = TRUE)) # Prior bWA
	  
	  nll <- nll - sum(dnorm(Alpha0, 0.6, sd = 0.45, log = TRUE)) # Prior (rM)
	  if(lhdiston) nll <- nll - sum(dnorm(Alpha02, 0, sd = 31.6, log = TRUE)) # Prior (rD)
	  
	  for (i in 1:N_Stk){
		nll <- nll - dnorm(logSMAX_re[i], 0, sd = 1, log = TRUE)
		nll <- nll - dnorm(Alpha_re[i], 0, sd = 1, log = TRUE)
		
		logSMAX[i] <- b0[1] + b0[2]*type[i] + (bWA[1] + bWA[2]*type[i]) * WAbase$logWAshifted[i] + logSMAX_re[i]*logSMAX_sd + biaslogSMAX
		SMAX[i] <- exp(logSMAX[i])
			
		if(lhdiston) logAlpha[i] <- Alpha0 + Alpha02*type[i] + Alpha_re[i]*Alpha_sd + biaslogAlpha
		else logAlpha[i] <- Alpha0 + Alpha_re[i]*Alpha_sd + biaslogAlpha

		nll <- nll - dgamma(tauobs[i], shape = 0.0001, scale = 1/0.0001, log = TRUE)
	  }

	  for (i in 1:N_Obs){
		Alpha_pred <- exp(logAlpha)
		logRS_pred[i] <- Alpha_pred[stk[i]] - S[i]/SMAX[stk[i]] + biaslogRS[stk[i]]

		if(!prioronly){ # If prioronly is 1, then likelihood is not calculated, if 0 then it is
		  nll <- nll - dnorm(logRS[i], logRS_pred[i], sd = sqrt(1/tauobs[stk[i]]), log = TRUE) 
		} 
	  }
	  
	  SMSY_r = numeric(nrow(WAbase))
	  BETA_r = numeric(nrow(WAbase))
	  SREP_r = numeric(nrow(WAbase))
	  
	  for (i in 1:N_Stk){
		BETA_r[i] <- 1/SMAX[i] 
		SMSY_r[i] <- (1 - LambertW0(exp(1 - Alpha_pred[i]))) / BETA_r[i] 
		SREP_r[i] <- Alpha_pred[i]/BETA_r[i] 
	  }

	  BETA = numeric(nrow(WAin))
	  SMSY = numeric(nrow(WAin))
	  SGEN = numeric(nrow(WAin))
	  SREP = numeric(nrow(WAin))

	  for (i in 1:N_Pred){
		
		if(lhdiston) logAlpha_tar[i] <- Alpha0 + Alpha02*type_tar[i] + biaslogAlpha # + biaslogAlpha + logAlpha_sd^2/2 
		else logAlpha_tar[i] <- Alpha0 + biaslogAlpha # + biaslogAlpha + logAlpha_sd^2/2 
		
		logSMAX_tar[i] <- b0[1] + b0[2]*type_tar[i] + (bWA[1] + bWA[2]*type_tar[i])*WAin$logWAshifted_t[i] + biaslogSMAX
		SMAX_tar[i] <- exp(logSMAX_tar[i]) 
		
		Alpha_tar <- exp(logAlpha_tar)
		BETA[i] <- 1/SMAX_tar[i] 
		SMSY[i] <- (1 - LambertW0(exp(1 - Alpha_tar[i])))/BETA[i] 
		SGEN[i] <- -1/BETA[i]*LambertW0(-BETA[i]*SMSY[i]/(exp(Alpha_tar[i]))) 
		SREP[i] <- Alpha_tar[i]/BETA[i] 
	  }
	  
	  for (i in 1:line){
		logSMAX_line_ocean[i] <- b0[1] + b0[2] + (bWA[1] + bWA[2])*lineWA[i] + biaslogSMAX
		SMAX_line_ocean[i] <- exp(logSMAX_line_ocean[i])
		
		logSMAX_line_stream[i] <- b0[1] + (bWA[1])*lineWA[i] + biaslogSMAX 
		SMAX_line_stream[i] <- exp(logSMAX_line_stream[i])
	  }
	  
	  REPORT(b0) # Testing simulate()
	  REPORT(bWA) # Testing simulate()

	  REPORT(logRS_pred)
	  
	  REPORT(logSMAX_re)
	  REPORT(logSMAX_sd)
	  REPORT(SMAX) # E (Srep) for all synoptic data set rivers (25)
	  REPORT(logSMAX)
	  REPORT(logAlpha) # model logAlpha (25)
	  REPORT(Alpha0)
	  REPORT(Alpha02)
	  REPORT(Alpha_re) # random effect parameter for resampling
	  REPORT(Alpha_sd)
	  REPORT(Alpha_pred)  
	  REPORT(SMSY_r)
	  REPORT(BETA_r)
	  REPORT(SREP_r) 
	  REPORT(tauobs) # Necessary to add back in observation error?
	  
	  REPORT(SMAX_tar)
	  REPORT(logSMAX_tar)
	  REPORT(logAlpha_tar) 
	  REPORT(Alpha_tar) 
	  
	  REPORT(BETA)
	  REPORT(SMSY)
	  REPORT(SGEN)
	  REPORT(SREP) 
	  
	  REPORT(SMAX_line_stream) 
	  REPORT(logSMAX_line_stream) 
	  REPORT(SMAX_line_ocean) 
	  REPORT(logSMAX_line_ocean)
	  
	  nll
	}

	obj <- RTMB::MakeADFun(f_smax, par,
		random = c("Alpha_re", "logSMAX_re"), 
		silent=TRUE)

	upper <- numeric(length(obj$par)) + Inf
	lower <- numeric(length(obj$par)) + -Inf
	lower[names(obj$par) == "tauobs"] <- 0
	upper[names(obj$par) == "logSMAX_sd"] <- 100 # Turn off for half dists. 
	lower[names(obj$par) == "logSMAX_sd"] <- 0 
	upper[names(obj$par) == "Alpha_sd"] <- 100 # Turn off for half dists.
	lower[names(obj$par) == "Alpha_sd"] <- 0

	init <- function() {
		listinit <- list(b0 = c(rnorm(1, 10, 1), rnorm(1, 0, 1)), # Contains negatives
			bWA = c(rnorm(1, 0, 1), rnorm(1, 0 ,1)), # Contains negatives
			logSMAX_re = rnorm(N_Stk, 0, 1), # Contains negatives 
			Alpha0 = rnorm(1, 0.6, 1), # Contains negatives
			Alpha_re = rnorm(nrow(dat$WAbase), 0, 1), # Contains negatives
			tauobs = runif(N_Stk, min = 0.005, max = 0.015), # Uniform to REMAIN positive
			logSMAX_sd = runif(1, 0.01, 3), # Positive 
			Alpha_sd = runif(1, 0.01, 3) # Positive
		)
	  
		if (lhdiston) {
			listinit$Alpha02 <- rnorm(1, 0, 1) # alpha0 prior for LH specific dists.
		}
	  
		return(listinit)
	}

	set.seed(1) ; fitstan <- tmbstan(obj, iter = 5000, warmup = 2500, #
					   init = init,
					   seed = 1, 
					   lower = lower, upper = upper,
					   chains = 4, open_progress = FALSE, silent = TRUE); beep(2)

	# End of IWAM simple function
	return(fitstan = fitstan)
}