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
library(here)

here::i_am("R/LambertWs.R") # For non-RStudio functionality
source(here::here("R/LambertWs.R")) # Lambert W function
source(here::here("R/helperFunctions.R")) # For bootstrapping

options(scipen = 999)

## Turn off byte compiler ####
	# This the only way to get RTMB functions to run within wrapper functions
compiler::enableJIT(0)

# parIWAM <- function(){} # pars for RTMB function
# priorIWAM_1 <- function(prior = c("b01", "b02", "bWA1", "bWA2", 
								# "Alpha0", "Alpha02", "logSMAX_re", "Alpha_re", "tauobs"),
					# mub01 = 10, sdb01 = 31.6,
					# mub02 = 0, sdb02 = 31.6,
					# mubWA1 = 0, sdbWA1 = 31.6,
					# mubWA2 = 0, sdbWA2 = 31.6,
					# muAlpha0 = 0.6, sdAlpha0 = 0.45,
					# muAlpha02 = 0, sdAlpha02 = 31.6,
					# mulogSMAX_re = 0, sdlogSMAX_re = 1,
					# muAlpha_re = 0, sdAlpha_re = 1,
					# shapetauobs = 0.0001, scaletauobs = 1/0.0001
	# ){
	
	# prior <- match.arg(prior)
	
	# switch(prior,
		# b01 = dnorm(par$b0[1], mub01, sd = sdb01, log = TRUE),
		# b02 = dnorm(par$b0[2], mub02, sd = sdb02, log = TRUE),
		# bWA1 = dnorm(par$bWA[1], mubWA1, sd = sdbWA1, log = TRUE),
		# bWA2 = dnorm(par$bWA[2], mubWA2, sd = sdbWA2, log = TRUE),
		
		# Alpha0 = dnorm(par$Alpha0, muAlpha0, sd = sdAlpha0, log = TRUE),
		# Alpha02 = dnorm(par$Alpha02, muAlpha02, sd = sdAlpha02, log = TRUE),
		
		# logSMAX_re = dnorm(par$logSMAX_re, mulogSMAX_re, sd = sdlogSMAX_re, log = TRUE),
		# Alpha_re = dnorm(par$Alpha_re, muAlpha_re, sd = sdAlpha_re, log = TRUE),
		
		# tauobs = dgamma(par$tauobs, shape = shapetauobs, scale = scaletauobs, log = TRUE)
	# )
# }

priorIWAM <- function(mub01 = 10,        sdb01 = 31.6,
					mub02 = 0,         sdb02 = 31.6,
					mubWA1 = 0,        sdbWA1 = 31.6,
					mubWA2 = 0,        sdbWA2 = 31.6,
					muAlpha0 = 0.6,    sdAlpha0 = 0.45,
					muAlpha02 = 0,     sdAlpha02 = 31.6,
					mulogSMAX_re = 0,  sdlogSMAX_re = 1,
					muAlpha_re = 0,    sdAlpha_re = 1,
					shapetauobs = 0.0001, scaletauobs = 1/0.0001
	) {
	
	list(b0 = c(mub01, mub02, sdb01, sdb02),
	bWA = c(mubWA1, mubWA2, sdbWA1, sdbWA2),
	Alpha0 = c(muAlpha0, sdAlpha0),
	Alpha02 = c(muAlpha02, sdAlpha02),
	logSMAX_re = c(mulogSMAX_re, sdlogSMAX_re),
	Alpha_re = c(muAlpha_re, sdAlpha_re),
	tauobs = c(shapetauobs, scaletauobs)
    )
	
	# instead do the dnorms for all the priors
	# and then save them as a vectored list as above so they can be called
	# then add a prior object to IWAM() so that it takes an object created
	# by priorIWAM
}

check_prior <- function(prior) {
  
  expected_lengths <- c(b0 = 4, bWA = 4, Alpha0 = 2, Alpha02 = 2, 
                        logSMAX_re = 2, Alpha_re = 2, tauobs = 2)
  
  # Check that prior is a list
  if (!is.list(prior)) {
    stop("'prior' must be a list.")
  }
  
  # Check lengths of each element
  for (nm in names(expected_lengths)) {
    if (length(prior[[nm]]) != expected_lengths[nm]) {
      stop(paste0("'prior$", nm, "' must have length ", expected_lengths[nm], 
                  ", but has length ", length(prior[[nm]]), "."))
    }
  }
  
  message("'prior' passed all checks.")
}

# initIWAM <- function(){} # Random inits matching priors for sampling 
# coreIWAM <- function(){} # Only core data processing, function - draws from above

IWAM <- function(WAin = c('DataIn/Parken_evalstocks.csv'),
					# srdat = c(''), # User inputted path for private SR data
					prioronly = 0, # 1 = ...
					lhdiston = TRUE,
					bias.cor = FALSE,
					
					# prior = priorIWAM,
					# parb0 = c(10,0),
					# parbWA = c(0,0),
					# parlogSMAX_re = numeric(), # is a numeric of length
					# parAlpha0 = 0.6,
					# parAlpha_re = numeric(), # is a numeric of length
					# partauobs = numeric(), # is a numeric of length
					# parlogSMAX_sd = 1, # try 0.01
					# parAlpha_sd = 1, # try 0.01
					# parAlpha02 = 0,
					
					# biaslogSMAX = 0,
					# biaslogAlpha = 0,
					# biaslogRS = numeric(), # is a numeric of length
					
					# priorb0 = c(10, 31.6, 0, 31.6),
					# priorbWA = c(0, 31.6, 0, 31.6),
					# priorAlpha0 = c(0.6, 0.45),
					# priorAlpha02 = c(0, 31.6),
					# priorlogSMAX_re = c(0, 1),
					# priorAlpha_re = c(0, 1),
					# priortauobs = c(0.0001, 1/0.0001),
					
					random = c("Alpha_re", "logSMAX_re"),
					
					lowerlimit = 0, # lower limit of uniforms
					upperlimit = 100, # upper limit of uniforms
					
					# initb0 = c(rnorm(1, 10, 1), rnorm(1, 0, 1)),
					# initbWA = c(rnorm(1, 0, 1), rnorm(1, 0 ,1)),
					# initlogSMAX_re = rnorm(N_Stk, 0, 1), # is a numeric of length
					# initAlpha0 = rnorm(1, 0.6, 1),
					# initAlpha_re = rnorm(nrow(dat$WAbase), 0, 1), # is a numeric of length
					# inittauobs =runif(N_Stk, min = 0.005, max = 0.015), # is a numeric of length
					# initlogSMAX_sd = runif(1, 0.01, 3),
					# initAlpha_sd = runif(1, 0.01, 3),
					# initAlpha02 = rnorm(1, 0, 1),
					
					myiter = 5000,
					mywarmup = 2500,
					myseed = 1,
					mychains = 4)
{

	options(scipen = 999)

	# New LambertW0 see: https://github.com/pbs-assess/renewassess/blob/main/code/RTMB/PosteriorPredictiveSample.r
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
	# WAin <- c("DataIn/Parken_evalstocks.csv")
	# c("DataIn/WCVIStocks.csv") or c("DataIn/Parken_evalstocks.csv") 
	# WAin <- c("DataIn/Ordered_backcalculated_noagg.csv")

	# Data Manipulations ####
	srdatwna <- read.csv(here::here("inst/extdata/SRinputfile.csv")) # Consider _TK **** as alternative with longer dataset
	WAbase <- read.csv(here::here("inst/extdata/WatershedArea.csv"))
	WAin <- read.csv(here::here(WAin))

	# Data setup ####
	# Remove Hoko and Hoh stocks - consider removing Skagit OR Chehalis
		# Due to sampling reasons explained in Parken 2006.
	srdatwna <- srdatwna %>% 
		filter(!Name %in% c("Hoko","Hoh")) # |>
		# (!(Name == "Skagit")) # |> # Skagit is #22
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

	# Shift log WA for the mean - base - makes estimation easier
	mean_logWA <- mean(WAbase$logWA)
	WAbase$logWAshifted <- WAbase$logWA - mean_logWA

	WAin$logWA <- log(WAin$WA)
	WAin$logWAshifted_t <- WAin$logWA - mean_logWA

	lifehist <- srdat %>% dplyr::select(Stocknumber, Name, Stream) %>% 
		group_by(Stocknumber) %>% 
		summarize(lh=max(Stream))



	## RTMB dat and par setup ####
	dat <- list(srdat = srdat,
				WAbase = WAbase,
				WAin = WAin,
				lineWA =  seq(min(WAbase$logWAshifted), 
							  max(WAbase$logWAshifted), 0.1),
				mean_logWA = mean_logWA,
				logRS = log(srdat$Rec) - log(srdat$Sp),
				prioronly = prioronly) # 0-run with data, 1-prior prediction mode

	N_Stk <- max(srdat$Stocknumber + 1)
	N_Obs <- nrow(srdat)
	stk = srdat$Stocknumber + 1

	# lhdiston <- T # T = LH Specific
	# bias.cor <- F # T = subtract bias correction terms from expontiated mean terms

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

	check_prior(prior) ## CHECK ##
	prior <- priorIWAM()

	# Please note all references to alpha are defined explicitly on the log-scale e.g. alpha = log(alpha)
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
	  
	  # Why is logRS_pred not a parameter or vector input here?
	  logRS_pred <- numeric(N_Obs) # Does this still report if not a vector?

	  SMAX_tar <- numeric(N_Pred)
	  logSMAX_tar <- numeric(N_Pred) 
	  logAlpha_tar <- numeric(N_Pred)
	  
	  # Simulated line vectors
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
	  
	  # Can I remove the sum() from these arguments?
	  nll <- nll - sum(dnorm(b0[1], prior$b0[1], sd = prior$b0[3], log = TRUE)) # Prior b0
	  nll <- nll - sum(dnorm(b0[2], prior$b0[2], sd = prior$b0[4], log = TRUE)) # Prior b0
	  nll <- nll - sum(dnorm(bWA[1], prior$bWA[1], sd = prior$bWA[3], log = TRUE)) # Prior bWA
	  nll <- nll - sum(dnorm(bWA[2], prior$bWA[2], sd = prior$bWA[4], log = TRUE)) # Prior bWA
	  
	  nll <- nll - sum(dnorm(Alpha0, prior$Alpha0[1], sd = prior$Alpha0[2], log = TRUE)) # Prior (rM)
	  if(lhdiston) nll <- nll - sum(dnorm(Alpha02, prior$Alpha02[1], sd = prior$Alpha02[2], log = TRUE)) # Prior (rD)
	  
	  ## Second level of hierarchy - Ricker parameters:
	  for (i in 1:N_Stk){
		nll <- nll - dnorm(logSMAX_re[i], prior$logSMAX_re[1], sd = prior$logSMAX_re[2], log = TRUE)
		nll <- nll - dnorm(Alpha_re[i], prior$Alpha_re[1], sd = prior$Alpha_re[2], log = TRUE)

		# nll <- nll - dmvnorm(c(logSMAX_re[i],Alpha_re[i]), c(0,0), Sig_RE, log = TRUE) # make logSMAX_re and Alpha_re correlated
		# Sig_RE = cbind(c(1,corr_par),c(corr_par,1)) or do it without correlation to test same results, where corr_par = 0
			## corr_par has to be from -1 to 1 - so I would set a limit specifically for corr_par
			## corr_par would need a prior e.g. Uniform -1 to 1
		
		logSMAX[i] <- b0[1] + b0[2]*type[i] + (bWA[1] + bWA[2]*type[i]) * WAbase$logWAshifted[i] + logSMAX_re[i]*logSMAX_sd + biaslogSMAX
		SMAX[i] <- exp(logSMAX[i])
			
		if(lhdiston) logAlpha[i] <- Alpha0 + Alpha02*type[i] + Alpha_re[i]*Alpha_sd + biaslogAlpha
		else logAlpha[i] <- Alpha0 + Alpha_re[i]*Alpha_sd + biaslogAlpha

		nll <- nll - dgamma(tauobs[i], shape = prior$tauobs[1], scale = prior$tauobs[2], log = TRUE)
	  }

	  ## First level of hierarchy: Ricker model:
	  for (i in 1:N_Obs){
		Alpha_pred <- exp(logAlpha)
		logRS_pred[i] <- Alpha_pred[stk[i]] - S[i]/SMAX[stk[i]] + biaslogRS[stk[i]]
		# logRS_pred[i] <- logAlpha[stk[i]] - S[i]/SMAX[stk[i]] + biaslogRS[stk[i]]
		# logRS_pred[i] <- logAlpha[stk[i]]*(1 - S[i]/SREP[stk[i]]) + biaslogRS[stk[i]] # Old Ricker parameterization

		if(!prioronly){ # If prioronly is 1, then likelihood is not calculated, if 0 then it is
		  nll <- nll - dnorm(logRS[i], logRS_pred[i], sd = sqrt(1/tauobs[stk[i]]), log = TRUE)
		} 
	  }
	  
	  ## Calculate SMSY for Synoptic set - for plotting
	  SMSY_r = numeric(nrow(WAbase))
	  BETA_r = numeric(nrow(WAbase))
	  SREP_r = numeric(nrow(WAbase))
	  
	  for (i in 1:N_Stk){
		BETA_r[i] <- 1/SMAX[i] 
		SMSY_r[i] <- (1 - LambertW0(exp(1 - Alpha_pred[i]))) / BETA_r[i] 
		SREP_r[i] <- Alpha_pred[i]/BETA_r[i] 
	  }

	  ## PREDICTIONS
	  BETA = numeric(nrow(WAin))
	  SMSY = numeric(nrow(WAin))
	  SGEN = numeric(nrow(WAin))
	  SREP = numeric(nrow(WAin))

	  # NOTE: For simplicity with RTMB and TMBstan, predictions were calculated as part of the posterior and excluded from the negative log-likelihood.
	  for (i in 1:N_Pred){
		# Impose a new alpha here ...
		
		if(lhdiston) logAlpha_tar[i] <- Alpha0 + Alpha02*type_tar[i] + biaslogAlpha # + biaslogAlpha + logAlpha_sd^2/2 
		else logAlpha_tar[i] <- Alpha0 + biaslogAlpha # + biaslogAlpha + logAlpha_sd^2/2 
		
		# posterior predictive alpha?
		# if(lhdiston) logAlpha_tar[i] <- Alpha0 + Alpha02*type_tar[i] + biaslogAlpha
		# else logAlpha_tar[i] <- Alpha0 + biaslogAlpha

		logSMAX_tar[i] <- b0[1] + b0[2]*type_tar[i] + (bWA[1] + bWA[2]*type_tar[i])*WAin$logWAshifted_t[i] + biaslogSMAX
		SMAX_tar[i] <- exp(logSMAX_tar[i]) 
		
		Alpha_tar <- exp(logAlpha_tar)
		# Predict BETA
		# BETA[i] <- logAlpha_tar[i]/SREP_tar[i] 
		BETA[i] <- 1/SMAX_tar[i] 
		# Predict SMSY
		SMSY[i] <- (1 - LambertW0(exp(1 - Alpha_tar[i])))/BETA[i] 
		# Predict SGEN
		SGEN[i] <- -1/BETA[i]*LambertW0(-BETA[i]*SMSY[i]/(exp(Alpha_tar[i]))) 
		# Predict SREP 
		SREP[i] <- Alpha_tar[i]/BETA[i] 
	  }
	  
	  # Create predictions on an simulated line
	  for (i in 1:line){
		logSMAX_line_ocean[i] <- b0[1] + b0[2] + (bWA[1] + bWA[2])*lineWA[i] + biaslogSMAX
		SMAX_line_ocean[i] <- exp(logSMAX_line_ocean[i])
		
		logSMAX_line_stream[i] <- b0[1] + (bWA[1])*lineWA[i] + biaslogSMAX 
		SMAX_line_stream[i] <- exp(logSMAX_line_stream[i])
	  }
	  
	  REPORT(b0) # Testing simulate()
	  REPORT(bWA) # Testing simulate()

	  REPORT(logRS_pred)
	  
	  # REPORT(logRS) # logRS for all 501 data points
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
	  
	  # Simulated line values for plotting
	  REPORT(SMAX_line_stream) 
	  REPORT(logSMAX_line_stream) 
	  REPORT(SMAX_line_ocean) 
	  REPORT(logSMAX_line_ocean)
	  
	  nll
	}

	## MakeADFun ####
	obj <- RTMB::MakeADFun(f_smax, par,
		random = c(random[1], random[2]), 
		silent=TRUE)



	# LIMITS: Set Upper and Lower Limits ####
	upper <- numeric(length(obj$par)) + Inf
	lower <- numeric(length(obj$par)) + -Inf
	lower[names(obj$par) == "tauobs"] <- lowerlimit
	upper[names(obj$par) == "logSMAX_sd"] <- upperlimit # Turn off for half dists. 
	lower[names(obj$par) == "logSMAX_sd"] <- lowerlimit 
	upper[names(obj$par) == "Alpha_sd"] <- upperlimit # Turn off for half dists.
	lower[names(obj$par) == "Alpha_sd"] <- lowerlimit



	# MCMC ####
	# RANDOM INIT - MATCHING PRIORS
	init <- function() {
		listinit <- list(b0 = c(rnorm(1, 10, 1), rnorm(1, 0, 1)), # Contains negatives
			bWA = c(rnorm(1, 0, 1), rnorm(1, 0 ,1)), # Contains negatives
		   
			# logRS_pred = rnorm(N_Obs, 0, 1),
			logSMAX_re = rnorm(N_Stk, 0, 1), # Contains negatives 
			Alpha0 = rnorm(1, 0.6, 1), # Contains negatives
			# logAlpha02 = rnorm(1, 0, 1) , # NEW: alpha0 prior for LH specific dists.
			Alpha_re = rnorm(nrow(dat$WAbase), 0, 1), # Contains negatives

			tauobs = runif(N_Stk, min = 0.005, max = 0.015), # Uniform to REMAIN positive
		   
			# Should these be 0.01 to 100s? Would that be more accurate?
			logSMAX_sd = runif(1, 0.01, 3), # Positive 
			Alpha_sd = runif(1, 0.01, 3) # Positive
		)
	
		if (lhdiston) {
			listinit$Alpha02 <- rnorm(1, 0, 1) # alpha0 prior for LH specific dists.
		}
	
		return(listinit)
	}

	# SAMPLE MCMC ####
	# The set.seed precedes this in order to have fixed "identical" runs with initfixed
	set.seed(1) ; fitstan <- tmbstan(obj, iter = myiter, warmup = mywarmup, # default iter/2 for warmup 
																	  # - Typically 5000 and 1000
					   init = init, # init = init function or "random" or initfixed for fixed points
					   seed = myseed, # set seed or leave out for random - now set at 1 above
						  # but not technically within sampling - unsure if each chain takes seed                  
					   # control = list(adapt_delta = 0.99, max_treedepth = 15),
					   lower = lower, upper = upper,
					   chains = mychains, open_progress = FALSE, silent = TRUE); beep(2)
	# tmbstan operates by default with NUTS MCMC sampler see: 
		# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0197954


	return(fitstan)
}