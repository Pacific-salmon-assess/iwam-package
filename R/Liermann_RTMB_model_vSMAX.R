# Liermann SMAX (B) RTMB Model with MCMC Sampling ####

# Libaries ####
library(RTMB)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyverse) 
library(progress) # Progress bar
library(tmbstan) # MCMC tmb model sampling
# library(TMB) # Original TMB
# library(tidybayes) # bayesian visualization
library(bayesplot) # bayesian visualization
library(coda) # bayesian package
library(beepr) # Sounds
library(viridis) # Colours
library(ggridges) # Ridge plots
library(data.table) # Create data tables for pivoting

here::i_am("R/LambertWs.R") # For non-RStudio functionality
source(here::here("R/LambertWs.R")) # Lambert W function
source(here::here("R/helperFunctions.R")) # For bootstrapping
source(here::here("R/derived_post.R")) # For posterior extraction


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
WAin <- c("DataIn/UpperSoGChinook.csv")
	# c("DataIn/WCVIStocks.csv") or 
	# c("DataIn/Parken_evalstocks.csv") or 
	# c("DataIn/Ordered_backcalculated_noagg.csv")
	# c("DataIn/UpperSoGChinook.csv")

# Data Manipulations ####
srdatwna <- read.csv(here::here("DataIn/SRinputfile.csv")) 
	# Consider _TK **** as alternative with longer dataset
WAbase <- read.csv(here::here("DataIn/WatershedArea.csv"))
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

# Shift log WA for the (mean - base), makes estimation easier
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
            prioronly = 0) # 0-run with data, 1-prior prediction mode

if(dat$prioronly == 1){print("Prior Prediction Mode")} else 
	{print("Posterior Prediction Mode")}

N_Stk <- max(srdat$Stocknumber + 1)
N_Obs <- nrow(srdat)
stk = srdat$Stocknumber + 1

lhdiston <- T # T = LH Specific
bias.cor <- F # T = subtract bias correction terms from expontiated mean terms

# Please note all references to alpha are defined explicitly on the log-scale e.g. alpha = log(alpha)
par <- list(b0 = c(10, 0), # Initial values for WA regression intercepts
			bWA = c(0, 0), # Inital values for WA regression slopes
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

# Please search for '@@' for the code cross-walk when comparing Tech Report equations.
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
  
  # IWAM EQUATION(S): 12 - 15 @@
	# Explanation: Prior distributions of the watershed area regression slope and intercept parameters.
  nll <- nll - sum(dnorm(b0[1], 10, sd = 31.6, log = TRUE)) 
  nll <- nll - sum(dnorm(b0[2], 0, sd = 31.6, log = TRUE)) 
  nll <- nll - sum(dnorm(bWA[1], 0, sd = 31.6, log = TRUE)) 
  nll <- nll - sum(dnorm(bWA[2], 0, sd = 31.6, log = TRUE)) 
	# Can I remove the sum() from these arguments?

  # IWAM EQUATION(S): 6 and 7 @@
	# Explanation: The prior for mean log(productivity) for stream-type and adjustment for ocean-type.
  nll <- nll - sum(dnorm(Alpha0, 0.6, sd = 0.45, log = TRUE)) # Prior (rM) in Liermann et al. (2010)
  if(lhdiston) nll <- nll - sum(dnorm(Alpha02, 0, sd = 31.6, log = TRUE)) # Prior (rD) in Liermann et al. (2010)
  
  ## Second level of hierarchy - Ricker parameters:
  for (i in 1:N_Stk){
	# IWAM EQUATION(S): 5 (8), and 10 (11) @@
	  # Explanation: Eq. 5 is the population level log productivity, with Uniform prior (8) for standard deviation.
	  # Eq. 10 is the random effect term accounting for process error, with a Uniform prior (11) for standard deviation.
    nll <- nll - dnorm(logSMAX_re[i], 0, sd = 1, log = TRUE)
	nll <- nll - dnorm(Alpha_re[i], 0, sd = 1, log = TRUE)

	# cov_par <- corxy*logSMAX_sd*Alpha_sd
	# Sig_RE = cbind(c(logSMAX_sd^2, cov_par), c(cov_par, Alpha_sd^2)) # or do it without correlation to test same results, where corr_par = 0
	# nll <- nll - RTMB::dmvnorm(x = c(logSMAX_re[i], Alpha_re[i]), mu = c(0, 0), Sigma = Sig_RE, log = TRUE) # make logSMAX_re and Alpha_re correlated
		## corr_par has to be from -1 to 1 - so I would set a limit specifically for corr_par
		## corr_par would need a prior e.g. Uniform -1 to 1
	
	# IWAM EQUATION(S): 9 @@
	  # Explanation: The watershed area regression equation, parameterized for SMAX.
    logSMAX[i] <- b0[1] + b0[2]*type[i] + (bWA[1] + bWA[2]*type[i]) * WAbase$logWAshifted[i] + logSMAX_re[i]*logSMAX_sd + biaslogSMAX
    
	SMAX[i] <- exp(logSMAX[i])
    
	# IWAM EQUATION(S): 4 @@
	  # Explanation: The Ricker stock-specific (i) productivity.
    if(lhdiston) logAlpha[i] <- Alpha0 + Alpha02*type[i] + Alpha_re[i]*Alpha_sd + biaslogAlpha
    else logAlpha[i] <- Alpha0 + Alpha_re[i]*Alpha_sd + biaslogAlpha

	# IWAM EQUATION(S): 3 @@
	  # Explanation: The prior for the stock-specific (i) precision for the spawner-recruit model.
    nll <- nll - dgamma(tauobs[i], shape = 0.0001, scale = 1/0.0001, log = TRUE)
  }

  ## First level of hierarchy: Ricker model:
  for (i in 1:N_Obs){
	Alpha_pred <- exp(logAlpha)
	
	# IWAM EQUATION(S): 1 @@
	  # Explanation: The Ricker spawner-recruit model, parameterized for SMAX.
	logRS_pred[i] <- Alpha_pred[stk[i]] - S[i]/SMAX[stk[i]] + biaslogRS[stk[i]]
	# logRS_pred[i] <- logAlpha[stk[i]] - S[i]/SMAX[stk[i]] + biaslogRS[stk[i]]
	# logRS_pred[i] <- logAlpha[stk[i]]*(1 - S[i]/SREP[stk[i]]) + biaslogRS[stk[i]] # Old Ricker parameterization

	# IWAM EQUATION(S): 2 @@
	  # Explanation: The deviation in recruits per spawners.
    if(!prioronly){ # If prioronly is 1, then likelihood is not calculated, if 0 then it is
      nll <- nll - dnorm(logRS[i], logRS_pred[i], sd = sqrt(1/tauobs[stk[i]]), log = TRUE) 
    } 
  }
  
  ## Calculate SMSY for Synoptic set - for plotting
  SMSY_r = numeric(nrow(WAbase))
  BETA_r = numeric(nrow(WAbase))
  SREP_r = numeric(nrow(WAbase))
  
  # IWAM EQUATION(S): 17 and 19 @@
	# Explanation: The calculation of the population management benchmarks SMSY and SREP, and Ricker
	# parameter Beta.
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
	# IWAM EQUATION(S): 20 @@
		# Explanation: The calculation of the population management benchmark SGEN using Lambert's W function
		# as described in Scheuerell (2006).
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
  
  REPORT(b0) 
  REPORT(bWA) 
  REPORT(logRS_pred)
  # REPORT(logRS) # logRS for all 501 data points
  REPORT(logSMAX_re)
  REPORT(logSMAX_sd)
  REPORT(SMAX) 
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
    random = c("Alpha_re", "logSMAX_re"), 
    silent=TRUE)



# LIMITS: Set Upper and Lower Limits ####
upper <- numeric(length(obj$par)) + Inf
lower <- numeric(length(obj$par)) + -Inf
lower[names(obj$par) == "tauobs"] <- 0
upper[names(obj$par) == "logSMAX_sd"] <- 100 # Turn off for half dists. 
lower[names(obj$par) == "logSMAX_sd"] <- 0 
upper[names(obj$par) == "Alpha_sd"] <- 100 # Turn off for half dists.
lower[names(obj$par) == "Alpha_sd"] <- 0



# RANDOM INIT - MATCHING PRIORS ####
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
set.seed(1) ; fitstan <- tmbstan(obj, iter = 5000, warmup = 2500, # default iter/2 for warmup 
                                                                  # - Typically 5000 and 1000
                   init = init, # init = init function or "random" or initfixed for fixed points
                   seed = 1, # set seed or leave out for random - now set at 1 above
                      # but not technically within sampling - unsure if each chain takes seed                  
                   # control = list(adapt_delta = 0.99, max_treedepth = 15),
                   lower = lower, upper = upper,
                   chains = 4, open_progress = FALSE, silent = TRUE); beep(2)
# tmbstan operates by default with NUTS MCMC sampler see: 
	# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0197954

# Save stan fit object
# if(dat$prioronly == 1) {save(fitstan, file = "fitstan_prioronly.RData")} else 
	# {save(fitstan, file = "fitstan.RData")}

# Acquire outputs of MCMC ####
derived_obj <- derived_post(fitstan, model = 'SMAX'); beep(2)
	# add stocknames - see extra code from _Plots.R
dsmax <- derived_obj
dsmaxs <- derived_obj$deripost_summary
dsmaxf <- derived_obj$deripost_full
fitsmax <- fitstan

outpp <- data.frame(
	Stock = WAin$Stock, WA = WAin$WA, lh = WAin$lh,
	SREP_median = dsmaxs$SREP_adj$Median, SREP_mean = dsmaxs$SREP_adj$Mean, SREP_lwr5 = dsmaxs$SREP_adj$LQ_5, SREP_upr95 = dsmaxs$SREP_adj$UQ_95, 
	SMSY_median = dsmaxs$SMSY_adj$Median, SMSY_mean = dsmaxs$SMSY_adj$Mean, SMSY_lwr5 = dsmaxs$SMSY_adj$LQ_5, SMSY_upr95 = dsmaxs$SMSY_adj$UQ_95,
	SGEN_median = dsmaxs$SGEN_adj$Median, SGEN_mean = dsmaxs$SGEN_adj$Mean, SGEN_lwr5 = dsmaxs$SGEN_adj$LQ_5, SGEN_upr95 = dsmaxs$SGEN_adj$UQ_95,
	SMAX_median = dsmaxs$SMAX_tar_adj$Median, SMAX_mean = dsmaxs$SMAX_tar_adj$Mean, SMAX_lwr5 = dsmaxs$SMAX_tar_adj$LQ_5, SMAX_upr95 = dsmaxs$SMAX_tar_adj$UQ_95
)
outpp <- outpp %>% mutate(across(where(is.numeric), round))
# outname <- paste("DataOut/", "UpperSoGChinook", "_out_posteriorpredictive.csv", sep = "")
# write.csv(outpp, here::here(outname), row.names = FALSE)


# Simulate alternative priors
# source(here::here("R/Liermann_RTMB_model_Bootstrap.R")) # Bootstrapping simulations of alternative Ricker alpha priors
# BS.smax.og <- dobootstrap(bsiters = 20000, # 20,000 for full iterations
						# adj = TRUE,
						# bias.cor = FALSE,
						# prod = c("LifeStageModel"),
						# MCMC = TRUE,
						# model = c("SMAX"),
						# Ricprior = c(1, 0.3),
						# prior_rho = c(-0.4),
						# round = FALSE,
						# WAinname = c("DataIn/Parken_evalstocks.csv")); beep(2)
						# c("DataIn/WCVIStocks.csv") or c("DataIn/Parken_evalstocks.csv") or c("DataIn/Nanaimo_test.csv") or c("DataIn/UpperSoGChinook.csv")
# BS.smax.og <- BS.smax.og$BS.dfout
# BS.bpar <- BS.smax$bpar

source(here::here("R/simalpha.r"))
BS.smax <- simalpha(bsiters = 20000, 
			newalpha = c(1, 0.3),
			# prior_rho = c(-0.4),
			WAinname = c("DataIn/Parken_evalstocks.csv")); beep(2)
BS.smax <- BS.smax$BS.dfout

# head(BS.smax$upr - BS.smax$lwr)
# head(BS.smax.og$upr - BS.smax.og$lwr)



# SAVING R OBJECTS: ####
# save(derived_obj, file = "derived_obj.RData")
# if(dat$prioronly == 1) {save(derived_obj, file = "derived_obj_prioronly.RData")} else 
	# {save(derived_obj, file = "derived_obj.RData")}