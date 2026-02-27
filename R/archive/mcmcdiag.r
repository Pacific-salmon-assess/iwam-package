# Example MCMC Diagnostics for stanfit objects

	# Most of these plots are available in the extended figures section:
	# https://pacific-salmon-assess.github.io/iwam-package/TechReportExtFigures.html

# Libaries ####
library(RTMB)
library(ggplot2)
library(dplyr)
library(tidyverse) 
library(progress) # Progress bar
library(tmbstan) # MCMC tmb model sampling
library(bayesplot) # bayesian visualization
library(coda) # bayesian package

# Test and diagnostic plots ####
names(fitstan) # for complete list of parameters from stan object

  # TRACE PLOTS
# mcmc_trace(as.array(fitstan)) # ALL PARAMETER TRACEPLOT
# mcmc_trace(as.array(fitstan), regex_pars = "b0") # Single regex parameter traceplot e.g. b0

  # PAIRS PLOTS
# pairs_pars <- c("b0", "bWA", "logAlpha0") # "logAlpha_sd" "logSREP_sd"
# pairs(fitstan, pars = pairs_pars) # for specific par names from above
# bayesplot::mcmc_pairs(fitstan, regex_pars = pairs_pars)
# bayesplot::mcmc_pairs(fitstan, regex_pars = c("tauobs"))
# bayesplot::mcmc_pairs(fitstan, regex_pars = c("logSREP_re"))
# bayesplot::mcmc_pairs(fitstan, regex_pars = c("logAlpha_re"))

  # ACF
# bayesplot::mcmc_acf(fitstan, regex_pars = "b0") # ACF plot for b0
# bayesplot::mcmc_acf(fitstan) # For every and all - should be a faster way to save them

  # RHAT
# fitstan |> rhat() |> mcmc_rhat() + yaxis_text() # rhat plot for assessing rhat of each parameter

  # LOO/WAIC etc.
# Need a log-likelihood matrix to do this - would have to add a generated
  # quantity in the model to achieve this. Not generated otherwise.

  # Matching Liermann et al. statistics
# mcmc_chains <- As.mcmc.list(fitstan)
# autocorr.plot(mcmc_chains)
# geweke_results <- geweke.diag(mcmc_chains)
# print(geweke_results)
# par(mfrow = c(2, 2))  # adjust depending on number of chains
# for (i in 1:length(mcmc_chains)) {
  # geweke.plot(mcmc_chains[[i]], main = paste("Chain", i))
# }
# heidel_results <- heidel.diag(mcmc_chains)
# print(heidel_results)
# effectiveSize(mcmc_chains)
# gelman.diag(mcmc_chains)
