# PLOTTING: Liermann Srep (E) RTMB Model with MCMC Sampling ####

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
library(HDInterval)

library(scales) 
library(ggnewscale)

source(here::here("R/LambertWs.R")) # Lambert W function
source(here::here("R/helperFunctions.R")) # For bootstrapping
source(here::here("R/derived_post.R")) # For posterior extraction

# Run Liermann_RTMB_model_vC.R to load in all necessary global objects
    # e.g. derived_obj, WAin, WAbase, srdat, fitstan, pars, names

# SAVING R OBJECTS: ####
# save(my_object, file = "my_object.RData") # and then load:
# load("my_object.RData")

# add stocknames - could add them in to each object of derived_obj?
# Stocknames <- WAin$Stock
# Srep_example <- cbind(derived_obj$deripost_summary$SREP_tar_adj, Stocknames) |> 
#   Srep_example[c(1, 7, 2, 3, 4, 5, 6)]



# PRIOR PUSHFORWARD/PREDICTIVE CHECKS ####
# APPROACH - Prior prediction with data exclusion
	# Make sure dat$prioronly <- 1 to not include data in the nll
	# Run all the way down to derived_obj
# Pushforward tests:
  # Take the Prior Prediction Mode (dat$prioronly == 1)
  # Then plot desired priors/parameters E.g. logAlpha
# Predictive tests: 
  # Exactly the same as for Posterior - just depends on whether data has been included in the likelihood or not.
  # Note: prior predictive checks (if taking the method of excluding data)
  # Can be repeated using the below posterior predictive checks.
  # Specifically: logRS, logAlpha



#### POSTERIOR PREDICTIVE CHECKS ####################################################################################################################
  # Run full model with dat$prioronly <- 0 for data included in nll
  # Run derived_post() to extract posterior from chains
  # Randomly sample for logRS
  # Compare logRS with logRS_pred (including observation error)
  # e.g. rnorm(1, derived_obj$deripost_summary$logRS_pred$Median, 
             # sd  = sqrt(1/derived_obj$deripost_summary$tauobs$Median))
slogRS_pred <- derived_obj$deripost_full$logRS_pred
stauobs <- derived_obj$deripost_full$tauobs
simlogRS <- matrix(NA, nrow = dim(slogRS_pred)[1], ncol = dim(slogRS_pred)[2])
for (i in 1:dim(slogRS_pred)[1]){
  simlogRS[i, ] <- rnorm(dim(slogRS_pred)[2], # 501
                        mean = slogRS_pred[i, ], 
                        sd = sqrt(1/stauobs[i, stk]))
} 

# Retrieve 9 samples out of the 10000 iterations
nsim <- 9
draws <- sample(1:dim(slogRS_pred)[1], nsim)
savedsims <- simlogRS[draws,] # [1:9,]

ppcheckdata0 <- data.frame(index = seq_along(dat$logRS), value = dat$logRS) # data.frame(dat$logRS)
ppcheckplot0 <- ggplot(ppcheckdata0, aes(x = index, y = value)) + 
	geom_point(alpha = 0.4) + 
	ylab("") + 
	xlab("") + 
	theme_classic() + 
	theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ppcheckdata <- lapply(1:9, function(i) data.frame())
ppcheckplot <- c()
for (i in 1:nsim){
	ppcheckdata[[i]] <- data.frame(index = seq_along(simlogRS[i, ]), value = simlogRS[i, ]) # or just data.frame(simlogRS[i,])
	ppcheckplot[[i]] <- ggplot(ppcheckdata[[i]], aes(x = index, y = value)) + 
		geom_point(alpha = 0.4) + 
		ylab("") + 
		xlab("") + 
		theme_classic() + 
		theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}
do.call(grid.arrange, c(list(ppcheckplot0), ppcheckplot[1:nsim], ncol = 2, 
		bottom = "Observation Order",
		left = "log(R/S)"))



#### POSTERIOR P-VALUES #############################################################################################################################
    # Take the mean and sd of each of the iterations e.g. do more than 9 in this case
    # Then plot the histogram of these means against the mean of
    # the observed data
    # Then the p-value is the proportion of simulated means that are
    # ABOVE the mean of the observations.
    # Report this value.
mean_simlogRS <- apply(simlogRS, 1, mean) # Mean of all simulated logRS
mean_simlogRS <- data.frame(mean_simlogRS)

hist(data.frame(mean_simlogRS), xlab = "Mean of Log(R/S)", 
     main = "Visualization of Posterior Predictive P-Value")
abline(v = mean(dat$logRS), col = "red", lwd = 3) # Mean of observed logRS

ggplot(mean_simlogRS, aes(x = mean_simlogRS)) + 
	geom_density(fill = "grey") +
	geom_vline(xintercept = mean(dat$logRS), color = "red", linetype = "dashed") + 
	xlab("Mean of Simulated log(R/S)") + 
	ylab("Density") + 
	theme_classic()

pvalpp <- sum(mean_simlogRS > mean(dat$logRS))/length(mean_simlogRS) # Proportion of simulated means above observed mean
print(paste0("Posterior Predictive P-Value = ", pvalpp))
  # A slightly higher bias than observed predictions



#### Posterior Predictive Distribution: logAlpha ####################################################################################################
slogAlpha0 <- derived_obj$deripost_full$logAlpha0 # dim(10000, 1)
slogAlpha02 <- derived_obj$deripost_full$logAlpha02 # dim(10000, 1)
slogAlpha_sd <- derived_obj$deripost_full$logAlpha_sd # dim(10000, 1)
simlogAlpha_s <- matrix(NA, nrow = dim(slogAlpha0)[1], ncol = dim(slogAlpha0)[2])
simlogAlpha_o <- matrix(NA, nrow = dim(slogAlpha0)[1], ncol = dim(slogAlpha0)[2])
for (i in 1:dim(slogAlpha0)[1]){
  simlogAlpha_s[i] <- slogAlpha0[i] + rnorm(1, mean = 0, sd = slogAlpha_sd[i])
  simlogAlpha_o[i] <- simlogAlpha_s[i] + slogAlpha02[i]
}

ggplot() +
  geom_density(aes(x = simlogAlpha_s), # For a new STREAM observation
               fill = "forestgreen", alpha = 0.4, color = "forestgreen", linewidth = 1.2) +
  geom_density(aes(x = simlogAlpha_o), # For a new OCEAN observation
              fill = "skyblue", alpha = 0.4, color = "skyblue", linewidth = 1.2) +
  theme_classic() +
  labs(x = "Mean of uncentered logAlpha Posterior Predictive Distribution (Stream and Ocean)", 
         y = "Density")
  # labs(x = "Mean of uncentered logAlpha Prior Predictive Distribution (Stream and Ocean)", 
       # y = "Density")



#### Pushforward (under prior predictive): logAlpha #################################################################################################
ggplot() +
  geom_density(aes(x = slogAlpha02), # For a pushforward alpha?
               fill = "grey", alpha = 0.4, color = "grey", linewidth = 1.2) +
  theme_classic() +
  labs(x = "Hierarchical Mean logAlpha_LH Prior Pushforward Distribution", y = "Density")

# Posterior OR Prior Distribution ridge plot: logAlpha ####
  # IF Prior - then its pushforward - its a per stock value - NOT a new observ.
dfalpharidge <- derived_obj$deripost_full$logAlpha[, 1:25]
Stocknames <- WAbase$Name
colnames(dfalpharidge) <- Stocknames
alpharidgetable <- as.data.table(dfalpharidge)
alpharidgetable_long <- melt(alpharidgetable, measure.vars = Stocknames,
                variable.name = "Stock", value.name = "Value")
alpharidgetable_long <- as.data.table(alpharidgetable_long)
TypeLabels <- ifelse(lifehist$lh == 0, "S", "O")
	alpharidgetable_long[, Type := TypeLabels[match(Stock, Stocknames)]]
n_S <- sum(TypeLabels == "S")
n_O <- sum(TypeLabels == "O")

ggplot(alpharidgetable_long, aes(x = Value, y = Stock, fill = interaction(Type, Stock))) +
  geom_density_ridges(color = "gray20", alpha = 0.8, scale = 1.2) +
  theme_classic() +
  labs(x = "logAlpha", y = "Stock") +
  theme(legend.position = "none") +
  scale_fill_manual(
    values = c(
      setNames(colorRampPalette(c("#a8e6a3", "#0b6e0b"))(n_S),
               paste("S.", Stocknames[TypeLabels == "S"], sep = "")),
      setNames(colorRampPalette(c("#a3d5ff", "#084a8b"))(n_O),
               paste("O.", Stocknames[TypeLabels == "O"], sep = ""))
    )
 )

# ggplot(alpharidge_long, aes(x = Value, y = factor(Index), fill = factor(Index))) +
#   geom_density_ridges(alpha = 0.5, scale = 1.2, color = "forestgreen") +
#   theme_classic() +
#   labs(x = "logAlpha Posterior", y = "Stock ID") +
#   theme(legend.position = "none")
# Change ordering of stocks e.g. by WA size as with point-wise comp. plots
# Can I add in ocean/stream-type identifiers and colour via?

# beta?

# Posterior Predictive Calculation for: SREP (E)
	# logE[i] <- b0[1] + b0[2]*type[i] + (bWA[1] + bWA[2]*type[i]) * WAbase$logWAshifted[i] + logSREP_re[i]*logSREP_sd
	# Need b0, bWA, WA, and random effects
	# Will get one for stream and one for ocean



#### Posterior Distributions of b0 and bWA ##########################################################################################################
ggplot() +
  geom_density(aes(x = derived_obj$deripost_full$b0[,1]), 
               fill = "forestgreen", alpha = 0.4, color = "forestgreen", linewidth = 1.2) +
  geom_density(aes(x = (derived_obj$deripost_full$b0[,2] + derived_obj$deripost_full$b0[,1])),
              fill = "skyblue", alpha = 0.4, color = "skyblue", linewidth = 1.2) +
  theme_classic() +
  labs(x = "b0", y = "Density")

ggplot() +
  geom_density(aes(x = derived_obj$deripost_full$bWA[,1]), 
               fill = "forestgreen", alpha = 0.4, color = "forestgreen", size = 1.2) +
  geom_density(aes(x = (derived_obj$deripost_full$bWA[,2] + derived_obj$deripost_full$bWA[,1])),
              fill = "skyblue", alpha = 0.4, color = "skyblue", size = 1.2) +
  theme_classic() +
  labs(x = "bWA", y = "Density")



#### Plot SR Relationship Curves ####################################################################################################################
lineSREPdraws <- derived_obj$deripost_full$SREP # 10000, 25
lineAlphadraws <- exp(derived_obj$deripost_full$logAlpha) # 10000, 25
SSdraws <- matrix(NA, nrow = 10000, ncol = 100) # this needs to be 10,000 by 100?
RRdraws <- matrix(NA, nrow = 10000, ncol = 100)
RRmedian <- lineAlphamedian <- lineSREPmedian <- SSmedian <- NA

rowsample <- sample(1:10000, 1000)

lineSREPdraws <- derived_obj$deripost_full$SREP # 10000, 25
lineAlphadraws <- exp(derived_obj$deripost_full$logAlpha) # 10000, 25
rowsample <- sample(1:10000, 100) # 1:10000, 1000

Smsylines <- derived_obj$deripost_summary$SMSY_r$Median

par(mfrow = c(5, 5), mar = c(2, 2, 1, 0.1) + 0.1, oma = c(3, 3, 1, 1))

for (i in 1:25) {
  stock_name <- unique(srdat$Name)[i]
  spawners   <- srdat$Sp[srdat$Name == stock_name]
  recruits   <- srdat$Rec[srdat$Name == stock_name]
  Smax       <- max(spawners)

  SSseq <- seq(Smax/100, Smax, length.out = 100)
  alpha_draws <- lineAlphadraws[, i]
  srep_draws  <- lineSREPdraws[, i]

  SSmat <- matrix(SSseq, nrow = 10000, ncol = 100, byrow = TRUE)
  RRmat <- SSmat * alpha_draws^(1 - SSmat / srep_draws)

  RRmed <- SSseq * median(alpha_draws)^(1 - SSseq / median(srep_draws))

  plot(spawners, recruits, xlim = c(0, Smax + Smax/10), ylim = c(0, max(recruits)))
  abline(v = Smsylines[i], col = 'red', lty = 'dashed')
  mtext(stock_name, side = 3, cex = 0.8)
  matlines(t(SSmat[rowsample, ]), t(RRmat[rowsample, ]),
           col = rgb(0, 0, 0, 0.1), lty = 1)
  lines(SSseq, RRmed, col = "red", lwd = 2)
}

mtext("Spawners", side = 1, line = 1, outer = TRUE, cex = 1.3)
mtext("Recruitment", side = 2, line = 1, outer = TRUE, cex = 1.3)



#### RESIDUALS ######################################################################################################################################
#### For the relationship between capacity and watershed size 
savelogSREP_re <- as.data.frame(derived_obj$deripost_full$logSREP_re) # dim 10000 25 (per stock)
savelogSREP_sd <- as.data.frame(derived_obj$deripost_full$logSREP_sd) # dim 10000 1 (single)

saveresiduals <- as.data.frame(derived_obj$deripost_full$logSREP_re * derived_obj$deripost_full$logSREP_sd[,1])

dflogSREP_re <- savelogSREP_re %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(
    cols = starts_with("V"),
	names_prefix = "V", 	
    names_to = "chain",
    values_to = "value"
  ) %>% 
  mutate(chain = as.numeric(chain))
  
dfresiduals <- saveresiduals %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(
    cols = starts_with("V"),
	names_prefix = "V", 	
    names_to = "chain",
    values_to = "value"
  ) %>% 
  mutate(chain = as.numeric(chain))

inter <- WAbase %>%
  mutate(chain = row_number())

dflogSREP_re <- dflogSREP_re %>%
  left_join(inter, by = "chain")
  
dfresiduals <- dfresiduals %>%
  left_join(inter, by = "chain")

# were plotted against log-centered watershed size and latitude
# to investigate potential patterns as in a standard regression (i.e. patterns in the mean or variance of the residuals).
res1 <- ggplot(data = dflogSREP_re, aes(x = logWAshifted, y = value)) + 
	geom_point(alpha = 0.1) + 
	theme_classic() + 
	labs(y = "Residual of WA Relationship", x = "Log Centered WA")
	
res2 <- ggplot(data = dflogSREP_re, aes(x = Latitude, y = value)) + 
	geom_point(alpha = 0.1) + 
	theme_classic() + 
	labs(y = "Residual of WA Relationship", x = "Latitude")

grid.arrange(res1, res2, nrow = 1)

#### The residuals from the spawner-recruit relationship (the Ricker model)
	# tauobs - watch transform out of precision?
	# or residuals of logRS? calculated manually?
		# e.g. observed logRS - posterior
residlogRS <- matrix(0, nrow = 10000, ncol = 501)
for (i in 1:501){
	residlogRS[,i] <- dat$logRS[i] - derived_obj$deripost_full$logRS_pred[,i]
}

df_residlogRS <- as.data.frame(residlogRS)

dfresidlogRS <- df_residlogRS %>%
  mutate(iteration = row_number()) %>%
  pivot_longer(
    cols = starts_with("V"),
	names_prefix = "V", 	
    names_to = "chain",
    values_to = "value"
  ) %>% 
  mutate(chain = as.numeric(chain))
  
# were plotted against year and spawners.
	# By year or observation
	# logRS is by observation
# connect to dat$srdat? In order to plot by year?
tempnumbersrdat <- dat$srdat %>%
  mutate(chain = row_number())
residplot <- dfresidlogRS %>%
  left_join(tempnumbersrdat, by = "chain")

residplot_median <- residplot %>%
  group_by(Name, chain, Stocknumber, Yr, Sp, Rec, Stream, yr_num, Comments, lh) %>%
  summarise(median_value = median(value), .groups = "drop")

plot_listSp = list()
for (i in 1:25) {
	plot_listSp[[i]] <- ggplot(data = residplot_median[residplot_median$Stocknumber == i - 1,], aes(x = Sp, y = median_value)) + 
		geom_point(alpha = 0.5) + 
		theme_classic() + 
		ggtitle(unique(residplot_median$Name[residplot_median$Stocknumber == i - 1])) + 
		labs(x = "Raw Observed Spawners", y = "Residual of logRS")
}

grid.arrange(grobs = plot_listSp, nrow = 5) # for a 5 x 5 grid of 25 synoptic stocks

plot_listYr = list()
for (i in 1:25) {
	plot_listYr[[i]] <- ggplot(data = residplot_median[residplot_median$Stocknumber == i - 1,], aes(x = Yr, y = median_value)) + 
		geom_point(alpha = 0.5) + 
		theme_classic() + 
		ggtitle(unique(residplot_median$Name[residplot_median$Stocknumber == i - 1])) + 
		labs(x = "Observation Year", y = "Residual of logRS")
}

grid.arrange(grobs = plot_listYr, nrow = 5) # for a 5 x 5 grid of 25 synoptic stocks



#### ACF PLOT PER STOCK OF RESIDUALS ################################################################################################################
# EXAMPLE line: acf(residplot$value[residplot$Name == 'Harrison'])
	# median of posterior - and then subtract to calculate the residual
	# calculate residual for each of 25 stocks
par(mfrow = c(5, 5), mar = c(2, 2, 1, 0.1) + 0.1, oma = c(3, 3, 1, 1))

for (i in 1:25) {
	acf(residplot_median$median_value[residplot$Stocknumber == i - 1])  
}
mtext("Lag", side = 1, line = 1, outer = TRUE, cex = 1.3)
mtext("Autocorrelation", side = 2, line = 1, outer = TRUE, cex = 1.3)


# Option 2: Residual from posterior median prediction
# Calculate median prediction across iterations for each observation
median_logRS_pred <- apply(derived_obj$deripost_full$logRS_pred, 2, median)

# Calculate residuals
residlogRS_opt2 <- dat$logRS - median_logRS_pred

# Create dataframe matching your structure - ASSUMING ORDER HASN'T CHANGED
residplot_median_opt2 <- dat$srdat %>%
  mutate(
    chain = row_number(),
    median_value = residlogRS_opt2
  )

# Plot ACF per stock
par(mfrow = c(5, 5), mar = c(2, 2, 1, 0.1) + 0.1, oma = c(3, 3, 1, 1))
for (i in 1:25) {
  acf(residplot_median_opt2$median_value[residplot_median_opt2$Stocknumber == i - 1])  
}
mtext("Lag", side = 1, line = 1, outer = TRUE, cex = 1.3)
mtext("Autocorrelation", side = 2, line = 1, outer = TRUE, cex = 1.3)



#### PLOT: POINT-WISE BENCHMARK COMPARISONS #########################################################################################################
# Prepare/load datasets for plotting #
ParkenCaseStudy <- read.csv(here::here("DataIn/Parken_evalstocks.csv")) # Case study stocks

dpars_srep <- dsrep$deripost_summary
dpars_smax <- dsmax$deripost_summary

targets <- WAin |>
  rename("Stock_name" = Stock) # This name can change depending on what data sets are being run

#### SREP MODEL RESULTS ####
  # SMSY Estimate for TARGET STOCKS
targets1_srep <- cbind(targets, dsrep$deripost_summary$SMSY) |> 
  rename("SMSY_mean" = Mean, "SMSY_median" = Median,
    "SMSY_LQ_5" = LQ_5, "SMSY_UQ_95" = UQ_95, "SMSY_Stocknum" = Stock,
    "SMSY_Mode" = PosteriorMode)
  # SGEN Estimate for TARGET STOCKS
targets2_srep <- cbind(targets1_srep, dsrep$deripost_summary$SGEN) |> 
  rename("SGEN_mean" = Mean, "SGEN_median" = Median,
    "SGEN_LQ_5" = LQ_5, "SGEN_UQ_95" = UQ_95, "SGEN_Stocknum" = Stock,
    "SGEN_Mode" = PosteriorMode)
  # Marginal Mean for SREP (E)
targets3_srep <- cbind(targets2_srep, dsrep$deripost_summary$SREP_tar_adj) |> 
  rename("SREP_tar_adj_mean" = Mean, "SREP_tar_adj_median" = Median,
    "SREP_tar_adj_LQ_5" = LQ_5, "SREP_tar_adj_UQ_95" = UQ_95, "SREP_tar_adj_Stocknum" = Stock,
    "SREP_tar_adj_Mode" = PosteriorMode)
  # SREP ESTIMATE FOR TARGET STOCKS
targetsAll_srep <- cbind(targets3_srep, dsrep$deripost_summary$SREP_tar) |> 
  rename("SREP_tar_mean" = Mean, "SREP_tar_median" = Median,
    "SREP_tar_LQ_5" = LQ_5, "SREP_tar_UQ_95" = UQ_95, "SREP_tar_Stocknum" = Stock,
    "SREP_tar_Mode" = PosteriorMode)

#### SMAX MODEL RESULTS ####
  # SMSY Estimate for TARGET STOCKS
targets1_smax <- cbind(targets, dsmax$deripost_summary$SMSY) |> 
  rename("SMSY_mean" = Mean, "SMSY_median" = Median,
    "SMSY_LQ_5" = LQ_5, "SMSY_UQ_95" = UQ_95, "SMSY_Stocknum" = Stock,
    "SMSY_Mode" = PosteriorMode)
  # SGEN Estimate for TARGET STOCKS
targets2_smax <- cbind(targets1_smax, dsmax$deripost_summary$SGEN) |> 
  rename("SGEN_mean" = Mean, "SGEN_median" = Median,
    "SGEN_LQ_5" = LQ_5, "SGEN_UQ_95" = UQ_95, "SGEN_Stocknum" = Stock,
    "SGEN_Mode" = PosteriorMode)
  # Marginal Mean for SREP (E)
targets3_smax <- cbind(targets2_smax, dsmax$deripost_summary$SREP_tar_adj) |> 
  rename("SREP_tar_adj_mean" = Mean, "SREP_tar_adj_median" = Median,
    "SREP_tar_adj_LQ_5" = LQ_5, "SREP_tar_adj_UQ_95" = UQ_95, "SREP_tar_adj_Stocknum" = Stock,
    "SREP_tar_adj_Mode" = PosteriorMode)
  # SREP ESTIMATE FOR TARGET STOCKS
targetsAll_smax <- cbind(targets3_smax, dsmax$deripost_summary$SREP_tar) |> 
  rename("SREP_tar_mean" = Mean, "SREP_tar_median" = Median,
    "SREP_tar_LQ_5" = LQ_5, "SREP_tar_UQ_95" = UQ_95, "SREP_tar_Stocknum" = Stock,
    "SREP_tar_Mode" = PosteriorMode)

	# SMSY, SGEN, and SREP from bootstrapping
# bstargets1 <- cbind(targetsAll), # bstargets2 <- cbind(), # pbenchmarks <- cbind()
BS_wide <- BS.dfout %>%
  pivot_wider(
    id_cols = c(Stock, WA, lh), 
    names_from = RP, 
    values_from = c(Value, lwr, upr),
    names_sep = "_"
  )
targetsAll_srep <- targetsAll_srep %>%
	left_join(BS_wide, by = c("Stock_name" = "Stock", "WA", "lh"))

	# SREP from IWAM
# dfout
dfout_wide <- dfout %>%
  pivot_wider(
    id_cols = c(Stock, WA, lh), 
    names_from = RP, 
    values_from = c(Value, lwr, upr),
    names_sep = "_",
    names_prefix = "IWAM_"
  )
targetsAll <- targetsAll %>%
  left_join(dfout_wide, by = c("Stock_name" = "Stock", "WA", "lh"))

# Ricker sigma for SYNOPTIC STOCKS 
  # Re-order Stocks to be in order of Ricker variance
# targetsAll <- cbind(targetsAll, derived_obj$deripost_summary$tauobs) |> 
#   rename("tauobs_mean" = Mean, "tauobs_median" = Median,
#     "tauobs_LQ_5" = LQ_5, "tauobs_UQ_95" = UQ_95, "tauobs_Stocknum" = Stock)
    # tauobs is based on the synoptic sets and will now have different lengths

parken <- ParkenCaseStudy |> 
  rename("SMSYp" = SMSY, "SMSYp_5" = SMSY_5, "SMSYp_95" = SMSY_95) |> 
  rename("SREPp" = SREP, "SREPp_5" = SREP_5, "SREPp_95" = SREP_95)

cols <- viridis(8, alpha=0.9, option = "mako", direction = -1)
options(scipen = 999)

#### Point-wise Benchmark Comparison SREP - BY LOG WA ####
ggplot() +
  
  geom_errorbar(data = parken, aes(x = fct_reorder(Stock, log(WA)), y = SREPp, ymax = SREPp_95, ymin = SREPp_5,
                                       color = "Parken",
                                       width=.1),
) +
  geom_point(data = parken,
             aes(x = fct_reorder(Stock, log(WA)), y = SREPp, color = "Parken")) +
  
  geom_errorbar(data = targetsAll_srep, aes(x = fct_reorder(Stock_name, log(WA)),
                                     y = SREP_tar_median,
                                     ymax = SREP_tar_UQ_95, 
                                     ymin = SREP_tar_LQ_5,
                                 color = "Liermann MCMC Cond.",
                                 width=.1),
                position = position_nudge(+0.1)) +
  geom_point(data = targetsAll_srep,
             position = position_nudge(+0.1),
             aes(x = fct_reorder(Stock_name, log(WA)), y = SREP_tar_median, color = "Liermann MCMC Cond. SREP")) +

  geom_errorbar(data = targetsAll_smax, aes(x = fct_reorder(Stock_name, log(WA)),
                                     y = SREP_tar_median,
                                     ymax = SREP_tar_UQ_95, 
                                     ymin = SREP_tar_LQ_5,
                                 color = "Liermann MCMC Cond. SMAX",
                                 width=.1),
                position = position_nudge(+0.1)) +
  geom_point(data = targetsAll_smax,
             position = position_nudge(+0.1),
             aes(x = fct_reorder(Stock_name, log(WA)), y = SREP_tar_median, color = "Liermann MCMC Cond. SMAX")) +
  
  # geom_errorbar(data = targetsAll, aes(x = fct_reorder(Stock_name, logWA),
                                     # y = SREP_tar_adj_mean,
                                     # ymax = SREP_tar_adj_UQ_95, 
                                     # ymin = SREP_tar_adj_LQ_5,
                                 # color = "Liermann MCMC Marg.",
                                 # width=.1),
                # position = position_nudge(+0.2)) +
  # geom_point(data = targetsAll,
             # position = position_nudge(+0.2),
             # aes(x = fct_reorder(Stock_name, logWA), y = SREP_tar_adj_mean, color = "Liermann MCMC Marg.")) +
  
  # geom_errorbar(data = targetsAll, aes(x = fct_reorder(Stock_name, log(WA)), # Bootstrap
                                     # y = Value_SREP,
                                     # ymax = lwr_SREP, 
                                     # ymin = upr_SREP,
                                 # color = "Liermann Bootstrap",
                                 # width=.1),
                # position = position_nudge(+0.2)) +
  # geom_point(data = targetsAll,
             # position = position_nudge(+0.2),
             # aes(x = fct_reorder(Stock_name, log(WA)), y = Value_SREP, color = "Liermann Bootstrap")) +
  
  geom_errorbar(data = targetsAll, aes(x = fct_reorder(Stock_name, log(WA)), # Bootstrap
                                     y = Value_IWAM_SREP,
                                     ymax = lwr_IWAM_SREP, 
                                     ymin = upr_IWAM_SREP,
                                 color = "IWAM",
                                 width=.1),
                position = position_nudge(+0.2)) +
  geom_point(data = targetsAll,
             position = position_nudge(+0.2),
             aes(x = fct_reorder(Stock_name, log(WA)), y = Value_IWAM_SREP, color = "IWAM")) +
  
  theme_classic() +
  scale_y_continuous(transform = "log", 
                     breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
  ylab(TeX("$S_{REP}$ Estimate")) +
  xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust = 1)) +
  scale_color_manual(name='Model',
                     breaks=c('Parken',
                              # 'RTMB MLE',
							  'IWAM',
                              'Liermann MCMC Cond.',
                              'Liermann MCMC Marg.',
							  'Liermann Bootstrap'),
                     values=c('Parken' = "black",
                              # 'RTMB MLE' = "orange",
							  'IWAM' = 'orange',
                              'Liermann MCMC Cond.' = "skyblue",
                              'Liermann MCMC Marg.' = "darkblue",
							  'Liermann Bootstrap' = 'forestgreen')
							  )
							  
# save
# ggsave("pointwise_example.png", width = 4, height = 3, dpi = 300) # NEED TO CHANGE WIDTH AND HEIGHT for new save

#### Point-wise comparison - SMSY - by logWA ####
ggplot() +
  
  geom_errorbar(data = parken, aes(x = fct_reorder(Stock, log(WA)), y = SMSYp, ymax = SMSYp_95, ymin = SMSYp_5,
                                       color = "Parken",
                                       width=.1),
) +
  geom_point(data = parken,
             aes(x = fct_reorder(Stock, log(WA)), y = SMSYp, color = "Parken")) +
  
  geom_errorbar(data = targetsAll, aes(x = fct_reorder(Stock_name, log(WA)),
                                     y = SMSY_median,
                                     ymax = SMSY_UQ_95, 
                                     ymin = SMSY_LQ_5,
                                 color = "Liermann MCMC Cond.",
                                 width=.1),
                position = position_nudge(+0.1)) +
  geom_point(data = targetsAll,
             position = position_nudge(+0.1),
             aes(x = fct_reorder(Stock_name, log(WA)), y = SMSY_median, color = "Liermann MCMC Cond.")) +
  
  geom_errorbar(data = targetsAll, aes(x = fct_reorder(Stock_name, log(WA)),
                                     y = Value_SMSY,
                                     ymax = lwr_SMSY, 
                                     ymin = upr_SMSY,
                                 color = "Liermann Bootstrap",
                                 width=.1),
                position = position_nudge(+0.2)) +
  geom_point(data = targetsAll,
             position = position_nudge(+0.2),
             aes(x = fct_reorder(Stock_name, log(WA)), y = Value_SMSY, color = "Liermann Bootstrap")) +
  
  theme_classic() +
  scale_y_continuous(transform = "log", 
                     breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
  ylab(TeX("$S_{MSY}$ Estimate")) +
  xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust = 1)) +
  scale_color_manual(name='Model',
                     breaks=c('Parken',
                              # 'RTMB MLE',
                              'Liermann MCMC Cond.',
                              'Liermann MCMC Marg.',
							  'Liermann Bootstrap'),
                     values=c('Parken' = "black",
                              # 'RTMB MLE' = "orange",
                              'Liermann MCMC Cond.' = "skyblue",
                              'Liermann MCMC Marg.' = "darkblue",
							  'Liermann Bootstrap' = 'forestgreen'))
							 


#### Point-wise comparison of SYNOPTIC POPULATIONS ####
Parkentable1 <- read.csv(here::here("DataIn/Parken_Table1n2.csv")) # Test stocks e.g. WCVI stocks
	# Double check ordering of pops
	# Double check if WA is available

synoptic <- WAbase |>
  rename("Stock_name" = Name) # This name can change depending on what data sets are being run
  # SREP ESTIMATE FOR SYNOPTIC POPULATIONS
targetsyn_srep <- cbind(synoptic, dsrep$deripost_summary$SREP) |> 
  rename("SREP_mean" = Mean, "SREP_median" = Median,
    "SREP_LQ_5" = LQ_5, "SREP_UQ_95" = UQ_95, "SREP_Stocknum" = Stock,
    "SREP_Mode" = PosteriorMode)
targetsyn_srep <- cbind(targetsyn_srep, dsrep$deripost_summary$SMSY_r) |> 
  rename("SMSY_mean" = Mean, "SMSY_median" = Median,
    "SMSY_LQ_5" = LQ_5, "SMSY_UQ_95" = UQ_95, "SMSY_Stocknum" = Stock,
    "SMSY_Mode" = PosteriorMode)
	
targetsyn_smax <- cbind(synoptic, dsmax$deripost_summary$SREP_r) |>
	 rename("SREP_mean" = Mean, "SREP_median" = Median,
    "SREP_LQ_5" = LQ_5, "SREP_UQ_95" = UQ_95, "SREP_Stocknum" = Stock,
    "SREP_Mode" = PosteriorMode)
targetsyn_smax <- cbind(targetsyn_smax, dsmax$deripost_summary$SMSY_r) |>
	 rename("SMSY_mean" = Mean, "SMSY_median" = Median,
    "SMSY_LQ_5" = LQ_5, "SMSY_UQ_95" = UQ_95, "SMSY_Stocknum" = Stock,
    "SMSY_Mode" = PosteriorMode)

#### SREP ####
ggplot() +
  
  # geom_errorbar(data = Parkentable1, aes(x = fct_reorder(Stock, log(WA)), y = SREP, ymax = SREPp_95, ymin = SREPp_5,
                                       # color = "Parken",
                                       # width=.1),
# ) +
  geom_point(data = Parkentable1,
             aes(x = fct_reorder(Stock, log(WA)), y = Srep, color = "Parken")) +
  
  geom_errorbar(data = targetsyn_srep, aes(x = fct_reorder(Stock_name, log(WA)),
                                     y = SREP_median,
                                     ymax = SREP_UQ_95, 
                                     ymin = SREP_LQ_5,
                                 color = "Liermann MCMC Cond.",
                                 width=.1),
                position = position_nudge(+0.2)) +
  geom_point(data = targetsyn_srep,
             position = position_nudge(+0.2),
             aes(x = fct_reorder(Stock_name, log(WA)), y = SREP_median, color = "Liermann MCMC Cond. (SREP)")) +
			 
  geom_errorbar(data = targetsyn_smax, aes(x = fct_reorder(Stock_name, log(WA)),
                                     y = SREP_median,
                                     ymax = SREP_UQ_95, 
                                     ymin = SREP_LQ_5,
                                 color = "Liermann MCMC Cond.",
                                 width=.1),
                position = position_nudge(-0.2)) +
  geom_point(data = targetsyn_smax,
             position = position_nudge(-0.2),
             aes(x = fct_reorder(Stock_name, log(WA)), y = SREP_median, color = "Liermann MCMC Cond. (SMAX)")) +
  
  # geom_errorbar(data = targetsyn, aes(x = fct_reorder(Stock_name, log(WA)), # Bootstrap
                                     # y = Value_IWAM_SREP,
                                     # ymax = lwr_IWAM_SREP, 
                                     # ymin = upr_IWAM_SREP,
                                 # color = "IWAM",
                                 # width=.1),
                # position = position_nudge(+0.2)) +
  # geom_point(data = targetsAll,
             # position = position_nudge(+0.2),
             # aes(x = fct_reorder(Stock_name, log(WA)), y = Value_IWAM_SREP, color = "IWAM")) +
  
  theme_classic() +
  scale_y_continuous(transform = "log", 
                     breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
  ylab(TeX("$S_{REP}$ Estimate")) +
  xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust = 1)) +
  scale_color_manual(name='Model',
                     breaks=c('Parken',
                              # 'RTMB MLE',
							  'IWAM',
                              'Liermann MCMC Cond. (SREP)',
							  'Liermann MCMC Cond. (SMAX)',
                              'Liermann MCMC Marg.',
							  'Liermann Bootstrap'),
                     values=c('Parken' = "black",
                              # 'RTMB MLE' = "orange",
							  'IWAM' = 'orange',
                              'Liermann MCMC Cond. (SREP)' = "skyblue",
							  'Liermann MCMC Cond. (SMAX)' = "forestgreen",
                              'Liermann MCMC Marg.' = "darkblue",
							  'Liermann Bootstrap' = 'forestgreen')
							  )



#### SMSY ####
ggplot() +
  
  # geom_errorbar(data = Parkentable1, aes(x = fct_reorder(Stock, log(WA)), y = SREP, ymax = SREPp_95, ymin = SREPp_5,
                                       # color = "Parken",
                                       # width=.1),
# ) +
  geom_point(data = Parkentable1,
             aes(x = fct_reorder(Stock, log(WA)), y = Smsy, color = "Parken")) +
  
  geom_errorbar(data = targetsyn_srep, aes(x = fct_reorder(Stock_name, log(WA)),
                                     y = SMSY_median,
                                     ymax = SMSY_UQ_95, 
                                     ymin = SMSY_LQ_5,
                                 color = "Liermann MCMC Cond. (SREP)",
                                 width=.1),
                position = position_nudge(+0.2)) +
  geom_point(data = targetsyn_srep,
             position = position_nudge(+0.2),
             aes(x = fct_reorder(Stock_name, log(WA)), y = SMSY_median, color = "Liermann MCMC Cond. (SREP)")) +
  
  geom_errorbar(data = targetsyn_smax, aes(x = fct_reorder(Stock_name, log(WA)),
                                     y = SMSY_median,
                                     ymax = SMSY_UQ_95, 
                                     ymin = SMSY_LQ_5,
                                 color = "Liermann MCMC Cond.",
                                 width=.1),
                position = position_nudge(-0.2)) +
  geom_point(data = targetsyn_smax,
             position = position_nudge(-0.2),
             aes(x = fct_reorder(Stock_name, log(WA)), y = SMSY_median, color = "Liermann MCMC Cond. (SMAX)")) +
  
  # geom_errorbar(data = targetsyn, aes(x = fct_reorder(Stock_name, log(WA)), # Bootstrap
                                     # y = Value_IWAM_SREP,
                                     # ymax = lwr_IWAM_SREP, 
                                     # ymin = upr_IWAM_SREP,
                                 # color = "IWAM",
                                 # width=.1),
                # position = position_nudge(+0.2)) +
  # geom_point(data = targetsAll,
             # position = position_nudge(+0.2),
             # aes(x = fct_reorder(Stock_name, log(WA)), y = Value_IWAM_SREP, color = "IWAM")) +
  
  theme_classic() +
  scale_y_continuous(transform = "log", 
                     breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
  ylab(TeX("$S_{MSY}$ Estimate")) +
  xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust = 1)) +
  scale_color_manual(name='Model',
                     breaks=c('Parken',
                              # 'RTMB MLE',
							  'IWAM',
                              'Liermann MCMC Cond. (SREP)',
							  'Liermann MCMC Cond. (SMAX)',
                              'Liermann MCMC Marg.',
							  'Liermann Bootstrap'),
                     values=c('Parken' = "black",
                              # 'RTMB MLE' = "orange",
							  'IWAM' = 'orange',
                              'Liermann MCMC Cond. (SREP)' = "skyblue",
							  'Liermann MCMC Cond. (SMAX)' = "forestgreen",
                              'Liermann MCMC Marg.' = "darkblue",
							  'Liermann Bootstrap' = 'forestgreen')
							  )



#### Linear Regression: Liermann vs. Parken model ###################################################################################################

# NEW VERSION: GGPLOT2
	# - Make a dataframe with:
		# - simulated x values of WA
		# - slopes (2), intercepts (2) (watch out for shifted WA)
		# - calculate the line as y = mx + b
		# - confidence e.g. quantiles?
	# - Plot scatter plot of points
	# - geom_ribbon with upper and lowers
	# - geom line
# Need posterior predictive means of the linear regression - of logSREP_tar w/ quantiles
# Plotted against simulated WA
# derived_obj$deripost_summary$SREP_line_ocean and SREP_line_ocean including their marginal variants
# Plot lineWA against SREP 
baselinescatter <- data.frame(
	WAbaseshifted = WAbase$logWAshifted,
	WAreal = WAbase$WA,
	WAtarget = dat$WAin$WA, # ONLY WORKS IF SAME NUMBER ****
	WAlh = WAbase$lh,
	WAlhtarget = dat$WAin$lh, # 0 and 1's # ONLY WORKS IF SAME NUMBER ****
	SREP = dpars_srep$SREP$Median, # ?
	SREP_tar = dpars_srep$SREP_tar$Median # ONLY WORKS IF SAME NUMBER ****
) # dat$mean_logWA

posteriorline <- data.frame(
	lineWA = dat$lineWA, # vector of 72
	median_stream_line = derived_obj$deripost_summary$logSREP_line_stream$Median,
	median_ocean_line = derived_obj$deripost_summary$logSREP_line_ocean$Median,
	lower_stream_line = derived_obj$deripost_summary$logSREP_line_stream$LQ_5 ,
	upper_stream_line = derived_obj$deripost_summary$logSREP_line_stream$UQ_95,
	lower_ocean_line = derived_obj$deripost_summary$logSREP_line_ocean$LQ_5,
	upper_ocean_line = derived_obj$deripost_summary$logSREP_line_ocean$UQ_95
)

# On Logged scales
ggplot(data = baselinescatter, aes(x = WAbaseshifted, y = log(SREP), color = WAlh)) + # base observation scatter plot ???????????????????????????????
	geom_point(alpha = 0.5) + # colours, etc. for observation scatter plot
	scale_color_manual(values = c('stream' = 'forestgreen', 'ocean' = 'skyblue'), guide = "none") +
	# geom_point(data = baselinescatter, aes(x = , y = )inherit.aes = FALSE) +
	# geom_point(inherit.aes = FALSE, data = targetscatter, aes(x = WAtarget, y = SREPtarget)) + # Plot posterior of SREP or posterior predictive?
	geom_line(data = posteriorline, aes(x = lineWA, y = median_stream_line), color = "forestgreen", size = 1) + 
	geom_line(data = posteriorline, aes(x = lineWA, y = median_ocean_line), color = "skyblue", size = 1) +
	geom_ribbon(data = posteriorline, aes(x = lineWA, ymin = lower_stream_line, ymax = upper_stream_line), fill = "forestgreen", alpha = 0.3, inherit.aes = FALSE) + 
	geom_ribbon(data = posteriorline, aes(x = lineWA, ymin = lower_ocean_line, ymax = upper_ocean_line), fill = "skyblue", alpha = 0.3, inherit.aes = FALSE) + 
	theme_classic() + 
	# scale_x_log10() +
	# scale_y_log10() + 
	labs(x = "Log Mean Centered Accessible Watershed Area", y = "LogSREP (Spawners at Replacement)")

# Processing
# Do I need to add back in mean_logWA? for the shifted WA's?
	# If I convert to real-scale for WA then exp(lineWA + mean_logWA)
eb01 <- derived_obj$deripost_full$b0[,1]
eb02 <- derived_obj$deripost_full$b0[,2]
ebWA1 <- derived_obj$deripost_full$bWA[,1]
ebWA2 <- derived_obj$deripost_full$bWA[,2]

realWAline <- exp(seq(log(min(dat$WAbase$WA, na.rm=TRUE)), log(max(dat$WAbase$WA, na.rm=TRUE)), length.out = 72))
posteriorline$logWAline <- log(realWAline)

unc1 <- eb01 -  ebWA1 * dat$mean_logWA
unc2 <- (eb01 + eb02) - (ebWA1 + ebWA2) * dat$mean_logWA

elogSREP1 <- outer(unc1, rep(1, length(posteriorline$logWAline))) + outer(ebWA1, posteriorline$logWAline)
elogSREP2 <- outer(unc2, rep(1, length(posteriorline$logWAline))) + outer(ebWA1 + ebWA2, posteriorline$logWAline)
eSREP1 <- exp(elogSREP1)
eSREP2 <- exp(elogSREP2)

posteriorline$eSREP1median <- apply(eSREP1, 2, median)
posteriorline$eSREP1_5 <- apply(eSREP1, 2, quantile, probs = c(0.05, 0.95))[1,]
posteriorline$eSREP1_95 <- apply(eSREP1, 2, quantile, probs = c(0.05, 0.95))[2,]
posteriorline$eSREP2median <- apply(eSREP2, 2, median)
posteriorline$eSREP2_5 <- apply(eSREP2, 2, quantile, probs = c(0.05, 0.95))[1,]
posteriorline$eSREP2_95 <- apply(eSREP2, 2, quantile, probs = c(0.05, 0.95))[2,]

hdi_list1 <- apply(eSREP1, 2, function(x) {
  out <- hdi(x, credMass = 0.95)
  c(lower = out[1], upper = out[2])
})
hdi_list2 <- apply(eSREP2, 2, function(x) {
  out <- hdi(x, credMass = 0.95)
  c(lower = out[1], upper = out[2])
})
posteriorline$hdi_lo1 <- hdi_list1[1, ]
posteriorline$hdi_hi1 <- hdi_list1[2, ]
posteriorline$hdi_lo2 <- hdi_list2[1, ]
posteriorline$hdi_hi2 <- hdi_list2[2, ]

ggplot(data = baselinescatter, aes(x = WAreal, y = SREP, color = WAlh)) + # base observation scatter plot ???????????????????????????????
	geom_point(alpha = 0.8, size = 3) + # colours, etc. for observation scatter plot
	scale_color_manual(name = "Life History", 
		values = c('stream' = 'forestgreen', 'ocean' = 'skyblue'), 
		labels = c('stream' = 'Stream Type', 'ocean' = 'Ocean Type'),
		guide = guide_legend(override.aes = list(size = 3))) + # "none"
	# geom_point(inherit.aes = FALSE, data = targetscatter, aes(x = WAtarget, y = SREPtarget)) + # Plot posterior of SREP or posterior predictive?
	# new_scale_color() +
	
	# geom_point(data = baselinescatter, aes(x = WAtarget, y = SREP_tar, alpha = 0.2), inherit.aes = FALSE) +
	# scale_color_manual(values = c('0' = 'forestgreen', '1' = 'skyblue'), guide = "none") +
	# THESE ARE MEDIAN VALUES e.g. ON THE LINE
	# COULD ALSO PLOT THE POSTERIOR PREDICTIVE SREP_tar VALUES?
	
	geom_line(data = posteriorline, aes(x = exp(logWAline), y = eSREP1median), 
		color = "forestgreen", size = 1, inherit.aes = FALSE) + 
	geom_ribbon(data = posteriorline, aes(x = exp(logWAline), ymin = eSREP1_5, ymax = eSREP1_95), 
		fill = "forestgreen", alpha = 0.3, inherit.aes = FALSE) + 
	# geom_ribbon(data = posteriorline, aes(x = exp(logWAline), ymin = hdi_lo1, ymax = hdi_hi1), 
		# fill = "forestgreen", alpha = 0.2, inherit.aes = FALSE) + # HDPI
	
	geom_line(data = posteriorline, aes(x = exp(logWAline), y = eSREP2median), 
		color = "skyblue", size = 1, inherit.aes = FALSE) + 
	geom_ribbon(data = posteriorline, aes(x = exp(logWAline), ymin = eSREP2_5, ymax = eSREP2_95), 
		fill = "skyblue", alpha = 0.3, inherit.aes = FALSE) +
	# geom_ribbon(data = posteriorline, aes(x = exp(logWAline), ymin = hdi_lo2, ymax = hdi_hi2), 
		# fill = "skyblue", alpha = 0.2, inherit.aes = FALSE) + # HDPI
	
	theme_classic() + 
	theme(legend.position = "none") + 
	
	theme(legend.position = c(0.02, 0.98),  # x, y coordinates (0-1 scale)
		legend.justification = c(0, 1),    # anchor at top-left of legend box
		# legend.background = element_rect(fill = "white", color = "black", size = 0.5),
		legend.title = element_text(size = 12, face = "bold"),
		legend.text = element_text(size = 10)) +
	
	scale_x_log10(labels = function(x) format(x, scientific = FALSE, trim = TRUE, big.mark = ",")) +
	scale_y_log10(labels = function(x) format(x, scientific = FALSE, trim = TRUE, big.mark = ",")) + 
	
	labs(x = "Accessible Watershed Area", y = "SREP (Spawners at Replacement)")

ggsave("figure.png", width = 4, height = 3, dpi = 300) # REALLY IMPORTANT VISUAL



#### Bar plot comparison of SYNOPTIC values of SREP #################################################################################################
  # Compare Parken and Liermann estimates of SREP for the SYNOPTIC STOCKS
tempSREPpars <- dpars_srep$SREP |> 
  rename("SREP_stock_temp" = Stock)
bardf <- cbind(Parkentable1, tempSREPpars)
bardf_long <- bardf %>%
  pivot_longer(cols = c(Srep, Mean), names_to = "Category", values_to = "Value")

# Now plot using a single geom_bar()
bcompare <- ggplot(bardf_long, aes(x = Stock, y = Value, fill = Category)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  
  theme_classic() +
  scale_fill_manual(
    name = "Model",  # Custom legend title
    values = c("Srep" = "black", "Mean" = "skyblue"),  # Custom colors
    labels = c("Srep" = "Parken", "Mean" = "Liermann")  # Custom category names
  ) +
  labs(x = "Stock",
       y = expression(S[REP])) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
bcompare



#####################################################################################################################################################
#### OLD PLOTTING CODE - DEPRECIATED ################################################################################################################
#####################################################################################################################################################

# OLD VERSION in BASE R PLOT - LINEAR REGRESSION PLOT
    # Step 1. Is the data prepared.
options(scipen = 5) # for ticks without sci. notation
col.use <- NA
for(i in 1:length(WAbase$lh)) {
  if (WAbase$lh[i]== 'stream') col.use[i] <- "forestgreen"  # stream
  else col.use[i] <- "dodgerblue3" # ocean
}

    # Step 2. PLOT base points
plot(y = dpars$SREP$Mean, x = WAbase$WA, pch = 20, 
     # col = ifelse(WAbase$Name == "Chehalis" | WAbase$Name == "Skagit", 'red', col.use), 
     col = col.use, cex = 1.5,
     xlab = expression("Accessible Watershed Area, km"^2), 
     ylab = expression(S[REP]), log = 'xy',
     xlim = c(50,200000) , ylim = c(200,2000000)
  )

# ADD Parken points of SYNOPTIC SET from Table 1 (Parken et al. 2006)
points(y = Parkentable1$Srep, x = Parkentable1$WA, pch = 20, col = alpha('black', 0.5))
    # Shown here because they have changed slighly in WA between model versions?

# Liermann points for TARGET ESTIMATES
    # Can change here betwen Median, Mean, and to Marginal vs. Conditional
points(y = dpars$SREP_tar$Median, x = WAin$WA, pch = 20, 
  col = alpha(ifelse(WAin$Stock == "Coldwater"| WAin$Stock == "Deadman", 'red', 'skyblue'), 0.5))
  # Can I write this such that for the 'red' points, the alpha is also higher?

    # Step 3. LINES
sum_pars <- summary(fitstan)
bWA1 <- sum_pars$summary[,1][3] 
bWA2 <- sum_pars$summary[,1][4] + bWA1
b01 <- sum_pars$summary[,1][1] 
b01 <- b01 - mean_logWA*bWA1 # To deal with centered/shifted watershed areas
b02 <- sum_pars$summary[,1][2] + sum_pars$summary[,1][1] - mean_logWA*bWA2

simWA <- seq(2, 13, 0.5)
Preds <- b01 + simWA*bWA1
Predso <- b02 + simWA*bWA2
lines(x = exp(simWA), y = exp(Preds), col = alpha("forestgreen", 0.5), lwd = 2, lty = 1)
lines(x = exp(simWA), y = exp(Predso), col = alpha("dodgerblue3", 0.5), lwd = 2, lty = 1)

    # Step 4. Error polygons
# Calculate quantile - uppers and lowers
  # pred_lnSMSY IS THE TARGET's
  # predlnWA should be using: WAin$logWAshifted_t
SREPline_stream <- derived_obj$deripost_summary$SREP_line_stream |> 
  rename("line_stocks" = Stock,
    "s_mean" = Mean,
    "s_median" = Median,
    "s_LQ_5" = LQ_5,
    "s_UQ_95" = UQ_95)
SREPline_ocean <- derived_obj$deripost_summary$SREP_line_ocean |> 
  rename("line_stocko" = Stock,
    "o_mean" = Mean,
    "o_median" = Median,
    "o_LQ_5" = LQ_5,
    "o_UQ_95" = UQ_95)
lineWA <- cbind(dat$lineWA, SREPline_stream, SREPline_ocean) # NEED NEW LINE VALUES

up_S <- lineWA$s_UQ_95
lo_S <- lineWA$s_LQ_5
up_O <- lineWA$o_UQ_95
lo_O <- lineWA$o_LQ_5

polygon(x=c(exp(lineWA$`dat$lineWA` + mean_logWA), exp(rev(lineWA$`dat$lineWA` + mean_logWA))), 
        y=c(up_S, rev(lo_S)), 
        col=rgb(0,0.4,0, alpha=0.2), border=NA)

polygon(x=c(exp(lineWA$`dat$lineWA` + mean_logWA), exp(rev(lineWA$`dat$lineWA` + mean_logWA))), 
        y=c(up_O, rev(lo_O)), 
        col=rgb(0,0.2,0.4, alpha=0.2), border=NA) 

    # Step 5. Grab Parken estimates for the line and add as y = mx + b
    # From Table 4. Srep Habitat Model
Preds_Parken <- 3.89 + simWA*0.693 + (0.240/2) # Stream-type
Predso_Parken <- 3.52 + simWA*0.878 + (0.133/2) # Ocean-type
lines(x = exp(simWA), y = exp(Preds_Parken), col = alpha("forestgreen", 0.5), lwd = 2, lty = 2)
lines(x = exp(simWA), y = exp(Predso_Parken), col = alpha("dodgerblue3", 0.5), lwd = 2, lty = 2)

    # Step 6. Add text to describe model equations
# q: Based on the plot coded above, how can I add text labels to each point to state the name associated with them?
# a: Use geom_text() to add text labels to the plot.
# geom_text()



#### SR Curves for individual stocks - DEPRECIATED ##################################################################################################
Stks <- unique(srdat$Stocknumber)
NStks <- length(Stks)
par(mfrow=c(5,5), mar=c(2, 2, 1, 0.1) + 0.1, oma=c(3,3,1,1)) # To fit all the plots on one grid
# par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1, oma=c(0,0,0,0)) # To plot each individual stock

# ADDING IWAM estimates - will produce errors if it is not the standard iwam_obj name

Parken_ab <- read.csv(here::here("DataIn/Parken_Table1n2.csv")) 
  # KSR, Andrew, Lewis, Columbia - just different names - but ORDER is the same
  # Needs correction for removal of Skagit and Chehalis
if (!'Skagit' %in% WAbase$Name) Parken_ab <- Parken_ab |> filter( !(Stock == "Skagit")) 
if (!'Chehalis' %in% WAbase$Name) Parken_ab <- Parken_ab |> filter( !(Stock == "Chehalis")) 
Parken_ab <- Parken_ab |> 
  mutate(Stocknumber = as.integer(factor(Stocknumber)) - 1) # Re-order Stocknumber to be ascending

for (i in Stks){
  # Get stocknames and numbers
      # names <- pars %>% dplyr::select ("Name", "Stocknumber") %>% distinct()
      # name <- pars %>% filter (Stocknumber==i) %>% 
      #   dplyr::select ("Name") %>% distinct()
  name <- names$Name[names$Stocknumber == i]
  
  R <- srdat %>% filter (Stocknumber==i) %>% 
    dplyr::select(Rec) 
  S <- srdat %>% filter (Stocknumber==i) %>% 
    dplyr::select(Sp) 
  
  if(name != "Skagit" & name != "KSR" & name != "Humptulips" & name != "Chehalis") 
    plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)+(max(S$Sp)/3)), ylim=c(0,max(R$Rec) ) )
  if(name == "Skagit") 
    plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)*3+(max(S$Sp)*2)), ylim=c(0,max(R$Rec) ) )
  if(name == "KSR") 
    plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,1150), ylim=c(0,max(R$Rec) ) ) # xlim was 500 or 1000
  if(name == "Humptulips") # Requires extra larger xlim
    plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)+10000), ylim=c(0,max(R$Rec) ) )
  if(name == "Chehalis") # Requires extra larger xlim
    plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)+50000), ylim=c(0,max(R$Rec) ) )

  # Get alpha and beta parameters
      # NOTE: can substitute $alpha$Mean if you want to use Mean
      # Median is resistant to monotonic transformations e.g. ln, exp, ... cube
  a <- as.numeric(derived_obj$deripost_summary$logAlpha$Median[derived_obj$deripost_summary$logAlpha$Stock - 1 == i])
  b <- as.numeric(derived_obj$deripost_summary$BETA_r$Median[derived_obj$deripost_summary$BETA_r$Stock - 1 == i])
  
    # ADD IN IWAM ESTIMATES
      # Require running the IWAM model function externally - see IWAMfunctionrun.R
      # Save as a global variable
  if(exists("iwamobj")) {
    aiwam <- iwamobj$modelpars |> filter (Stocknumber==i) |>  
      filter(Param=="logA") %>%
      summarise(A=exp(Estimate)) %>%
      as.numeric()
    Sc <- iwamobj$srdat %>% filter (Stocknumber==i) %>% 
       dplyr::select(scale) %>% distinct() %>% 
       as.numeric()
    biwam <- iwamobj$modelpars %>% filter (Stocknumber==i) %>%
      filter(Param=="logB") %>%
      summarise(B=exp(Estimate)/Sc) %>%
      as.numeric()
    biwam_se <- iwamobj$modelpars %>%
      filter(Stocknumber == i, Param == "logB") %>%
      mutate(B_lower = exp(Estimate - 1.96 * Std..Error) / Sc,
        B_upper = exp(Estimate + 1.96 * Std..Error) / Sc )
  }
  
  # Parken values for skagit - from Parken et al. 2006 Table 2 (Ocean life-histories)
  skagit_alpha <- 7.74
  skagit_beta <- 0.0000657
  RR_skagit <- NA
  SS <- RR <- RR_parken <- RRiwam <- NA
  ap <- Parken_ab$Alpha # vector
  bp <- Parken_ab$Beta # vector
  
  for (j in 1:250){ # Creates a step-wise sample line by which to create a line on
    if ("Skagit" %in% WAbase$Name) if (i!=22 & i!=7) SS[j] <- j*(max(S$Sp)/100) # IF NOT SKAGIT OR KSR - 
      # When Skagit exists
    if ("Skagit" %in% WAbase$Name) if (i==22) SS[j] <- j*(max(S$Sp*3)/100) # If Skagit exists
    
    if(!"Skagit" %in% WAbase$Name) if (i!=7) SS[j] <- j*(max(S$Sp)/100) # If not KSR - When Skagit doesn't exist
    if (i==7) SS[j] <- j*(500/100) # if KSR - could also add: if ("King Salmon" %in% WAbase$Name) 
    
    RR[j] <- exp(a) * SS[j] * exp(-b * SS[j])
    
    RR_parken[j] <- ap[i+1] * SS[j] * exp(-bp[i+1] * SS[j])
    
    if ("Skagit" %in% WAbase$Name) if(i==22) {RR_skagit[j] <- skagit_alpha * SS[j] * exp(-skagit_beta * SS[j])} 
      # Skagit Line based on alpha and beta from Table 1 and 2 from Parken et al. 2006
    
    if(exists("iwamobj")) RRiwam[j] <- aiwam * SS[j] * exp(-biwam * SS[j])
  }
  
  mtext(name, side=3, cex=0.8)
  
  col.use <- "black"
  lines(x=SS, y=RR, col='black') 
  
  # For Skagit, add Parken et al. 2006 model curve
  if ("Skagit" %in% WAbase$Name) if(i==22) lines(x=SS, y=RR_skagit, lty="dashed") # }
  
  # For all stocks, added in Parken et al. 2006 model curve
  lines(x=SS, y=RR_parken, lty="dashed", col="red")

  # For all stocks, add IWAM model curve
  if(exists("iwamobj")) lines(x=SS, y=RRiwam, lty="dashed", col="forestgreen")
  
  # Calculate and Plot VERTICAL LINES FOR SMSY, SMAX, OR SREP
  # SMSY <- derived_obj$deripost_summary$SMSY_r$Median[derived_obj$deripost_summary$SMSY_r$Stock - 1 == i]
  # SMSY_ul <- derived_obj$deripost_summary$SMSY_r$UQ_95[derived_obj$deripost_summary$SMSY_r$Stock - 1 == i]
  # SMSY_ll <- derived_obj$deripost_summary$SMSY_r$LQ_5[derived_obj$deripost_summary$SMSY_r$Stock - 1 == i]
  SMAX <- 1/derived_obj$deripost_summary$BETA_r$Median[derived_obj$deripost_summary$BETA_r$Stock - 1 == i]
  SMAX_ul <- 1/derived_obj$deripost_summary$BETA_r$UQ_95[derived_obj$deripost_summary$BETA_r$Stock - 1 == i] # Rev.
  SMAX_ll <- 1/derived_obj$deripost_summary$BETA_r$LQ_5[derived_obj$deripost_summary$BETA_r$Stock - 1 == i] # Rev.
  # SREP <- derived_obj$deripost_summary$SREP$Median[derived_obj$deripost_summary$SREP$Stock - 1 == i]
  # SREP_ul <- derived_obj$deripost_summary$SREP$UQ_95[derived_obj$deripost_summary$SREP$Stock - 1 == i]
  # SREP_ll <- derived_obj$deripost_summary$SREP$LQ_5[derived_obj$deripost_summary$SREP$Stock - 1 == i]
  if(exists("iwamobj")) SMAXiwam <- 1/biwam
  if(exists("iwamobj")) SMAXiwam_ll <- 1/(biwam_se$B_lower)
  if(exists("iwamobj")) SMAXiwam_ul <- 1/(biwam_se$B_upper)
  
  # abline(v = SMSY, col=col.use, lty='dotted')
  abline(v = SMAX, col=col.use, lty='dotted')
  # abline(v = SREP, col=col.use, ly='dotted')
  
  # CI' shaded polygons - Repeat for SREP or SMSY if desired
  # IWAM SMAX CI
  if(exists("iwamobj")) polygon(x=c(SMAXiwam_ul, SMAXiwam_ll, SMAXiwam_ll, SMAXiwam_ul),
        y=c(-10000,-10000,10000+max(R$Rec),10000+max(R$Rec)),
        col=rgb(0,0.4,0, alpha=0.1), border=NA )
  
  # polygon(x=c(SMSY_ul, SMSY_ll, SMSY_ll, SMSY_ul), 
  #         y=c(-10000,-10000,10000+max(R$Rec),10000+max(R$Rec)), 
  #         col=grey(0.8, alpha=0.4), border=NA )

  polygon(x=c(SMAX_ul, SMAX_ll, SMAX_ll, SMAX_ul), 
        y=c(-10000,-10000,10000+max(R$Rec),10000+max(R$Rec)), 
        col=grey(0.8, alpha=0.4), border=NA )
  
  # Parken Smsy Estimate (vert. line) from Table 1/2 Parken et al. 2006
  # Parken_smsy <- Parken_ab$Smsy[Parken_ab$Stocknumber == Parken_ab$Stocknumber[i+1]]
  # abline(v = Parken_smsy, col="red", lty='dotted')

  Parken_smax <- 1 / Parken_ab$Beta[Parken_ab$Stocknumber == Parken_ab$Stocknumber[i+1]]
  abline(v = Parken_smax, col="red", lty='dotted')

  # Parken_srep <- Parken_ab$Srep[Parken_ab$Stocknumber == Parken_ab$Stocknumber[i+1]]
  # abline(v = Parken_srep, col="red", lty='dotted')
  
  if(exists("iwamobj")) IWAM_smax <- 1 / biwam
  if(exists("iwamobj")) abline(v = IWAM_smax, col="forestgreen", lty='dotted')
}

# Add an GLOBAL figure axis label across par()
  # x = Spawners, y = Recruitment
mtext("Spawners", side = 1, line = 1, outer = TRUE, cex = 1.3)
mtext("Recruitment", side = 2, line  = 1, outer = TRUE, cex = 1.3, las = 0)



# Compare Deviation of Marginal and Conditional Means ####
edeviation <- data.frame(re = derived_obj$deripost_summary$logE_tar_adj,
                      tar = derived_obj$deripost_summary$logE_tar,
                      dev = derived_obj$deripost_summary$logE_tar_adj$Median - derived_obj$deripost_summary$logE_tar$Median,
                      og_re = derived_obj$deripost_summary$logE_re$Median,
                      name = WAin$Stock)
  # Note: logE_re is not used - unsure what to plot this against

ggplot(edeviation, aes(x = name, y = dev)) +
  geom_hline(yintercept = 0, color = "gray") +
  geom_point(size = 4, color = "darkred") +
  labs(y = "Deviation in predicted E (More positive is\n more random effect)", x = "Stock Name") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#### Point-wise Benchmark Comparison - BY RICKER VARIANCE - ONLY for re-predicting SYNOPTIC STOCKS ##################################################
# full <- cbind(targetsAll, parken$Stock, parken$SREPp, parken$SREPp_5, parken$SREPp_95)
# 
# ggplot() +
#   
#   # Re-written as the Parken model is included in WAin - bound to rtmb MLE version
#   geom_errorbar(data = full, aes(x = fct_reorder(parken$Stock, tauobs_mean), 
#                                  y = parken$SREPp, ymax = parken$SREPp_5, ymin = parken$SREPp_95,
#                                  color = "Parken", width=.1)) +
#   geom_point(data = full, aes(x = fct_reorder(parken$Stock, tauobs_mean), 
#                               y = parken$SREPp, color = "Parken")) +
#   
#   # Add in LIERMANN from Liermann_RTMB_model.R as a global object
#   geom_errorbar(data = full, aes(x = fct_reorder(Stock_name, tauobs_mean),
#                                  y = E_tar_median,
#                                  ymax = E_tar_UQ_95, 
#                                  ymin = E_tar_LQ_5,
#                                  color = "Liermann MCMC",
#                                  width=.1),
#                 position = position_nudge(+0.2)) +
#   geom_point(data = full,
#              position = position_nudge(+0.2),
#              aes(x = fct_reorder(Stock_name, tauobs_mean), 
#                  y = E_tar_mean, color = "Liermann MCMC")) +
#   
#   theme_classic() +
#   scale_y_continuous(transform = "log", 
#                      breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
#   ylab(TeX("$S_{REP}$ Estimate")) +
#   xlab("Stock Name") + 
#   theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust = 1)) +
#   scale_color_manual(name='Model',
#                      breaks=c('Parken',
#                               # 'RTMB MLE',
#                               'Liermann MCMC'),
#                      values=c('Parken' = "black",
#                               # 'RTMB MLE' = "orange",
#                               'Liermann MCMC' = "skyblue"))



#### Testing deviations when re-predicting SYNOPTIC STOCKS ##########################################################################################
smsy_deviation <- derived_obj$deripost_summary$SMSY_r$Median - derived_obj$deripost_summary$SMSY$Median
smax_deviation <- (1/derived_obj$deripost_summary$BETA_r$Median) - (1/derived_obj$deripost_summary$BETA$Median)
testdf <- data.frame(name = WAin$Name, smsy_deviation = smsy_deviation, smax_deviation = smax_deviation)

# Plot of deviations in SMSY
ggplot(testdf, aes(x = name, y = smsy_deviation)) +
  geom_hline(yintercept = 0, color = "gray") +
  geom_point(size = 4, color = "darkred") +
  labs(y = "Deviation in SMSY", x = "Name") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(testdf, aes(x = name, y = smax_deviation)) +
  geom_hline(yintercept = 0, color = "gray") +
  geom_point(size = 4, color = "darkred") +
  labs(y = "Deviation in SMAX", x = "Name") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# As a regression plot
smax_r <- (1/derived_obj$deripost_summary$BETA_r$Median)
smax <- (1/derived_obj$deripost_summary$BETA$Median)
testdf2 <- data.frame(smax_r = smax_r, smax = smax)

ggplot(data = testdf2, aes(x = smax_r, y = smax)) +
  geom_point() + 
  labs(y = "Ricker Smax", x = "Predicted Smax") + 
  theme_classic() + 
  geom_abline(intercept = 0, slope = 1, color = "gray")

