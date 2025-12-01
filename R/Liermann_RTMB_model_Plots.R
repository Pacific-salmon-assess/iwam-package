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
# In script_A.R
# save(my_object, file = "my_object.RData")
# In script_B.R
# load("my_object.RData")

# add stocknames - could add them in to each object of derived_obj?
# Stocknames <- WAin$Stock
# Srep_example <- cbind(derived_obj$deripost_summary$SREP_tar_adj, Stocknames) |> 
#   Srep_example[c(1, 7, 2, 3, 4, 5, 6)]

# Prepare/load datasets for plotting ####
Parkentable1 <- read.csv(here::here("DataIn/Parken_Table1n2.csv")) # Test stocks e.g. WCVI stocks
ParkenCaseStudy <- read.csv(here::here("DataIn/Parken_evalstocks.csv")) # Case study stocks

dpars <- derived_obj$deripost_summary

targets <- WAin |>
  rename("Stock_name" = Stock) # This name can change depending on what data sets are being run

  # SMSY Estimate for TARGET STOCKS
targets1 <- cbind(targets, derived_obj$deripost_summary$SMSY) |> 
  rename("SMSY_mean" = Mean, "SMSY_median" = Median,
    "SMSY_LQ_5" = LQ_5, "SMSY_UQ_95" = UQ_95, "SMSY_Stocknum" = Stock,
    "SMSY_Mode" = PosteriorMode)
  # SGEN Estimate for TARGET STOCKS
targets2 <- cbind(targets1, derived_obj$deripost_summary$SGEN) |> 
  rename("SGEN_mean" = Mean, "SGEN_median" = Median,
    "SGEN_LQ_5" = LQ_5, "SGEN_UQ_95" = UQ_95, "SGEN_Stocknum" = Stock,
    "SGEN_Mode" = PosteriorMode)
  # Marginal Mean for SREP (E)
targets3 <- cbind(targets2, derived_obj$deripost_summary$SREP_tar_adj) |> 
  rename("SREP_tar_adj_mean" = Mean, "SREP_tar_adj_median" = Median,
    "SREP_tar_adj_LQ_5" = LQ_5, "SREP_tar_adj_UQ_95" = UQ_95, "SREP_tar_adj_Stocknum" = Stock,
    "SREP_tar_adj_Mode" = PosteriorMode)
  # SREP ESTIMATE FOR TARGET STOCKS
targetsAll <- cbind(targets3, derived_obj$deripost_summary$SREP_tar) |> 
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
targetsAll <- targetsAll %>%
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

# Point-wise Benchmark Comparison SREP - BY LOG WA ####
ggplot() +
  
  geom_errorbar(data = parken, aes(x = fct_reorder(Stock, log(WA)), y = SREPp, ymax = SREPp_95, ymin = SREPp_5,
                                       color = "Parken",
                                       width=.1),
) +
  geom_point(data = parken,
             aes(x = fct_reorder(Stock, log(WA)), y = SREPp, color = "Parken")) +
  
  geom_errorbar(data = targetsAll, aes(x = fct_reorder(Stock_name, log(WA)),
                                     y = SREP_tar_median,
                                     ymax = SREP_tar_UQ_95, 
                                     ymin = SREP_tar_LQ_5,
                                 color = "Liermann MCMC Cond.",
                                 width=.1),
                position = position_nudge(+0.1)) +
  geom_point(data = targetsAll,
             position = position_nudge(+0.1),
             aes(x = fct_reorder(Stock_name, log(WA)), y = SREP_tar_median, color = "Liermann MCMC Cond.")) +
  
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

# Point-wise comparison - SMSY - by logWA ####
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
							  
# Point-wise comparison - SGEN ####



# Point-wise Benchmark Comparison - UNSORTED ####
# ggplot() +
#   
#   # Re-written as the Parken model is included in WAin - bound to rtmb MLE version
#   geom_errorbar(data = parken, aes(x = Stock, y = SREPp, ymax = SREPp_95, ymin = SREPp_5,
#                                        color = "Parken",
#                                        width=.1),
#                 # position = position_nudge(-0.4)
#     ) +
#   geom_point(data = parken,
#              # position = position_nudge(-0.4),
#              aes(x = Stock, y = SREPp, color = "Parken")) +
#   
#   # Add in LIERMANN from Liermann_RTMB_model.R as a global object
#   geom_errorbar(data = targetsAll, aes(x = Stock_name,
#                                      y = E_tar_median,
#                                      ymax = E_tar_UQ_95, 
#                                      ymin = E_tar_LQ_5,
#                                  color = "Liermann MCMC",
#                                  width=.1),
#                 position = position_nudge(+0.2)) +
#   geom_point(data = targetsAll,
#              position = position_nudge(+0.2),
#              aes(x = Stock_name, y = E_tar_median, color = "Liermann MCMC")) +
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



# Linear Regression: Liermann vs. Parken model ####

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
	SREP = dpars$SREP$Median, # ?
	SREP_tar = dpars$SREP_tar$Median # ONLY WORKS IF SAME NUMBER ****
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


# OLD VERSION in BASE R PLOT
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



# Bar plot comparison of SYNOPTIC values of SREP ####
  # Compare Parken and Liermann estimates of SREP for the SYNOPTIC STOCKS
tempSREPpars <- dpars$SREP |> 
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



# SR Curves for individual stocks - NOT FINISHED ####
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
# Point-wise Benchmark Comparison - BY RICKER VARIANCE - ONLY for re-predicting SYNOPTIC STOCKS ####
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



# Testing deviations when re-predicting SYNOPTIC STOCKS ####
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

