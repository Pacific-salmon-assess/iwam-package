# Library
library(tidyverse)
library(ggplot2)
library(patchwork)
library(latex2exp)
# set.seed(1)
# To create the sequence and dnorm for plotting illustrative distributions
	# This is repeated alot
make_dnorm_dist <- function(mu, sigma, n = 500) {
  x <- seq(mu - 4*sigma, mu + 4*sigma, length.out = n)
  list(seq = x, dnorm = dnorm(x, mean = mu, sd = sigma))
}

# Create a distribution seq and dnorm df for New Alpha
Ricprior = c(1, 0.3)
# RicAseq <- seq(Ricprior[1] - 4*Ricprior[2], Ricprior[1] + 4*Ricprior[2], length.out = 500)
# RicAdnorm <- dnorm(RicAseq, mean = Ricprior[1], sd = Ricprior[2])
RicA <- make_dnorm_dist(mu = Ricprior[1], sigma = Ricprior[2])
RicAdf <- data.frame(seq = RicA$seq, dnorm = RicA$dnorm)


# Get draws of posterior and calculate change in b0
# poplen <- ncol(dsmax$deripost_full$Alpha_tar_adj) # 25 stock predictions to be made
# logSMAX_sd <- dsmax$deripost_summary$logSMAX_sd$Mean
b0 <- dsmax$deripost_full$b0
bWA <- dsmax$deripost_full$bWA
Alpha0 <- exp(dsmax$deripost_full$Alpha0) # logAlpha - NOT log(log(alpha))

post <- cbind(Alpha0, b0[,1], bWA[,1]) # Full posteriors for stream-type and global mean alpha
post2 <- cbind(exp(dsmax$deripost_full$Alpha0 + dsmax$deripost_full$Alpha02), b0[,1] + b0[,2], bWA[,1] + bWA[,2])
covmatrix <- cov(post)
covmatrix2 <- cov(post2)
mu <- apply(post, 2, mean)
mu2 <- apply(post2, 2, mean)

# re-calculate covmatrix for a new value of rho
# prior_rho <- -0.4
# covmatrix[1,2] <- covmatrix[2,1] <- prior_rho * sqrt(covmatrix[1,1])*sqrt(covmatrix[2,2]) 

# Create a distribution seq and dnorm df for Original b0
b0seq <- seq(mean(post[,2]) - 4*sd(post[,2]), mean(post[,2]) + 4*sd(post[,2]), length.out = 500)
b0dnorm <- dnorm(b0seq, mean = mean(post[,2]), sd = sd(post[,2]))
b0df <- data.frame(x = b0seq, y = b0dnorm)

b02seq <- seq(mean(post2[,2]) - 4*sd(post2[,2]), mean(post2[,2]) + 4*sd(post2[,2]), length.out = 500)
b02dnorm <- dnorm(b02seq, mean = mean(post2[,2]), sd = sd(post2[,2]))
b02df <- data.frame(x = b02seq, y = b02dnorm)

# Create a distribution seq and dnorm df for Original Alpha0
alpha0seq <-  seq(mean(post[,1]) - 4*sd(post[,1]), mean(post[,1]) + 4*sd(post[,1]), length.out = 500)
alpha0dnorm <- dnorm(alpha0seq, mean = mean(post[,1]), sd = sd(post[,1]))
alpha0df <- data.frame(x = alpha0seq, y = alpha0dnorm)

alpha02seq <-  seq(mean(post2[,1]) - 4*sd(post2[,1]), mean(post2[,1]) + 4*sd(post2[,1]), length.out = 500)
alpha02dnorm <- dnorm(alpha02seq, mean = mean(post2[,1]), sd = sd(post2[,1]))
alpha02df <- data.frame(x = alpha02seq, y = alpha02dnorm)

# logalphak <- rnorm(1, Ricprior[1], Ricprior[2])
logalphak <- Ricprior[1]
# stream type (base case)
mu_new <- mu[-1] + covmatrix[-1,1]/covmatrix[1,1]*(logalphak - mu[1])
var_new <- covmatrix[-1,-1] - (covmatrix[-1,1, drop = FALSE] / covmatrix[1,1]) %*% covmatrix[1,-1, drop = FALSE]
# bnew <- MASS::mvrnorm(1, mu = mu_new, var_new)

# ocean type (+ additive case)
mu2_new <- mu2[-1] + covmatrix2[-1,1]/covmatrix2[1,1]*(logalphak - mu2[1])
var2_new <- covmatrix2[-1,-1] - (covmatrix2[-1,1, drop = FALSE] / covmatrix2[1,1]) %*% covmatrix2[1,-1, drop = FALSE]
# b2new <- MASS::mvrnorm(1, mu = mu2_new, var2_new)

# Create a distribution seq and dnorm df for New b0
b0newseq <-  seq(mu_new[1] - 4*sqrt(var_new[1,1]), mu_new[1] + 4*sqrt(var_new[1,1]), length.out = 500)
b0newdnorm <- dnorm(b0newseq, mean = mu_new[1], sd = sqrt(var_new[1,1]))
b0newdf <- data.frame(x = b0newseq, y = b0newdnorm)

b02newseq <-  seq(mu2_new[1] - 4*sqrt(var2_new[1,1]), mu2_new[1] + 4*sqrt(var2_new[1,1]), length.out = 500)
b02newdnorm <- dnorm(b02newseq, mean = mu2_new[1], sd = sqrt(var2_new[1,1]))
b02newdf <- data.frame(x = b02newseq, y = b02newdnorm)

# Loop to create new random draws of b0 and Alpha
logalphai <- c()
bnewi <- c()
for (i in 1:100){
	logalphai[i] <- rnorm(1, Ricprior[1], Ricprior[2])
	mu_newi <- mu[-1] + covmatrix[-1,1]/covmatrix[1,1]*(logalphai[i] - mu[1])
	var_newi <- covmatrix[-1,-1] - (covmatrix[-1,1, drop = FALSE] / covmatrix[1,1]) %*% covmatrix[1,-1, drop = FALSE]
	bnewi[i] <- (MASS::mvrnorm(1, mu = mu_newi, var_newi))[1]
}

logalphai2 <- c()
bnewi2 <- c()
for (i in 1:100){
	logalphai2[i] <- rnorm(1, Ricprior[1], Ricprior[2])
	mu_newi2 <- mu2[-1] + covmatrix2[-1,1]/covmatrix2[1,1]*(logalphai2[i] - mu2[1])
	var_newi2 <- covmatrix2[-1,-1] - (covmatrix2[-1,1, drop = FALSE] / covmatrix2[1,1]) %*% covmatrix2[1,-1, drop = FALSE]
	bnewi2[i] <- (MASS::mvrnorm(1, mu = mu_newi2, var_newi2))[1]
}


# SET LIMITS TO MATCH ACROSS FACETS
# x_lim <- range(c(exp(dsmax$deripost_full$Alpha0), RicA$seq))
# y_lim <- range(c(dsmax$deripost_full$b0[,1], b0newseq))

# x_lim <- range(c(exp(dsmax$deripost_full$Alpha0 + dsmax$deripost_full$Alpha02), RicA$seq))
# y_lim <- range(c(dsmax$deripost_full$b0[,1] + dsmax$deripost_full$b0[,2], b02newseq))

x_lim <- c(-0.2, 2.3)
y_lim <- c(7.5, 10.6)

# p_scatter: scatter w/o lines
p_scatter <- ggplot(data.frame(x = exp(dsmax$deripost_full$Alpha0), y = dsmax$deripost_full$b0[,1]),
	aes(x = x, y = y)) +
	geom_point(alpha = 0.2, show.legend = FALSE) + 
	# point for mean of new distribution
	annotate("point", x = mean(RicA$seq), y = mean(b0newseq), color = 'red', size = 4) + 
	geom_point(data = data.frame(x = logalphai, y = bnewi), aes(x = logalphai, y = bnewi), color = "pink", size = 3, alpha = 0.4) +
	xlab(TeX("$log(\\alpha)$")) + 
	ylab(TeX("$b0$ (Intercept of the Regression) for stream-type")) + 
	theme_classic() + 
	theme(legend.position="none") + 
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim, ylim = y_lim) + # to map all 3 figures to the same coordinate grid
	labs(title = "A")

# p_top: logAlpha0 vs. RicAdf
p_top <- ggplot() + 
	geom_line(data = alpha0df, aes(x = alpha0seq, y = alpha0dnorm), linewidth = 1.5) + 
	geom_line(data = RicAdf, aes(x = seq, y = dnorm), color = "red", linewidth = 1.5) + 
	geom_vline(xintercept = mean(alpha0seq), alpha = 0.4, linetype = "dashed", linewidth = 1.2) + # mean
	geom_vline(xintercept = mean(RicA$seq), alpha = 0.4, color = "red", linetype = "dashed", linewidth = 1.2) + # mean
	geom_segment(aes(x = mean(alpha0seq), xend = mean(RicA$seq), y = 2, yend = 2),
		arrow = arrow(length = unit(0.3, "cm")), color = "red", linewidth = 1.2) + # change in mean
	xlab(TeX("$log(\\alpha)$")) + 
	ylab("Density") + 
	theme_classic() + 
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim) +
	labs(title = "B")

# p_right: b0
p_right <- ggplot() + 
	geom_line(data = b0df, aes(x = b0seq, y = b0dnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(data = b0newdf, aes(x = b0newseq, y = b0newdnorm), color = "red", linewidth = 1.2) + # smoothed new smax
	
	geom_vline(xintercept = mean(b0seq), alpha = 0.2, linetype = "dashed", linewidth = 1.2) + 
	geom_vline(xintercept = mean(b0newseq), alpha = 0.2, linetype = "dashed", color = "red", linewidth = 1.2) +
	
	geom_segment(aes(x = mean(b0seq), xend = mean(b0newseq), y = 0.6, yend = 0.6),
		arrow = arrow(length = unit(0.3, "cm")), color = "red", linewidth = 1.2) + 
	
	xlab(TeX("$b0$ (Intercept of the Regression) for stream-type")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_flip(xlim = y_lim) +
	labs(title = "C")
	
	

o_scatter <- ggplot(data.frame(x = exp(dsmax$deripost_full$Alpha0 + dsmax$deripost_full$Alpha02), y = dsmax$deripost_full$b0[,1] + dsmax$deripost_full$b0[,2]),
	aes(x = x, y = y)) +
	geom_point(alpha = 0.2, show.legend = FALSE) + 
	annotate("point", x = mean(RicA$seq), y = mean(b02newseq), color = 'red', size = 4) + 
	geom_point(data = data.frame(x = logalphai2, y = bnewi2), aes(x = logalphai2, y = bnewi2), color = "pink", size = 3, alpha = 0.4) +
	xlab(TeX("$log(\\alpha)$")) + 
	ylab(TeX("$b0$ (Intercept of the Regression) for ocean-type")) + 
	theme_classic() + 
	theme(legend.position="none") + 
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim, ylim = y_lim) + # to map all 3 figures to the same coordinate grid
	labs(title = "A")


# p_top: logAlpha0 vs. RicAdf
o_top <- ggplot() + 
	geom_line(data = alpha02df, aes(x = alpha02seq, y = alpha02dnorm), linewidth = 1.5) + 
	geom_line(data = RicAdf, aes(x = seq, y = dnorm), color = "red", linewidth = 1.5) + 
	geom_vline(xintercept = mean(alpha02seq), alpha = 0.4, linetype = "dashed", linewidth = 1.2) + # mean
	geom_vline(xintercept = mean(RicA$seq), alpha = 0.4, color = "red", linetype = "dashed", linewidth = 1.2) + # mean
	geom_segment(aes(x = mean(alpha02seq), xend = mean(RicA$seq), y = 2, yend = 2),
		arrow = arrow(length = unit(0.3, "cm")), color = "red", linewidth = 1.2) + # change in mean
	xlab(TeX("$log(\\alpha)$")) + 
	ylab("Density") + 
	theme_classic() + 
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim) +
	labs(title = "B")


# p_right: b0
o_right <- ggplot() + 
	geom_line(data = b02df, aes(x = b02seq, y = b02dnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(data = b02newdf, aes(x = b02newseq, y = b02newdnorm), color = "red", linewidth = 1.2) + # smoothed new smax
	geom_vline(xintercept = mean(b02seq), alpha = 0.2, linetype = "dashed", linewidth = 1.2) + 
	geom_vline(xintercept = mean(b02newseq), alpha = 0.2, linetype = "dashed", color = "red", linewidth = 1.2) +
	geom_segment(aes(x = mean(b02seq), xend = mean(b02newseq), y = 0.6, yend = 0.6),
		arrow = arrow(length = unit(0.3, "cm")), color = "red", linewidth = 1.2) + 
	xlab(TeX("$b0$ (Intercept of the Regression) for ocean-type")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_flip(xlim = y_lim) +
	labs(title = "C")

# patchwork setup
(p_top + plot_spacer() + p_scatter + p_right) + 
	plot_layout(ncol = 2, widths  = c(3, 1), heights = c(1, 3)) + 
	plot_annotation(title = "Relationship of logAlpha and the Regression Intercept for Stream-type Populations")

(o_top + plot_spacer() + o_scatter + o_right) + 
	plot_layout(ncol = 2, widths  = c(3, 1), heights = c(1, 3)) + 
	plot_annotation(title = "Relationship of logAlpha and the Regression Intercept for Ocean-type Populations")

# Second patchwork of Benchmarks 
# logSMAX_sd <- dsmax$deripost_summary$logSMAX_sd$Mean
# logWAshifted <- log(WAin$WA) - mean(WAbase$logWA) # turn this into a sequence of 500?

# Stream type example - Louis (Middle sized WA)
# newlogSMAXLouis <- mu_new[1] + mu_new[2] * logWAshifted[10] # + rnorm(1, 0, sd = logSMAX_sd)
# newlogSMAXseqLouis <- seq(newlogSMAXLouis - 4*logSMAX_sd, newlogSMAXLouis + 4*logSMAX_sd, length.out = 500)
# newlogSMAXdnormLouis <- dnorm(newlogSMAXseqLouis, mean = newlogSMAXLouis, sd = logSMAX_sd)

# Ocean type example - Maria (Smallest WA)
# newlogSMAXMaria <- mu2_new[1] + mu2_new[2] * logWAshifted[25] # + rnorm(1, 0, sd = logSMAX_sd)
# newlogSMAXseqMaria <- seq(newlogSMAXMaria - 4*logSMAX_sd, newlogSMAXMaria + 4*logSMAX_sd, length.out = 500)
# newlogSMAXdnormMaria <- dnorm(newlogSMAXseqMaria, mean = newlogSMAXMaria, sd = logSMAX_sd)

# original logSMAX distribution - posterior predictive
SMAXmu <- log(dsmax$deripost_summary$SMAX_tar_adj$Mean) # mean of the medians 
SMAXsd <- apply(log(dsmax$deripost_full$SMAX_tar_adj), 2, sd)

# Then calculate SMSY and SGEN based on new SMAX and alpha?
	# draw rnorms for new SMAX distribution
# newSMAXrnormLouis <- rnorm(newlogSMAXseqLouis, mean = newlogSMAXLouis, logSMAX_sd)	
# RicArnorm <- rnorm(RicA$seq, mean = Ricprior[1], sd = Ricprior[2])
# BETALouis <- 1/exp(newSMAXrnormLouis)

# newSMAXrnormMaria <- rnorm(newlogSMAXseqMaria, mean = newlogSMAXMaria, logSMAX_sd)	
# BETAMaria <- 1/exp(newSMAXrnormMaria)

# NEW SMSY
# SMSYLouis <- (1 - LambertW0(exp(1 - RicArnorm))) / BETALouis
# logSMSYLouis <- log(SMSYLouis)
# newSMSYseqLouis <- seq(mean(logSMSYLouis) - 4*sd(logSMSYLouis), mean(logSMSYLouis) + 4*sd(logSMSYLouis), length.out = 500)
# newSMSYdnormLouis <- dnorm(newSMSYseqLouis, mean = mean(logSMSYLouis), sd = sd(logSMSYLouis))

# SMSYMaria <- (1 - LambertW0(exp(1 - RicArnorm))) / BETAMaria
# logSMSYMaria <- log(SMSYMaria)
# newSMSYseqMaria <- seq(mean(logSMSYMaria) - 4*sd(logSMSYMaria), mean(logSMSYMaria) + 4*sd(logSMSYMaria), length.out = 500)
# newSMSYdnormMaria <- dnorm(newSMSYseqMaria, mean = mean(logSMSYMaria), sd = sd(logSMSYMaria))

# ORIGINAL SMSY
SMSYmu <- log(dsmax$deripost_summary$SMSY_adj$Median) # mean of the medians 
SMSYsd <- apply(log(dsmax$deripost_full$SMSY_adj), 2, sd)

# NEW SGEN
# SGENLouis <- -1/BETALouis * LambertW0(-BETALouis * SMSYLouis/(exp(RicArnorm)))
# logSGENLouis <- log(SGENLouis)
# newSGENseqLouis <- seq(mean(logSGENLouis) - 4*sd(logSGENLouis), mean(logSGENLouis) + 4*sd(logSGENLouis), length.out = 500)
# newSGENdnormLouis <- dnorm(newSGENseqLouis, mean = mean(logSGENLouis), sd = sd(logSGENLouis))

# SGENMaria <- -1/BETAMaria * LambertW0(-BETAMaria * SMSYMaria/(exp(RicArnorm)))
# logSGENMaria <- log(SGENMaria)
# newSGENseqMaria <- seq(mean(logSGENMaria) - 4*sd(logSGENMaria), mean(logSGENMaria) + 4*sd(logSGENMaria), length.out = 500)
# newSGENdnormMaria <- dnorm(newSGENseqMaria, mean = mean(logSGENMaria), sd = sd(logSGENMaria))

# ORIGINAL SGEN
SGENmu <- log(dsmax$deripost_summary$SGEN_adj$Median) # mean of the medians 
SGENsd <- apply(log(dsmax$deripost_full$SGEN_adj), 2, sd)

# NEW SREP
# logSREPLouis <- log(RicArnorm/BETALouis)
# newSREPseqLouis <- seq(mean(logSREPLouis) - 4*sd(logSREPLouis), mean(logSREPLouis) + 4*sd(logSREPLouis), length.out = 500)
# newSREPdnormLouis <- dnorm(newSREPseqLouis, mean = mean(logSREPLouis), sd = sd(logSREPLouis))

# logSREPMaria <- log(RicArnorm/BETAMaria)
# newSREPseqMaria <- seq(mean(logSREPMaria) - 4*sd(logSREPMaria), mean(logSREPMaria) + 4*sd(logSREPMaria), length.out = 500)
# newSREPdnormMaria <- dnorm(newSREPseqMaria, mean = mean(logSREPMaria), sd = sd(logSREPMaria))

# ORIGINAL SREP
SREPmu <- log(dsmax$deripost_summary$SREP_adj$Median) # mean of the medians 
SREPsd <- apply(log(dsmax$deripost_full$SREP_adj), 2, sd)

# ORIGINAL POSTERIORS:
originalparam <- list(
	SMAX = list(mu = SMAXmu, sigma = SMAXsd),
	SMSY = list(mu = SMSYmu, sigma = SMSYsd),
	SGEN = list(mu = SGENmu, sigma = SGENsd),
	SREP = list(mu = SREPmu, sigma = SREPsd)
)
Louis <- lapply(originalparam, function(p) make_dnorm_dist(p$mu[10], p$sigma[10]))
Maria <- lapply(originalparam, function(p) make_dnorm_dist(p$mu[25], p$sigma[25]))
Thompson <- lapply(originalparam, function(p) make_dnorm_dist(p$mu[24], p$sigma[24]))
	# Now can access Louis$SMAX$seq and Louis$SMAX$dnorm

# What if we just extract mean and sd of each benchmark from bootstrapped simulations?
newSMAXmu <- log(BS.smax$Value[BS.smax$Stock == "Louis" & BS.smax$RP == "SMAX"]) # SMAX Median
newSMAXsd <- (log(BS.smax$Value[BS.smax$Stock == "Louis" & BS.smax$RP == "SMAX"]) - log(BS.smax$lwr[BS.smax$Stock == "Louis" & BS.smax$RP == "SMAX"]))/1.96
newSMAXseqLouis <- seq(newSMAXmu - 4*newSMAXsd, newSMAXmu + 4*newSMAXsd, length.out = 500)
newSMAXdnormLouis <- dnorm(newSMAXseqLouis, mean = newSMAXmu, sd = newSMAXsd)

newSMSYmu <- log(BS.smax$Value[BS.smax$Stock == "Louis" & BS.smax$RP == "SMSY"]) # SMAX Median
newSMSYsd <- (log(BS.smax$Value[BS.smax$Stock == "Louis" & BS.smax$RP == "SMSY"]) - log(BS.smax$lwr[BS.smax$Stock == "Louis" & BS.smax$RP == "SMSY"]))/1.96
newSMSYseqLouis <- seq(newSMSYmu - 4*newSMSYsd, newSMSYmu + 4*newSMSYsd, length.out = 500)
newSMSYdnormLouis <- dnorm(newSMSYseqLouis, mean = newSMSYmu, sd = newSMSYsd)

newSGEMmu <- log(BS.smax$Value[BS.smax$Stock == "Louis" & BS.smax$RP == "SGEN"]) # SMAX Median
newSGENsd <- (log(BS.smax$Value[BS.smax$Stock == "Louis" & BS.smax$RP == "SGEN"]) - log(BS.smax$lwr[BS.smax$Stock == "Louis" & BS.smax$RP == "SGEN"]))/1.96
newSGENseqLouis <- seq(newSGEMmu - 4*newSGENsd, newSGEMmu + 4*newSGENsd, length.out = 500)
newSGENdnormLouis <- dnorm(newSGENseqLouis, mean = newSGEMmu, sd = newSGENsd)

newSREPmu <- log(BS.smax$Value[BS.smax$Stock == "Louis" & BS.smax$RP == "SREP"]) # SMAX Median
newSREPsd <- (log(BS.smax$Value[BS.smax$Stock == "Louis" & BS.smax$RP == "SREP"]) - log(BS.smax$lwr[BS.smax$Stock == "Louis" & BS.smax$RP == "SREP"]))/1.96
newSREPseqLouis <- seq(newSREPmu - 4*newSREPsd, newSREPmu + 4*newSREPsd, length.out = 500)
newSREPdnormLouis <- dnorm(newSREPseqLouis, mean = newSREPmu, sd = newSREPsd)

# Maria TEST
newSMAXmuMaria <- log(BS.smax$Value[BS.smax$Stock == "Maria" & BS.smax$RP == "SMAX"]) # SMAX Median
newSMAXsdMaria <- (log(BS.smax$Value[BS.smax$Stock == "Maria" & BS.smax$RP == "SMAX"]) - log(BS.smax$lwr[BS.smax$Stock == "Maria" & BS.smax$RP == "SMAX"]))/1.96
newSMAXseqMaria <- seq(newSMAXmuMaria - 4*newSMAXsdMaria, newSMAXmuMaria + 4*newSMAXsdMaria, length.out = 500)
newSMAXdnormMaria <- dnorm(newSMAXseqMaria, mean = newSMAXmuMaria, sd = newSMAXsdMaria)

newSMSYmuMaria <- log(BS.smax$Value[BS.smax$Stock == "Maria" & BS.smax$RP == "SMSY"]) # SMAX Median
newSMSYsdMaria <- (log(BS.smax$Value[BS.smax$Stock == "Maria" & BS.smax$RP == "SMSY"]) - log(BS.smax$lwr[BS.smax$Stock == "Maria" & BS.smax$RP == "SMSY"]))/1.96
newSMSYseqMaria <- seq(newSMSYmuMaria - 4*newSMSYsdMaria, newSMSYmuMaria + 4*newSMSYsdMaria, length.out = 500)
newSMSYdnormMaria <- dnorm(newSMSYseqMaria, mean = newSMSYmuMaria, sd = newSMSYsdMaria)

newSGEMmuMaria <- log(BS.smax$Value[BS.smax$Stock == "Maria" & BS.smax$RP == "SGEN"]) # SMAX Median
newSGENsdMaria <- (log(BS.smax$Value[BS.smax$Stock == "Maria" & BS.smax$RP == "SGEN"]) - log(BS.smax$lwr[BS.smax$Stock == "Maria" & BS.smax$RP == "SGEN"]))/1.96
newSGENseqMaria <- seq(newSGEMmuMaria - 4*newSGENsdMaria, newSGEMmuMaria + 4*newSGENsdMaria, length.out = 500)
newSGENdnormMaria <- dnorm(newSGENseqMaria, mean = newSGEMmuMaria, sd = newSGENsdMaria)

newSREPmuMaria <- log(BS.smax$Value[BS.smax$Stock == "Maria" & BS.smax$RP == "SREP"]) # SMAX Median
newSREPsdMaria <- (log(BS.smax$Value[BS.smax$Stock == "Maria" & BS.smax$RP == "SREP"]) - log(BS.smax$lwr[BS.smax$Stock == "Maria" & BS.smax$RP == "SREP"]))/1.96
newSREPseqMaria <- seq(newSREPmuMaria - 4*newSREPsdMaria, newSREPmuMaria + 4*newSREPsdMaria, length.out = 500)
newSREPdnormMaria <- dnorm(newSREPseqMaria, mean = newSREPmuMaria, sd = newSREPsdMaria)

# S. Thompson TEST
newSMAXmuThompson <- log(BS.smax$Value[BS.smax$Stock == "S. Thompson" & BS.smax$RP == "SMAX"]) # SMAX Median
newSMAXsdThompson <- (log(BS.smax$Value[BS.smax$Stock == "S. Thompson" & BS.smax$RP == "SMAX"]) - log(BS.smax$lwr[BS.smax$Stock == "S. Thompson" & BS.smax$RP == "SMAX"]))/1.96
newSMAXseqThompson <- seq(newSMAXmuThompson - 4*newSMAXsdThompson, newSMAXmuThompson + 4*newSMAXsdThompson, length.out = 500)
newSMAXdnormThompson <- dnorm(newSMAXseqThompson, mean = newSMAXmuThompson, sd = newSMAXsdThompson)

newSMSYmuThompson <- log(BS.smax$Value[BS.smax$Stock == "S. Thompson" & BS.smax$RP == "SMSY"]) # SMAX Median
newSMSYsdThompson <- (log(BS.smax$Value[BS.smax$Stock == "S. Thompson" & BS.smax$RP == "SMSY"]) - log(BS.smax$lwr[BS.smax$Stock == "S. Thompson" & BS.smax$RP == "SMSY"]))/1.96
newSMSYseqThompson <- seq(newSMSYmuThompson - 4*newSMSYsdThompson, newSMSYmuThompson + 4*newSMSYsdThompson, length.out = 500)
newSMSYdnormThompson <- dnorm(newSMSYseqThompson, mean = newSMSYmuThompson, sd = newSMSYsdThompson)

newSGEMmuThompson <- log(BS.smax$Value[BS.smax$Stock == "S. Thompson" & BS.smax$RP == "SGEN"]) # SMAX Median
newSGENsdThompson <- (log(BS.smax$Value[BS.smax$Stock == "S. Thompson" & BS.smax$RP == "SGEN"]) - log(BS.smax$lwr[BS.smax$Stock == "S. Thompson" & BS.smax$RP == "SGEN"]))/1.96
newSGENseqThompson <- seq(newSGEMmuThompson - 4*newSGENsdThompson, newSGEMmuThompson + 4*newSGENsdThompson, length.out = 500)
newSGENdnormThompson <- dnorm(newSGENseqThompson, mean = newSGEMmuThompson, sd = newSGENsdThompson)

newSREPmuThompson <- log(BS.smax$Value[BS.smax$Stock == "S. Thompson" & BS.smax$RP == "SREP"]) # SMAX Median
newSREPsdThompson <- (log(BS.smax$Value[BS.smax$Stock == "S. Thompson" & BS.smax$RP == "SREP"]) - log(BS.smax$lwr[BS.smax$Stock == "S. Thompson" & BS.smax$RP == "SREP"]))/1.96
newSREPseqThompson <- seq(newSREPmuThompson - 4*newSREPsdThompson, newSREPmuThompson + 4*newSREPsdThompson, length.out = 500)
newSREPdnormThompson <- dnorm(newSREPseqThompson, mean = newSREPmuThompson, sd = newSREPsdThompson)



x_lim2 <- c(2, 14)

Louis_bottom1 <- ggplot() + # SMAX
	geom_line(aes(x = Louis$SMAX$seq, y = Louis$SMAX$dnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(aes(x = newSMAXseqLouis, y = newSMAXdnormLouis), 
		color = "red", linetype = "dashed", linewidth = 1.2, alpha = 0.5) +
	xlab(TeX("$log(S_{MAX})$")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim2) +
	labs(title = "A")
	
Louis_bottom2 <- ggplot() + # SMSY
	geom_line(aes(x = Louis$SMSY$seq, y = Louis$SMSY$dnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(aes(x = newSMSYseqLouis, y = newSMSYdnormLouis),
		color = "red", linetype = "dashed", linewidth = 1.2, alpha = 0.5) +
	xlab(TeX("$log(S_{MSY})$")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim2) +
	labs(title = "B")
	
Louis_bottom3 <- ggplot() + # SGEN
	geom_line(aes(x = Louis$SGEN$seq, y = Louis$SGEN$dnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(aes(x = newSGENseqLouis, y = newSGENdnormLouis),
		color = "red", linetype = "dashed", linewidth = 1.2, alpha = 0.5) +
	xlab(TeX("$log(S_{GEN})$")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim2) +
	labs(title = "C")
	
Louis_bottom4 <- ggplot() + # SREP
		geom_line(aes(x = Louis$SREP$seq, y = Louis$SREP$dnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(aes(x = newSREPseqLouis, y = newSREPdnormLouis),
		color = "red", linetype = "dashed", linewidth = 1.2, alpha = 0.5) +
	xlab(TeX("$log(S_{REP})$")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim2) +
	labs(title = "D")
	
# (Louis_bottom1 + Louis_bottom2 + Louis_bottom3 + Louis_bottom4) + 
	# plot_layout(ncol = 1) +
	# plot_annotation(title = "Example Benchmarks from Louis")
	

# FOR Maria
Maria_bottom1 <- ggplot() + # SMAX
	geom_line(aes(x = Maria$SMAX$seq, y = Maria$SMAX$dnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(aes(x = newSMAXseqMaria, y = newSMAXdnormMaria), 
		color = "red", linetype = "dashed", linewidth = 1.2, alpha = 0.5) +
	xlab(TeX("$log(S_{MAX})$")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim2) +
	labs(title = "A")
	
Maria_bottom2 <- ggplot() + # SMSY
	geom_line(aes(x = Maria$SMSY$seq, y = Maria$SMSY$dnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(aes(x = newSMSYseqMaria, y = newSMSYdnormMaria),
		color = "red", linetype = "dashed", linewidth = 1.2, alpha = 0.5) +
	xlab(TeX("$log(S_{MSY})$")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim2) +
	labs(title = "B")
	
Maria_bottom3 <- ggplot() + # SGEN
	geom_line(aes(x = Maria$SGEN$seq, y = Maria$SGEN$dnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(aes(x = newSGENseqMaria, y = newSGENdnormMaria),
		color = "red", linetype = "dashed", linewidth = 1.2, alpha = 0.5) +
	xlab(TeX("$log(S_{GEN})$")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim2) +
	labs(title = "C")
	
Maria_bottom4 <- ggplot() + # SREP
		geom_line(aes(x = Maria$SREP$seq, y = Maria$SREP$dnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(aes(x = newSREPseqMaria, y = newSREPdnormMaria),
		color = "red", linetype = "dashed", linewidth = 1.2, alpha = 0.5) +
	xlab(TeX("$log(S_{REP})$")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim2) +
	labs(title = "D")
	
# (Maria_bottom1 + Maria_bottom2 + Maria_bottom3 + Maria_bottom4) + 
	# plot_layout(ncol = 1) +
	# plot_annotation(title = "Example Benchmarks from Maria")

# FOR Thompson
Thomp_bottom1 <- ggplot() + # SMAX
	geom_line(aes(x = Thompson$SMAX$seq, y = Thompson$SMAX$dnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(aes(x = newSMAXseqThompson, y = newSMAXdnormThompson), 
		color = "red", linetype = "dashed", linewidth = 1.2, alpha = 0.5) +
	xlab(TeX("$log(S_{MAX})$")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim2) +
	labs(title = "A")
	
Thomp_bottom2 <- ggplot() + # SMSY
	geom_line(aes(x = Thompson$SMSY$seq, y = Thompson$SMSY$dnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(aes(x = newSMSYseqThompson, y = newSMSYdnormThompson),
		color = "red", linetype = "dashed", linewidth = 1.2, alpha = 0.5) +
	xlab(TeX("$log(S_{MSY})$")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim2) +
	labs(title = "B")
	
Thomp_bottom3 <- ggplot() + # SGEN
	geom_line(aes(x = Thompson$SGEN$seq, y = Thompson$SGEN$dnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(aes(x = newSGENseqThompson, y = newSGENdnormThompson),
		color = "red", linetype = "dashed", linewidth = 1.2, alpha = 0.5) +
	xlab(TeX("$log(S_{GEN})$")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim2) +
	labs(title = "C")
	
Thomp_bottom4 <- ggplot() + # SREP
		geom_line(aes(x = Thompson$SREP$seq, y = Thompson$SREP$dnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(aes(x = newSREPseqThompson, y = newSREPdnormThompson),
		color = "red", linetype = "dashed", linewidth = 1.2, alpha = 0.5) +
	xlab(TeX("$log(S_{REP})$")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim2) +
	labs(title = "D")
	
# two columns
p_col <- wrap_plots(Louis_bottom1, Louis_bottom2, Louis_bottom3, Louis_bottom4, ncol = 1) + 
  plot_annotation(title = "Louis (Stream)")

m_col <- wrap_plots(Maria_bottom1, Maria_bottom2, Maria_bottom3, Maria_bottom4, ncol = 1) + 
  plot_annotation(title = "Maria (Ocean)")
  
t_col <- wrap_plots(Thomp_bottom1, Thomp_bottom2, Thomp_bottom3, Thomp_bottom4, ncol = 1) + 
  plot_annotation(title = "S. Thompson (Ocean)")

wrap_elements(patchworkGrob(m_col)) | wrap_elements(patchworkGrob(p_col)) | wrap_elements(patchworkGrob(t_col))