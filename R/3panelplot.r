# Library
library(tidyverse)
library(ggplot2)
library(patchwork)
library(latex2exp)



# Create a distribution seq and dnorm df for New Alpha
Ricprior = c(1, 0.3)
RicAseq <- seq(Ricprior[1] - 4*Ricprior[2], Ricprior[1] + 4*Ricprior[2], length.out = 500)
RicAdnorm <- dnorm(RicAseq, mean = Ricprior[1], sd = Ricprior[2])
RicAdf <- data.frame(x = RicAseq, y = RicAdnorm)

# Get draws of posterior and calculate change in b0
# poplen <- ncol(dsmax$deripost_full$Alpha_tar_adj) # 25 stock predictions to be made
# logSMAX_sd <- dsmax$deripost_summary$logSMAX_sd$Mean
b0 <- dsmax$deripost_full$b0
bWA <- dsmax$deripost_full$bWA
Alpha0 <- exp(dsmax$deripost_full$Alpha0) # logAlpha - NOT log(log(alpha))

post <- cbind(Alpha0, b0[,1], bWA[,1]) # Full posteriors for stream-type and global mean alpha
covmatrix <- cov(post)
mu <- apply(post, 2, mean)

# Create a distribution seq and dnorm df for Original b0
b0seq <- seq(mean(post[,2]) - 4*sd(post[,2]), mean(post[,2]) + 4*sd(post[,2]), length.out = 500)
b0dnorm <- dnorm(b0seq, mean = mean(post[,2]), sd = sd(post[,2]))
b0df <- data.frame(x = b0seq, y = b0dnorm)
# Create a distribution seq and dnorm df for Original Alpha0
alpha0seq <-  seq(mean(post[,1]) - 4*sd(post[,1]), mean(post[,1]) + 4*sd(post[,1]), length.out = 500)
alpha0dnorm <- dnorm(alpha0seq, mean = mean(post[,1]), sd = sd(post[,1]))
alpha0df <- data.frame(x = alpha0seq, y = alpha0dnorm)

logalphak <- rnorm(1, Ricprior[1], Ricprior[2])
# stream type (base case)
mu_new <- mu[-1] + covmatrix[-1,1]/covmatrix[1,1]*(logalphak - mu[1])
var_new <- covmatrix[-1,-1] - (covmatrix[-1,1, drop = FALSE] / covmatrix[1,1]) %*% covmatrix[1,-1, drop = FALSE]
# bnew <- MASS::mvrnorm(1, mu = mu_new, var_new)

# Create a distribution seq and dnorm df for New b0
b0newseq <-  seq(mu_new[1] - 4*sqrt(var_new[1,1]), mu_new[1] + 4*sqrt(var_new[1,1]), length.out = 500)
b0newdnorm <- dnorm(b0newseq, mean = mu_new[1], sd = sqrt(var_new[1,1]))
b0newdf <- data.frame(x = b0newseq, y = b0newdnorm)

# Loop to create new random draws of b0 and Alpha
logalphai <- c()
# mu_newi <- c()
bnewi <- c()
for (i in 1:100){
	logalphai[i] <- rnorm(1, Ricprior[1], Ricprior[2])
	mu_newi <- mu[-1] + covmatrix[-1,1]/covmatrix[1,1]*(logalphai[i] - mu[1])
	var_newi <- covmatrix[-1,-1] - (covmatrix[-1,1, drop = FALSE] / covmatrix[1,1]) %*% covmatrix[1,-1, drop = FALSE]
	bnewi[i] <- (MASS::mvrnorm(1, mu = mu_newi, var_newi))[1]
}


# SET LIMITS TO MATCH ACROSS FACETS
x_lim <- range(c(exp(dsmax$deripost_full$Alpha0), RicAseq))
y_lim <- range(c(dsmax$deripost_full$b0[,1], b0newseq))

# p_scatter: scatter w/o lines
p_scatter <- ggplot(data.frame(x = exp(dsmax$deripost_full$Alpha0), y = dsmax$deripost_full$b0[,1]),
	aes(x = x, y = y)) +
	geom_point(alpha = 0.2, show.legend = FALSE) + 
	# point for mean of new distribution
	geom_point(aes(x = mean(RicAseq), y = mean(b0newseq)), color = "red", size = 4) + # MEAN LINE POINT
	geom_point(data = data.frame(x = logalphai, y = bnewi), aes(x = logalphai, y = bnewi), color = "royalblue", size = 3, alpha = 0.4) +
	xlab(TeX("global mean $log(\\alpha)$")) + 
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
	geom_line(data = RicAdf, aes(x = RicAseq, y = RicAdnorm), color = "red", linewidth = 1.5) + 
	geom_vline(xintercept = mean(alpha0seq), alpha = 0.4, linetype = "dashed", linewidth = 1.2) + # mean
	geom_vline(xintercept = mean(RicAseq), alpha = 0.4, color = "red", linetype = "dashed", linewidth = 1.2) + # mean
	geom_segment(aes(x = mean(alpha0seq), xend = mean(RicAseq), y = 2, yend = 2),
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

# patchwork setup
(p_top + plot_spacer() + p_scatter + p_right) + 
	plot_layout(ncol    = 2, widths  = c(3, 1), heights = c(1, 3)) # + 
	# plot_annotation(title = "From the posterior - conditional on random effects being equal to zero")
	# plot_annotation(title = "From the posterior predictive (random effects marginalized)")
		# using _adj posterior values
		


# patchwork 2 setup
# bottom lowest plot
	# original logSMAX distribution - posterior predictive
SMAXmu <- log(dsmax$deripost_summary$SMAX_tar_adj$Mean) # mean of the medians 
SMAXsd <- apply(log(dsmax$deripost_full$SMAX_tar_adj), 2, sd)
SMAXseq <- seq(SMAXmu[10] - 4*SMAXsd[10], SMAXmu[10] + 4*SMAXsd[10], length.out = 500)
SMAXdnorm <- dnorm(SMAXseq, mean = SMAXmu[10], sd = SMAXsd[10])
SMAXdf <- data.frame(x = SMAXseq, y = SMAXdnorm)

# Should this be PER population or mean of all populations?
logSMAX_sd <- dsmax$deripost_summary$logSMAX_sd$Mean
logWAshifted <- log(WAin$WA) - mean(WAbase$logWA) # turn this into a sequence of 500?
mu_new <- mu[-1] + covmatrix[-1,1]/covmatrix[1,1]*(logalphak - mu[1])
var_new <- covmatrix[-1,-1] - (covmatrix[-1,1, drop = FALSE] / covmatrix[1,1]) %*% covmatrix[1,-1, drop = FALSE]
bnew <- MASS::mvrnorm(1, mu = mu_new, var_new)
newlogSMAX <- bnew[1] + bnew[2] * logWAshifted[10] + rnorm(1, 0, sd = logSMAX_sd)

newlogSMAXseq <- seq(newlogSMAX - 4*logSMAX_sd, newlogSMAX + 4*logSMAX_sd, length.out = 500)
newlogSMAXdnorm <- dnorm(newlogSMAXseq, mean = newlogSMAX, sd = logSMAX_sd)
newlogSMAXdf <- data.frame(x = newlogSMAXseq, y = newlogSMAXdnorm)

# Then calculate SMSY and SGEN based on new SMAX and alpha?
	# draw rnorms for new SMAX distribution
newSMAXrnorm <- rnorm(newlogSMAXseq, mean = newlogSMAX, logSMAX_sd)	
RicArnorm <- rnorm(RicAseq, mean = Ricprior[1], sd = Ricprior[2])
BETA <- 1/exp(newSMAXrnorm)
SMSY <- (1 - LambertW0(exp(1 - RicArnorm))) / BETA
logSMSY <- log(SMSY)
newSMSYseq <- seq(mean(logSMSY) - 4*sd(logSMSY), mean(logSMSY) + 4*sd(logSMSY), length.out = 500)
newSMSYdnorm <- dnorm(newSMSYseq, mean = mean(logSMSY), sd = sd(logSMSY))

# Original SMSY
SMSYmu <- log(dsmax$deripost_summary$SMSY_adj$Median) # mean of the medians 
SMSYsd <- apply(log(dsmax$deripost_full$SMSY_adj), 2, sd)
SMSYseq <- seq(SMSYmu[10] - 4*SMSYsd[10], SMSYmu[10] + 4*SMSYsd[10], length.out = 500)
SMSYdnorm <- dnorm(SMSYseq, mean = SMSYmu[10], sd = SMSYsd[10])

SGEN <- -1/BETA * LambertW0(-BETA * SMSY/(exp(RicArnorm)))
logSGEN <- log(SGEN)
newSGENseq <- seq(mean(logSGEN) - 4*sd(logSGEN), mean(logSGEN) + 4*sd(logSGEN), length.out = 500)
newSGENdnorm <- dnorm(newSGENseq, mean = mean(logSGEN), sd = sd(logSGEN))

# ORIGINAL SGEN
SGENmu <- log(dsmax$deripost_summary$SGEN_adj$Median) # mean of the medians 
SGENsd <- apply(log(dsmax$deripost_full$SGEN_adj), 2, sd)
SGENseq <- seq(SGENmu[10] - 4*SGENsd[10], SGENmu[10] + 4*SGENsd[10], length.out = 500)
SGENdnorm <- dnorm(SGENseq, mean = SGENmu[10], sd = SGENsd[10])

# SET LIMITS TO MATCH ACROSS FACETS
x_lim2 <- c(4, 11)
# y_lim2 <- c(3, 12)

p_bottom1 <- ggplot() + # SMAX
	geom_line(aes(x = SMAXseq, y = SMAXdnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(aes(x = newlogSMAXseq, y = newlogSMAXdnorm), color = "red", linewidth = 1.2) +
	xlab(TeX("$log(S_{MAX})$")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim2) +
	labs(title = "A")
	
p_bottom2 <- ggplot() + # SMSY
	geom_line(aes(x = SMSYseq, y = SMSYdnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(aes(x = newSMSYseq, y = newSMSYdnorm), color = "red", linewidth = 1.2) +
	xlab(TeX("$log(S_{MSY})$")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim2) +
	labs(title = "B")
	
p_bottom3 <- ggplot() + # SGEN
	geom_line(aes(x = SGENseq, y = SGENdnorm), linewidth = 1.2) + # smoothed SMAX original
	geom_line(aes(x = newSGENseq, y = newSGENdnorm), color = "red", linewidth = 1.2) +
	xlab(TeX("$log(S_{GEN})$")) + 
	ylab("Density") + 
	theme_classic() +
	theme(axis.text = element_text(size = 18),
		axis.title = element_text(size = 18)) +
	coord_cartesian(xlim = x_lim2) +
	labs(title = "C")
	
(p_bottom1 + p_bottom2 + p_bottom3) + 
	plot_layout(ncol = 1)