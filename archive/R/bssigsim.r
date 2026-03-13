# Simulating different SD for Ricker Alpha prior for Bootstrapping Simulations

# library() 

# Source
source(here::here("R/Liermann_RTMB_model_Bootstrap.R")) 

bssig0.62 <- dobootstrap(bsiters = 20000, Ricprior = c(1,0.62))
bssig0.51 <- dobootstrap(bsiters = 20000, Ricprior = c(1, 0.51)) # Default
bssig0.45 <- dobootstrap(bsiters = 20000, Ricprior = c(1, 0.45))
bssig0.4 <- dobootstrap(bsiters = 20000, Ricprior = c(1, 0.4))
bssig0.3 <- dobootstrap(bsiters = 20000, Ricprior = c(1, 0.3))

bspivot <- function(bsout) {
	wide_bsout <- bsout %>%
		pivot_wider(
		id_cols = c(Stock, WA, lh), 
		names_from = RP, 
		values_from = c(Value, lwr, upr),
		names_sep = "_"
	)
	return(wide_bsout)
}

wide_bssig0.62 <- bspivot(bssig0.62$BS.dfout)
wide_bssig0.51 <- bspivot(bssig0.51$BS.dfout)
wide_bssig0.45 <- bspivot(bssig0.45$BS.dfout)
wide_bssig0.4 <- bspivot(bssig0.4$BS.dfout)
wide_bssig0.3 <- bspivot(bssig0.3$BS.dfout)

ggplot() + 

	# wide_bssig0.62
	geom_errorbar(data = wide_bssig0.62, aes(x = fct_reorder(Stock, log(WA)), 
											y = Value_SMSY, 
											ymax = upr_SMSY, 
											ymin = pmax(lwr_SMSY, 1),
											color = "sig0.62",
											width=.1),
					position = position_nudge(-0.2)) +
	geom_point(data = wide_bssig0.62, aes(x = fct_reorder(Stock, log(WA)), 
										y = Value_SMSY, 
										color = "sig0.62"),
				position = position_nudge(-0.2)) +

	# wide_bssig0.51
	geom_errorbar(data = wide_bssig0.51, aes(x = fct_reorder(Stock, log(WA)), 
											y = Value_SMSY, 
											ymax = upr_SMSY, 
											ymin = pmax(lwr_SMSY, 1),
											color = "sig0.51",
											width=.1),
					position = position_nudge(-0.1)) +
	geom_point(data = wide_bssig0.51, aes(x = fct_reorder(Stock, log(WA)), 
										y = Value_SMSY, 
										color = "sig0.51"),
				position = position_nudge(-0.1)) +
								
	# wide_bssig0.45
	geom_errorbar(data = wide_bssig0.45, aes(x = fct_reorder(Stock, log(WA)), 
											y = Value_SMSY, 
											ymax = upr_SMSY, 
											ymin = lwr_SMSY, # pmax(lwr_SREP, 1)
											color = "sig0.45",
											width=.1)) +
	geom_point(data = wide_bssig0.45, aes(x = fct_reorder(Stock, log(WA)), 
										y = Value_SMSY, 
										color = "sig0.45")) +
										
	# wide_bssig0.4
	geom_errorbar(data = wide_bssig0.4, aes(x = fct_reorder(Stock, log(WA)), 
											y = Value_SMSY, 
											ymax = upr_SMSY, 
											ymin = lwr_SMSY, # pmax(lwr_SREP, 1)
											color = "sig0.4",
											width=.1),
					position = position_nudge(+0.1)) +
	geom_point(data = wide_bssig0.4, aes(x = fct_reorder(Stock, log(WA)), 
										y = Value_SMSY, 
										color = "sig0.4"),
				position = position_nudge(+0.1)) +
	
	# wide_bssig0.3
	geom_errorbar(data = wide_bssig0.3, aes(x = fct_reorder(Stock, log(WA)), 
											y = Value_SMSY, 
											ymax = upr_SMSY, 
											ymin = lwr_SMSY, # pmax(lwr_SREP, 1)
											color = "sig0.3",
											width=.1),
					position = position_nudge(+0.2)) +
	geom_point(data = wide_bssig0.3, aes(x = fct_reorder(Stock, log(WA)), 
										y = Value_SMSY, 
										color = "sig0.3"),
				position = position_nudge(+0.2)) +
										
	geom_hline(yintercept = 1, lty = 'dashed', colour = 'grey') +

	theme_classic() + 
	scale_y_continuous(transform = "log", 
                    breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
	ylab(TeX("Bootstrapped $S_{MSY}$ Estimate")) +
	xlab("") + 
	theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust = 1)) +
	scale_color_manual(name='Model',
                    breaks=c('sig0.51',
							'sig0.62',
                            'sig0.45',
                            'sig0.4',
							'sig0.3'),
                    values=c('sig0.51' = "black",
							'sig0.62' = 'grey',
                            'sig0.45' = "skyblue",
                            'sig0.4' = "royalblue",
							'sig0.3' = 'forestgreen'))
