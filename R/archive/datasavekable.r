# Saving IWAM Model Posteriors for Sharing

# Libaries ####
library(dplyr)
library(tidyverse) 
library(beepr) # Sounds
library(data.table) # Create data tables for pivoting

# CREATE CSV/TABLES ####
targetnames <- WAin |>
  rename("Stock_name" = Stock)
  
  # SREP ESTIMATE FOR SYNOPTIC POPULATIONS
# dtable <- cbind(targetnames, dsmax$deripost_summary$SREP) |> 
	# rename("SREP_mean" = Mean, "SREP_median" = Median,
    # "SREP_LQ_5" = LQ_5, "SREP_UQ_95" = UQ_95, "SREP_Stocknum" = Stock,
    # "SREP_Mode" = PosteriorMode)
# dtable <- cbind(dtable, dsmax$deripost_summary$SMSY) |> 
	# rename("SMSY_mean" = Mean, "SMSY_median" = Median,
    # "SMSY_LQ_5" = LQ_5, "SMSY_UQ_95" = UQ_95, "SMSY_Stocknum" = Stock,
    # "SMSY_Mode" = PosteriorMode)
# dtable <- cbind(dtable, dsmax$deripost_summary$SGEN) |> 
	# rename("SGEN_mean" = Mean, "SGEN_median" = Median,
    # "SGEN_LQ_5" = LQ_5, "SGEN_UQ_95" = UQ_95, "SGEN_Stocknum" = Stock,
    # "SGEN_Mode" = PosteriorMode)

# dtable <- dtable |> 
	# select(
		# CU_INDEX, Stock_name, WA, lh,
		# SGEN_mean, SGEN_median, SGEN_LQ_5, SGEN_UQ_95,
		# SREP_mean, SREP_median, SREP_LQ_5, SREP_UQ_95,
		# SMSY_mean, SMSY_median, SMSY_LQ_5, SMSY_UQ_95
	# ) |> # INSERTION METHOD
	# mutate(
		# WA = round(WA),
		# SGEN_mean = round(SGEN_mean), SGEN_median = round(SGEN_median),
		# SGEN_LQ_5 = round(SGEN_LQ_5), SGEN_UQ_95 = round(SGEN_UQ_95),
		# SREP_mean = round(SREP_mean), SREP_median = round(SREP_median),
		# SREP_LQ_5 = round(SREP_LQ_5), SREP_UQ_95 = round(SREP_UQ_95),
		# SMSY_mean = round(SMSY_mean), SMSY_median = round(SMSY_median),
		# SMSY_LQ_5 = round(SMSY_LQ_5), SMSY_UQ_95 = round(SMSY_UQ_95)
	# )
	
# Sample table code w/ kable
# kable(complete, caption = "IWAM Smax Model: SGEN, SREP, and SMSY Estimates for All Stocks",
      # format.args = list(big.mark = ","))
# write.csv(dtable, here::here("docs/TableofEstimates_SEPdraftbenchmarks_20260127.csv"), row.names = FALSE)


# dtable_adj <- cbind(targetnames, dsmax$deripost_summary$SREP_adj) |> 
	# rename("SREP_mean" = Mean, "SREP_median" = Median,
    # "SREP_LQ_5" = LQ_5, "SREP_UQ_95" = UQ_95, "SREP_Stocknum" = Stock,
    # "SREP_Mode" = PosteriorMode)
# dtable_adj <- cbind(dtable_adj, dsmax$deripost_summary$SMSY_adj) |> 
	# rename("SMSY_mean" = Mean, "SMSY_median" = Median,
    # "SMSY_LQ_5" = LQ_5, "SMSY_UQ_95" = UQ_95, "SMSY_Stocknum" = Stock,
    # "SMSY_Mode" = PosteriorMode)
# dtable_adj <- cbind(dtable_adj, dsmax$deripost_summary$SGEN_adj) |> 
	# rename("SGEN_mean" = Mean, "SGEN_median" = Median,
    # "SGEN_LQ_5" = UQ_95, "SGEN_UQ_95" = LQ_5, "SGEN_Stocknum" = Stock,
    # "SGEN_Mode" = PosteriorMode)
	
# dtable_adj <- dtable_adj |> 
	# select(CU_INDEX, Stock_name, WA, lh,
	# SGEN_mean, SGEN_median, SGEN_LQ_5, SGEN_UQ_95,
	# SREP_mean, SREP_median, SREP_LQ_5, SREP_UQ_95,
	# SMSY_mean, SMSY_median, SMSY_LQ_5, SMSY_UQ_95
	# ) |> # INSERTION METHOD
	# mutate(
		# WA = round(WA),
		# SGEN_mean = round(SGEN_mean), SGEN_median = round(SGEN_median),
		# SGEN_LQ_5 = round(SGEN_LQ_5), SGEN_UQ_95 = round(SGEN_UQ_95),
		# SREP_mean = round(SREP_mean), SREP_median = round(SREP_median),
		# SREP_LQ_5 = round(SREP_LQ_5), SREP_UQ_95 = round(SREP_UQ_95),
		# SMSY_mean = round(SMSY_mean), SMSY_median = round(SMSY_median),
		# SMSY_LQ_5 = round(SMSY_LQ_5), SMSY_UQ_95 = round(SMSY_UQ_95)
	# )
	
# write.csv(dtable_adj, here::here("docs/TableofEstimates_marginaladjusted_SEPdraftbenchmarks_20262201.csv"), row.names = FALSE)
