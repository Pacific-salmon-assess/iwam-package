# IWAM Model Run for RPA on Nanaimo and Puntledge stocks

# Model steps available in the following files:
  # Model function available in IWAM_model.R
  # Complete TMB model in IWAM_Liermann.cpp
  # Bootstrapping method in Get_LRP_bs.R

# Libraries and Sources ####
# library(RTMB)
# library(ggplot2)
# library(dplyr)
# library(tidyverse)
# library(progress)
# library(diffr)

# source
# source(here::here("R/helperFunctions.R"))
# source(here::here("R/PlotFunctions.R"))
source(here::here("R/IWAM_model.R"))

# Function
    # 15-30 minutes of run time for 20,000 bootstrap iterations
iwamnanpunt <- IWAM_func(WAinraw = c("DataIn/WCVIStocks_NanPunt.csv"),
                          targetname = "WCVI_NanPunt_test",
                          bs_seed = 1,
                          bs_nBS = 20000,
                          bias.cor = TRUE,
                          # Remaining function parameters are default listed in IWAM_model.R
                          prod = "Parken") # prod = "LifeStageModel")

# There are two versions of PredInt functions within
  # helperFunctions.r
  # This is part of the source of the discrepancy.
  # There remains a small difference in the upr confidence interval estimate.

# Other known differences

# Bias correction
  # Bias correction terms can be turned on or off in the model.
  # But must be physically commented out within the bootstrapping stage.

# Penalty terms:
  # Utilizing original penalty terms - most notably dgamma on precision for 
  # Ricker and alpha terms.

# For reference:
  # You can get explicit Parken estimates of loga and beta
  # if you change RPs to RPs_e within the Get_LRP_bs.R.
