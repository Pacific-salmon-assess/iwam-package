
# Libs ####
library(RTMB)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(progress)
library(diffr)

# source
source(here::here("R/helperFunctions.R"))
source(here::here("R/PlotFunctions.R"))
source(here::here("R/IWAM_model.R"))
  # Had to add TMB::sdreport to all_pars function
  # otherwise can't find sdreport through rtmb
# source(here::here("R/IWAMsrep_RTMB_model.R"))
source(here::here("R/IWAMsmax_RTMB_model.R"))

# Test Runs ####
# smaxrun <- IWAMsmax_rtmb(bs_nBS = 20000,
#                          bs_seed = 1)
smaxruntest2 <- IWAMsmax_rtmb(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
                             bs_nBS = 20000,
                             lamberts = "tmb")
# smaxruntest <- IWAMsmax_rtmb()
# iwamrun <- IWAM_func(bs_nBS = 20000,
#                      bs_seed = 1)
iwamruntest2 <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
                         bs_nBS = 20000,
                         plot = TRUE)
# iwamruntest <- IWAM_func()
# model opt's converge
  # tmb vs. rtmb model by itself is the same

# IWAM PRIOR TESTING ####
  # Add in the save RDS files for the plotRicA boxplot
# Ricker Log A prior tests
iwam_r_default <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
                            bs_nBS = 20000,
                            plot = TRUE,
                            SigRicPrior = c(F,T,F), # Ric
                            SigDeltaPrior = c(F,T,F), # WA default
                            TauPrior = c(0.1, 1)) # [1] is Ric, [2] is WA
iwam_r_halfnorm <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
                          bs_nBS = 20000,
                          plot = TRUE,
                          SigRicPrior = c(T,F,F), # Ric half normal
                          SigDeltaPrior = c(F,T,F),
                          TauPrior = c(0.1, 1)) 
iwam_r_halfcauchy <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
                             bs_nBS = 20000,
                             plot = TRUE,
                             SigRicPrior = c(F,F,T), # Ric half cauchy
                             SigDeltaPrior = c(F,T,F), # WA Default
                             TauPrior = c(0.1, 1))  # default

# RIC Gamma Priors x 3
iwam_r_gamma0.1 <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
                           bs_nBS = 20000,
                           plot = TRUE,
                           SigRicPrior = c(F,T,F), # Ric half cauchy
                           SigDeltaPrior = c(F,T,F), # WA Default
                           TauPrior = c(0.1, 1))  # default
iwam_r_gamma0.01 <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
                            bs_nBS = 20000,
                            plot = TRUE,
                            SigRicPrior = c(F,T,F), # Ric half cauchy
                            SigDeltaPrior = c(F,T,F), # WA Default
                            TauPrior = c(0.01, 1))  # default
iwam_r_gamma0.001 <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
                             bs_nBS = 20000,
                             plot = TRUE,
                             SigRicPrior = c(F,T,F), # Ric half cauchy
                             SigDeltaPrior = c(F,T,F), # WA Default
                             TauPrior = c(0.001, 1))  # default

# WA Priors x 3
  # half normal and half cauchy models do not run
iwam_wa_default <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
                          bs_nBS = 10,
                          plot = TRUE,
                          SigRicPrior = c(F,T,F), # Ric
                          SigDeltaPrior = c(F,T,F), # WA default
                          TauPrior = c(0.1, 1)) # [1] is Ric, [2] is WA
iwam_wa_halfnorm <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
                          bs_nBS = 10,
                          plot = TRUE,
                          SigRicPrior = c(F,T,F), # Ric
                          SigDeltaPrior = c(T,F,F), # WA default
                          TauPrior = c(0.1, 1)) # [1] is Ric, [2] is WA
  # NaNs produced in sqrt
iwam_wa_cauchy <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
                          bs_nBS = 10,
                          plot = TRUE,
                          SigRicPrior = c(F,T,F), # Ric
                          SigDeltaPrior = c(F,F,T), # WA default
                          TauPrior = c(0.1, 1)) # [1] is Ric, [2] is WA
  # NAs produced - missing value where TRUE/FALSE needed

# WA Gamma Priors x 3
  # All models compile, run, and converge
iwam_wa_gamma1 <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
                          bs_nBS = 10,
                          plot = TRUE,
                          SigRicPrior = c(F,T,F), # Ric
                          SigDeltaPrior = c(F,T,F), # WA default
                          TauPrior = c(0.1, 1)) # [1] is Ric, [2] is WA
iwam_wa_gamma0.1 <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
                            bs_nBS = 10,
                            plot = TRUE,
                            SigRicPrior = c(F,T,F), # Ric
                            SigDeltaPrior = c(F,T,F), # WA default
                            TauPrior = c(0.1, 0.1)) # [1] is Ric, [2] is WA
iwam_wa_gamma0.01 <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
                            bs_nBS = 10,
                            plot = TRUE,
                            SigRicPrior = c(F,T,F), # Ric
                            SigDeltaPrior = c(F,T,F), # WA default
                            TauPrior = c(0.1, 0.01)) # [1] is Ric, [2] is WA


# bootstrap convergence: ####
  # Tor:Seems like there is not exact convergence for bootstraps
  # Plot the differences

iwamrundel <- iwamrun$dfout %>% 
  filter(!Stock %in% c("Barkley", "Clayoquot", "Kyuquot", "Nootka/Esperanza",
                      "Quatsino", "WCVI Nootka & Kyuquot", "WCVI North",
                      "WCVI South"))

plot(iwamrundel$Value[iwamrundel$RP == 'SGEN'], col = 'red')
plot(smaxrun$dfout$Value[smaxrun$dfout$RP == 'SGEN'], col = 'black')

temp1 <- ggplot(data = iwamrundel[iwamrundel$RP == 'SGEN', ], aes(x = Stock, 
                                       y = Value)) + 
  geom_errorbar(mapping = aes(ymax = upr, 
                    ymin = lwr, color = 'red')) +
  geom_point(aes(color = 'red')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1))

temp2 <- temp1 + 
  geom_point(data = smaxrun$dfout[smaxrun$dfout$RP == 'SGEN', ], mapping = aes(x = Stock, y = Value)) +
  geom_errorbar(aes(x = Stock, ymax = upr, ymin = lwr))

temp2

# example
SMSY_parken <- SMSY_pl + 
  geom_point(eval_data, mapping = aes(x = Stock, y = PA_SMSY), 
             position = position_nudge(-0.3)) + 
  geom_errorbar(aes(x = Stock, ymax = PA_SMSY + (1.96 * PA_SE_SMSY) , ymin = PA_SMSY - (1.96 * PA_SE_SMSY)), 
                width = 0.2,                     
                position = position_nudge(-0.3), 
                inherit.aes = FALSE) 



# Lambert Testing ####
compile(here::here("TMB_Files/lambert.cpp"))
dyn.load(dynlib(here::here("TMB_Files/lambert")))

objlw <- MakeADFun(data=list(), parameters=list(x=1), DLL="lambert")
objlw$fn(7 * exp(7))
lambert_W0(7 * exp(7)) # eg

