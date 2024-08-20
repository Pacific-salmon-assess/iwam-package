
# Libs
library(RTMB)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(progress)

# source
source(here::here("R/helperFunctions.R"))
source(here::here("R/PlotFunctions.R"))
source(here::here("R/IWAM_model.R"))

# source rtmb
# source(here::here("R/IWAMsrep_RTMB_model.R"))
source(here::here("R/IWAMsmax_RTMB_model.R"))

# Run's
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



# bootstrap convergence:
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

