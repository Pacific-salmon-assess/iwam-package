
# Libs ####
library(RTMB)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(progress)
library(diffr)
library(knitr)

# source
source(here::here("R/helperFunctions.R"))
source(here::here("R/PlotFunctions.R"))
source(here::here("R/PDF_plot.r"))
source(here::here("R/IWAM_model.R"))
  # Had to add TMB::sdreport to all_pars function
  # otherwise can't find sdreport through rtmb
# source(here::here("R/IWAMsrep_RTMB_model.R"))
source(here::here("R/IWAMsmax_RTMB_model.R"))

# Test Runs ####
  # RTMB SMAX Run test
  # Using tmb lambertW version
# smaxruntest <- IWAMsmax_rtmb(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
#                              bs_nBS = 20000,
#                              lamberts = "tmb")
# 
# # IWAM test with no aggregation
# iwamruntest <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
#                          bs_nBS = 10,
#                          plot = TRUE)
# 
# # IWAM test with no inlet aggregation
# iwamtest_noinlet <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoInlet.csv"),
#                           bs_nBS = 10,
#                           plot = FALSE)
# 
# # IWAM test with ALL aggregations (CU and Inlet)
# iwamtest_full <- IWAM_func(WAin = c("DataIn/WCVIStocks.csv"),
#                           bs_nBS = 10,
#                           plot = FALSE)
# 
# # IWAM test with Parken alpha bootstrapping
# iwamtestparken <- IWAM_func(WAinraw = c("DataIn/WCVIStocks_NoAgg.csv"),
#                            bs_nBS = 10,
#                            plot = TRUE,
#                            prod = "Parken")

# iwamruntest <- IWAM_func()
# model opt's converge
  # tmb vs. rtmb model by itself is the same

# IWAM PRIOR TESTING ####
    # only 10 iter's as no bootstrapping occurs before the plotting or RDS creation
    # 20000 iterations for boostrapped values
# Ricker Log A prior tests
  # "DataOut/WAregSREP_target_ricgamma_0.1_wagamma_1_IWAM_Liermann"
iwam_d <- IWAM_func(WAinraw = c("DataIn/WCVIStocks.csv"),
                          targetname = "biascor_OFF",
                    bias.cor = FALSE) # [1] is Ric, [2] is WA

iwam_default <- IWAM_func(WAinraw = c("DataIn/WCVIStocks_NoAgg.csv"),
                          targetname = "techreport",
                          run.predict = TRUE,
                          run.bootstraps = TRUE,
                            bs_nBS = 10,
                            plot = TRUE,
                            SigRicPrior = c(F,T,F), # Ric
                            SigDeltaPrior = c(F,T,F), # WA default
                            TauPrior = c(0.1, 1)) # [1] is Ric, [2] is WA
  # "DataOut/WAregSREP_target_richalfnorm_0.1_wagamma_1_IWAM_Liermann"
iwam_r_halfnorm <- IWAM_func(WAinraw = c("DataIn/WCVIStocks_NoAgg.csv"),
                             targetname = "techreport",
                          bs_nBS = 10,
                          plot = TRUE,
                          SigRicPrior = c(T,F,F), # Ric half normal
                          SigDeltaPrior = c(F,T,F),
                          TauPrior = c(0.1, 1)) 
  # "DataOut/WAregSREP_target_riccauchy_0.1_wagamma_1_IWAM_Liermann"
iwam_r_halfcauchy <- IWAM_func(WAinraw = c("DataIn/WCVIStocks_NoAgg.csv"),
                               targetname = "techreport",
                             bs_nBS = 10,
                             plot = TRUE,
                             SigRicPrior = c(F,F,T), # Ric half cauchy
                             SigDeltaPrior = c(F,T,F), # WA Default
                             TauPrior = c(0.1, 1))  # default

# RIC Gamma Priors x 3
  # iwam_default to compare against
iwam_r_gamma0.01 <- IWAM_func(WAinraw = c("DataIn/WCVIStocks_NoAgg.csv"),
                              targetname = "techreport",
                            bs_nBS = 10,
                            plot = TRUE,
                            SigRicPrior = c(F,T,F), 
                            SigDeltaPrior = c(F,T,F), # WA Default
                            TauPrior = c(0.01, 1))  # default
iwam_r_gamma0.001 <- IWAM_func(WAinraw = c("DataIn/WCVIStocks_NoAgg.csv"),
                               targetname = "techreport",
                             bs_nBS = 10,
                             plot = TRUE,
                             SigRicPrior = c(F,T,F), 
                             SigDeltaPrior = c(F,T,F), # WA Default
                             TauPrior = c(0.001, 1))  # default

# WA Priors x 3
  # half normal and half cauchy models do not run
  # iwam_default to compare against
iwam_wa_halfnorm <- IWAM_func(WAinraw = c("DataIn/WCVIStocks_NoAgg.csv"),
                              targetname = "techreport",
                          bs_nBS = 10,
                          plot = TRUE,
                          SigRicPrior = c(F,T,F), # Ric
                          SigDeltaPrior = c(T,F,F), # WA default
                          TauPrior = c(0.1, 1)) # [1] is Ric, [2] is WA
  # NaNs produced in sqrt
iwam_wa_cauchy <- IWAM_func(WAinraw = c("DataIn/WCVIStocks_NoAgg.csv"),
                            targetname = "techreport",
                          bs_nBS = 10,
                          plot = TRUE,
                          SigRicPrior = c(F,T,F), # Ric
                          SigDeltaPrior = c(F,F,T), # WA default
                          TauPrior = c(0.1, 1)) # [1] is Ric, [2] is WA
  # NAs produced - missing value where TRUE/FALSE needed

# WA Gamma Priors x 3
  # All models compile, run, and converge
  # iwam_default to compare against
iwam_wa_gamma0.1 <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
                              targetname = "techreport",
                            bs_nBS = 10,
                            plot = TRUE,
                            SigRicPrior = c(F,T,F), # Ric
                            SigDeltaPrior = c(F,T,F), # WA default
                            TauPrior = c(0.1, 0.1)) # [1] is Ric, [2] is WA
iwam_wa_gamma0.01 <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
                               targetname = "techreport",
                            bs_nBS = 10,
                            plot = TRUE,
                            SigRicPrior = c(F,T,F), # Ric
                            SigDeltaPrior = c(F,T,F), # WA default
                            TauPrior = c(0.1, 0.01)) # [1] is Ric, [2] is WA

# A total of 12 runs
  # Removed 3 runs due to being the same as default - 9 runs in total
  # 9 RDS files to find
iwam_fixedeff <- IWAM_func(WAinraw = c("DataIn/WCVIStocks_NoAgg.csv"),
                           targetname = "techreport",
                          predict = FALSE,
                          bs_nBS = 10,
                          plot = TRUE,
                          run.bootstraps = FALSE,
                          mod = "IWAM_FixedEffects")
  # PRODUCING ERROR: object datain not found

# SR_curves and WA Linear Reg plots ####
# Basic plots will be created provided plot = TRUE for desired penalty variants.

# RicA Boxplot call ####
png(paste("DataOut/RicADist_ComparePriors_wBC.png", sep=""), width=7, height=7, units="in", res=500)
plotRicA_reduc()
dev.off()
  # Only 1 variant

# PDF Prior Plotting ####
   # WA sigma priors - InvGamma
# png(paste("DataOut/DeltaPriors_InvGamma.png", sep=""), width=7, height=7, units="in", res=500)
png(paste("DataOut/PDF_sigWA_wBC.png", sep=""), width=7, height=7, units="in", res=500)
plotPriors(plot_inv_gamma_only=TRUE, Delta=TRUE, modelobject = iwam_default$all_Deltas)
  # FALSE, TRUE is also WA related
# plotPriors(plot_inv_gamma_only=TRUE, Delta=TRUE, modelobject = iwam_wa_gamma0.1$all_Deltas) # prior comparison
dev.off()

  # Ricker sigma priors
# See TWG ppt - slide 18 and 19
  # plotPriors(plot_inv_gamma_only=TRUE, Delta=FALSE) or
  # plotPriors(plot_inv_gamma_only=FALSE, Delta=FALSE)
# png statement
png(paste("DataOut/PDF_sigRicker_wBC.png", sep=""), width=7, height=7, units="in", res=500)
plotPriors(plot_inv_gamma_only=TRUE, Delta=FALSE, modelobject = iwam_default$all_Deltas)
dev.off()
# plotPriors(plot_inv_gamma_only=FALSE, Delta=FALSE, modelobject = iwam_wa_gamma0.1$all_Deltas)

# SREP Point Estimate comparison with CI's - IWAM vs. Parken ####
  # What file am I pulling from?
    # IWAM_modelcompare.Rmd and SEP_BackCalcRefPoints.RMD (maybe)

# Prepare estimates
SMSY_tmb <- iwam_default$dfout %>% # return: "dfout"
  filter(RP=='SMSY')
SREP_tmb <- iwam_default$dfout %>%
  filter(RP=='SREP')

Parken_WCVI <- read.csv(here::here("DataIn/WCVI_Parken.csv"))

eval_dat <- Parken_WCVI %>% left_join(SMSY_tmb, by=join_by(Stock)) %>% 
  rename("SMSY" = Value, "UL SMSY" = upr, "LL SMSY" = lwr) %>%
  select(-RP) %>% 
  left_join(SREP_tmb, by=join_by(Stock)) %>%
  rename("SREP" = Value, "UL SREP" = upr, "LL SREP" = lwr) %>%
  select(-RP) %>% 
  mutate(PA_UL_SMSY = PA_SMSY + (1.96 * PA_SE_SMSY)) %>% 
  mutate(PA_LL_SMSY = PA_SMSY - (1.96 * PA_SE_SMSY)) %>% 
  mutate(PA_UL_SREP = PA_SREP + (1.96 * PA_SE_SREP)) %>% 
  mutate(PA_LL_SREP = PA_SREP - (1.96 * PA_SE_SREP))
benchmarks <- eval_dat

# Plot
png(paste("DataOut/PWC_WCVI_IWAM_PARKEN_wBC.png", sep=""), width=7, height=7, units="in", res=500)
ggplot(benchmarks, aes(x=Stock, y = SREP)) +
  geom_errorbar(aes(ymax = `UL SREP`, ymin = `LL SREP`, color='TMB',), width = 0.2, position = position_nudge(0.1)) +
  geom_point(aes(color = 'TMB'), position = position_nudge(0.1)) +
  geom_point(aes(x = Stock, y = PA_SREP, color='Parken'), position = position_nudge(-0.1)) +
  geom_errorbar(aes(x = Stock, ymax = PA_UL_SREP, ymin = PA_LL_SREP, color='Parken'), width = 0.2, position = position_nudge(-0.1)) +
  theme_classic() +
  ylab(TeX("$S_{REP}$ Estimate")) +
  theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust = 1)) +
  scale_color_manual(name='Model',
                     breaks=c('Parken', 'TMB'),
                     values=c('Parken'='black', 'TMB'='red'))
dev.off()

# Point-wise-comparison - IWAM vs. Parken - Synoptic stocks ####
# Parken_est <- read.csv(here::here("DataIn/Parkenest.csv"))

# Plot

# Table of SGEN, SMSY, and SREP for target WCVI ####
  # What file can I pull code from?
# Prepare Data 
SREP <- iwam_default$dfout %>%
  filter(RP=='SREP') %>%
  rename('Lower Quantile'=lwr, 'SREP'=Value, 'Upper Quantile'=upr) %>%
  mutate(RP = NULL) %>%
  relocate(Stock)
SMSY <- iwam_default$dfout %>%
  filter(RP=='SMSY') %>%
  rename('Lower Quantile'=lwr, 'SMSY'=Value, 'Upper Quantile'=upr) %>%
  mutate(RP = NULL) %>%
  relocate(Stock)
SGEN <- iwam_default$dfout %>%
  filter(RP=='SGEN') %>%
  rename('Lower Quantile'=lwr, 'SGEN'=Value, 'Upper Quantile'=upr) %>%
  mutate(RP = NULL) %>%
  relocate(Stock)

# WAin <- read.csv(here::here("DataIn/Backcalc_targetstocks_NoAgg.csv")) %>%
WAin <- read.csv(here::here("DataIn/WCVIStocks_NoAgg.csv")) %>%
  # filter(Stock != "Cypre") %>% # remove Cypre - see Get_LRP_bs.R
  mutate(WA = round(WA,0))

complete <- data.frame(SGEN, SREP, SMSY) %>%
  # confirm that all stocks line-up
  select(-Stock.1, -Stock.2, -WA.1, -WA.2) %>%
  rename("SGEN LQ" = Lower.Quantile, "SGEN UQ" = Upper.Quantile, "SREP LQ" = Lower.Quantile.1, "SREP UQ" = Upper.Quantile.1,
         "SMSY LQ" = Lower.Quantile.2, "SMSY UQ" = Upper.Quantile.2)

complete <- complete %>%
  relocate(WA, .after=Stock) %>%
  mutate_at(vars(SGEN, "SGEN LQ", "SGEN UQ", SMSY, "SMSY UQ", "SMSY LQ"), round, 0)
# rename("LH" = type)
# arrange(Stock)

# kable table
kable(complete, caption = "IWAM Smax Model: SGEN, SREP, and SMSY Estimates for All Stocks",
      format.args = list(big.mark = ","))

# export dataframe to csv - and then transfer excel table into word
write.csv(complete, here::here("DataOut/TableofEstimates_kable_IWAM.csv"), row.names = FALSE)

# Bootstrap Convergence: ####
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



# LambertW function Testing ####
compile(here::here("TMB_Files/lambert.cpp"))
dyn.load(dynlib(here::here("TMB_Files/lambert")))

objlw <- MakeADFun(data=list(), parameters=list(x=1), DLL="lambert")
objlw$fn(7 * exp(7))
lambert_W0(7 * exp(7)) # eg

#
temp <- IWAM_func(WAin = c("DataIn/WCVIStocks_NoAgg.csv"),
          bs_nBS = 10,
          plot = FALSE,
          mod = "IWAM_Liermann_srep")
