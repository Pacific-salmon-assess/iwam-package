
# Testing multiple distributions with several IWAM runs looped

# Libaries ####
library(RTMB)
library(TMB)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(progress)
library(diffr)
library(knitr)
library(viridis)

source(here::here("R/PDF_plot.r"))
source(here::here("R/IWAM_model.R"))

# source(here::here("R/Liermann_RTMB_model.R")) # RTMB Liermann model function
  # Contains: derived_post.R source
# source(here::here("R/IWAMsrep_RTMB_model.R")) # Original SREP RTMB model for comparison

# Functions
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

# IWAM Initial run ####
  # Default was -0.412 (log-scale) = 0.66
  # Therefore try 0.4, 0.6, 0.8, 1.0 <- c(-0.9162907,  -0.5108256, -0.2231436, 0)

# DataIn/WCVIStocks_NoAgg.csv

iwam_big <- IWAM_func(WAinraw = c("DataIn/Parken_evalstocks.csv"),
                    targetname = "static_ver",
                    bs_nBS = 1000,
                    static_nusigma = c(-5), # test value of nu sigma)
                    mod = "IWAM_static"
)

# # 
iwam_big_logSigmaA <- logSigmaA_SAVED
# # plot logSigmaA against static values plotted and default
# # 

SMSY_tmb_b <- iwam_big$dfout %>%
  filter(RP=='SMSY') %>% 
  rename("SMSY_i" = Value)
SREP_tmb_b <- iwam_big$dfout %>%
  filter(RP=='SREP') %>% 
  rename("SREP_i" = Value)

iwam_v <- IWAM_func(WAinraw = c("DataIn/Parken_evalstocks.csv"),
                     targetname = "static_ver",
                     bs_nBS = 1000,
                     static_nusigma = c(-0.412), # test value of nu sigma)
                     mod = "IWAM_static"
)

iwam_v_logSigmaA <- logSigmaA_SAVED

SMSY_tmb_v <- iwam_v$dfout %>%
  filter(RP=='SMSY') %>% 
  rename("SMSY_i" = Value)
SREP_tmb_v <- iwam_v$dfout %>%
  filter(RP=='SREP') %>% 
  rename("SREP_i" = Value)

iwam_v.4 <- IWAM_func(WAinraw = c("DataIn/Parken_evalstocks.csv"),
                     targetname = "static_ver",
                     bs_nBS = 1000,
                     static_nusigma = c(-0.9162907), # test value of nu sigma)
                     mod = "IWAM_static"
)

iwam_v.4_logSigmaA <- logSigmaA_SAVED

SMSY_tmb_v.4 <- iwam_v.4$dfout %>%
  filter(RP=='SMSY') %>% 
  rename("SMSY_i" = Value)
SREP_tmb_v.4 <- iwam_v.4$dfout %>%
  filter(RP=='SREP') %>% 
  rename("SREP_i" = Value)

iwam_v.6 <- IWAM_func(WAinraw = c("DataIn/Parken_evalstocks.csv"),
                     targetname = "static_ver",
                     bs_nBS = 1000,
                     static_nusigma = c(-0.5108256), # test value of nu sigma)
                     mod = "IWAM_static"
)

iwam_v.6_logSigmaA <- logSigmaA_SAVED

SMSY_tmb_v.6 <- iwam_v.6$dfout %>%
  filter(RP=='SMSY') %>% 
  rename("SMSY_i" = Value)
SREP_tmb_v.6 <- iwam_v.6$dfout %>%
  filter(RP=='SREP') %>% 
  rename("SREP_i" = Value)

iwam_v.8 <- IWAM_func(WAinraw = c("DataIn/Parken_evalstocks.csv"),
                     targetname = "static_ver",
                     bs_nBS = 1000,
                     static_nusigma = c(-0.2231436), # test value of nu sigma)
                     mod = "IWAM_static"
)

iwam_v.8_logSigmaA <- logSigmaA_SAVED

SMSY_tmb_v.8 <- iwam_v.8$dfout %>%
  filter(RP=='SMSY') %>% 
  rename("SMSY_i" = Value)
SREP_tmb_v.8 <- iwam_v.8$dfout %>%
  filter(RP=='SREP') %>% 
  rename("SREP_i" = Value)

  # change to 1 instead of 0
iwam_v1 <- IWAM_func(WAinraw = c("DataIn/Parken_evalstocks.csv"),
                     targetname = "static_ver",
                     bs_nBS = 1000,
                     static_nusigma = c(0), # test value of nu sigma)
                     mod = "IWAM_static"
)

iwam_v1_logSigmaA <- logSigmaA_SAVED

SMSY_tmb_v1 <- iwam_v1$dfout %>%
  filter(RP=='SMSY') %>% 
  rename("SMSY_i" = Value)
SREP_tmb_v1 <- iwam_v1$dfout %>%
  filter(RP=='SREP') %>% 
  rename("SREP_i" = Value)

iwam_default <- IWAM_func(WAinraw = c("DataIn/Parken_evalstocks.csv"),
                          targetname = "test",
                          bs_nBS = 1000, 
                          SigRicPrior = c(F, F, F, T),
                          SigDeltaPrior = c(F, F, F, T, F),
                          # TauPrior = c(7.5, 0.1, 3, 1, 0.75), # Currently set in hard
                          TauDist = c(0.1, 1), # [1] is Ricker, [2] is WA
                          mod = "IWAM_Liermann"
)
  
SMSY_tmb_d <- iwam_default$dfout %>%
  filter(RP=='SMSY') %>% 
  rename("SMSY_i" = Value)
SREP_tmb_d <- iwam_default$dfout %>%
  filter(RP=='SREP') %>% 
  rename("SREP_i" = Value)

# Plotting ####
  # Run all of the above with: "DataIn/Parken_evalstocks.csv"

## Plotting ####
parken <- read.csv(here::here("DataIn/Parken_evalstocks.csv"))

parken <- parken |> 
  rename("SMSY_i" = SMSY, "lwr" = SMSY_5, "upr" = SMSY_95)

cols <- viridis(8, alpha=0.9, option = "mako", direction = -1)

ggplot() + # SMSY_tmb_v, aes(x = Stock, y = SMSY_i)
  
  # geom_errorbar(data = SMSY_tmb_v, aes(x = Stock, y = SMSY_i, ymax = lwr, ymin = upr,
  #                                      color = 'Static 0.66',
  #                                      width=.1)) +
  # geom_point(data = SMSY_tmb_v,
  #            aes(x = Stock, y = SMSY_i, color = 'Static 0.66')) +
  
  
  # geom_errorbar(data = SMSY_tmb_v.4, aes(x = Stock, y = SMSY_i, ymax = lwr, ymin = upr,
  #                                        color = "Static 0.4",
  #                                        width=.1),
  #               position = position_nudge(0.1)) +
  # geom_point(data = SMSY_tmb_v.4,
  #            position = position_nudge(0.1),
  #            aes(x = Stock, y = SMSY_i, color = "Static 0.4")) +
  
  # geom_errorbar(data = SMSY_tmb_v.8, aes(x = Stock, y = SMSY_i, ymax = lwr, ymin = upr,
  #                                        color = "Static 0.8",
  #                                        width=.1),
  #               position = position_nudge(-0.1)) +
  # geom_point(data = SMSY_tmb_v.8,
  #            position = position_nudge(-0.1),
  #            aes(x = Stock, y = SMSY_i, color = "Static 0.8")) + 
  
  # geom_errorbar(data = SMSY_tmb_v1, aes(x = Stock, y = SMSY_i, ymax = lwr, ymin = upr,
  #                                       color = "Static 1",
  #                                       width=.1),
  #               position = position_nudge(-0.2)) +
  # geom_point(data = SMSY_tmb_v1,
  #            position = position_nudge(-0.2),
  #            aes(x = Stock, y = SMSY_i, color = "Static 1")) +
  
  # geom_errorbar(data = SMSY_tmb_d, aes(x = Stock, y = SMSY_i, ymax = lwr, ymin = upr,
  #                                       color = "IWAM Default",
  #                                      width=.1),
  #               position = position_nudge(-0.3)) +
  # geom_point(data = SMSY_tmb_d,
  #            position = position_nudge(-0.3),
  #            aes(x = Stock, y = SMSY_i, color = "IWAM Default")) +
  
  # geom_errorbar(data = SMSY_tmb_b, aes(x = Stock, y = SMSY_i, ymax = lwr, ymin = upr,
  #                                      color = "close to zero",
  #                                      width=.1),
  #               position = position_nudge(.2)) +
  # geom_point(data = SMSY_tmb_b,
  #            position = position_nudge(.2),
  #            aes(x = Stock, y = SMSY_i, color = "close to zero")) +
  
  # Re-written as the Parken model is included in WAin - bound to rtmb MLE version
  geom_errorbar(data = rtmb_full, aes(x = Stock_name, y = SMSY, ymax = SMSY_95, ymin = SMSY_5,
                                       color = "Parken",
                                       width=.1),
                # position = position_nudge(-0.4)
    ) +
  geom_point(data = rtmb_full,
             # position = position_nudge(-0.4),
             aes(x = Stock_name, y = SMSY, color = "Parken")) +
  
  # Add in RTMB from IWAMsrep_RTMB_model.R as a global object (run internal - function is broken)
    # otherwise comment this one out
  # rtmb <- no nll with predictions
  # rtmb_full <- original
  geom_errorbar(data = testfull, aes(x = Stock_name, y = SMSY.Mean, ymax = SMSY.UQ, ymin = SMSY.LQ,
                                   color = "RTMB MLE",
                                   width=.1),
                position = position_nudge(-0.2)) +
  geom_point(data = testfull,
             position = position_nudge(-0.2),
             aes(x = Stock_name, y = SMSY.Mean, color = "RTMB MLE")) +
  
  # Add in LIERMANN from Liermann_RTMB_model.R as a global object
    # otherwise comment this one out
    # Global object transform: derived_obj --> derived_obj$deripost_summary$SMSY
  geom_errorbar(data = testfull, aes(x = Stock_name,
                                     y = Mean,
                                     ymax = UQ_95, 
                                     ymin = LQ_5,
                                 color = "Liermann MCMC",
                                 width=.1),
                position = position_nudge(+0.2)) +
  geom_point(data = testfull,
             position = position_nudge(+0.2),
             aes(x = Stock_name, y = Mean, color = "Liermann MCMC")) +
  
  theme_classic() +
  scale_y_continuous(transform = "log", 
                     breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000, 10000000)) +
  ylab(TeX("$S_{MSY}$ Estimate")) +
  xlab("Stock Name") + 
  theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust = 1)) +
  scale_color_manual(name='Model',
                     breaks=c(# 'Static 0.4',
                              # 'Static 0.66',
                              # 'Static 0.8',
                              # 'Static 1',
                              # 'IWAM Default',
                              # 'close to zero',
                              'Parken',
                              'RTMB MLE',
                              'Liermann MCMC'),
                     values=c(# 'Static 0.4' = cols[1],
                              # 'Static 0.66' = cols[2],
                              # 'Static 0.8' = cols[3],
                              # 'Static 1' = cols[4],
                              # 'close to zero' = cols[5],
                              # 'IWAM Default' = "orange",
                              'Parken' = "black",
                              'RTMB MLE' = "orange",
                              'Liermann MCMC' = "skyblue"))




# ggplot(benchmarks, aes(x=Stock, y = SMSY_i)) +
#   # error bars for default IWAM model
#   geom_errorbar(aes(ymax = `UL SMSY`, ymin = `LL SMSY`, 
#                     color='Integrated Watershed Area Model',), 
#                 width = 0.2, position = position_nudge(0.1)) +
#   # points for default iwam model
#   geom_point(aes(color = 'Integrated Watershed Area Model'), 
#              position = position_nudge(0.1)) +
#   # points for parken est.
#   geom_point(aes(x = Stock, y = SMSY, color='Parken et al. 2006'), 
#              position = position_nudge(-0.1)) +
#   # errorbars for parken est.
#   geom_errorbar(aes(x = Stock, ymax = SMSY_95, ymin = SMSY_5, 
#                     color='Parken et al. 2006'), width = 0.2, 
#                 position = position_nudge(-0.1)) +
#   
#   theme_classic() +
#   
#   scale_y_continuous(transform = "log", 
#                      breaks = c(0, 10, 100, 1000, 10000, 100000)) +
#   
#   ylab(TeX("$S_{MSY}$ Estimate")) +
#   
#   xlab("Stock Name") + 
#   
#   theme(axis.text.x = element_text(angle = 90, vjust=0.3, hjust = 1)) +
#   
#   scale_color_manual(name='Model',
#                      breaks=c('Parken et al. 2006', 
#                               'Integrated Watershed Area Model'),
#                      values=c('Parken et al. 2006'='black', 
#                               'Integrated Watershed Area Model'='red'))