# Summary ----------------------------------------------------------------------
# Code to estimate LRPs for WCVI CK from watershed-area based Sgen by
# bootstrapping from SREP estimates from the watershed-area model and Ricker 
# a values from a plausible range derived from expert opinion and a 
# life-history model
# ---------------------------------------------------------------------------- #
# Libraries and Functions ------------------------------------------------------
# ---------------------------------------------------------------------------- #
library(tidyverse)
library(ggplot2)
library(gsl)
#library(TMB) # Required if using: run_logReg
library(viridis)

# library(future)
# library(furrr)
# future::plan(multisession, workers = availableCores() - 1)
# future::plan(sequential)
# plan(multisession, workers = 8)

# Functions
source (here::here("R/helperFunctions.R"))

# Function to estimate LRPs for WCVI CK ----------------------------------------
# Arguments; 
# remove.EnhStocks = A logical reflecting if enhanced stock are to be 
# included
# prod = character specifying which assumption about productivity is made,
# either "LifeStageModel" (default), where productivity is derived life-
# stage model with expert opinion (W. LUedke pers. comm.) or from a run
# reconstruction assuming same harvest rates across WCVI Chinook stocks 
# estimated from Robertson Creek Hatchery fish (D. Dobson, pers. comm.) 
# Bern_logistic = logical (TRUE/FALSE), indicating if a Bernoulli logistic 
# regression is used to estimate LRPs based on aggregate abundances 
# (TRUE, default) or if binomial logistic regression is used (FALSE)
# LOO = numeric for leave-one-out cross validation of the logistic regression
# This number is the index of the time-series of ppn of CUs and aggregate 
# abundances that are removed prior to implementing the logistic regression 
# in TMB. Set to NA as default (no values removed). Note, the outputted 
# time-series ('out') contain all the data, but parameter estimates are 
# derived from time-series without LOO index value
# The code requires that the IWAM model has been run and 
# "WCVI_SMSY_noEnh.csv" or "WCVI_SMSY_wEnh.csv" exist
# Returns:
# csv file, DataOut/wcviRPs_noEnh.csv or DataOut/wcviRPs_wEnh.csv of stock, 
# inlet, and CU level Sgen, adjusted SMSY values(adjusted for expert derived 
# productivity) and SREP values from integrated watershed-area model
# Dataframe $out
# Dataframe $WCVIEsc
# Dataframe $SMU_Esc
# Dataframe $CU_Status
# Dataframe $SMU_ppn

# Current wcvi case study example uses 
# remove.EnhStocks <- FALSE
# Bern_logistic <- FALSE
# prod <- "LifeStageModel"
# LOO <- NA
# run_logReg <- FALSE
# run.bootstraps <- TRUE

# datain <- c("DataOut/dataout_target_ocean_noEnh.csv") # RUNS CLEAN because NO CU INFORMATION to overlap
# datain <- c("DataOut/dataout_target_wEnh.csv") # ERRORS BECAUSE CU DUPLICATES CAUSED BY TWO LH's

Get.LRP.bs <- function(datain = "DataOut/dataout_target_ocean_noEnh.csv", # file name/path of output of IWAM Model
                       dataraw = WAinraw, # should be whatever the original input in IWAM_function is e.g. WAinraw
                       remove.EnhStocks = TRUE,  
                       Bern_logistic = FALSE, 
                       bias.cor = TRUE, # add in bias.cor on/off from IWAM_model inputs
                       prod = "LifeStageModel", 
                       LOO = NA, 
                       run_logReg = FALSE){
  
  #--------------------------------------------------------------------------- #
  # Read in watershed area-based reference points (SREP and SMSY) --------------
  #--------------------------------------------------------------------------- #
  # *TOR*: updated all object naming to be non-specific with wcvi
    # E.g.: wcviRPs_long to RPs_long, wcviRPs to RPs
  
  # Core data: 
  if(remove.EnhStocks) RPs_long <- read.csv(here::here(datain))
  if(!remove.EnhStocks) RPs_long <- read.csv(here::here(datain))
  #if(!remove.EnhStocks) RPs_long<- read.csv(here::here("DataOut/dataout_target_ocean_wEnh.csv"))
  # Old files:
  # if (remove.EnhStocks) RPs_long <- read.csv("DataOut/WCVI_SMSY_noEnh_wBC.csv")
  # if (!remove.EnhStocks) RPs_long <- read.csv("DataOut/WCVI_SMSY_wEnh_wBC.csv")
  
  # Remove Cypre as it's not a core indicator (Diana McHugh, 22 Oct 2020)
  stock_SMSY <- RPs_long %>% filter(Stock != "Cypre") %>%  # ********************************************** CYPRE FLAG
    filter (Param == "SMSY") %>% 
    rename(SMSY=Estimate, SMSYLL=LL, SMSYUL=UL) %>% 
    dplyr::select (-Param, -X) #, -CU)
  stock_SREP <- RPs_long %>% filter(Stock != "Cypre") %>% # ********************************************** CYPRE FLAG
    filter (Param == "SREP") %>% 
    rename(SREP=Estimate, SREPLL=LL, SREPUL=UL) %>% 
    dplyr::select (-Param, -X)
  
  # Different cleaning tech - Is this still needed?
    # *TOR*: I think the issue here was fixed through a change to IWAM_model
    # to eliminate duplicate CU's with both life histories.
  # RPs_short <- RPs_long %>% 
  #   # filter(Stock != "Cypre") %>% 
  #   mutate(LH = case_when(grepl("ocean", X) ~ "ocean", grepl("stream", X) ~ "stream", TRUE ~ "other")) %>% 
  #   pivot_wider(id_cols = c(Stock, LH), names_from = c(Param), values_from = c(Estimate, LL, UL)) %>%
  #   rename(SMSY=Estimate_SMSY, SREP=Estimate_SREP, SMSYLL=LL_SMSY, SMSYUL=UL_SMSY, SREPLL=LL_SREP, SREPUL=UL_SREP)
  
  # RPs <- RPs_short - for testing
  
  # **********************************************************************************************
  RPs <- stock_SMSY %>% left_join(stock_SREP, by="Stock") # ERROR CAUSED BY DUPLICATES AFTER SPLIT
  # **********************************************************************************************
  
  # Calculate scale for each stock
  digits <- count.dig(stock_SMSY$SMSY)
  # Scale <- 10^(digits) # Not the same as what is used for the IWAM model?
  Scale <- 10^(digits-1)
    # IWAM USES: 10^(maxdigits-1)
  
  #SREP_SE <- RPs %>% mutate(SE = ((RPs$SREP) - (RPs$SREPLL)) / 1.96) **********
  SREP_logSE <- RPs %>% mutate(SE = (log(RPs$SREP) - log(RPs$SREPLL)) / 1.96)
  # print(SREP_logSE$SE)
  # The UpperLimit-MLE gives same answer
  #SREP_logSE <- RPs %>% mutate(SE = (log(RPs$SREPUL) - log(RPs$SREP)) / 1.96)
  SREP_logSE <- SREP_logSE %>% dplyr::select(Stock, SE)
  
  #--------------------------------------------------------------------------- #
  # Calculate Sgen 2 ways ------------------------------------------------------
  #--------------------------------------------------------------------------- #
  
  # 1.Use Watershed-area SMSY and SREP to estimate Sgen (assuming productivity
  # SGENcalcs <- purrr::map2_dfr (RPs$SMSY/Scale,RPs$SREP/Scale, Sgen.fn) 
  
  # 2. Assume independent estimate of productivity and watershed-area 
  # estimate of SREP
  
  # For base case assume no variablity in RicA among stocks. 
  # Add variability in Ric.A when drawing MC samples from prediction intervals 
  # of WA model so that Ric.A is drawn multiple times for each stock from 
  # rnorm distribution
  
  # There are two assumptions about productivity, (1) from life-stage model with 
  # expert opinion (W. Luedke pers. comm.), and (2) from "RunReconstruction
  # which assumes same harvest rates across wCVI Chinook stocks (D. Dobson 
  # pers. comm.)  The default is the life-stage model with expert opinion (1)
  
  # Lower estimate of Ricker a derived from life-stage model (Luedke pers.
  # comm.) 
  # DEFAULT
  # prod <-  "LifeStageModel"
  if(prod == "LifeStageModel"){
    # print(paste("Productivity Assumption running: LifeStageModel"))
    Mean.Ric.A <- 1 # Derived from life-history model (Luedke pers.comm.) and 
    # WCVI CK run reconstruction SR analysis (Dobson pers. comm.)
    Ric.A <- exp(rnorm(length(Scale), Mean.Ric.A, 0))
    
    # When incorporating uncertainty in Ricker A:
    Sig.Ric.A <- 0.51 #0.255 #0.51 for a wider plausible bound
    # Sig.Ric.A derived from 95% CL of lower and upper plausible limits = 
    # 0.5 logA - 1.5 logA (Luedke pers. comm. Dec 2020)
    # See distribution below:
    # test <- seq(0,4, len=40)
    # plot(x=test, y=dnorm(test, 1,0.255), type="l", xlab="LogA", 
    # ylab="Probability Density", ylim=c(0,5))
    # # With this sigma, 95% of probablity density is within bounds mean 
    # +/- 0.50 (assuming range 0.5-1.5, mean=1). 0.255*1.96 = 0.50
    # lines(x=test, y=dnorm(test, 1,0.51))# With this sigma, 95% of probablity 
    # density is within bounds mean +/- 1.0 
    # (assuming range 0-2.0, mean=1). 0.510*1.96 = 1.0
    
    Ric.A <- exp(rnorm(length(Scale), Mean.Ric.A, Sig.Ric.A))
    if(min(Ric.A)<=0) Ric.A <- exp(rnorm(length(Scale), Mean.Ric.A, Sig.Ric.A))
    
    if (bias.cor == TRUE) {
      # print("Bias correction added in bootstrapping")
      sREP <- exp(rnorm(length(Scale), log(RPs$SREP) - 0.5*SREP_logSE$SE^2, SREP_logSE$SE)) # **
      if(min(sREP)<=0)   sREP <- exp(rnorm(length(Scale), log(RPs$SREP) - 0.5*SREP_logSE$SE^2,
        SREP_logSE$SE))
    } else {
      sREP <- exp(rnorm(length(Scale), log(RPs$SREP), SREP_logSE$SE))
      if(min(sREP)<=0)   sREP <- exp(rnorm(length(Scale), RPs$SREP, SREP_SE$SE)) # SREP_SE doesn't exist ***********
    }
    # sREP <- exp(rnorm(length(Scale), log(RPs$SREP), SREP_logSE$SE)) # **
    # sREP <- exp(rnorm(length(Scale), log(RPs$SREP) - 0.5*SREP_logSE$SE^2, SREP_logSE$SE)) # **
    
    # if(min(sREP)<=0)   sREP <- exp(rnorm(length(Scale), RPs$SREP,
    #                                      SREP_SE$SE))
    # if(min(sREP)<=0)   sREP <- exp(rnorm(length(Scale), log(RPs$SREP) - 0.5*SREP_logSE$SE^2,
    #                                      SREP_SE$SE))

    SGENcalcs <- purrr::map2_dfr (Ric.A, sREP/Scale, Sgen.fn2)
      # what are the default parameters of sgen.fn2?
      # explicit=TRUE by default
    # SGENcalcs <- future_map2_dfr(Ric.A, sREP/Scale, Sgen.fn2)


    RPs <- RPs %>% mutate (SGEN = SGENcalcs$SGEN) %>%
      mutate(SGEN=round(SGEN*Scale,0))
    RPs <- RPs %>% mutate (a.par = SGENcalcs$apar) %>% 
      mutate(a.par=round(a.par,2))
    RPs <- RPs %>% mutate (SMSY = SGENcalcs$SMSY) %>% 
      mutate(SMSY=round(SMSY*Scale,0))
    
    RPs <- RPs[c("Stock", "SGEN", "SMSY", "SMSYLL", "SMSYUL", "SREP", 
                         "SREPLL", "SREPUL", "a.par")] #"CU"
    
  }#End of if(prod == "LifeStageModel6")
  
  # Ricker a's from Diana Dobson's Run Reconstruction (pers.comm) coded in TMB
  # Higher estimate of Ricker a (lower Sgen)  
  
  if(prod == "RunReconstruction"){
    # print(paste("Productivity Assumption running: RunReconstruction"))
    lnalpha_inlet <- read.csv("DataIn/CUPars_wBC.csv") %>% 
      select(alpha,stkName) %>% rename(inlets=stkName, lnalpha=alpha)
    lnalpha_nBC_inlet <- read.csv("DataIn/CUPars_nBC.csv") %>% 
      select(alpha,stkName) %>% rename(inlets=stkName, lnalpha_nBC=alpha)
    targetstocks <- read.csv("DataIn/WCVIStocks.csv") %>% # Previously WCVIStocks - should this be the same as "datain"?
        # now named WAinraw
      # filter (Stock != "Cypre") %>% 
      rename(inlets=Inlet)
    Ric.A <- lnalpha_inlet %>% left_join(targetstocks, by="inlets") %>% select(c(lnalpha,inlets,CU,Stock))
    
    RPs <- RPs %>% left_join(Ric.A) %>% mutate(a.RR=exp(lnalpha))
    RPs[RPs$Stock=="Nitinat",]$a.RR <- exp(1)
    RPs[RPs$Stock=="San Juan",]$a.RR <- exp(1)
    RPs[RPs$Stock=="Nitinat",]$a.RR <- exp(1)
    
    RPs[RPs$Stock=="Barkley",]$a.RR <- 
      RPs[RPs$inlets=="Barkley",]$a.RR[1]
    RPs[RPs$Stock=="Clayoquot",]$a.RR <- 
      RPs[RPs$inlets=="Clayoquot",]$a.RR[1]
    RPs[RPs$Stock=="Kyuquot",]$a.RR <- 
      RPs[RPs$inlets=="Kyuquot",]$a.RR[1]
    RPs[RPs$Stock=="Nootka/Esperanza",]$a.RR <- 
      RPs[RPs$inlets=="Nootka/Esperanza",]$a.RR[1]
    RPs[RPs$Stock=="Quatsino",]$a.RR <- 
      RPs[RPs$inlets=="Quatsino",]$a.RR[1]
    RPs[RPs$Stock=="WCVI South",]$a.RR <- 
      RPs[RPs$inlets=="Barkley",]$a.RR[1]
    RPs[RPs$Stock=="WCVI Nootka & Kyuquot",]$a.RR <- 
      RPs[RPs$inlets=="Nootka/Esperanza",]$a.RR[1]
    RPs[RPs$Stock=="WCVI North",]$a.RR <- 
      RPs[RPs$inlets=="Quatsino",]$a.RR[1]
    
    RPs <- RPs %>% select(-c(inlets, CU, lnalpha)) %>% rename(a.par=a.RR)
    
    # When incorporating uncertainty in Ricker A:
    Sig.Ric.A <- 0.51 #0.255 for a narrower plausible bound
    # Sig.Ric.A derived from 95% CL of lower and upper plausible limits = 
    # 0.5 logA - 1.5 logA (Luedke pers. comm. Dec 2020)
    # See distribution below:
    # test <- seq(0,4, len=40)
    # plot(x=test, y=dnorm(test, 1,0.255), type="l", xlab="LogA", 
    # ylab="Probability Density", ylim=c(0,5))
    # # With this sigma, 95% of probablity density is within bounds mean 
    # +/- 0.50 (assuming range 0.5-1.5, mean=1). 0.255*1.96 = 0.50
    # lines(x=test, y=dnorm(test, 1,0.51))# With this sigma, 95% of probablity 
    # density is within bounds mean +/- 1.0 
    # (assuming range 0-2.0, mean=1). 0.510*1.96 = 1.0
    
    
    
    Ric.A.hi <- exp(rnorm(length(Scale), log(RPs$a.par), Sig.Ric.A))
    if(min(Ric.A.hi)<0) Ric.A <- exp(rnorm(length(Scale), RPs$a.RR, Sig.Ric.A))
    
    sREP <- exp(rnorm(length(Scale), log(RPs$SREP), SREP_logSE$SE))
    if(min(sREP)<0)   sREP <- exp(rnorm(length(Scale), RPs$SREP, 
                                        SREP_SE$SE))
    
    SGENcalcs <- purrr::map2_dfr (Ric.A.hi, sREP/Scale, Sgen.fn2)
    
    RPs <- RPs %>% mutate (SGEN = SGENcalcs$SGEN) %>% # POSSIBLE ERROR
      mutate(SGEN=round(SGEN*Scale,0))
    RPs <- RPs %>% mutate (a.par = SGENcalcs$apar) %>% 
      mutate(a.par=round(a.par,2))
    RPs <- RPs %>% mutate (SMSY = SGENcalcs$SMSY) %>% 
      mutate(SMSY=round(SMSY*Scale,0))
    
    RPs <- RPs[c("Stock", "SGEN", "SMSY", "SMSYLL", "SMSYUL", "SREP", 
                         "SREPLL", "SREPUL", "a.par")] #"CU"
    
  } # if(prod == "RunReconstruction"){
  
  # PARKEN Method ####
  if(prod == "Parken"){
    # RPs_e <- RPs # save a version of RPs for explicit calculation
    
    # print(paste("Productivity Assumption running: Parken"))
    est_loga <- function(SMSY, SREP, shortloga=FALSE){
      
      loga <- nlminb(start = (0.5 - SMSY/SREP) / 0.07, 
                     objective = calc_loga, # Try to fix with a Scheurel version of LW if possible
                     SMSY= SMSY, 
                     SREP=SREP)$par
      if(shortloga) loga <- (0.5 - SMSY/SREP) / 0.07
      beta <- loga/SREP
      return( list( loga = loga , beta = beta, SMSY = SMSY, SREP = SREP) )
    }
    
    # Re-label csv's to make sure they are pulling the right ones
      # See above
    # if(!ExtInd) {
    # WCVIStocks <- read.csv(here::here(dataraw)) %>% 
    #     filter (Stock != "Cypre") # %>% rename(inlets=Inlet)
          # removed Inlet naming - as most of this code is done without CU's and Inlets
      # if (remove.EnhStocks) wcviRPs_long <- 
      #     read.csv("DataOut/WCVI_SMSY_noEnh_wBC.csv")
      # if (!remove.EnhStocks) wcviRPs_long <- 
      #     read.csv("DataOut/WCVI_SMSY_wEnh_wBC.csv")
    if(remove.EnhStocks) RPs_long <- read.csv(here::here(datain))
    if(!remove.EnhStocks) RPs_long <- read.csv(here::here(datain))
    # }
    
    # TK: I think ExtInd is FALSE for these runs
    # if(ExtInd) {
    #   WCVIStocks <- read.csv("DataIn/WCVIStocks_ExtInd.csv") %>% 
    #     rename(inlets=Inlet)
    #   wcviRPs_long <- read.csv("DataOut/WCVI_SMSY_ExtInd.csv")
    # }
    
    SMSY <- RPs_long %>% filter(Param == "SMSY") %>% select(Estimate) 
    SREP <- RPs_long %>% filter(Param == "SREP") %>% select(Estimate) 
    lnalpha_Parkin <- purrr::map2_dfr (SMSY, SREP, shortloga=FALSE, 
                                       est_loga)
    
    # Explicit solution for logA and beta
    # SREP_e <- as.numeric(SREP$Estimate)
    # SMSY_e <- as.numeric(SMSY$Estimate)
    # loga_e <- SREP_e*(SMSY_e*gsl::lambert_W0(-exp(1-SREP_e/SMSY_e)*(SREP_e-SMSY_e)/SMSY_e) + 
                # SREP_e - SMSY_e)/(SMSY_e*(SREP_e-SMSY_e))
    # beta_e <- loga_e/SREP_e
    
    if (bias.cor == TRUE) {
      # print("Bias correction added in bootstrapping")
          sREP <- exp(rnorm(length(Scale), log(RPs$SREP) - 0.5*SREP_logSE$SE^2,
                      SREP_logSE$SE))
      if(min(sREP)<0)    sREP <- exp(rnorm(length(Scale), log(RPs$SREP) -
                                           0.5*SREP_logSE$SE^2, SREP_logSE$SE))
    } else {
      sREP <- exp(rnorm(length(Scale), log(RPs$SREP), SREP_logSE$SE))
      if(min(sREP)<0)   sREP <- exp(rnorm(length(Scale), RPs$SREP, SREP_SE$SE)) # SREP_SE doesn't exist *************
    }
    # Without log-normal bias correction when sampling log-beta
    # sREP <- exp(rnorm(length(Scale), log(RPs$SREP), SREP_logSE$SE))
    # With a log-normal bias correction when sampling log-beta
    # sREP <- exp(rnorm(length(Scale), log(RPs$SREP) - 0.5*SREP_logSE$SE^2,
    #                   SREP_logSE$SE))
    # if(min(sREP)<0)   sREP <- exp(rnorm(length(Scale), RPs$SREP, SREP_SE$SE))
    # if(min(sREP)<0)    sREP <- exp(rnorm(length(Scale), log(RPs$SREP) -
    #                                        0.5*SREP_logSE$SE^2, SREP_logSE$SE))
    
    SGENcalcs <- purrr::map2_dfr (exp(median(lnalpha_Parkin$loga)), sREP/Scale, Sgen.fn2)
    
    RPs <- RPs %>% mutate (SGEN = SGENcalcs$SGEN) %>% 
      mutate(SGEN=round(SGEN*Scale,0))
    RPs <- RPs %>% mutate (a.par = SGENcalcs$apar) %>% 
      mutate(a.par=round(a.par,2))
    RPs <- RPs %>% mutate (SMSY = SGENcalcs$SMSY) %>% 
      mutate(SMSY=round(SMSY*Scale,0))
    
    # Explicit version
    # SGENcalcs_e <- purrr::map2_dfr (exp(median(loga_e)), sREP/Scale, Sgen.fn2)

    # RPs_e <- RPs_e %>% mutate (SGEN = SGENcalcs_e$SGEN) %>% 
    #   mutate(SGEN=round(SGEN*Scale,0))
    # RPs_e <- RPs_e %>% mutate (a.par = SGENcalcs_e$apar) %>% 
    #   mutate(a.par=round(a.par,2))
    # RPs_e <- RPs_e %>% mutate (SMSY = SGENcalcs_e$SMSY) %>% 
    #   mutate(SMSY=round(SMSY*Scale,0))
    
    RPs <- RPs[c("Stock", "SGEN", "SMSY", "SMSYLL", "SMSYUL", "SREP", 
                         "SREPLL", "SREPUL", "a.par")] #"CU"
    # RPs_e <- RPs_e[c("Stock", "SGEN", "SMSY", "SMSYLL", "SMSYUL", "SREP", 
    #                      "SREPLL", "SREPUL", "a.par")]
  }
  
  #--------------------------------------------------------------------------- #
  # Add Sgen and revised SMSY to RPs data frame ----------------------------
  #--------------------------------------------------------------------------- #
  
  
  
  RPs # Could also output RPs_e?
  # # Write this to a csv file so that it can be called in plotting functions
  # # write.csv(RPs, "DataOut/wcviRPs.csv") # or "targetRPs"
  # if (remove.EnhStocks) write.csv(wcviRPs, "DataOut/wcviRPs_noEnh.csv")
  # if (!remove.EnhStocks) write.csv(wcviRPs, "DataOut/wcviRPs_wEnh.csv")

  
  #--------------------------------------------------------------------------- #
  # Scenario if productivity is reducted by half, as in WSP SBC CK assessment ----
  # (DFO 2014). NOT NEEDED
  #--------------------------------------------------------------------------- #
  # 
  # SGENcalcsv2 <- map2_dfr (RPs$SMSY/Scale,RPs$SREP/Scale, 
  #                          Sgen.fn, half.a = TRUE, const.SMAX = FALSE)
  # RPs <- RPs %>% mutate (SGENha.cSREP = SGENcalcsv2$SGEN) %>% 
  #   mutate( SGENha.cSREP = round( SGENha.cSREP*Scale, 0 ) )
  # RPs <- RPs %>% mutate (SMSYha.cSREP = SGENcalcsv2$SMSY) %>% 
  #   mutate( SMSYha.cSREP = round( SMSYha.cSREP*Scale, 0 ) )
  # RPs <- RPs %>% mutate (SREPha.cSREP = SGENcalcsv2$SREP) %>% 
  #   mutate( SREPha.cSREP = round( SREPha.cSREP*Scale, 0 ) )
  # ###RPs <- RPs %>% mutate (SMAXrev = 1/SGENcalcs$bpar) %>% 
  # ###mutate(SMAXrev=round(SMAXrev,0))
  # 
  # SGENcalcsv3 <- map2_dfr (RPs$SMSY/Scale, RPs$SREP/Scale, 
  #                          Sgen.fn, half.a = TRUE, const.SMAX = TRUE)
  # RPs <- RPs %>% mutate (SGENha.cSMAX = SGENcalcsv3$SGEN) %>% 
  #   mutate( SGENha.cSMAX = round( SGENha.cSMAX*Scale, 0 ) )
  # RPs <- RPs %>% mutate (SMSYha.cSMAX = SGENcalcsv3$SMSY) %>% 
  #   mutate( SMSYha.cSMAX = round( SMSYha.cSMAX*Scale, 0 ) )
  # RPs <- RPs %>% mutate (SREPha.cSMAX = SGENcalcsv3$SREP) %>% 
  #   mutate( SREPha.cSMAX = round( SREPha.cSMAX*Scale, 0 ) )
  
  run_logReg <- FALSE
  if(run_logReg==FALSE){
    return(list(bench= select(SGENcalcs, -apar, -bpar)*Scale))
  }
  #--------------------------------------------------------------------------- #
  # Sum escapements across indicators within inlets ----------------------------
  #--------------------------------------------------------------------------- #
  # Removing run_logReg for the purpose of the IWAM function
  # if(run_logReg==TRUE){
  #   WCVIEsc <- data.frame(read.csv("DataIn/WCVIEsc.csv", row.names="Yr")) %>% 
  #     dplyr::select (-"Little.Zeballos")
  #   
  #   # Take "." out of name as in escapement data
  #   WCVIEsc_names <- sapply(colnames(WCVIEsc), 
  #                           function(x) (gsub(".", " ", x, fixed=TRUE) ) )
  #   WCVIEsc_names <- sapply(WCVIEsc_names, function(x) 
  #     (gsub("Bedwell Ursus", "Bedwell/Ursus", x, fixed=TRUE) ) )
  #   WCVIEsc_names <- sapply(WCVIEsc_names, function(x) 
  #     (gsub("Nootka Esperanza", "Nootka/Esperanza", x, fixed=TRUE) ) )
  #   colnames(WCVIEsc) <- WCVIEsc_names 
  #   
  #   EnhStocks <- data.frame(read.csv("DataIn/WCVIstocks.csv")) %>% filter (Enh==1) %>%
  #     pull(Stock)
  #   EnhStocks <- as.character(EnhStocks)
  #   
  #   #EnhStocks <- c("Burman",  "Conuma", "Leiner", "Nitinat", "Sarita",  
  #   #               "Somass",  "Zeballos", "San Juan", "Tranquil")
  #   # Artlish removed from Enhanced stocks 23 Dec. 2020
  #   # Tranquil added 18 Jan 2021
  #   
  #   
  #   if (remove.EnhStocks) {WCVIEsc <- WCVIEsc %>% dplyr::select(-EnhStocks) }
  #   
  #   Years <- rownames(WCVIEsc)
  #   
  #   # Get stock information for WCVI Chinook & Remove Cypre as it's not an 
  #   # indicator stocks
  #   targetstocks <- read.csv("DataIn/WCVIStocks.csv") %>% 
  #     filter (Stock != "Cypre")
  #   if (remove.EnhStocks) targetstocks <- targetstocks %>% 
  #     filter(Stock %not in% EnhStocks)
  #   
  #   Inlet_Names <- unique(targetstocks$Inlet)
  #   Inlet_Sum <- matrix(NA, nrow=length(Years), ncol=length(Inlet_Names))
  #   colnames(Inlet_Sum) <- Inlet_Names
  #   CU_Names <- unique(targetstocks$CU)
  #   
  #   
  #   # Sum escapements across stocks within inlets
  #   for (i in 1:length(Inlet_Names)) {
  #     # For each inlet, which are the component indicator stocks
  #     Ins <- targetstocks %>% filter(Inlet==Inlet_Names[i]) %>% pull(Stock)
  #     WCVIEsc_Inlets <- matrix(NA, nrow= length(Years), ncol= length(Ins))
  #     
  #     #  Make a matrix of escapements of component indicators
  #     for (j in 1:length(Ins)){
  #       WCVIEsc_Inlets[,j] <- WCVIEsc %>% 
  #         dplyr::select(as.character(Ins[j])) %>% pull()
  #       
  #     }
  #     
  #     # Sum the escapement of component indicators, setting sum=NA for years 
  #     # where there are any NAs
  #     Inlet_Sum[,i] <- apply(WCVIEsc_Inlets, 1, sum, na.rm=F)
  #   }
  #   
  #   
  #   #------------------------------------------------------------------------- #
  #   # Sum escapements across indicators within CUs -----------------------------
  #   #------------------------------------------------------------------------- #
  #   
  #   nCU <- length(unique(targetstocks$CU))
  #   CU_Sum <- matrix(NA, nrow=length(Years), ncol=nCU)
  #   colnames(CU_Sum) <- CU_Names
  #   
  #   for (k in 1:length(CU_Names)) {
  #     # For each CU, which are the component indicators
  #     CUss <- targetstocks %>% filter(CU==CU_Names[k]) %>% pull(Stock)
  #     WCVIEsc_CUs <- matrix(NA, nrow= length(Years), ncol= length(CUss))
  #     
  #     #  Make a matrix of escapements of component indicators
  #     for (j in 1:length(CUss)){
  #       WCVIEsc_CUs[,j] <- WCVIEsc %>% dplyr::select(as.character(CUss[j])) %>% 
  #         pull()
  #     }
  #     
  #     # Sum the escapement of component indicators, setting sum=NA for years 
  #     # where there are any NAs
  #     CU_Sum[,k] <- apply(WCVIEsc_CUs, 1, sum, na.rm=F)
  #   }
  #   # Remove double incidence of Nitinat and San Juan (= stock and an inlet) 
  #   # when enhancement is included
  #   if(!remove.EnhStocks) WCVIEsc <- WCVIEsc %>% 
  #     dplyr::select(-Nitinat, -'San Juan')
  #   
  #   WCVIEsc <- cbind(WCVIEsc, Inlet_Sum, CU_Sum) 
  #   
  #   #------------------------------------------------------------------------- #
  #   # Assess status for each inlet relative to inlet-level SGEN for each year ----
  #   #   Is Inlet level escapement above inlet-level Sgen: T/F?
  #   #   Inlet_Status = FALSE if summed escapement is below Sgen
  #   #   Inlet_Status = TRUE if summed escapement is below Sgen
  #   #------------------------------------------------------------------------- #
  #   
  #   Inlet_Status <- matrix(NA, nrow=length(Years), ncol=length(Inlet_Names) )
  #   colnames(Inlet_Status) <- Inlet_Names
  #   
  #   for (i in 1:length(Inlet_Names)) {
  #     Inlet_Status[,i] <- (Inlet_Sum[,i] > 
  #                            (RPs %>% 
  #                               filter(Stock == as.character(Inlet_Names[i])) %>% 
  #                               pull(SGEN)) )
  #   }
  #   
  #   Inlet_Status <- as.data.frame(Inlet_Status)
  #   
  #   #------------------------------------------------------------------------- #
  #   # Assess status for each CU for each year of the time-series ---------------
  #   #   (floor of summed CU-level numeric statuses)
  #   #   CU_Status = below LB if any inlet within the CU below their Sgen = 0 
  #   #   CU_Status = above LB if all inlets within the CU above their Sgen = 1
  #   #------------------------------------------------------------------------- #
  #   
  #   CU_Status <- matrix(NA, nrow=length(Years), ncol=length(CU_Names))
  #   colnames(CU_Status) <- CU_Names
  #   
  #   stock.LRP <- TRUE
  #   if(stock.LRP){
  #     for (k in 1:length(CU_Names)) {
  #       # For each CU, which are the component indicators
  #       CU_ins <- unique( targetstocks %>% filter(CU==CU_Names[k]) %>% pull(Inlet))
  #       
  #       isAbove <- matrix(NA, nrow= length(Years), ncol= length(CU_ins))
  #       
  #       # Make a matrix of status of component inlets. Is each inlet > 
  #       # Sgen values, T/F?
  #       for (i in 1:length(CU_ins)){
  #         isAbove[,i] <- Inlet_Status %>% 
  #           dplyr::select(as.character(CU_ins[i])) %>% pull()
  #       }
  #       
  #       # CU-level status: are ALL inlets above their Sgen values?
  #       # Sum the "true"
  #       isAboveFun <- function(x){ 
  #         floor(sum( as.numeric(x), na.rm=F) / length(x) ) }
  #       CU_Status[,k] <- apply(X= isAbove, MARGIN = 1, FUN=isAboveFun)
  #     }
  #     CU_Status <- as.data.frame(CU_Status)
  #   }
  #   
  #   # Alternatively, are there CU-level Sgen values to derive CU-level status?
  #   CU.LRP <- FALSE
  #   if(CU.LRP){
  #     for (k in 1:length(CU_Names)){
  #       CU_Status[,k] <- (CU_Sum[,k] > 
  #                           (RPs %>% filter(Stock == 
  #                                                 as.character(CU_Names[k])) %>% 
  #                              pull(SGEN)) )
  #     }
  #     CU_Status <- as.data.frame(CU_Status)
  #     
  #   }
  #   
  #   #------------------------------------------------------------------------- #
  #   # Proportion of CUs that are not in the red zone ---------------------------
  #   #------------------------------------------------------------------------- #
  #   
  #   ppnAboveFun <- function(x) {sum( as.numeric(x), na.rm=F) / length(x) }
  #   SMU_ppn <- apply(X=CU_Status, MARGIN=1, FUN=ppnAboveFun)
  #   
  #   #------------------------------------------------------------------------- #
  #   # Logistic regression ------------------------------------------------------
  #   #------------------------------------------------------------------------- #
  #   # Get SMU-level escapement time-series
  #   SMU_Esc <- apply(Inlet_Sum, 1, sum, na.rm=F)
  #   
  #   SMUlogisticData <- data.frame(SMU_Esc) %>% 
  #     add_column(ppn=SMU_ppn, Years=as.numeric(Years)) %>% 
  #     filter(SMU_Esc != "NA")
  #   
  #   data <- list()
  #   data$N_Stks <- length(CU_Names)
  #   digits <- count.dig(SMU_Esc)
  #   ScaleSMU <- min(10^(digits -1 ), na.rm=T)
  #   
  #   data$LM_Agg_Abund <- SMUlogisticData$SMU_Esc/ScaleSMU
  #   data$N_Above_BM <- SMUlogisticData$ppn * data$N_Stks
  #   
  #   if(!is.na(LOO)) { #If applying leave-one-out cross validation, remove that
  #     #year
  #     data$LM_Agg_Abund <- data$LM_Agg_Abund[-LOO]
  #     data$N_Above_BM <- data$N_Above_BM[-LOO]
  #   }
  #   data$Pred_Abund <- seq(0, max(data$LM_Agg_Abund)*1.1, 0.1)
  #   if(remove.EnhStocks) data$Pred_Abund <- 
  #     seq(0, max(data$LM_Agg_Abund)*1.5, 0.1)
  #   data$p <- 0.95#0.67
  #   
  #   if(Bern_logistic==FALSE) data$Penalty <- as.numeric(TRUE)
  #   if(Bern_logistic==TRUE) data$Penalty <- as.numeric(FALSE)
  #   data$Bern_logistic <- as.numeric(Bern_logistic)
  #   
  #   # Add a normally distributed penalty on aggregate abundances 
  #   # when p is very small (0.01) 
  #   # Lower 95% CL = abundance of the smallest CU in its lowest abundance 
  #   # year (most CUs lost, only 1 remains below LB)
  #   # Upper 95% CL = abundance of the ave annual abundance of the sum of 
  #   # across CUs. 
  #   # If CU-level benchmarks exist can sum those benchmarks for Upper 95%CL
  #   min <- min(apply(CU_Sum, 2, min, na.rm=T), na.rm=T)
  #   max <- mean(apply(CU_Sum, 1, sum, na.rm=F), na.rm=T)
  #   # Parameters for normal penalty (mu, sig):
  #   # mu  = ave of min and max values. 
  #   # sig = SD which allows the density = 0.05 at min and max values
  #   # sum(dnorm(seq(min,max,1), mean=mean(c(min,max)), sd=22400))# Should 
  #   # give 95% density
  #   # plot(x=seq(min,max,100), y=dnorm(seq(min,max,100), 
  #   # mean=mean(c(min,max)), sd=22400),type="l")
  #   data$B_penalty_mu <- mean(c(min,max))/ScaleSMU
  #   if (!remove.EnhStocks) data$B_penalty_sig <- 22400/ScaleSMU
  #   # Should give 95% density:
  #   # sum(dnorm(seq(min,max,1), mean=mean(c(min,max)), sd=2700))
  #   if (remove.EnhStocks) data$B_penalty_sig <- 2700/ScaleSMU
  #   
  #   param <- list()
  #   param$B_0 <- -2
  #   param$B_1 <- 0.1
  #   
  #   
  #   # dyn.unload(dynlib(paste("TMB_Files/Logistic_LRPs", sep="")))
  #   # compile(paste("TMB_Files/Logistic_LRPs.cpp", sep=""))
  #   
  #   dyn.load(dynlib(paste("TMB_Files/Logistic_LRPs", sep=""))) 
  #   
  #   obj <- MakeADFun(data, param, DLL="Logistic_LRPs", silent=TRUE) 
  #   
  #   opt <- nlminb(obj$par, obj$fn, obj$gr, control = 
  #                   list(eval.max = 1e5, iter.max = 1e5)) 
  #   pl <- obj$env$parList(opt$par) 
  #   #summary(sdreport(obj), p.value=TRUE)
  #   
  #   # Get parameter estimates and logit predicted values for CIs
  #   All_Ests <- data.frame(summary(sdreport(obj), p.value=TRUE))
  #   All_Ests$Param <- row.names(All_Ests)
  #   All_Ests$Param <- sapply(All_Ests$Param, function(x) 
  #     (unlist(strsplit(x, "[.]"))[[1]]))
  #   Preds <- All_Ests %>% filter(Param == "Logit_Preds")
  #   All_Ests <- All_Ests %>% filter(!(Param %in% c( "Logit_Preds"))) 
  #   
  #   
  #   
  #   
  #   out <- list()
  #   out$All_Ests <- All_Ests
  #   
  #   
  #   Logistic_Data <- data.frame(yr = SMUlogisticData$Years, 
  #                               yy = SMUlogisticData$ppn, 
  #                               xx = SMUlogisticData$SMU_Esc)
  #   
  #   out$Logistic_Data <- Logistic_Data
  #   
  #   Logistic_Fits <- data.frame(xx = data$Pred_Abund*ScaleSMU, 
  #                               fit = inv_logit(Preds$Estimate),
  #                               lwr = inv_logit(Preds$Estimate - 
  #                                                 1.96*Preds$Std..Error),
  #                               upr = inv_logit(Preds$Estimate + 
  #                                                 1.96*Preds$Std..Error))
  #   
  #   out$Preds <- Logistic_Fits
  #   
  #   out$LRP <- data.frame(fit = (All_Ests %>% filter(Param == "Agg_LRP") %>% 
  #                                  pull(Estimate))*ScaleSMU, 
  #                         lwr = (All_Ests %>% filter(Param == "Agg_LRP") %>% 
  #                                  mutate(xx =Estimate - 1.96*Std..Error) %>% 
  #                                  pull(xx) ) * ScaleSMU,
  #                         upr = (All_Ests %>% filter(Param == "Agg_LRP") %>% 
  #                                  mutate(xx =Estimate + 1.96*Std..Error) %>% 
  #                                  pull(xx) ) * ScaleSMU)
  #   
  #   return(list(out=out, # C
  #               WCVIEsc=WCVIEsc, # C
  #               SMU_Esc=SMU_Esc, # C
  #               CU_Status=CU_Status, # C
  #               SMU_ppn=SMU_ppn, # C
  #               LRPppn= data$p, # C
  #               nLL=obj$report()$ans, # within
  #               LOO=LOO, # Outside
  #               bench=select(SGENcalcs, -apar, -bpar)*Scale))
  # }
  
  # FUNCTION STALLING HERE
  # return(out=out, WCVIEsc=WCVIEsc, SMU_Esc=SMU_Esc,
  #             CU_Status=CU_Status, SMU_ppn=SMU_ppn,
  #             LRPppn= data$p, nLL=obj$report()$ans, LOO=LOO,
  #             bench= select(SGENcalcs,-apar, -bpar)*Scale)

} 

# Return statements are not working
# When removed function runs completely. 

# **********************************************************************************************************************
# End ----
# **********************************************************************************************************************

