#-------------------------------------------------------------------------
# Integrated Watershed Area Model
# Steps:
# 1. Read in stock-recruitment data, life-history type, and watershed areas
#   for synoptic survey data, and watershed area and life-history type for 
#   additional stocks
# 2. Create data and parameter lists for TMB
# 3. Run TMB model to estimate Ricker parameters and SMSY & SREP for synoptic 
#   data sets, estimate paraemeters of watershed-area regression, and 
#   estimate SMSY and SREP for additional stocks 
# 4. Compile model outputs
# 5. Calculate diagnostics for SR models in synoptic data set and plot SR 
#   curves, etc.
# 6. Calculate prediction intervals for SMSY and SREP estimates for additional 
#   "test" stocks. These are written to a *.csv file
#---------------------------------------------------------

#-------------------------------------------------------------------------------
# Libaries
#-------------------------------------------------------------------------------

library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
library(zoo)
library(viridis)
library(hrbrthemes)

source ("R/helperFunctions.R")
source ("R/PlotFunctions.r")

#-------------------------------------------------------------------------------
# Function to run integrated watershed-area model

# Arguments:
# remove.EnhStocks <- TRUE # A logical representing if enhanced stocks should 
  # be removed for WCVI CK case study
# removeSkagit <- FALSE # A logical representing if Skagit should be removed
 mod <- "Liermann_PriorRicSig_PriorDeltaSig" # A character for the TMB model 
  # to be used. Optoins include:
  #"Liermann_HalfNormRicVar_FixedDelta"
  #"Ricker_AllMod"#"Liermann"#""Ricker_AllMod"
  #"IWAM_FixedSep_RicStd"##"IWAM_FixedSep_Constm"
  #"IWAM_FixedSep_Constyi"#"IWAM_FixedSep_RicStd"
  #"IWAM_FixedSep"#"IWAM_FixedCombined"
# plot <- FALSE # A logical representing if plots should be produced

# Returns:
# csv file of SMSY and SREP values from watershed-are model
  # DataOut/WCVI_SMSY_noEnh,csv or 
  # DataOutWCVI_SMSY_wEnh.csv
# plots of SR diagnostics and Watershed-area regression

runIWAM <- function(remove.EnhStocks = TRUE, removeSkagit = FALSE, 
                    mod = "Liermann_PriorRicSig_PriorDeltaSig", plot = FALSE){
  
  if( plot== TRUE) {
    # For SR plotting purposes below, need to estimate std Ricker SMSY for AR1 
    # stocks, "SMSY_std":
    source ("R/CheckAR1.r")
  }
  
  
  #-----------------------------------------------------------------------------
  # 1. Read in data
  #-----------------------------------------------------------------------------
  
  SRDatwNA <- read.csv("DataIn/SRinputfile.csv")
  SRDatwNA <- SRDatwNA %>% filter(Name != "Hoko" & Name != "Hoh") 
  #remove two stocks not used in Parken et al, and not documented in Liermann 
  if (removeSkagit==TRUE) {
    SRDatwNA <- SRDatwNA %>% filter(Name != "Skagit")
    # Skagit is Stocknumber=22. Need to re-align stock numbers of last 2 
    # stocks: 23 & 24
    SRDatwNA [which(SRDatwNA$Stocknumber==23),2] = 22
    SRDatwNA [which(SRDatwNA$Stocknumber==24),2] = 23
  }
  
  # Which stocks have NAs?
  stockwNA <- SRDatwNA %>% filter (is.na(Rec) == TRUE) %>% 
    dplyr::select (Stocknumber) %>%  unique() %>% unlist() 
  #Do not use AR(1) model on  stocks with NAs, Humptulips and Queets (20 & 21)
  
  # Remove years with NAs
  SRDat <- SRDatwNA %>% filter(Rec != "NA") 
  
  # Revise yr_num list where NAs have been removed
  test <- SRDat %>% filter(Stocknumber == stockwNA[1]|Stocknumber == stockwNA[2])
  if( max(SRDat$Stocknumber) >= stockwNA[1]) {
    for (i in 1:length(stockwNA)) {
      len <- length (SRDat [which (SRDat$Stocknumber == stockwNA[i]), ]$yr_num) - 1
      SRDat [which (SRDat$Stocknumber == stockwNA[i]), ]$yr_num <- c (0:len)
    }
  }
  test <- SRDat %>% filter(Stocknumber == stockwNA[1]| Stocknumber == stockwNA[2])
    # Removed years - wanted to have the time series consistent - re-indexing the stocks - so that there are 
    # no gaps in the. E.g. 0, 1, 2, 3, remove 2 - 0, 1, 3 (now has a break-point)
    # mask?
    # simulate NAN's in future - COSEWIC example
  
  
  # Calculate scale for each stock
  digits <- SRDat %>% group_by(Stocknumber) %>% 
    summarize(maxDigits = count.dig(max(Sp)))
    # requires source loads per run - stored in environment **TOR**
  SRDat <- left_join(SRDat, digits)
  SRDat <- SRDat %>% mutate(Scale = 10^(maxDigits-1))
  
  
  stks_ar <- c("Chikamin", "Keta", "Blossom", "Situk", "Siletz", "Columbia Sp")
  # Cowichan, stk-23, not included here becuase modelled as per Tompkins with a 
  # surival covariate
  stksNum_ar <- c(4,5,6,10,11,16)
  
  # Ricker with a surival co-variate. 
  stksNum_surv <- c(0,23)
  stks_surv <- c("Harrison", "Cowichan")
  if (removeSkagit==TRUE) {stksNum_surv <- c(0,22)}
  
  len_stk <- length(unique(SRDat$Stocknumber))
  stksNum_std <- which(0:(len_stk-1) %not in%c(stksNum_ar, stksNum_surv)==TRUE)-1 
  # Assuming there are only 25 stocks (0:24 StockNumber)
  
  # When aggregated standard, ar1, surv, "ModelOrder" is the order of stocks for 
  # aligning with WA and life-history data
  stksOrder <- data.frame(Stocknumber =  c(stksNum_std, stksNum_ar, stksNum_surv), 
                          ModelOrder = 0:(len_stk-1))
  
  # What is the scale of S,R and SMSY,SREP data,
  # ordered when aggregated by std, AR1, surv?
  SRDat_Scale <- SRDat %>% dplyr::select(Stocknumber, Scale) %>% distinct() 
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedCombined"|mod=="IWAM_FixedSep_Constyi"|
     mod=="IWAM_FixedSep_Constm") {
    SRDat_Scale <- SRDat_Scale %>% left_join(stksOrder) %>% arrange(ModelOrder)
  }
  if(mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|
     mod=="Liermann_PriorRicSig_PriorDeltaSig"|
     mod=="Liermann_HalfNormRicVar_FixedDelta") {
    SRDat_Scale <- SRDat_Scale %>% left_join(stksOrder) %>% arrange(Stocknumber)
  }
  SRDat_Scale <- SRDat_Scale$Scale
  
  # When 3 model forms are used, separate stocks accordingly
  SRDat_std <- SRDat %>% filter(Stocknumber %not in% c(stksNum_ar,stksNum_surv)) 
  SRDat_ar <- SRDat %>% filter(Stocknumber %in% stksNum_ar) 
  SRDat_surv <- SRDat %>% filter(Stocknumber %in% stksNum_surv) 
  # Othewise, use standard Ricker model
  if(mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|
     mod=="Liermann_PriorRicSig_PriorDeltaSig"|
     mod=="Liermann_HalfNormRicVar_FixedDelta") SRDat_std <- SRDat
  
  
  # Assign new stock numbers to each stock so that they are sequential within each 
  # model form. These are used in TMB When only one model form is use (std Ricker)
  # then ind_std = Stocknumber
  ind_std <- tibble(ind_std= 0:(length(unique(SRDat_std$Name))-1))
  ind_std <- add_column(ind_std, Stocknumber = (unique(SRDat_std$Stocknumber)))
  SRDat_std <- SRDat_std %>% left_join(ind_std)
  
  ind_ar <- tibble(ind_ar= 0:(length(unique(SRDat_ar$Name))-1))
  ind_ar <- add_column(ind_ar, Stocknumber = (unique(SRDat_ar$Stocknumber)))
  SRDat_ar <- SRDat_ar %>% left_join(ind_ar)
  
  ind_surv <- tibble(ind_surv= 0:(length(unique(SRDat_surv$Name))-1))
  ind_surv <- add_column(ind_surv, Name = (unique(SRDat_surv$Name)))
  Surv <- as.data.frame(read.csv("DataIn/Surv.csv")) 
  SRDat_surv <- SRDat_surv %>% left_join(ind_surv, by="Name") %>% left_join(Surv)
  
  
  # Remove years 1981-1984, 1986-1987  from Cowichan (Stocknumber 23) as per 
  # Tompkins et al. 2005
  SRDat_surv_Cow <- SRDat_surv %>% filter(Name == "Cowichan" & 
                                            Yr >= 1985 & Yr !=1986 &Yr != 1987) 
  n_surv_Cow <- length(SRDat_surv_Cow$Yr)
  SRDat_surv_Cow$yr_num <- 0:(n_surv_Cow-1)
  if("Cowichan" %in% stks_surv) SRDat_surv <- SRDat_surv %>% filter(Name != "Cowichan") %>% 
    bind_rows(SRDat_surv_Cow)
  if(mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=="Liermann_PriorRicSig_PriorDeltaSig"|
     mod=="Liermann_HalfNormRicVar_FixedDelta") {
    SRDat_Cow <- SRDat_std %>% filter(Name == "Cowichan" & Yr >= 1985 & Yr !=1986 & Yr != 1987) 
    n_Cow <- length(SRDat_Cow$Yr)
    SRDat_Cow$yr_num <- 0:(n_surv_Cow-1)
    SRDat_std <- SRDat_std %>%  filter(Name != "Cowichan") %>% bind_rows(SRDat_Cow) %>%
      arrange(ind_std) 
  }
  
  # Read in watershed area data and life-history type (stream vs ocean)
  WA <- read.csv("DataIn/WatershedArea.csv")
  names <- SRDat %>% dplyr::select (Stocknumber, Name) %>% distinct() 
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedCombined"|mod=="IWAM_FixedSep_Constm"|
     mod=="IWAM_FixedSep_Constyi") {
    WA <- WA %>% full_join(names, by="Name") %>% full_join (stksOrder, by="Stocknumber") %>% 
      arrange(ModelOrder)
  }
  if(mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=="Liermann_PriorRicSig_PriorDeltaSig"|
     mod=="Liermann_HalfNormRicVar_FixedDelta") {
    WA <- WA %>% full_join(names, by="Name") %>% full_join (stksOrder, by="Stocknumber") %>% 
      arrange(Stocknumber)
  }
  
  if (removeSkagit==TRUE) {WA <- WA %>% filter(Name !="Skagit")}
  Stream <- SRDat %>% dplyr::select(Stocknumber, Name, Stream) %>% group_by(Stocknumber) %>% 
    summarize(lh=max(Stream))
  Stream <- Stream %>% full_join(stksOrder, by="Stocknumber") %>% arrange(ModelOrder)
  if(mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=="Liermann_PriorRicSig_PriorDeltaSig"|
     mod=="Liermann_HalfNormRicVar_FixedDelta") Stream <- Stream  %>% arrange(Stocknumber)
  
  
  # ---------------------------------------------------------------------------------------
  # 2. Create data and parameter lists for TMB
  # ---------------------------------------------------------------------------------------
  
  # Only rho_Start used. Disregard other TMB_Input
  TMB_Inputs <- list(rho_Start = 0.0, logDelta1_start=3.00, logDelta2_start =log(0.72), 
                     logDeltaSigma_start = -0.412, logMuDelta1_mean= 5, logMuDelta1_sig= 10, 
                     logMuDelta2_mean=-0.5, logMuDelta2_sig= 10, Tau_Delta1_dist= 0.1, 
                     Tau_Delta2_dist= 0.1, Tau_sigma = 0.01) 
  
  
  # Data 
  data <- list()
  Scale_std <- SRDat_std$Scale 
  data$S_std <- SRDat_std$Sp/Scale_std 
  data$logRS_std <- log( (SRDat_std$Rec/Scale_std) / (SRDat_std$Sp/Scale_std) )
  data$stk_std <- as.numeric(SRDat_std$ind_std)
  N_Stocks_std <- length(unique(SRDat_std$Name))
  data$yr_std <- SRDat_std$yr_num
  
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedCombined"|mod=="IWAM_FixedSep_Constm"|
     mod=="IWAM_FixedSep_Constyi"|mod=="Ricker_AllMod"){
    Scale_ar <- SRDat_ar$Scale 
    data$S_ar <- SRDat_ar$Sp/Scale_ar 
    data$logRS_ar <- log( (SRDat_ar$Rec/Scale_ar) / (SRDat_ar$Sp/Scale_ar) ) 
    data$stk_ar <- as.numeric(SRDat_ar$ind_ar)
    N_Stocks_ar <- length(unique(SRDat_ar$Name))
    data$yr_ar <- SRDat_ar$yr_num
    
    Scale_surv <- SRDat_surv$Scale 
    data$S_surv <- SRDat_surv$Sp/Scale_surv
    data$logRS_surv <- log( (SRDat_surv$Rec/Scale_surv) / (SRDat_surv$Sp/Scale_surv) )
    data$stk_surv <- as.numeric(SRDat_surv$ind_surv)
    N_Stocks_surv <- length(unique(SRDat_surv$Name))
    data$yr_surv <- SRDat_surv$yr_num
    data$Surv_surv <- log(SRDat_surv$Surv) #Tompkins et al. used Ln(Surv+1)
    meanLogSurv <- SRDat_surv %>% group_by(Stocknumber) %>% 
      summarize(meanLogSurv = mean(log(Surv))) %>% dplyr::select(meanLogSurv) 
    data$MeanLogSurv_surv <- meanLogSurv$meanLogSurv
    
  }
  
  if(mod=="IWAM_FixedSep_RicStd"){
    data$biasCor <- as.numeric(TRUE)
  }
  if (mod=="Liermann"){
    data$Tau_dist <- 0.01#TMB_Inputs$Tau_sigma
    data$logMuA_mean <- 1.5
    data$logMuA_sig <- 2
    data$Tau_A_dist <- 0.01#TMB_Inputs$Tau_sigma
  }
  
  if (mod=="Liermann_PriorRicSig_PriorDeltaSig"){
    data$logMuAs_mean <- 1.5
    data$logMuAs_sig <- 2
    data$logMuAo_mean <- 0#1.5
    data$logMuAo_sig <- 2
    data$HalfNormMean <- 0#TMB_Inputs$Tau_sigma
    data$HalfNormSig <- 1#TMB_Inputs$Tau_sigma
    data$HalfNormMeanA <- 0#0.44#TMB_Inputs$Tau_sigma
    data$HalfNormSigA <- 1#0.5#TMB_Inputs$Tau_sigma
    data$SigRicPriorNorm <- as.numeric(F)
    data$SigRicPriorGamma <- as.numeric(T)
    data$SigRicPriorCauchy <- as.numeric(F)
    data$biasCor <- as.numeric(TRUE)
    data$Tau_dist <- 0.1
    
    data$sigDelta_mean <- 0.80# See KFrun.R, #For half-normal use N(0,1)
    data$sigDelta_sig <- 0.28# See KFrun.R,
    data$sigNu_mean <- 0.84# See KFrun.R,
    data$sigNu_sig <- 0.275# See KFrun.R,
    data$SigDeltaPriorNorm <- as.numeric(F)
    data$SigDeltaPriorGamma <- as.numeric(T)
    data$SigDeltaPriorCauchy <- as.numeric(F)
    data$Tau_D_dist <- 1
    data$TestlnWAo <- read.csv("DataIn/WCVIStocks.csv") %>% mutate (lnWA=log(WA)) %>%
      filter(lh==1) %>% pull(lnWA)
    # Add aggregated WAs at inlet level
    InletlnWA <- data.frame(read.csv("DataIn/WCVIStocks.csv")) %>% 
      filter(Stock != "Cypre") %>% group_by(Inlet) %>%
      summarize(InletlnWA = log(sum(WA))) %>% filter(Inlet != "San Juan") %>%
      filter(Inlet !="Nitinat")
    InletlnWAnoEnh <- data.frame(read.csv("DataIn/WCVIStocks.csv")) %>% 
      filter(Stock != "Cypre") %>% filter(Enh==0) %>%
      group_by(Inlet) %>% summarize(InletlnWA = log(sum(WA))) %>% 
      filter(Inlet != "San Juan") %>%
      filter(Inlet !="Nitinat")
    CUlnWA <- data.frame(read.csv("DataIn/WCVIStocks.csv")) %>% 
      filter(Stock != "Cypre") %>% group_by(CU) %>%
      summarize(CUlnWA = log(sum(WA)))
    CUlnWAnoEnh <- data.frame(read.csv("DataIn/WCVIStocks.csv")) %>% 
      filter(Stock != "Cypre") %>% filter(Enh==0) %>%
      group_by(CU) %>% summarize(CUlnWA = log(sum(WA)))
    if(remove.EnhStocks) data$TestlnWAo <- c(data$TestlnWAo, InletlnWAnoEnh$InletlnWA,
                                             CUlnWAnoEnh$CUlnWA)
    if(!remove.EnhStocks) data$TestlnWAo <- c(data$TestlnWAo, InletlnWA$InletlnWA,
                                              CUlnWA$CUlnWA )

    ## Code for using Parken et al. test stocks:
    # data$TestlnWAs <- read.csv("DataIn/ParkenTestStocks.csv") %>% mutate (lnWA=log(WA)) %>%
    # filter(lh==0) %>% pull(lnWA)
    # data$TestlnWAo <- read.csv("DataIn/ParkenTestStocks.csv") %>% mutate (lnWA=log(WA)) %>%
    # filter(lh==1) %>% pull(lnWA)
    
    
  }
  
  if (mod=="Liermann_HalfNormRicVar_FixedDelta"){
    data$logMuAs_mean <- 1.5
    data$logMuAs_sig <- 2
    data$logMuAo_mean <- 0#1.5
    data$logMuAo_sig <- 2
    data$HalfNormMean <- 0#TMB_Inputs$Tau_sigma
    data$HalfNormSig <- 1#TMB_Inputs$Tau_sigma
    data$HalfNormMeanA <- 0#0.44#TMB_Inputs$Tau_sigma
    data$HalfNormSigA <- 1#0.5#TMB_Inputs$Tau_sigma
    
    data$logDeltaSigma <- log(0.21)# See KFrun.R, "medSDlogSmsy" = log(0.21)
    ## Or take Parken et al. (2006) WA sigmas, averaging var of straeam and ocean type models, 
    ## and taking sqrt = log(sqrt((0.293+0.146)/2)=0.47)
    data$logNuSigma <- log(0.29)# See KFrun.R, "medSDlogSrep" = log(0.29)
    ## Or take Parken et al. (2006) WA sigmas, averaging var of straeam and ocean type models, 
    ## and taking sqrt = log(sqrt((0.240+0.133)/2)=0.43)
    # data$TestlnWAs <- read.csv("DataIn/ParkenTestStocks.csv") %>% mutate (lnWA=log(WA)) %>% 
    # filter(lh==0) %>% pull(lnWA)
    # data$TestlnWAo <- read.csv("DataIn/ParkenTestStocks.csv") %>% mutate (lnWA=log(WA)) %>% 
    # filter(lh==1) %>% pull(lnWA)
    
    data$TestlnWAo <- read.csv("DataIn/WCVIStocks.csv") %>% mutate (lnWA=log(WA)) %>% 
      filter(lh==1) %>% pull(lnWA)
    
  }
  
  
  # Read in wateshed area data and life-history type and scale
  if (mod!="Ricker_AllMod") data$WA <- WA$WA
  if (mod!="Ricker_AllMod") data$Stream <- Stream$lh
  data$Scale <- SRDat_Scale #ordered by std, AR1, surv, if all 3 Ricker models uses. 
  # Otherwise ordered by Stocknumber
  if (mod=="IWAM_FixedSep") data$order_noChick <- c(0:23)##c(0:15,19)#, 17:23)
  
  
  # Read in log(watershed area) for additional stocks
  # Predicated lnWA for plottig CIs:
  if (mod!="Ricker_AllMod") data$PredlnWA <- seq(min(log(WA$WA)), max(log(WA$WA)), 0.1) 
  
  # Parameters
  param <- list()
  Scale.stock_std <- (SRDat %>% group_by(Stocknumber) %>% 
                        filter(Stocknumber %not in% c(stksNum_ar,stksNum_surv)) %>% 
                        summarize(Scale.stock_std = max(Scale)))$Scale.stock_std
  Scale.stock_ar <- (SRDat %>% group_by(Stocknumber) %>% 
                       filter(Stocknumber %in% stksNum_ar) %>% 
                       summarize(Scale.stock_ar = max(Scale)))$Scale.stock_ar
  Scale.stock_surv <- (SRDat %>% group_by(Stocknumber) %>% 
                         filter(Stocknumber %in% stksNum_surv) %>% 
                         summarize(Scale.stock_surv = max(Scale)))$Scale.stock_surv
  
  if(mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=="Liermann_PriorRicSig_PriorDeltaSig"|
     mod=="Liermann_HalfNormRicVar_FixedDelta"){
    Scale.stock_std <- (SRDat %>% group_by(Stocknumber) %>% 
                          summarize(Scale.stock_std = max(Scale)))$Scale.stock_std
  }
  
  # Parameters for stocks without AR1
  param$logA_std <- ( SRDat_std %>% group_by (Stocknumber) %>% 
                        summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) )$yi
  B_std <- SRDat_std %>% group_by(Stocknumber) %>% 
    summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
  param$logB_std <- log ( 1/ ( (1/B_std$m)/Scale.stock_std ))#log(B_std$m/Scale.stock)
  param$logSigma_std <- rep(-2, N_Stocks_std)
  
  
  if(mod=="Liermann"){
    param$logMuA <- 1.5
    param$logSigmaA <- 5#1.5#1
  }
  
  if(mod=="Liermann_PriorRicSig_PriorDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta"){
    param$logMuAs <- 1.5
    #param$logSigmaAs <- 1
    param$logMuAo <- 0#1.5
    param$logSigmaA <- -2#5
    #param$logSigmaAo <- 1
    
  }
  
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedCombined"|mod=="IWAM_FixedSep_Constm"|
     mod=="IWAM_FixedSep_Constyi"|mod=="Ricker_AllMod"){
    # Parameters for stocks with AR1
    param$logA_ar <- ( SRDat_ar %>% group_by(Stocknumber) %>% 
                         summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) )$yi
    B_ar <- SRDat_ar %>% group_by(Stocknumber) %>% summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
    # Take inverse of B (=Smax and apply scale), the take the inverse again and log to 
    # get logB of scaled Smax
    param$logB_ar <- log ( 1/ ( (1/B_ar$m)/Scale.stock_ar ))
    param$rho <- rep(TMB_Inputs$rho_Start, N_Stocks_ar)
    param$logSigma_ar <- rep (-2, N_Stocks_ar)
    
    # Parameters for stock with survival covariate
    param$logA_surv <- ( SRDat_surv %>% group_by(Stocknumber) %>% 
                           summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) )$yi
    B_surv <- SRDat_surv %>% group_by(Stocknumber) %>% 
      summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
    # Take inverse of B (=Smax and apply scale), the take the inverse again and log to 
    # get logB of scaled Smax
    param$logB_surv <- log ( 1/ ( (1/B_surv$m)/Scale.stock_surv ))
    param$logSigma_surv <- rep (-2, N_Stocks_surv)
    param$gamma <- rep (0, N_Stocks_surv)
    
    #param$logSgen <- log((SRDat %>% group_by(CU_Name) %>%  
    # summarise(x=quantile(Spawners, 0.5)))$x/Scale) 
    
  }
  
  
  if (mod=="IWAM_FixedCombined"){
    param$logDelta1 <- 3.00# with skagit 2.881
    param$logDelta2 <- log(0.72)#log(0.72/(1-0.72)) 
    #param$Delta2Ocean <- 0
    #param$Delta2 <- log(0.72/(1-0.72)) 
    param$logDeltaSigma <- -0.412 #from Parken et al. 2006 where sig=0.662
    # without Skagit lnDelta1_start <- 2.999911
    # without Skagit lnDelta2_start <- -0.3238648, or Delta2 = 0.723348
  }
  
  ## Separate Stream and Ocean type models
  #param$slogDelta1 <- 2.744 #best estimates from run of stream-specific WAreg TMB
  #param$sDelta2 <- 0.857 
  #param$slogDeltaSigma <- -0.709 
  
  #param$ologDelta1 <- 3.00#1.519 #best estimates from run of ocean-specific WAreg TMB 
  #param$ologDelta2 <- log(0.94)#0#21.2 
  #param$ologDeltaSigma <-  -0.412#-0.94 
  
  if (mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|
      mod=="Liermann_PriorRicSig_PriorDeltaSig"){
    ## Lierman model
    param$logDelta1 <- 3#10# with skagit 2.881
    param$logDelta1ocean <- 0# with skagit 2.881
    param$logDelta2 <- log(0.72)#log(0.72/(1-0.72)) 
    param$Delta2ocean <- 0#log(0.72/(1-0.72)) 
    param$logDeltaSigma <- -0.412 #from Parken et al. 2006 where sig=0.662
  }
 if (mod=="Liermann"| mod=="Liermann_PriorRicSig_PriorDeltaSig"){ 
    param$logNu1 <- 3#10# with skagit 2.881
    param$logNu1ocean <- 0# with skagit 2.881
    param$logNu2 <- log(0.72)#log(0.72/(1-0.72)) 
    param$Nu2ocean <- 0#log(0.72/(1-0.72)) 
    param$logNuSigma <- -0.412 #from Parken et al. 2006 where sig=0.66
  }
  
  if (mod=="Liermann_HalfNormRicVar_FixedDelta"){
    ## Lierman model
    param$logDelta1 <- 3#10# with skagit 2.881
    param$logDelta1ocean <- 0# with skagit 2.881
    param$logDelta2 <- log(0.72)#log(0.72/(1-0.72)) 
    param$Delta2ocean <- 0#log(0.72/(1-0.72)) 
    #param$logDeltaSigma <- -0.412 #from Parken et al. 2006 where sig=0.662
    
    param$logNu1 <- 3#10# with skagit 2.881
    param$logNu1ocean <- 0# with skagit 2.881
    param$logNu2 <- log(0.72)#log(0.72/(1-0.72)) 
    param$Nu2ocean <- 0#log(0.72/(1-0.72)) 
    #param$logNuSigma <- -0.412 #from Parken et al. 2006 where sig=0.662
    
    
  }
  
  if (mod=="IWAM_FixedSep_Constm"){
    ## Lierman model
    param$logDelta1 <- 3#10# with skagit 2.881
    param$logDelta1ocean <- 0# with skagit 2.881
    param$logDelta2 <- log(0.72)#log(0.72/(1-0.72)) 
    param$logDeltaSigma <- -0.412 #from Parken et al. 2006 where sig=0.662
    
  }
  if (mod=="IWAM_FixedSep_Constyi"){
    ## Lierman model
    param$logDelta1 <- 3#10# with skagit 2.881
    param$logDelta2 <- log(0.72)#log(0.72/(1-0.72)) 
    param$Delta2ocean <- 0#log(0.72/(1-0.72)) 
    param$logDeltaSigma <- -0.412 #from Parken et al. 2006 where sig=0.662
    
  }
  ## Hierarchcial model
  
  #param$logDelta1 <- TMB_Inputs$logDelta1_start
  #param$logDelta2 <- rep(TMB_Inputs$logDelta2_start, 2)
  #param$logDeltaSigma <-TMB_Inputs$logDeltaSigma_start 
  #param$logMuDelta1 <- TMB_Inputs$logDelta1_start
  #param$SigmaDelta1 <- 10
  #param$logMuDelta2 <- TMB_Inputs$logDelta2_start
  #param$SigmaDelta2 <- 1
  
  
  # -----------------------------------------------------------------------
  # 3. Estimate SR parameters from synoptic data set and SMSY and SREPs
  # -----------------------------------------------------------------------
  
  # Compile model if changed:
  #dyn.unload(dynlib(paste("TMB_Files/", mod, sep="")))
  #compile(paste("TMB_Files/", mod, ".cpp", sep=""))
  dyn.load(dynlib(paste("TMB_Files/", mod, sep="")))
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedCombined"|mod=="IWAM_FixedSep_Constyi"|
     mod=="IWAM_FixedSep_Constm"|mod=="IWAM_FixedSep_RicStd"|mod=="Ricker_AllMod"){
    obj <- MakeADFun(data, param, DLL=mod, silent=TRUE)#random = c( "logDelta2"), 
  }
  
  if(mod=="Liermann"|mod=="Liermann_PriorRicSig_PriorDeltaSig"|
     mod=="Liermann_HalfNormRicVar_FixedDelta"){
    obj <- MakeADFun(data, param, DLL=mod, silent=TRUE, random = c("logA_std"))
    
  }
  #
  
  # For phasing, (not needed)
  # map <- list(logDelta1=factor(NA), Delta2=factor(NA), logDeltaSigma=factor(NA)) 
  # obj <- MakeADFun(data, param, DLL="Ricker_AllMod", silent=TRUE, map=map)
  
  
  upper<-unlist(obj$par)
  upper[1:length(upper)]<- Inf
  
  lower<-unlist(obj$par)
  lower[1:length(lower)]<- -Inf
  
  if(mod=="Liermann_PriorRicSig_PriorDeltaSig"){
    if(data$SigDeltaPriorNorm == 1){
      upper[names(upper) == "logDeltaSigma"] <- log(1.39) # See KFrun.R, "SDlSMSYParken"
      upper[names(upper) == "logNuSigma"] <- log(1.38)# See KFrun.R, "SDlSREPParken"
      lower[names(lower) == "logDeltaSigma"] <- log(0.21) # See KFrun.R, "medSDlogSmsy"
      lower[names(lower) == "logNuSigma"] <- log(0.29) # See KFrun.R, "medSDlogSrep"
    }
  }
  # Don't need this any more - setting limits on the parameters
  
  # ## For uniform priors on Ricker var (set all Ricker Sig priors above to F)
  # upper[names(upper) == "logSigma_std"] <- log(2) 
  # lower[names(lower) == "logSigma_std"] <- log(0) #Inf!
  # upper[names(upper) == "logSigmaA"] <- log(2) 
  # lower[names(lower) == "logSigmaA"] <- log(0) #Inf!
  
  
  #### RUNNING THE MODEL ####
  opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5), 
                lower=lower, upper=upper)
  pl <- obj$env$parList(opt$par) # Gives the parameter estimates from the model
  #summary(sdreport(obj), p.value=TRUE)

  
  #library(tmbstan)
  # fitmcmc <- tmbstan(obj, chains=3, iter=1000, init=list(opt$par), 
  # control = list(adapt_delta = 0.95))
  # fitmcmc <- tmbstan(obj, chains=3, iter=1000, init=list("last.par.best"), 
  # control = list(adapt_delta = 0.95))
  
  # -----------------------------------------------------------------------
  # 4. Compile model outputs
  # -----------------------------------------------------------------------
  
  # Create Table of outputs
  All_Ests <- data.frame(summary(sdreport(obj)))
  All_Ests$Param <- row.names(All_Ests)
  
  
  # Put together readable data frame of values
  All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
  
  # By first spliting out stocks modelled with standard Ricker, RickerAR(1), 
  # and Ricker-survival models
  All_Ests_std <- data.frame()
  All_Ests_std <- All_Ests %>% filter (Param %in% c("logA_std", "logB_std", "logSigma_std",  
                                                    "SMSY_std", "SREP_std"))
  SN_std <- unique(SRDat_std[, c("Stocknumber")])
  All_Ests_std$Stocknumber <- rep(SN_std)
  All_Ests_std <- left_join(All_Ests_std, unique(SRDat_std[, c("Stocknumber", "Name")]))
  
  logDeltaSigma <- All_Ests %>% filter (Param %in% c("logDeltaSigma")) 
  DeltaSigmaUCL <- exp(logDeltaSigma$Estimate + logDeltaSigma$Std..Error*1.96)
  DeltaSigmaLCL <- exp(logDeltaSigma$Estimate - logDeltaSigma$Std..Error*1.96) 
  DeltaSigma <- exp(logDeltaSigma$Estimate)
  
  
  All_Ests_ar <- data.frame()
  All_Ests_surv <- data.frame()
  
  if (mod=="IWAM_FixedCombined"|mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_Constm"|
      mod=="IWAM_FixedSep_Constyi"|mod=="Ricker_AllMod"){
    All_Ests_ar<- All_Ests %>% filter (Param %in% c("logA_ar", "logB_ar", "rho",  
                                                    "logSigma_ar", "SMSY_ar", "SREP_ar" ))
    SN_ar <- unique(SRDat_ar[, c("Stocknumber")])
    All_Ests_ar$Stocknumber <- rep(SN_ar)
    All_Ests_ar <- left_join(All_Ests_ar, unique(SRDat_ar[, c("Stocknumber", "Name")]))
    
    All_Ests_surv<- All_Ests %>% filter (Param %in% c("logA_surv", "logB_surv", "gamma", 
                                                      "logSigma_surv", "SMSY_surv", 
                                                      "SREP_surv" ))#"gamma"
    SN_surv <- unique(SRDat_surv[, c("Stocknumber")])
    All_Ests_surv$Stocknumber <- rep(SN_surv)
    All_Ests_surv <- left_join(All_Ests_surv, unique(SRDat_surv[, c("Stocknumber", 
                                                                    "Name")]))
    
  }
  
  # Combine again
  All_Est <- bind_rows(All_Ests_std, All_Ests_ar, All_Ests_surv) 
  All_Est$ar <- All_Est$Stocknumber %in% stksNum_ar
  All_Est$surv <- All_Est$Stocknumber %in% stksNum_surv
  All_Est$Param <- sapply(All_Est$Param, function(x) (unlist(strsplit(x, "[_]"))[[1]]))
  All_Est <- All_Est %>%left_join(Stream, by="Stocknumber")
  
  All_Deltas <- data.frame()
  All_Deltas <- All_Ests %>% filter (Param %in% c("logDelta1", "logDelta2","sigma_delta", 
                                                  "Delta2_bounded", "logDelta1ocean", 
                                                  "logDelta2ocean", "Delta2ocean", "logNu1", 
                                                  "logNu2", "sigma_nu", "logNu1ocean", 
                                                  "Nu2ocean"))
  
  # -----------------------------------------------------------------------
  # 5. Calculate diagnostics and plot SR curves, etc.
  # -----------------------------------------------------------------------
  
  # Calculate AIC
  
  nLL_std <- data.frame(nLL_std=obj$report()$nLL_std) %>% 
    add_column(Stocknumber=SRDat_std$Stocknumber) %>% group_by(Stocknumber) %>% 
    summarize(CnLL=sum(nLL_std))
  aic_std <- nLL_std %>% mutate(aic = 2 * 3 + 2*CnLL) 
  if (mod=="IWAM_FixedCombined"|mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_Constm"|
      mod=="IWAM_FixedSep_Constyi"|mod=="Ricker_AllMod"){
    nLL_ar <- data.frame(nLL_ar=obj$report()$nLL_ar) %>% 
      add_column(Stocknumber=SRDat_ar$Stocknumber) %>% group_by(Stocknumber) %>% 
      summarize(CnLL=sum(nLL_ar))
    aic_ar <- nLL_ar %>% mutate(aic = 2 * 4 + 2*CnLL)
    nLL_surv <- data.frame(nLL_surv=obj$report()$nLL_surv) %>% 
      add_column(Stocknumber=SRDat_surv$Stocknumber) %>% group_by(Stocknumber) %>% 
      summarize(CnLL=sum(nLL_surv))
    aic_surv <- nLL_surv %>% mutate(aic = 2 * 4 + 2*CnLL)
    
  }
  
  # Get predicted values and calculate r2
  Pred_std <- data.frame()
  Pred_std <- All_Ests %>% filter (Param %in% c("LogRS_Pred_std"))
  Preds_std <- SRDat_std %>% dplyr::select("Stocknumber","yr_num", "Sp", "Rec", "Scale", "Name") %>% 
    add_column(Pred=Pred_std$Estimate)
  Preds_std <- Preds_std %>% mutate(ObsLogRS = log ( (Rec / Scale) / (Sp/Scale) ) )
  r2 <- Preds_std %>% group_by(Stocknumber) %>% summarize(r2=cor(ObsLogRS,Pred)^2)
  
  if (mod=="IWAM_FixedCombined"|mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_Constm"|
      mod=="IWAM_FixedSep_Constyi"|mod=="Ricker_AllMod"){
    Pred_ar <- data.frame()
    Pred_ar <- All_Ests %>% filter (Param %in% c("LogRS_Pred_ar"))
    Preds_ar <- SRDat_ar %>% dplyr::select("Stocknumber","yr_num", "Sp", "Rec", "Scale", "Name") %>% 
      add_column(Pred=Pred_ar$Estimate)
    Preds_ar <- Preds_ar %>% mutate(ObsLogRS = log ( (Rec / Scale) / (Sp / Scale) ) ) 
    r2_ar <- Preds_ar %>% group_by(Stocknumber) %>% summarize(r2=cor(ObsLogRS,Pred)^2)
    
    Pred_surv <- data.frame()
    Pred_surv <- All_Ests %>% filter (Param %in% c("LogRS_Pred_surv"))
    Preds_surv <- SRDat_surv %>% dplyr::select("Stocknumber","yr_num", "Sp", "Rec", "Scale", "Name") %>% 
      add_column(Pred=Pred_surv$Estimate)
    Preds_surv <- Preds_surv %>% mutate(ObsLogRS = log ( (Rec / Scale) / (Sp / Scale) ) ) 
    r2_surv <- Preds_surv %>% group_by(Stocknumber) %>% summarize(r2=cor(ObsLogRS,Pred)^2)
    r2 <- bind_rows(r2, r2_ar, r2_surv) %>% arrange(Stocknumber)
    
  }
  
  
  
  # Get predicted values and their SEs to plot CIs
  PredlnSMSY <- data.frame() 
  PredlnSMSY <- All_Ests %>% filter (Param %in% c("PredlnSMSY_S", "PredlnSMSY_O", "PredlnSMSY_CI", 
                                                  "PredlnSMSYs_CI", "PredlnSMSYo_CI"))
  PredlnSREP <- data.frame() 
  PredlnSREP <- All_Ests %>% filter (Param %in% c("PredlnSREP_S", "PredlnSREP_O", "PredlnSREP_CI", 
                                                  "PredlnSREPs_CI", "PredlnSREPo_CI"))
  
  
  # Calculate standardized residuals
  if (mod=="IWAM_FixedCombined"|mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_Constm"|
      mod=="IWAM_FixedSep_Constyi"|mod=="Ricker_AllMod"){
    SRes <- bind_rows(Preds_std, Preds_ar, Preds_surv) %>% arrange (Stocknumber)
  }
  if (mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=="Liermann_PriorRicSig_PriorDeltaSig"|
      mod=="Liermann_HalfNormRicVar_FixedDelta"){
    SRes <- Preds_std %>% arrange (Stocknumber)
  }
  SRes <- SRes %>% mutate ( Res = ObsLogRS- Pred) #%>% mutate (StdRes = Res/??)
  sigma <- All_Est %>% filter(Param=="logSigma") %>% dplyr::select(Stocknumber, Estimate, Name)
  SRes <- SRes %>% left_join(sigma) %>% rename(logSig = Estimate)
  SRes <- SRes %>% mutate (StdRes = Res/exp(logSig))
  
  
  #Plot SR curves. linearized model, standardized residuals, autocorrleation plots for synoptic data set
  # if using a Liermann model, use SRDat=SRDat_std; otherwise SRDat=SRDat
  
  if (plot==TRUE){
    png(paste("DataOut/SR_", mod, ".png", sep=""), width=7, height=7, units="in", res=500)
    PlotSRCurve(SRDat=SRDat_std, All_Est=All_Est, SMSY_std=SMSY_std, stksNum_ar=stksNum_ar, 
                stksNum_surv=stksNum_surv, r2=r2, removeSkagit=removeSkagit, mod=mod)
    dev.off()
    png(paste("DataOut/SRLin_", mod, ".png", sep=""), width=7, height=7, units="in", res=1000)
    PlotSRLinear(SRDat=SRDat_std, All_Est=All_Est, SMSY_std=SMSY_std, stksNum_ar=stksNum_ar,
                 stksNum_surv=stksNum_surv, r2=r2, removeSkagit=removeSkagit) 
    dev.off()
    png(paste("DataOut/StdResid_", mod, ".png", sep=""), width=7, height=7, units="in", res=1000)
    PlotStdResid(SRes)
    dev.off()
    png(paste("DataOut/ACF_", mod, ".png", sep=""), width=7, height=7, units="in", res=1000)
    Plotacf(SRes)
    dev.off()
    
  }
  
  #---------------------------------------------------------------------------------
  # Save parameter estimates to plot distribution of logA values (box plots)
  #saveRDS( All_Est, paste( "DataOut/All_Est_", mod, ".RDS", sep="") )
  #saveRDS( All_Est, paste( "DataOut/All_Est_Ricker_std_wBC.RDS", sep="") )# mod "IWAM_FixedSep_RicStd"
  
  #saveRDS( All_Est, paste( "DataOut/All_Est_", mod, "_invGamma0.1_wBC.RDS", sep="") )
  #saveRDS( All_Est, paste( "DataOut/All_Est_", mod, "_invGamma0.01_wBC.RDS", sep="") )
  #saveRDS( All_Est, paste( "DataOut/All_Est_", mod, "_invGamma0.001_wBC.RDS", sep="") )
  #saveRDS( All_Est, paste( "DataOut/All_Est_", mod, "_invGamma0.01_invGammaA0.001.RDS", sep="") )
  #saveRDS( All_Est, paste( "DataOut/All_Est_", mod, "_invGamma0.001_invGammaA0.01.RDS", sep="") )
  #saveRDS( All_Est, paste( "DataOut/All_Est_", mod, "_uniformSigmaPrior_wBC.RDS", sep="") )
  #saveRDS( All_Est, paste( "DataOut/All_Est_", mod, "_uniform1SigmaPrior.RDS", sep="") )
  #saveRDS( All_Est, paste( "DataOut/All_Est_", mod, "_noSigmaAPrior.RDS", sep="") )
  #saveRDS( All_Est, paste( "DataOut/All_Est_", mod, "_noSigmaPrior.RDS", sep="") )
  ## For Liermann_HalfCauchyRicVar and Liermann_HalfNormalRicVar
  ##saveRDS( All_Est, paste( "DataOut/All_Est_", mod, ".RDS", sep="") ) 
  #saveRDS( All_Est, paste( "DataOut/All_Est_", mod, "_HalfNormRicVar_wBC.RDS", sep="") ) 
  #saveRDS( All_Est, paste( "DataOut/All_Est_", mod, "_HalfCauchyRicVar_wBC.RDS", sep="") ) 
  #saveRDS( All_Est, paste( "DataOut/All_Est_", mod, "_noPriorRicVar_wBC.RDS", sep="") ) 
  
  #---------------------------------------------------------------------------------
  # Plot WA regression
  if(plot==TRUE){
    png(paste("DataOut/WAregSMSY_", mod, "_wBC.png", sep=""), width=7, height=7, units="in", res=500)
    par(mfrow=c(1,1), mar=c(4, 4, 4, 2) + 0.1)
    if (mod=="IWAM_FixedCombined") title_plot <- "Fixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
    if (mod=="IWAM_FixedSep") title_plot <- "Separate life-histories\nFixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
    if (mod=="IWAM_FixedSep_RicStd") title_plot <- "Separate life-histories: all Std Ricker\nFixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
    if (mod=="IWAM_FixedSep_Constm") title_plot <- "Separate life-histories\n Fixed-effect yi (logDelta1), \nConstant Fixed-effect slope (Delta2)"
    if (mod=="IWAM_FixedSep_Constyi") title_plot <- "Separate life-histories\n Constant Fixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
    if (mod=="Liermann") title_plot <- "Separate life-histories: Random Ricker a\n Fixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
    if (mod=="Liermann_PriorRicSig_PriorDeltaSig") title_plot <- "Prior Ricker sigma and prior WA regression sigma"
    if (mod=="Liermann_HalfNormRicVar_FixedDelta") title_plot <- "Half-normal Ricker sigma and fixed WA regression sigma"
    #title_plot <- "Separate life-histories: n=17\nFixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
    plotWAregressionSMSY (All_Est, All_Deltas, SRDat, Stream, WA, PredlnSMSY, 
                          PredlnWA = data$PredlnWA, title1=title_plot, mod)
    dev.off()
    
    png(paste("DataOut/WAregSREP_", mod, "_wBC.png", sep=""), width=7, height=7, units="in", res=500)
    #png(paste("DataOut/WAreg_Liermann_SepRicA_UniformSigmaAPrior.png", sep=""), width=7, height=7, units="in", res=500)
    par(mfrow=c(1,1), mar=c(4, 4, 4, 2) + 0.1)
    if (mod=="IWAM_FixedCombined") title_plot <- "Fixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
    if (mod=="IWAM_FixedSep") title_plot <- "Separate life-histories\nFixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
    if (mod=="IWAM_FixedSep_RicStd") title_plot <- "Separate life-histories: all Std Ricker\nFixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
    if (mod=="IWAM_FixedSep_Constm") title_plot <- "Separate life-histories\n Fixed-effect yi (logDelta1), \nConstant Fixed-effect slope (Delta2)"
    if (mod=="IWAM_FixedSep_Constyi") title_plot <- "Separate life-histories\n Constant Fixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
    if (mod=="Liermann") title_plot <- "Separate life-histories: Random Ricker a\n Fixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
    if (mod=="Liermann_PriorRicSig_PriorDeltaSig") title_plot <- "Prior Ricker sigmas and prior on WA regression sigma"
    if (mod=="Liermann_HalfNormRicVar_FixedDelta") title_plot <- "Half-normal Ricker sigma and fixed WA regression sigma"
    #title_plot <- "Separate life-histories: n=17\nFixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
    plotWAregressionSREP (All_Est, All_Deltas, SRDat, Stream, WA, PredlnSREP, 
                          PredlnWA = data$PredlnWA, title1=title_plot, mod)
    dev.off()
    #plotWAregression (All_Est, All_Deltas, SRDat, Stream, WA, PredlnSMSY, PredlnWA = data$PredlnWA, 
    # title1="Common, fixed yi (logDelta1), \nRandom slope (Delta2)")
    
  }
  
  
  
  # -------------------------------------------------------------------------
  # 6. Calculate prediction intervals for SMSY and SREP for additional stocks
  # -------------------------------------------------------------------------
  
  # Get predicted values to estimate prediction intervals
  PredlnSMSY_PI <- data.frame()
  PredlnSMSY_PI <- All_Ests %>% filter (Param %in% c("PredlnSMSY", "lnSMSY"))
  PredlnSREP_PI <- data.frame()
  PredlnSREP_PI <- All_Ests %>% filter (Param %in% c("PredlnSREP", "lnSREP"))
  
  PredlnSMSY_PI$Stocknumber <- rep(SN_std)
  PredlnSREP_PI$Stocknumber <- rep(SN_std)
  
  # To calculate prediction intervals, first get predicted and observed logSMSY 
  # and logSREP values for synoptic data set
  #   (actually only need observed logSMSY and logSREP values)
  
  #  First need to get the scale for each stock
  
  Scale_PI <- SRDat %>% dplyr::select(Stocknumber, Scale) %>% distinct()
  PredlnSMSY_PI <- PredlnSMSY_PI %>% left_join(unique(SRDat_std[, c("Stocknumber", "Name")])) %>% 
    left_join(Scale_PI)
  PredlnSREP_PI <- PredlnSREP_PI %>% left_join(unique(SRDat_std[, c("Stocknumber", "Name")])) %>% 
    left_join(Scale_PI)
  
  # Then need to separate observed stream vs ocean type
  
  # Plsmsys = predicted log SMSY for stream type
  # Plsmsyo = predicted log SMSY for ocean type
  # Olsmsys = observed log SMSY for stream type
  # Olsmsyo = observed log SMSY for ocean type
  
  Plsmsys <- PredlnSMSY_PI %>% filter(Param=="PredlnSMSY") %>% left_join(Stream) %>% 
    filter(lh == 0) %>% pull(Estimate) #Predicted lnSMSY from WA regression- Stream
  Plsmsyo <- PredlnSMSY_PI %>% filter(Param=="PredlnSMSY") %>% left_join(Stream) %>% 
    filter(lh == 1) %>% pull(Estimate) #Predicted lnSMSY from WA regression- ocean
  Olsmsys <- PredlnSMSY_PI %>% filter(Param=="lnSMSY") %>% left_join(Stream) %>% 
    filter( lh== 0) %>% pull(Estimate) # "observed" lnSMSY data output from SR models- stream
  Olsmsyo <- PredlnSMSY_PI %>% filter(Param=="lnSMSY") %>% left_join(Stream) %>% 
    filter(lh == 1) %>% pull(Estimate) # "observed" lnSMSY data output from SR models- ocean
  Plsreps <- PredlnSREP_PI %>% filter(Param=="PredlnSREP") %>% left_join(Stream) %>% 
    filter(lh == 0) %>% pull(Estimate) #Predicted lnSMSY from WA regression- Stream
  Plsrepo <- PredlnSREP_PI %>% filter(Param=="PredlnSREP") %>% left_join(Stream) %>% 
    filter(lh == 1) %>% pull(Estimate) #Predicted lnSMSY from WA regression- ocean
  Olsreps <- PredlnSREP_PI %>% filter(Param=="lnSREP") %>% left_join(Stream) %>% 
    filter(lh == 0) %>% pull(Estimate) # "observed" lnSMSY data output from SR models- stream
  Olsrepo <- PredlnSREP_PI %>% filter(Param=="lnSREP") %>% left_join(Stream) %>% 
    filter(lh == 1) %>% pull(Estimate) # "observed" lnSMSY data output from SR models- ocean
  
  
  # Get watershed areas for synoptic data set to calculate PIs (stream =s and ocean =o)
  WAs <- WA %>% left_join(Stream) %>% filter(lh == 0) %>% pull(WA)
  WAo <- WA %>% left_join(Stream) %>% filter(lh == 1) %>% pull(WA)
  
  # Get names of WCVI stocks
  sn <- read.csv("DataIn/WCVIStocks.csv")
  StockNames <- c(as.vector(sn$Stock), as.vector(InletlnWA$Inlet), as.vector(CUlnWA$CU))
  
  # Get Predicted SMSY and SREP values for new "test" WCVI stocks and their Prediction Intervals
  TestSMSY <- data.frame() 
  TestSMSYs <- data.frame() 
  TestSMSYo <- data.frame() 
  TestSMSY_SREP <- data.frame() 
  
  TestSMSY <- All_Ests %>% filter (Param %in% c("TestlnSMSYo")) %>% add_column(Stock = StockNames)
  TestSREP <- All_Ests %>% filter (Param %in% c("TestlnSREPo")) %>% add_column(Stock = StockNames)
  TestSMSYpull <- TestSMSY %>% pull(Estimate)
  TestSREPpull <- TestSREP %>% pull(Estimate)
  
  # Use custom function: PredInt() to estimate prediction intervals 
  TestSMSY_PI <- PredInt(x=log(WAo), y=Olsmsyo, Predy=TestSMSYpull, Newx= data$TestlnWAo)
  TestSREP_PI <- PredInt(x=log(WAo), y=Olsrepo, Predy=TestSREPpull, Newx= data$TestlnWAo)
  
  TestSMSY <- TestSMSY %>% add_column(LL=exp(TestSMSY_PI$lwr), UL=exp(TestSMSY_PI$upr))
  TestSREP <- TestSREP %>% add_column(LL=exp(TestSREP_PI$lwr), UL=exp(TestSREP_PI$upr))
  TestSMSY <- TestSMSY %>% mutate (Estimate = exp(Estimate)) %>% dplyr::select(-Std..Error, - Param) %>% 
    add_column(Param = "SMSY")
  TestSREP <- TestSREP %>% mutate (Estimate = exp(Estimate)) %>% dplyr::select(-Std..Error, - Param) %>% 
    add_column(Param = "SREP")
  WCVISMSY <- TestSMSY %>% mutate(Estimate=round(Estimate, 0), LL=round(LL,0), UL=round(UL,0))
  WCVISREP <- TestSREP %>% mutate(Estimate=round(Estimate, 0), LL=round(LL,0), UL=round(UL,0))
  WCVISMSY <- WCVISMSY %>% bind_rows(WCVISREP)
  
  # Write SMSY and SREP with PIs to file
  if(remove.EnhStocks) write.csv(WCVISMSY, "DataOut/WCVI_SMSY_noEnh_wBC.csv")
  if(!remove.EnhStocks) write.csv(WCVISMSY, "DataOut/WCVI_SMSY_wEnh_wBC.csv")
  
}# End of runIWAM() function


#---------------------------------------------------------------------------------------------------
# #Code for deriving SMSY and SREP for Parken et al. 2006 test stocks
#  StockNamess <- read.csv("DataIn/ParkenTestStocks.csv") %>% filter(lh == 0) %>% pull(Stock)
#  StockNameso <- read.csv("DataIn/ParkenTestStocks.csv") %>% filter(lh == 1) %>% pull(Stock)
# 
# #To get confidence intervals:
# ## TestSMSY <- TestSMSY %>% mutate (UL = exp(Estimate + 1.96*Std..Error), LL = exp(Estimate - 1.96*Std..Error)) %>%
# ## add_column(Source="IWAM")
#  #Split this out by stream and ocean type as they have different linear regerssions
#  TestSMSYs <- All_Ests %>% filter (Param %in% c("TestlnSMSYs"))  %>% add_column(Stock = StockNamess)
#  TestSMSYo <- All_Ests %>% filter (Param %in% c("TestlnSMSYo"))  %>% add_column(Stock = StockNameso)
#  WAs <- WA %>% left_join(Stream) %>% filter(lh == 0) %>% pull(WA)
#  WAo <- WA %>% left_join(Stream) %>% filter(lh == 1) %>% pull(WA)
#  TestSMSYs_PI <- PredInt(x=log(WAs), y=Olsmsys, Predy=TestSMSYs$Estimate, Newx= c(data$TestlnWAs))
# #Compare my calculations of prediction intervals to those from Stan: quite close, but STAN didn't fully converge..?strange
# #PredInt(x=log(WAs), y=Olsmsys, Predy=Plsmsys, Newx= log(WAs))
# #stantest <- stan("stan_Files/linreg.stan", data = list(x = log(WAs), y = Olsmsys, N = length(WAs)), iter = 8000)
#  TestSMSYo_PI <- PredInt(x=log(WAo), y=Olsmsyo, Predy=TestSMSYo$Estimate, Newx= c(data$TestlnWAo))
#  TestSMSYs <- TestSMSYs %>% add_column(LL=exp(TestSMSYs_PI$lwr), UL=exp(TestSMSYs_PI$upr))
#  TestSMSYo <- TestSMSYo %>% add_column(LL=exp(TestSMSYo_PI$lwr), UL=exp(TestSMSYo_PI$upr))
#  TestSMSY <- TestSMSYs %>% bind_rows(TestSMSYo)
#  TestSMSY <- TestSMSY %>% mutate (SMSY = exp(Estimate)) %>% dplyr::select(-Std..Error, - Param, -Estimate) %>%
#    add_column( Source="IWAM")
# #Compare IWAM estiamates of SMSY with those in Parken et al for test stocks (with UL and LL (these are 5th and 95th bootstrap estimates not CIs)
#  ParkenTestStocks <- read.csv("DataIn/ParkenTestStocks.csv") %>% rename(LL=SMSY5th, UL=SMSY95th) %>%
#    dplyr::select(-WA, -CV, -lh, -Area) %>% add_column (Source="Parken")
#  ParkenTestSMSY <- full_join(TestSMSY, ParkenTestStocks)
# ---------------------------------------------------------------------------------------------------



#---------------------------------------------------------------------------------
# Additional code; not currently needed 
#---------------------------------------------------------------------------------
# What initial values to use for WA model parameters?

getInits <- FALSE
if(getInits){
  SMSY <- All_Est %>% filter(Param=="SMSY") %>% mutate(ModelOrder=0:(length(unique(All_Est$Stocknumber))-1))
  # what is scale of SMSY?
  Sc <- SRDat %>% dplyr::select(Stocknumber, Scale) %>% distinct()
  SMSY <- SMSY %>% left_join(Sc) %>% mutate(rawSMSY=Estimate*Scale)
  SMSY <- SMSY %>% left_join(Stream, by=c("Stocknumber","ModelOrder"))
  lnSMSY <- log(SMSY$rawSMSY)
  lnWA <- log(WA$WA)
  order <- SMSY %>% dplyr::select(Stocknumber, ModelOrder)
  ParkenSMSY <- as.data.frame(read.csv("DataIn/ParkenSMSY.csv"))
  ParkenSMSY <- ParkenSMSY %>% left_join(order) %>% arrange(ModelOrder) %>% mutate(lnSMSY=log(SMSY))
  lnPSMSY <- ParkenSMSY$lnSMSY
  
  # lm(lnPSMSY ~ lnWA) #Get same coefficients as Parken et al. Table 4 for pooled data
  lnDelta1_start <- coef(lm(lnSMSY ~ lnWA))[1]
  lnDelta2_start <- log(coef(lm(lnSMSY ~ lnWA))[2])
  # without Skagit lnDelta1_start <- 2.999911
  # without Skagit lnDelta2_start <- -0.3238648, or Delta2 = 0.723348
  # With Skagit lnDelta1_start <- 2.881
  # with Skagit nDelta2_start <- -0.288
  
  #plot(y=exp(lnSMSY), x=exp(lnWA))
  #pdf("ParkenSMSYWA.pdf", width=4)
  #  par(mfcol=c(2,1))
  #  plot(y=exp(lnPSMSY), x=exp(lnWA), xlab="Watershed Area, km2", ylab="SMSY, Parken et al. 2006")
  #  plot(y=exp(lnPSMSY), x=exp(lnWA), xlim=c(0,2000), ylim=c(0,6000), xlab="Watershed Area, km2", ylab="SMSY, Parken et al. 2006")
  #dev.off()
}



#-----------------------------------------------------------------------------------------
test <- FALSE
if(test==TRUE){

#TEST: Watereshed-Area Regression with data inputs
SMSY <- read.csv("DataIn/SMSY_3mods.csv")
Sca <- SMSY %>% dplyr::select(Name, Scale)
SMSY <- SMSY %>% left_join(Stream) %>% left_join(WA)
SMSY_stream <- SMSY %>% filter(lh==0)
SMSY_ocean <- SMSY %>% filter(lh==1)

ParkenSMSY_stream <- read.csv("DataIn/ParkenSMSY.csv") %>% right_join(Stream) %>% 
  filter(lh==0) %>% mutate(Estimate=SMSY/SMSY_stream$Scale)


data <- list()
data$SMSY <- as.data.frame(read.csv("DataIn/WA_Parken.csv"))$SMSY / 
  as.data.frame(read.csv("DataIn/WA_Parken.csv"))$Scale
data$WA <- as.data.frame(read.csv("DataIn/WA_Parken.csv"))$WA 
data$Scale <- as.data.frame(read.csv("DataIn/WA_Parken.csv"))$Scale
data$Stream <- as.data.frame(read.csv("DataIn/WA_Parken.csv"))$lh
#data$Tau_dist <- 0.1

param <- list()
param$logDelta1 <- 3.00# with skagit 2.881
param$logDelta1ocean <- 0
param$logDelta2 <- log(0.72)
param$logDelta2ocean <- 0
#param$Delta2 <- log(0.72/(1-0.72)) #logit 0f 0.72 #with skagit logDelta2 = -0.288
param$logDeltaSigma <- -0.412 #from Parken et al. 2006 where sig=0.662

#dyn.unload(dynlib("TMB_Files/WAregression_sep"))
#compile("TMB_Files/WAregression_sep.cpp")
#dyn.unload(dynlib("TMB_Files/WAregression"))
#compile("TMB_Files/WAregression.cpp")
plot(x=log(data$WA), y=log(data$SMSY*data$Scale), xlab="ln(WA)", ylab="ln(SMSY)")
abline(lm(log(data$SMSY*data$Scale) ~ log(data$WA)))
plot(x=(data$WA), y=(data$SMSY*data$Scale), xlab="WA, km2", ylab="SMSY")

#dyn.load(dynlib("TMB_Files/WAregression"))
dyn.load(dynlib("TMB_Files/WAregression_sep"))
#obj <- MakeADFun(data, param, DLL="WAregression", silent=TRUE)
obj <- MakeADFun(data, param, DLL="WAregression_sep", silent=TRUE)

opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))#, upper=upper)
summary(sdreport(obj), p.value=TRUE)

All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)

# Put together readable data frame of values
All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))

All_Deltas <- data.frame()
All_Deltas <- All_Ests %>% filter (Param %in% c("logDelta1", "Delta2_bounded", "logDelta1ocean", 
                                                "logDelta2", "logDelta2ocean"))

png(paste("DataOut/Parken_comb.png"), width=7, height=7, units="in", res=500)
plotWAregression_Parken(data, All_Deltas)
dev.off()

png(paste("DataOut/Parken_sep.png"), width=7, height=7, units="in", res=500)
plotWAregression_ParkenSep(data, All_Deltas)
dev.off()

# I get same answers as Parken when I use his SMSY data, for streams (haven't checked ocean, 
# as Skagit is ocean type and numbering is wonky)
# For my SMSY_stream data slogDelta1 <- 2.744
# For my SMSY_stream data sDelta2 <- 0.857
# For my SMSY_stream data slogDeltaSigma <- -0.709

# make sure my output from the TMB model SMSY_stream and WA_stream match the values from:
#All_Est %>% filter(Param=="SMSY") %>% left_join(WA) %>% left_join(Stream) %>% filter(lh==1)
# Yes, they match

# WHen I run ocean and stream-type specific regressions, the SMSY values come out as a line! 
test <- All_Est %>% left_join(Sca) %>% filter(Param=="SMSY") %>% left_join(WA) %>% 
  left_join(Stream) %>% filter(lh==1)#switch to lh==0 for stream
lm(log(test$Estimate*test$Scale) ~ log(test$WA))
plot( y= log(test$Estimate*test$Scale), x=log(test$WA))

}
#-----------------------------------------------------------------------------------------



