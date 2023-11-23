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
library(TMB)
library(viridis)

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

# Current wcvi example uses 
remove.EnhStocks <- TRUE
Bern_logistic <- FALSE
prod <- "LifeStageModel"
LOO <- NA
run_logReg <- FALSE

Get.LRP.bs <- function(remove.EnhStocks=TRUE,  Bern_logistic=FALSE, 
                        prod="LifeStageModel", LOO = NA, run_logReg=FALSE){

  #--------------------------------------------------------------------------- #
  # Read in watershed area-based reference points (SREP and SMSY) --------------
  #--------------------------------------------------------------------------- #
  # Core data: 
    # New files:
  if(remove.EnhStocks) wcviRPs_long <- read.csv(here::here("DataOut/dataout_target_ocean_noEnh.csv"))
  #if(!remove.EnhStocks) wcviRPs_long<- read.csv(here::here("DataOut/dataout_target_ocean_wEnh.csv"))
    # Old files:
  # if (remove.EnhStocks) wcviRPs_long <- read.csv("DataOut/WCVI_SMSY_noEnh_wBC.csv")
  # if (!remove.EnhStocks) wcviRPs_long <- read.csv("DataOut/WCVI_SMSY_wEnh_wBC.csv")
  
  # Remove Cypre as it's not a core indicator (Diana McHugh, 22 Oct 2020)
  stock_SMSY <- wcviRPs_long %>% filter(Stock != "Cypre") %>% 
    filter (Param == "SMSY") %>% 
    rename(SMSY=Estimate, SMSYLL=LL, SMSYUL=UL) %>% 
    dplyr::select (-Param, -X)#, -CU)
  stock_SREP <- wcviRPs_long %>% filter(Stock != "Cypre") %>% 
    filter (Param == "SREP") %>% 
    rename(SREP=Estimate, SREPLL=LL, SREPUL=UL) %>% 
    dplyr::select (-Param, -X)
  wcviRPs <- stock_SMSY %>% left_join(stock_SREP, by="Stock")
  
  # Calculate scale for each stock
  digits <- count.dig(stock_SMSY$SMSY)
  Scale <- 10^(digits)
  
  #SREP_SE <- wcviRPs %>% mutate(SE = ((wcviRPs$SREP) - (wcviRPs$SREPLL)) / 1.96)
  SREP_logSE <- wcviRPs %>% mutate(SE = (log(wcviRPs$SREP) - log(wcviRPs$SREPLL)) / 1.96)
  # The UpperLimit-MLE gives same answer
  #SREP_logSE <- wcviRPs %>% mutate(SE = (log(wcviRPs$SREPUL) - log(wcviRPs$SREP)) / 1.96)
  SREP_logSE <- SREP_logSE %>% dplyr::select(Stock, SE)
  
  #--------------------------------------------------------------------------- #
  # Calculate Sgen 2 ways ------------------------------------------------------
  #--------------------------------------------------------------------------- #
  
  # 1.Use Watershed-area SMSY and SREP to estimate Sgen (assuming productivity
  # SGENcalcs <- purrr::map2_dfr (wcviRPs$SMSY/Scale,wcviRPs$SREP/Scale, Sgen.fn) 
  
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
    Mean.Ric.A <- 1 # Derived from life-history model (Luedke pers.comm.) and 
    # WCVI CK run reconstruction SR analysis (Dobson pers. comm.)
    Ric.A <- exp(rnorm(length(Scale), Mean.Ric.A, 0))
    
    
    # When incorporating uncertainty in Ricker A:
    Sig.Ric.A <- 0.51#0.255 #0.51 for a wider plausible bound
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
    if(min(Ric.A)<=0) Ric.A <- exp(rnorm(length(Scale), Mean.Ric.A, Sig.Ric.A))
    
    sREP <- exp(rnorm(length(Scale), log(wcviRPs$SREP), SREP_logSE$SE))
    if(min(sREP)<=0)   sREP <- exp(rnorm(length(Scale), wcviRPs$SREP, 
                                        SREP_SE$SE))
    if(min(sREP)<=0)   sREP <- exp(rnorm(length(Scale), wcviRPs$SREP, 
                                         SREP_SE$SE))
    
    SGENcalcs <- purrr::map2_dfr (Ric.A, sREP/Scale, Sgen.fn2)
   
     wcviRPs <- wcviRPs %>% mutate (SGEN = SGENcalcs$SGEN) %>% 
      mutate(SGEN=round(SGEN*Scale,0))
    wcviRPs <- wcviRPs %>% mutate (a.par = SGENcalcs$apar) %>% 
      mutate(a.par=round(a.par,2))
    wcviRPs <- wcviRPs %>% mutate (SMSY = SGENcalcs$SMSY) %>% 
      mutate(SMSY=round(SMSY*Scale,0))
    
    wcviRPs <- wcviRPs[c("Stock", "SGEN", "SMSY", "SMSYLL", "SMSYUL", "SREP", 
                         "SREPLL", "SREPUL", "a.par")]#"CU"
    
  }#End of if(prod == "LifeStageModel")
    
  # Ricker a's from Diana Dobson's Run Reconstruction (pers.comm) coded in TMB
  # Higher estimate of Ricker a (lower Sgen)  
  
  if(prod == "RunReconstruction"){
    lnalpha_inlet <- read.csv("DataIn/CUPars_wBC.csv") %>% 
      select(alpha,stkName) %>% rename(inlets=stkName, lnalpha=alpha)
    lnalpha_nBC_inlet <- read.csv("DataIn/CUPars_nBC.csv") %>% 
      select(alpha,stkName) %>% rename(inlets=stkName, lnalpha_nBC=alpha)
    WCVIStocks <- read.csv("DataIn/WCVIStocks.csv") %>% 
      filter (Stock != "Cypre") %>% rename(inlets=Inlet)
    Ric.A <- lnalpha_inlet %>% left_join(WCVIStocks, by="inlets") %>% select(c(lnalpha,inlets,CU,Stock))
    
    wcviRPs <- wcviRPs %>% left_join(Ric.A) %>% mutate(a.RR=exp(lnalpha))
    wcviRPs[wcviRPs$Stock=="Nitinat",]$a.RR <- exp(1)
    wcviRPs[wcviRPs$Stock=="San Juan",]$a.RR <- exp(1)
    wcviRPs[wcviRPs$Stock=="Nitinat",]$a.RR <- exp(1)
    
    wcviRPs[wcviRPs$Stock=="Barkley",]$a.RR <- 
      wcviRPs[wcviRPs$inlets=="Barkley",]$a.RR[1]
    wcviRPs[wcviRPs$Stock=="Clayoquot",]$a.RR <- 
      wcviRPs[wcviRPs$inlets=="Clayoquot",]$a.RR[1]
    wcviRPs[wcviRPs$Stock=="Kyuquot",]$a.RR <- 
      wcviRPs[wcviRPs$inlets=="Kyuquot",]$a.RR[1]
    wcviRPs[wcviRPs$Stock=="Nootka/Esperanza",]$a.RR <- 
      wcviRPs[wcviRPs$inlets=="Nootka/Esperanza",]$a.RR[1]
    wcviRPs[wcviRPs$Stock=="Quatsino",]$a.RR <- 
      wcviRPs[wcviRPs$inlets=="Quatsino",]$a.RR[1]
    wcviRPs[wcviRPs$Stock=="WCVI South",]$a.RR <- 
      wcviRPs[wcviRPs$inlets=="Barkley",]$a.RR[1]
    wcviRPs[wcviRPs$Stock=="WCVI Nootka & Kyuquot",]$a.RR <- 
      wcviRPs[wcviRPs$inlets=="Nootka/Esperanza",]$a.RR[1]
    wcviRPs[wcviRPs$Stock=="WCVI North",]$a.RR <- 
      wcviRPs[wcviRPs$inlets=="Quatsino",]$a.RR[1]
    
    wcviRPs <- wcviRPs %>% select(-c(inlets, CU, lnalpha)) %>% rename(a.par=a.RR)
    
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
    
    
    
    Ric.A.hi <- exp(rnorm(length(Scale), log(wcviRPs$a.par), Sig.Ric.A))
    if(min(Ric.A.hi)<0) Ric.A <- exp(rnorm(length(Scale), wcviRPs$a.RR, Sig.Ric.A))
    
    sREP <- exp(rnorm(length(Scale), log(wcviRPs$SREP), SREP_logSE$SE))
    if(min(sREP)<0)   sREP <- exp(rnorm(length(Scale), wcviRPs$SREP, 
                                        SREP_SE$SE))
    
    
    SGENcalcs <- purrr::map2_dfr (Ric.A.hi, sREP/Scale, Sgen.fn2)
    
    wcviRPs <- wcviRPs %>% mutate (SGEN = SGENcalcs$SGEN) %>% 
      mutate(SGEN=round(SGEN*Scale,0))
    wcviRPs <- wcviRPs %>% mutate (a.par = SGENcalcs$apar) %>% 
      mutate(a.par=round(a.par,2))
    wcviRPs <- wcviRPs %>% mutate (SMSY = SGENcalcs$SMSY) %>% 
      mutate(SMSY=round(SMSY*Scale,0))
    
    wcviRPs <- wcviRPs[c("Stock", "SGEN", "SMSY", "SMSYLL", "SMSYUL", "SREP", 
                         "SREPLL", "SREPUL", "a.par")]#"CU"
 
  } # if(prod == "RunReconstruction"){
  
  
  
  #--------------------------------------------------------------------------- #
  # Add Sgen and revised SMSY to wcviRPs data frame ----------------------------
  #--------------------------------------------------------------------------- #
  

  
  wcviRPs
  # # Write this to a csv file so that it can be called in plotting functions
  # # write.csv(wcviRPs, "DataOut/wcviRPs.csv")
  # if (remove.EnhStocks) write.csv(wcviRPs, "DataOut/wcviRPs_noEnh.csv")
  # if (!remove.EnhStocks) write.csv(wcviRPs, "DataOut/wcviRPs_wEnh.csv")
  
  
  #--------------------------------------------------------------------------- #
  # Scenario if productivity is reducted by half, as in WSP SBC CK assessment ----
  # (DFO 2014). NOT NEEDED
  #--------------------------------------------------------------------------- #
  # 
  # SGENcalcsv2 <- map2_dfr (wcviRPs$SMSY/Scale,wcviRPs$SREP/Scale, 
  #                          Sgen.fn, half.a = TRUE, const.SMAX = FALSE)
  # wcviRPs <- wcviRPs %>% mutate (SGENha.cSREP = SGENcalcsv2$SGEN) %>% 
  #   mutate( SGENha.cSREP = round( SGENha.cSREP*Scale, 0 ) )
  # wcviRPs <- wcviRPs %>% mutate (SMSYha.cSREP = SGENcalcsv2$SMSY) %>% 
  #   mutate( SMSYha.cSREP = round( SMSYha.cSREP*Scale, 0 ) )
  # wcviRPs <- wcviRPs %>% mutate (SREPha.cSREP = SGENcalcsv2$SREP) %>% 
  #   mutate( SREPha.cSREP = round( SREPha.cSREP*Scale, 0 ) )
  # ###wcviRPs <- wcviRPs %>% mutate (SMAXrev = 1/SGENcalcs$bpar) %>% 
  # ###mutate(SMAXrev=round(SMAXrev,0))
  # 
  # SGENcalcsv3 <- map2_dfr (wcviRPs$SMSY/Scale, wcviRPs$SREP/Scale, 
  #                          Sgen.fn, half.a = TRUE, const.SMAX = TRUE)
  # wcviRPs <- wcviRPs %>% mutate (SGENha.cSMAX = SGENcalcsv3$SGEN) %>% 
  #   mutate( SGENha.cSMAX = round( SGENha.cSMAX*Scale, 0 ) )
  # wcviRPs <- wcviRPs %>% mutate (SMSYha.cSMAX = SGENcalcsv3$SMSY) %>% 
  #   mutate( SMSYha.cSMAX = round( SMSYha.cSMAX*Scale, 0 ) )
  # wcviRPs <- wcviRPs %>% mutate (SREPha.cSMAX = SGENcalcsv3$SREP) %>% 
  #   mutate( SREPha.cSMAX = round( SREPha.cSMAX*Scale, 0 ) )

  # run_logReg <- TRUE
  if(run_logReg==FALSE){
    return(list(bench= select(SGENcalcs,-apar, -bpar)*Scale))
    
  }
  #--------------------------------------------------------------------------- #
  # Sum escapements across indicators within inlets ----------------------------
  #--------------------------------------------------------------------------- #
  if(run_logReg==TRUE){
    WCVIEsc <- data.frame(read.csv("DataIn/WCVIEsc.csv", row.names="Yr")) %>% 
      dplyr::select (-"Little.Zeballos")
    
    # Take "." out of name as in escapement data
    WCVIEsc_names <- sapply(colnames(WCVIEsc), 
                            function(x) (gsub(".", " ", x, fixed=TRUE) ) )
    WCVIEsc_names <- sapply(WCVIEsc_names, function(x) 
      (gsub("Bedwell Ursus", "Bedwell/Ursus", x, fixed=TRUE) ) )
    WCVIEsc_names <- sapply(WCVIEsc_names, function(x) 
      (gsub("Nootka Esperanza", "Nootka/Esperanza", x, fixed=TRUE) ) )
    colnames(WCVIEsc) <- WCVIEsc_names 
    
    EnhStocks <- data.frame(read.csv("DataIn/WCVIstocks.csv")) %>% filter (Enh==1) %>%
      pull(Stock)
    EnhStocks <- as.character(EnhStocks)
    
    #EnhStocks <- c("Burman",  "Conuma", "Leiner", "Nitinat", "Sarita",  
    #               "Somass",  "Zeballos", "San Juan", "Tranquil")
    # Artlish removed from Enhanced stocks 23 Dec. 2020
    # Tranquil added 18 Jan 2021
    
    
    if (remove.EnhStocks) {WCVIEsc <- WCVIEsc %>% dplyr::select(-EnhStocks) }
    
    Years <- rownames(WCVIEsc)
    
    # Get stock information for WCVI Chinook & Remove Cypre as it's not an 
    # indicator stocks
    WCVIStocks <- read.csv("DataIn/WCVIStocks.csv") %>% 
      filter (Stock != "Cypre")
    if (remove.EnhStocks) WCVIStocks <- WCVIStocks %>% 
      filter(Stock %not in% EnhStocks)
    
    Inlet_Names <- unique(WCVIStocks$Inlet)
    Inlet_Sum <- matrix(NA, nrow=length(Years), ncol=length(Inlet_Names))
    colnames(Inlet_Sum) <- Inlet_Names
    CU_Names <- unique(WCVIStocks$CU)
    
    
    # Sum escapements across stocks within inlets
    for (i in 1:length(Inlet_Names)) {
      # For each inlet, which are the component indicator stocks
      Ins <- WCVIStocks %>% filter(Inlet==Inlet_Names[i]) %>% pull(Stock)
      WCVIEsc_Inlets <- matrix(NA, nrow= length(Years), ncol= length(Ins))
      
      #  Make a matrix of escapements of component indicators
      for (j in 1:length(Ins)){
        WCVIEsc_Inlets[,j] <- WCVIEsc %>% 
          dplyr::select(as.character(Ins[j])) %>% pull()
        
      }
      
      # Sum the escapement of component indicators, setting sum=NA for years 
      # where there are any NAs
      Inlet_Sum[,i] <- apply(WCVIEsc_Inlets, 1, sum, na.rm=F)
    }
    
    
    #------------------------------------------------------------------------- #
    # Sum escapements across indicators within CUs -----------------------------
    #------------------------------------------------------------------------- #
    
    nCU <- length(unique(WCVIStocks$CU))
    CU_Sum <- matrix(NA, nrow=length(Years), ncol=nCU)
    colnames(CU_Sum) <- CU_Names
    
    for (k in 1:length(CU_Names)) {
      # For each CU, which are the component indicators
      CUss <- WCVIStocks %>% filter(CU==CU_Names[k]) %>% pull(Stock)
      WCVIEsc_CUs <- matrix(NA, nrow= length(Years), ncol= length(CUss))
      
      #  Make a matrix of escapements of component indicators
      for (j in 1:length(CUss)){
        WCVIEsc_CUs[,j] <- WCVIEsc %>% dplyr::select(as.character(CUss[j])) %>% 
          pull()
      }
      
      # Sum the escapement of component indicators, setting sum=NA for years 
      # where there are any NAs
      CU_Sum[,k] <- apply(WCVIEsc_CUs, 1, sum, na.rm=F)
    }
    # Remove double incidence of Nitinat and San Juan (= stock and an inlet) 
    # when enhancement is included
    if(!remove.EnhStocks) WCVIEsc <- WCVIEsc %>% 
      dplyr::select(-Nitinat, -'San Juan')
    
    WCVIEsc <- cbind(WCVIEsc, Inlet_Sum, CU_Sum) 
    
    #------------------------------------------------------------------------- #
    # Assess status for each inlet relative to inlet-level SGEN for each year ----
    #   Is Inlet level escapement above inlet-level Sgen: T/F?
    #   Inlet_Status = FALSE if summed escapement is below Sgen
    #   Inlet_Status = TRUE if summed escapement is below Sgen
    #------------------------------------------------------------------------- #
    
    Inlet_Status <- matrix(NA, nrow=length(Years), ncol=length(Inlet_Names) )
    colnames(Inlet_Status) <- Inlet_Names
    
    for (i in 1:length(Inlet_Names)) {
      Inlet_Status[,i] <- (Inlet_Sum[,i] > 
                             (wcviRPs %>% 
                                filter(Stock == as.character(Inlet_Names[i])) %>% 
                                pull(SGEN)) )
    }
    
    Inlet_Status <- as.data.frame(Inlet_Status)
    
    #------------------------------------------------------------------------- #
    # Assess status for each CU for each year of the time-series ---------------
    #   (floor of summed CU-level numeric statuses)
    #   CU_Status = below LB if any inlet within the CU below their Sgen = 0 
    #   CU_Status = above LB if all inlets within the CU above their Sgen = 1
    #------------------------------------------------------------------------- #
    
    CU_Status <- matrix(NA, nrow=length(Years), ncol=length(CU_Names))
    colnames(CU_Status) <- CU_Names
    
    stock.LRP <- TRUE
    if(stock.LRP){
      for (k in 1:length(CU_Names)) {
        # For each CU, which are the component indicators
        CU_ins <- unique( WCVIStocks %>% filter(CU==CU_Names[k]) %>% pull(Inlet))
        
        isAbove <- matrix(NA, nrow= length(Years), ncol= length(CU_ins))
        
        # Make a matrix of status of component inlets. Is each inlet > 
        # Sgen values, T/F?
        for (i in 1:length(CU_ins)){
          isAbove[,i] <- Inlet_Status %>% 
            dplyr::select(as.character(CU_ins[i])) %>% pull()
        }
        
        # CU-level status: are ALL inlets above their Sgen values?
        # Sum the "true"
        isAboveFun <- function(x){ 
          floor(sum( as.numeric(x), na.rm=F) / length(x) ) }
        CU_Status[,k] <- apply(X= isAbove, MARGIN = 1, FUN=isAboveFun)
      }
      CU_Status <- as.data.frame(CU_Status)
    }
    
    # Alternatively, are there CU-level Sgen values to derive CU-level status?
    CU.LRP <- FALSE
    if(CU.LRP){
      for (k in 1:length(CU_Names)){
        CU_Status[,k] <- (CU_Sum[,k] > 
                            (wcviRPs %>% filter(Stock == 
                                                  as.character(CU_Names[k])) %>% 
                               pull(SGEN)) )
      }
      CU_Status <- as.data.frame(CU_Status)
      
    }
    
    #------------------------------------------------------------------------- #
    # Proportion of CUs that are not in the red zone ---------------------------
    #------------------------------------------------------------------------- #
    
    ppnAboveFun <- function(x) {sum( as.numeric(x), na.rm=F) / length(x) }
    SMU_ppn <- apply(X=CU_Status, MARGIN=1, FUN=ppnAboveFun)
    
    #------------------------------------------------------------------------- #
    # Logistic regression ------------------------------------------------------
    #------------------------------------------------------------------------- #
    # Get SMU-level escapement time-series
    SMU_Esc <- apply(Inlet_Sum, 1, sum, na.rm=F)
    
    SMUlogisticData <- data.frame(SMU_Esc) %>% 
      add_column(ppn=SMU_ppn, Years=as.numeric(Years)) %>% 
      filter(SMU_Esc != "NA")
    
    data <- list()
    data$N_Stks <- length(CU_Names)
    digits <- count.dig(SMU_Esc)
    ScaleSMU <- min(10^(digits -1 ), na.rm=T)
    
    data$LM_Agg_Abund <- SMUlogisticData$SMU_Esc/ScaleSMU
    data$N_Above_BM <- SMUlogisticData$ppn * data$N_Stks
    
    if(!is.na(LOO)) { #If applying leave-one-out cross validation, remove that
      #year
      data$LM_Agg_Abund <- data$LM_Agg_Abund[-LOO]
      data$N_Above_BM <- data$N_Above_BM[-LOO]
    }
    data$Pred_Abund <- seq(0, max(data$LM_Agg_Abund)*1.1, 0.1)
    if(remove.EnhStocks) data$Pred_Abund <- 
      seq(0, max(data$LM_Agg_Abund)*1.5, 0.1)
    data$p <- 0.95#0.67
    
    if(Bern_logistic==FALSE) data$Penalty <- as.numeric(TRUE)
    if(Bern_logistic==TRUE) data$Penalty <- as.numeric(FALSE)
    data$Bern_logistic <- as.numeric(Bern_logistic)
    
    # Add a normally distributed penalty on aggregate abundances 
    # when p is very small (0.01) 
    # Lower 95% CL = abundance of the smallest CU in its lowest abundance 
    # year (most CUs lost, only 1 remains below LB)
    # Upper 95% CL = abundance of the ave annual abundance of the sum of 
    # across CUs. 
    # If CU-level benchmarks exist can sum those benchmarks for Upper 95%CL
    min <- min(apply(CU_Sum, 2, min, na.rm=T), na.rm=T)
    max <- mean(apply(CU_Sum, 1, sum, na.rm=F), na.rm=T)
    # Parameters for normal penalty (mu, sig):
    # mu  = ave of min and max values. 
    # sig = SD which allows the density = 0.05 at min and max values
    # sum(dnorm(seq(min,max,1), mean=mean(c(min,max)), sd=22400))# Should 
    # give 95% density
    # plot(x=seq(min,max,100), y=dnorm(seq(min,max,100), 
    # mean=mean(c(min,max)), sd=22400),type="l")
    data$B_penalty_mu <- mean(c(min,max))/ScaleSMU
    if (!remove.EnhStocks) data$B_penalty_sig <- 22400/ScaleSMU
    # Should give 95% density:
    # sum(dnorm(seq(min,max,1), mean=mean(c(min,max)), sd=2700))
    if (remove.EnhStocks) data$B_penalty_sig <- 2700/ScaleSMU
    
    param <- list()
    param$B_0 <- -2
    param$B_1 <- 0.1
    
    
    # dyn.unload(dynlib(paste("TMB_Files/Logistic_LRPs", sep="")))
    # compile(paste("TMB_Files/Logistic_LRPs.cpp", sep=""))
    
    dyn.load(dynlib(paste("TMB_Files/Logistic_LRPs", sep=""))) 
    
    obj <- MakeADFun(data, param, DLL="Logistic_LRPs", silent=TRUE) 
    
    opt <- nlminb(obj$par, obj$fn, obj$gr, control = 
                    list(eval.max = 1e5, iter.max = 1e5)) 
    pl <- obj$env$parList(opt$par) 
    #summary(sdreport(obj), p.value=TRUE)
    
    # Get parameter estimates and logit predicted values for CIs
    All_Ests <- data.frame(summary(sdreport(obj), p.value=TRUE))
    All_Ests$Param <- row.names(All_Ests)
    All_Ests$Param <- sapply(All_Ests$Param, function(x) 
      (unlist(strsplit(x, "[.]"))[[1]]))
    Preds <- All_Ests %>% filter(Param == "Logit_Preds")
    All_Ests <- All_Ests %>% filter(!(Param %in% c( "Logit_Preds"))) 
    
    
    
    
    out <- list()
    out$All_Ests <- All_Ests
    
    
    Logistic_Data <- data.frame(yr = SMUlogisticData$Years, 
                                yy = SMUlogisticData$ppn, 
                                xx = SMUlogisticData$SMU_Esc)
    
    out$Logistic_Data <- Logistic_Data
    
    Logistic_Fits <- data.frame(xx = data$Pred_Abund*ScaleSMU, 
                                fit = inv_logit(Preds$Estimate),
                                lwr = inv_logit(Preds$Estimate - 
                                                  1.96*Preds$Std..Error),
                                upr = inv_logit(Preds$Estimate + 
                                                  1.96*Preds$Std..Error))
    
    out$Preds <- Logistic_Fits
    
    out$LRP <- data.frame(fit = (All_Ests %>% filter(Param == "Agg_LRP") %>% 
                                   pull(Estimate))*ScaleSMU, 
                          lwr = (All_Ests %>% filter(Param == "Agg_LRP") %>% 
                                   mutate(xx =Estimate - 1.96*Std..Error) %>% 
                                   pull(xx) ) * ScaleSMU,
                          upr = (All_Ests %>% filter(Param == "Agg_LRP") %>% 
                                   mutate(xx =Estimate + 1.96*Std..Error) %>% 
                                   pull(xx) ) * ScaleSMU)
    
    return(list(out=out, # C
                WCVIEsc=WCVIEsc, # C
                SMU_Esc=SMU_Esc, # C
                CU_Status=CU_Status, # C
                SMU_ppn=SMU_ppn, # C
                LRPppn= data$p, # C
                nLL=obj$report()$ans, # within
                LOO=LOO, # Outside
                bench=select(SGENcalcs, -apar, -bpar)*Scale))
  }
  # FUNCTION STALLING HERE
  # return(out=out, WCVIEsc=WCVIEsc, SMU_Esc=SMU_Esc,
  #             CU_Status=CU_Status, SMU_ppn=SMU_ppn,
  #             LRPppn= data$p, nLL=obj$report()$ans, LOO=LOO,
  #             bench= select(SGENcalcs,-apar, -bpar)*Scale)
} 

# Return statements are not working
# When removed function runs completely. 

# **********************************************************************************************************************
# Eng of Get.LRP.bs() function ----
# **********************************************************************************************************************



#----------------------------------------------------------------------------- #
# Run bootstraps to derive benchmarks with uncertainty (see LRP code) ------
run.bootstraps <- FALSE

if (run.bootstraps){
  set.seed(100)#10#12#13(work for 1000)
  nBS <- 5000 # number trials for bootstrapping
  outBench <- list() 
  
  for (k in 1:nBS) {
    out <- Get.LRP.bs()
    outLRP <- as.data.frame(out$out$LRP) 
    if(k==1) LRP.bs <- data.frame(fit=outLRP$fit, upr=outLRP$upr, lwr=outLRP$lwr)
    if(k>1) LRP.bs <- add_row(LRP.bs, outLRP)
    
    outBench[[k]] <- out$bench
  }
  
  # # Is 200 enough trials? Yes
  # running.mean <- cumsum(LRP.bs$fit) / seq_along(LRP.bs$fit) 
  # plot(running.mean)
  
  # Calculate distribution of overall LRPs by integrating bootstrapped LRP 
  # values with uncertainty of each LRP value from TMB
  LRP.samples <- rnorm(nBS*10, LRP.bs$fit, (LRP.bs$fit - LRP.bs$lwr) / 1.96)
  hist(LRP.samples)
  LRP.boot <- quantile(LRP.samples, probs=c(0.025, 0.5, 0.975))
  names(LRP.boot) <- c("lwr", "LRP", "upr")
  
  # Compile bootstrapped estimates of Sgen, SMSY, and SREP, and identify 5th and 
  # 95th percentiles
  SGEN.bs <- select(as.data.frame(outBench), starts_with("SGEN"))
  stockNames <- read.csv("DataOut/WCVI_SMSY_noEnh_wBC.csv") %>% 
    filter(Stock != "Cypre") %>% pull(Stock)
  stockNames <- unique(stockNames)
  
  rownames(SGEN.bs) <- stockNames
  SGEN.boot <- data.frame(SGEN= apply(SGEN.bs, 1, quantile, 0.5), 
                          lwr=apply(SGEN.bs, 1, quantile, 0.025),
                          upr=apply(SGEN.bs, 1, quantile, 0.975) )
  
  SMSY.bs <- select(as.data.frame(outBench), starts_with("SMSY"))
  rownames(SMSY.bs) <- stockNames
  SMSY.boot <- data.frame(SMSY= apply(SMSY.bs, 1, quantile, 0.5), 
                          lwr=apply(SMSY.bs, 1, quantile, 0.025),
                          upr=apply(SMSY.bs, 1, quantile, 0.975) )
  
  SREP.bs <- select(as.data.frame(outBench), starts_with("SREP"))
  rownames(SREP.bs) <- stockNames
  SREP.boot <- data.frame(SREP= apply(SREP.bs, 1, quantile, 0.5), 
                          lwr=apply(SREP.bs, 1, quantile, 0.025),
                          upr=apply(SREP.bs, 1, quantile, 0.975) )
  
  boot <- list(LRP.boot=LRP.boot, SGEN.boot=SGEN.boot, SMSY.boot=SMSY.boot, 
               SREP.boot=SREP.boot)
  
  df1 <- data.frame(boot[["SGEN.boot"]], Stock=rownames(boot[["SGEN.boot"]]), RP="SGEN") 
  df1 <- df1 %>% rename(Value=SGEN)
  df2 <- data.frame(boot[["SREP.boot"]], Stock=rownames(boot[["SREP.boot"]]), RP="SREP")
  df2 <- df2 %>% rename(Value=SREP)
  df3 <- data.frame(boot[["SMSY.boot"]], Stock=rownames(boot[["SMSY.boot"]]), RP="SMSY")
  df3 <- df3 %>% rename(Value=SMSY)  
  dfout <- add_row(df1, df2)
  dfout <- add_row(dfout, df3)
  rownames(dfout) <- NULL
  write.csv(dfout, "DataOut/wcviCK-BootstrappedRPs1000v3.csv") 
}

#----------------------------------------------------------------------------- #
# Run bootstraps to derive LRPs with uncertainty in benchmarks -------------
  # ************* FORMAL SUBMISSION VALUES *******************************************************************************
# See implementation of this in WCVI_LRPs.Rmd

#run.bootstraps <- FALSE
run.bootstraps <- TRUE

if (run.bootstraps){
  set.seed(1) #10#12#13 (work for 1000), for 100, 200, 300, (for 5000trials), 1, 2, 3 (for 20000trials)
  nBS <- 10 # number trials for bootstrapping
  outBench <- list() 
  
  for (k in 1:nBS) {
    out <- Get.LRP.bs(run_logReg=FALSE)
    outBench[[k]] <- out$bench
  }

  # # Is 200 enough trials? Yes
  # running.mean <- cumsum(LRP.bs$fit) / seq_along(LRP.bs$fit) 
  # plot(running.mean)
  

  # Compile bootstrapped estimates of Sgen, SMSY, and SREP, and identify 5th and 
  # 95th percentiles
  SGEN.bs <- select(as.data.frame(outBench), starts_with("SGEN"))
  # stockNames <- read.csv("DataOut/WCVI_SMSY_noEnh_wBC.csv") %>% 
  stockNames <- read.csv(here::here("DataOut/dataout_target_ocean_noEnh.csv")) %>% 
    filter(Stock != "Cypre") %>% pull(Stock)
  stockNames <- unique(stockNames)

  rownames(SGEN.bs) <- stockNames
  SGEN.boot <- data.frame(SGEN= apply(SGEN.bs, 1, quantile, 0.5), 
                          lwr=apply(SGEN.bs, 1, quantile, 0.025),
                          upr=apply(SGEN.bs, 1, quantile, 0.975) )
  
  SMSY.bs <- select(as.data.frame(outBench), starts_with("SMSY"))
  rownames(SMSY.bs) <- stockNames
  SMSY.boot <- data.frame(SMSY= apply(SMSY.bs, 1, quantile, 0.5), 
                          lwr=apply(SMSY.bs, 1, quantile, 0.025),
                          upr=apply(SMSY.bs, 1, quantile, 0.975) )
  
  SREP.bs <- select(as.data.frame(outBench), starts_with("SREP"))
  rownames(SREP.bs) <- stockNames
  SREP.boot <- data.frame(SREP= apply(SREP.bs, 1, quantile, 0.5), 
                          lwr=apply(SREP.bs, 1, quantile, 0.025),
                          upr=apply(SREP.bs, 1, quantile, 0.975) )
  
  boot <- list(SGEN.boot=SGEN.boot, SMSY.boot=SMSY.boot, 
               SREP.boot=SREP.boot)
  
  df1 <- data.frame(boot[["SGEN.boot"]], Stock=rownames(boot[["SGEN.boot"]]), RP="SGEN") 
  df1 <- df1 %>% rename(Value=SGEN)
  df2 <- data.frame(boot[["SREP.boot"]], Stock=rownames(boot[["SREP.boot"]]), RP="SREP")
  df2 <- df2 %>% rename(Value=SREP)
  df3 <- data.frame(boot[["SMSY.boot"]], Stock=rownames(boot[["SMSY.boot"]]), RP="SMSY")
  df3 <- df3 %>% rename(Value=SMSY)  
  dfout <- add_row(df1, df2)
  dfout <- add_row(dfout, df3)
  rownames(dfout) <- NULL
  # now round to 2 signif digits
  dfout <- dfout %>% mutate(Value=signif(Value, 2)) %>% 
    mutate(lwr=signif(lwr,2)) %>% 
    mutate (upr=signif(upr,2))
  
  write.csv(dfout, here::here("DataOut/wcviCK-BootstrappedRPs.csv"))

}
#Found with 20000 the results have stabilized to two significant digits. Use 20000v1, 
# and recommend 2 signifcant digits to users
# To check the SE is within 2% with 5 trials of 20,000

# Commented out for the purpose of the wcvi_workedexample.RMD
# d1 <- read.csv("DataOut/wcviCK-BootstrappedRPs20000v1.csv")
# d2 <- read.csv("DataOut/wcviCK-BootstrappedRPs20000v2.csv")
# d3 <- read.csv("DataOut/wcviCK-BootstrappedRPs20000v3.csv")
# d4 <- read.csv("DataOut/wcviCK-BootstrappedRPs20000v4.csv")
# d5 <- read.csv("DataOut/wcviCK-BootstrappedRPs20000v5.csv")
# n<-length(d1$Value)
# sd.Value <- NA; sd.lwr <- NA; sd.upr <- NA
# mean.Value <- NA; mean.lwr <- NA; mean.upr <- NA
# se.Value <- NA; se.lwr <- NA; se.upr <- NA
# for (i in 1:n){
#  sd.Value[i] <- sd(c(d1$Value[i], d2$Value[i], d3$Value[i], d4$Value[i], d5$Value[i]))
#  mean.Value[i] <- mean(c(d1$Value[i], d2$Value[i], d3$Value[i], d4$Value[i], d5$Value[i]))
#  se.Value[i] <- sd.Value[i]/mean.Value[i]
#  sd.lwr[i] <- sd(c(d1$lwr[i], d2$lwr[i], d3$lwr[i], d4$lwr[i], d5$lwr[i]))
#  mean.lwr[i] <- mean(c(d1$lwr[i], d2$lwr[i], d3$lwr[i], d4$lwr[i], d5$lwr[i]))
#  se.lwr[i] <- sd.lwr[i]/mean.lwr[i]
#  sd.upr[i] <- sd(c(d1$upr[i], d2$upr[i], d3$upr[i], d4$upr[i], d5$upr[i]))
#  mean.upr[i] <- mean(c(d1$upr[i], d2$upr[i], d3$upr[i], d4$upr[i], d5$upr[i]))
#  se.upr[i] <- sd.upr[i]/mean.upr[i]
# }

#----------------------------------------------------------------------------- #
# R version of logistic regression ---------------------------------------------
#   This matches results from TMB code when penalty=FALSE
#----------------------------------------------------------------------------- #
R.logReg <- FALSE

if (R.logReg) {
  ModDat <- data.frame(xx=data$LM_Agg_Abund, yy=SMUlogisticData$ppn)
  #or family=binomial, which gives much larger SEs, and assumes var=1.
  Fit_Mod <- glm( yy ~ xx , family = quasibinomial, data=ModDat)
  summary(Fit_Mod)$coefficients
  LRP <- (log(data$p/(1-data$p)) - Fit_Mod$coefficients[[1]])/ 
    Fit_Mod$coefficients[[2]]
  # use MASS function to get "dose" 
  library(MASS)
  Dose <- dose.p(Fit_Mod, p=data$p)
  Dose
  
  #  - Make x vector to predict with this model, for plotting
  xx <- data.frame(xx = seq(0, max(data$LM_Agg_Abund*1.25), 
                            by=(max(data$LM_Agg_Abund*1.25)/1000)))
  
  # - Create model predictions that include standard error fit
  preds <- predict.glm(Fit_Mod, newdata = xx, type = "link", se.fit = TRUE)
  
  # Create a confidence interval (lwr, upr) on the link scale as the fitted 
  # value plus or minus 1.96 times the standard error:
  critval <- 1.96 ## approx 95% CI
  upr <- preds$fit + (critval * preds$se.fit)
  lwr <- preds$fit - (critval * preds$se.fit)
  fit <- preds$fit
  
  # Transform confidence interval using the inverse of the link function to 
  # map the fitted values and the upper and lower limits of the interval
  # back on to the response scale:
  fit2 <- Fit_Mod$family$linkinv(fit)
  upr2 <- Fit_Mod$family$linkinv(upr)
  lwr2 <- Fit_Mod$family$linkinv(lwr)
  
  # Load predicted data on response scale into dataframe
  preddata<-xx
  preddata$fit<-fit2
  preddata$lwr <- lwr2 
  preddata$upr <- upr2 
  
  ### Calculate confidence intervals for LRP 
  LRP_lwr <- LRP - critval*attr(Dose, "SE")[[1]]
  LRP_upr <- LRP + critval*attr(Dose, "SE")[[1]]
  # These seem much tighter than other method?
  
  # Outputs for logistic regression in R: not needed
  list.out<-list()
  list.out$Logistic_Data <- ModDat
  list.out$model <- Fit_Mod
  list.out$Preds <- preddata
  list.out$LRP<-data.frame(fit = LRP, lwr = LRP_lwr, upr = LRP_upr)
  #list.out
  
}


#----------------------------------------------------------------------------- #
# TMB version of code to estimate Sgen from SMSY and SREP (WA_Sgen.cpp): -------
# NOT WORKING
#----------------------------------------------------------------------------- #

# # Do not need TMB code given I need to run this over bootstraps PRIOR to 
# # logististic regression. 
# # I could input all boostrapped PI draws, and estimate LRP internally for 
# # each draw, however wrangling with data is difficult 
# # in TMB and since there is only one estimation step: logistic regression, 
# # could simply implement in R
# 
# 
# SMSY <- wcviRPs %>% pull(SMSY)
# SREP <- wcviRPs %>% pull(SREP)
# 
# # Calculate scale for each stock
# digits <- count.dig(SMSY)
# Scale <- 10^(digits)
# 
# SMSY <- SMSY/Scale
# SREP <- SREP/Scale
# 
# 
# data <- list()
# data$SMSY <- SMSY
# data$SREP <- SREP
# data$Inlets <- read.csv("DataIn/WCVIStocks.csv") %>% 
# filter (Stock != "Cypre") %>% pull(SA_ind)
# data$N_inlets <- length(unique(read.csv("DataIn/WCVIStocks.csv") %>% 
# filter (Stock != "Cypre") %>% pull(SA_ind)))
# 
# #data$Scale <- Scale
# 
# 
# param <- list()
# param$RicB <- 1/(data$SREP/2) #initialize SMAX at half SREP
# param$logSgen <- log(data$SMSY/2)#initialize SMAX at half SMSY
# 
# # Compile model if changed:
# #dyn.unload(dynlib("TMB_Files/WA_Sgen"))
# #compile("TMB_Files/WA_Sgen.cpp")
# dyn.load(dynlib("TMB_Files/WA_Sgen"))
# obj <- MakeADFun(data, param, DLL="WA_Sgen", silent=TRUE)
# 
# # b is bounded between 1/3 of SREP and SREP
# lower <- 1/data$SREP
# upper <- 3/data$SREP
# opt <- nlminb(obj$par, obj$fn, obj$gr, control = 
# list(eval.max = 1e5, iter.max = 1e5))#, lower=rep(0,20), upper=log(SMSY))
# pl <- obj$env$parList(opt$par) 
# #summary(sdreport(obj), p.value=TRUE)
# 
# 
# # The summation of Sgens across inlets is not working because Sgen's 
# # are each scaled differntly
# # Actually, do this in R as organizing data is easier there 
# # (and just almost as fast to run)
# 
# #exp(pl$logSgen)*Scale
# #1/((1/pl$RicB)*Scale)
