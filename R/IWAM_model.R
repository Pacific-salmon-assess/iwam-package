
#### Introduction -------------------------------------------------------------

# This model estimates biological benchmarks for Chinook populations based on 
# watershed area of the spawning habitat
# The underlying model uses the relationship between watershed area and 
# stock-recruitment parameters for a set of Chinook Populations across the NE 
# (synoptic data set) to derive stock-recruitment parameters (and associated 
# benchmarks) from novel populations from their watershed areas. 
# The model is adapted from Parken et al. (2007) and Liermann et al. (2012)
# This version is taken from IWAM.R (developed by C. Holt, adapted by T. 
# Kitching)

# Integrated Watershed Area Model
# Section Descriptions:
# 1. Read in stock-recruitment data, life-history type, and watershed areas
#   for synoptic data set, and watershed area and life-history type for 
#   additional stocks
#   1.a. Data cleaning
#   1.b. Scale calculation
#   1.c. Setup Watershed area sets
# 2. Create data and parameter lists for TMB
# 3. Run TMB model to estimate Ricker parameters and SMSY & SREP for synoptic 
#   data sets, estimate paraemeters of watershed-area regression, and 
#   estimate SMSY and SREP for additional stocks 
# 4. Compile model outputs
# 5. Calculate diagnostics for SR models in synoptic data set and plot SR 
#   curves, etc.
# 6. Calculate prediction intervals for SMSY and SREP estimates for additional 
#   "test" stocks. These are written to a *.csv file

# This reduced code snippet has removed the main wrapper function- To be added

# Future as a function:
#   - Section 1 e.g. data cleaning and setup not included
#   - Section 1 scale calculation should be embedded as either a separate 
#     function or as an input
#   - Section 2 requires a list of inputs both from the "core" data file
#     and a params list from TMB
#   - Section 3 and 4 are basic
#   - Section 5 could be made so that the ouputs are easily then used in 
#     the already created PlotFunctions.R OR to have an embedded plots=TRUE
#   - Section 6 could be worked in similarily to 5
#   - *Remember to add libraries/dependencies to DESCRIPTION and NAMESPACE*

#### Libraries -----------------------------------------------------------------

library(rsample)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(reshape2)
library(TMB)
library(zoo)
library(viridis)
library(hrbrthemes)

# Both helperFunctions and PlotFunctions are required.

# Tor- based on your experience with COSEWIC R package, can we remove here::here 
# to avoid problems when using in pkg? And in the meantime, if we're running 
# from *.Proj files, a simple source("R/helperFunctions.R") should work.Correct?
source (here::here("R/helperFunctions.R"))
source(here::here("R/PlotFunctions.R"))
# source(here::here("R/WCVILRPs_bootstrap.R"))
source(here::here("R/Get_LRP_bs.R"))

# Consider renaming all model.R scripts to include "mod" or some
  # other suffix/prefix
# To do - rename this to IWAM_model.R and rename IWAM.R to IWAMarchived_model.R


#### Remaining wrapper function objects ----------------------------------------

# Originally part of the main wrapper function stated outright
# remove.EnhStocks <- TRUE
# WA <- read.csv("DataIn/WatershedArea.csv") # PRIVATE
WAin <- c("DataIn/WCVIStocks.csv") 

#### New Wrapper function ------------------------------------------------------

# Testing WCVIStocks_NoAgg.csv
# Testing WCVIStocks_NoInlet.csv
# Original WCVIStocks.csv

IWAM_func <- function(WAin = "DataIn/WCVIStocks_NoAgg.csv", # insert Watershed areas file location within the base repository
                        # NEEDS RENAMING to _____
                      remove.EnhStocks = TRUE,
                      run.bootstraps = TRUE, # to turn on or off the bootstrap function added at the end
                      bs_seed = 1, # seed for bootstrapping
                      bs_nBS = 10, # trials for bootstrapping
                      plot = FALSE, # whether or not to create plots stored in DataOut/
                      est.table = TRUE # store kable tables as per wcvi_workedexample.RMD
)
  {
  
  #### 1. Read in data -------------------------------------------------
  
  # Our data includes: Stock name, stock number, year, spawners, recruits, 
  # stream num., and year num.
  # NA's are present in this sample data set and will be removed in the 
  # following sections.
  srdatwna <- read.csv(here::here("DataIn/SRinputfile.csv")) # PRIVATE
  
  # Added watersheds # PRIVATE
  # WA <- read.csv(here::here(WA)) # PRIVATE
  WA <- read.csv(here::here("DataIn/WatershedArea.csv"))
  
  # Add WCVIstocks # PUBLIC
  wcvistocks <- read.csv(here::here(WAin))
    # change wcvistocks to WAin (for now)
  # wcvistocks <- read.csv(here::here("DataIn/WCVIStocks.csv"))
  # sn <- read.csv(here::here("DataIn/WCVIStocks.csv")) # overwritten by wcvistocks
  
  #Naming issue:
    # All columns in the above base files are capcase
    # Insert lines below to convert all column names into lowercase
  
  # * Data Removals and Cleaning ----
  # First, remove any unused stocks using filter()
  # For e.g., two stocks not used in Parken et al, and not documented in Liermann
  srdatwna <- srdatwna %>% filter(Name != "Hoko" & Name != "Hoh") 
  
  
  # Determine Which stocks have NAs? Below filter returns only stock numbers.
  stockwna <- srdatwna %>% filter (is.na(Rec) == TRUE) %>% 
    dplyr::select (Stocknumber) %>%  unique() %>% unlist() 
  
  # Remove years with NAs
  srdat <- srdatwna %>% filter(Rec != "NA") 
  
  # Revise yr_num list where NAs have been removed to be continuous
  # Create a test df to check the order of stock numbers by yr_num
    # Within the subset of stocks with NA's identified earlier as stocks 20 and 21
    # test_1 is not a required object for the model. It is only for checking
  order_test_1 <- srdat %>% filter(Stocknumber == stockwna[1] | 
                               Stocknumber == stockwna[2])
  
  # if/for loop to adjust main df (srdat) to have a continuous year list
  if(max(srdat$Stocknumber) >= stockwna[1]) { # if the max stock number (24)
      # is greater or equal then the stock's identifed (20), then
    for (i in 1:length(stockwna)) { # for  stocks identified with NAs (2)
      len <- length (srdat[which (srdat$Stocknumber == stockwna[i]), ]$yr_num) - 1
        # Create a single value object based on the length of:
        # the number of year's - 1
      srdat [which (srdat$Stocknumber == stockwna[i]), ]$yr_num <- c (0:len)
        # re-write the year numbering for  srdat for the selected stock for a new
        # total length calculated in the line before
    }
  }
  
  # Check for the correction to yr_num - wanted to have the time series - 
  # consistent = re-indexing the stocks - so that there are 
  # no gaps in the. E.g. 0, 1, 2, 3, remove 2 - 0, 1, 3 (now has a break-point)
  # test_2 is not a required object for the model
  order_test_2 <- srdat %>% filter(Stocknumber == stockwna[1] | 
                               Stocknumber == stockwna[2])
  
  # **Future update: mask or simulate NAN's in future - COSEWIC example
  
  # At this point in the function:
  # - Data is read
  # - Undesired stocks removed
  # - Data has a continuous year list
  
  # * Scale Calculation ----------------------------------------------------------
  # Desired scale: 1000 - 0.1 to 100 - responsible for scaling the spawners
  # Points of scaling application and removal:
    # - This scaling is APPLIED in section 2. Create data and parameter lists
    # - This scaling is REMOVED when plotting within the plot functions defined
    # in the file PlotSR.R
    # - This scaling is REMOVED for the calculation of predicted values, R2,
    # and Standard Residuals
  
  # digit_scaling() is now a function within helperFunctions.R
  # Calculate scale for each stock as a tibble (tidyverse df)
  srdat <- digit_scaling(srdat)
  
  # What is the scale of S, R, SMSY, and SREP data,
  # Produces df with two columns: stock number, and scale
  srdat_scale <- srdat %>% dplyr::select(Stocknumber, scale) %>% distinct()
  # Creates the obj. srdat_scale into a vector of only the scales 
  srdat_scale <- srdat_scale$scale 
  
  # Remove years 1981-1984, 1986-1987  from Cowichan (Stocknumber 23) as per 
    # Tompkins et al. 2005
  srdat_cow <- srdat %>% filter(Name == "Cowichan" & 
                                  Yr >= 1985 & 
                                  Yr !=1986 & 
                                  Yr != 1987) # length = 10
  n_cow <- length(srdat_cow$Yr)
  srdat_cow$yr_num <- 0:(n_cow-1)
  srdat <- srdat %>%  filter(Name != "Cowichan") %>% bind_rows(srdat_cow) %>%
    arrange(Stocknumber)
  
  # * Watershed area data and life-history type --------------------------
    # (stream vs ocean)
  # Create a df of names and corresponding stock numbers to use in joining
  names <- srdat %>% dplyr::select (Stocknumber, Name) %>% distinct()
  
  WA <- WA %>% full_join(names, by="Name") %>% arrange(Stocknumber)
  
  # rename stream to --> "lh"
  stream <- srdat %>% dplyr::select(Stocknumber, Name, Stream) %>% 
    group_by(Stocknumber) %>% 
    summarize(lh=max(Stream)) %>% 
    arrange (Stocknumber)
  
  
  #### 2. Create data and parameter lists for TMB --------------------------------
  
  #### * DATA ####
  # Data list for TMB DATA and PARAMETER list - labelled as matches
    # *TOR*: Re-ordered to match TMB input organization
  data <- list()
  
  scale_TMB <- srdat$scale # scale enters the TMB data as: scale
  data$S <- srdat$Sp/scale_TMB # Spawners / scale 
  data$logRS <- log( (srdat$Rec/srdat$scale) / (srdat$Sp/srdat$scale) )
  # logged: scaled recruits / scaled spawners
  data$stk <- as.numeric(srdat$Stocknumber) # stock number
  data$yr <- srdat$yr_num
  N_Stocks <- length(unique(srdat$Name))
  
  data$logMuA_stream_mean <- 1.5 
  data$logMuA_stream_sig <- 2
  data$logMuA_ocean_mean <- 0 #1.5
  data$logMuA_ocean_sig <- 2
  data$HalfNormMean <- 0 #TMB_Inputs$Tau_sigma
  data$HalfNormSig <- 1 #TMB_Inputs$Tau_sigma
  data$HalfNormMeanA <- 0 #0.44 #TMB_Inputs$Tau_sigma
  data$HalfNormSigA <- 1 #0.5 #TMB_Inputs$Tau_sigma
  
  # Read in watershed area data and life-history type and scale
  data$WA <- WA$WA
  data$stream <- stream$lh
  data$scale <- srdat_scale # Ordered by Stocknumber
  
  data$SigRicPriorNorm <- as.numeric(F)
  data$SigRicPriorGamma <- as.numeric(T)
  data$SigRicPriorCauchy <- as.numeric(F)
  data$biasCor <- as.numeric(TRUE)
  data$SigDeltaPriorNorm <- as.numeric(F)
  data$SigDeltaPriorGamma <- as.numeric(T)
  data$SigDeltaPriorCauchy <- as.numeric(F)
  data$Tau_dist <- 0.1
  data$Tau_D_dist <- 1
  # logDeltaSigma # currently listed as param in R, but data_scalar in TMB
  # logNuSigma # currently listed as param in R, but data_scalar in TMB
  
  data$SigDelta_mean <- 0.80 # See KFrun.R, #For half-normal use N(0,1)
  data$SigDelta_sig <- 0.28 # See KFrun.R,
  data$SigNu_mean <- 0.84 # See KFrun.R,
  data$SigNu_sig <- 0.275 # See KFrun.R,
  
  # Read in log(watershed area) for additional stocks
  # Predicted lnWA for plottig CIs:
  data$pred_lnWA <- seq(min(log(WA$WA)), max(log(WA$WA)), 0.1)
  # TestlnWAo
  
  data$target_lnWA_ocean <- wcvistocks %>% # PUBLIC
    mutate (lnWA=log(WA)) %>%
    filter(lh==1) %>% 
    pull(lnWA)
  
  data$target_lnWA_stream <- wcvistocks %>% # PUBLIC
    mutate (lnWA=log(WA)) %>%
    filter(lh==0) %>% # not present in wcvistocks
    pull(lnWA)
  
  
  # Add aggregated wa_stream at inlet level
    # Not part of main function, not generic enough
    # Want a list of WA aggregation

  if (exists("Inlet", where = wcvistocks)){
# Make sure that this runs if BLANK
  InletlnWA <- data.frame(wcvistocks) %>%
    filter(Stock != "Cypre") %>% 
    group_by(Inlet) %>%
    summarize(InletlnWA = log(sum(WA))) %>% 
    filter(Inlet != "San Juan") %>%
    filter(Inlet !="Nitinat")
  
  InletlnWAnoEnh <- data.frame(wcvistocks) %>% 
    filter(Stock != "Cypre") %>% filter(Enh==0) %>%
    group_by(Inlet) %>% 
    summarize(InletlnWA = log(sum(WA))) %>% 
    filter(Inlet != "San Juan") %>%
    filter(Inlet !="Nitinat")
  }
  
  if (exists("CU", where = wcvistocks)){
# Caveat: CU's might be blank - meaning no aggregation
  CUlnWA <- data.frame(wcvistocks) %>% 
    filter(Stock != "Cypre") %>% 
    group_by(CU) %>%
    summarize(CUlnWA = log(sum(WA)))
  
  CUlnWAnoEnh <- data.frame(wcvistocks) %>% 
    filter(Stock != "Cypre") %>% 
    filter(Enh==0) %>%
    group_by(CU) %>% 
    summarize(CUlnWA = log(sum(WA)))
  }
  # Remove aggregation of populations into inlets
    # This code WILL BE REMOVED FROM MAIN REPOSITORY
    # Function: set watershed areas at various spatial scales included
  # remove.EnhStocks <- TRUE
  # To remain in main function
  # if(remove.EnhStocks) data$target_lnWA_ocean <- c(data$target_lnWA_ocean, 
  #                                                  InletlnWAnoEnh$InletlnWA,
  #                                                  CUlnWAnoEnh$CUlnWA)
  #  
  # if(!remove.EnhStocks) data$target_lnWA_ocean <- c(data$target_lnWA_ocean, 
  #                                                   InletlnWA$InletlnWA,
  #                                                   CUlnWA$CUlnWA)
  
  # new statement with the following conditions
    # IF CU and Inlet (above two statements)
    # IF ONLY CU ()
    # IF NOTHING
  # This will most likely work well as a IF, ELIF, ELSE statement
  
  if (all(sapply(c("Inlet","CU"), function(col) exists(col, where = wcvistocks)))) { # Complete aggregation
    if(remove.EnhStocks) data$target_lnWA_ocean <- c(data$target_lnWA_ocean, 
                                                     InletlnWAnoEnh$InletlnWA,
                                                     CUlnWAnoEnh$CUlnWA)
    if(!remove.EnhStocks) data$target_lnWA_ocean <- c(data$target_lnWA_ocean, 
                                                      InletlnWA$InletlnWA,
                                                      CUlnWA$CUlnWA)
    # (!exists("Inlet", where = wcvistocks)) 
  } else if (all(sapply("CU", function(col) exists(col, where = wcvistocks)) &
    !sapply("Inlet", function(col) exists(col, where = wcvistocks, inherits = FALSE)))) { # Just CU's
    if(remove.EnhStocks) data$target_lnWA_ocean <- c(data$target_lnWA_ocean, 
                                                     # InletlnWAnoEnh$InletlnWA,
                                                     CUlnWAnoEnh$CUlnWA)
    if(!remove.EnhStocks) data$target_lnWA_ocean <- c(data$target_lnWA_ocean, 
                                                      # InletlnWA$InletlnWA,
                                                      CUlnWA$CUlnWA)
  } else { # Just overwriting at this point - this is overkill
    if(remove.EnhStocks) data$target_lnWA_ocean <- c(data$target_lnWA_ocean # , 
                                                     # InletlnWAnoEnh$InletlnWA,
                                                     # CUlnWAnoEnh$CUlnWA
                                                     )
    if(!remove.EnhStocks) data$target_lnWA_ocean <- c(data$target_lnWA_ocean # , 
                                                      # InletlnWA$InletlnWA,
                                                      # CUlnWA$CUlnWA
                                                      )
  }
  
  # data$target_lnWA_ocean
  
  #### * PARAMETERS ####
  param <- list()
  
  # Parameters for stocks without AR1
  param$logA <- ( srdat %>% group_by (Stocknumber) %>% 
                        summarise(yi = lm(log( Rec / Sp) ~ Sp )$coef[1] ) )$yi
    # srdat_std: Rec and Sp are not scaled
  
  B <- srdat %>% group_by(Stocknumber) %>% 
    summarise( m = - lm(log( Rec / Sp) ~ Sp )$coef[2] )
    # *Tor* why the negative here?
  param$logB <- log ( 1/ ( (1/B$m)/data$scale ))
    # *Carrie* Need to apply the scale to the inverse of Beta, and then re-invert 
    # and log it. This way the initial parameter is log(scaled beta)
    # logB is scaled
    # in the TMB Ricker model - S is scaled
  param$logSigma <- rep(-2, N_Stocks)
  
  param$logMuA_stream <- 1.5
  param$logSigmaA <- -2
  param$logMuA_ocean <- 0
  
  ## Liermann model
  param$logDelta1 <- 3 
  param$logDelta1_ocean <- 0
  param$logDelta2 <- log(0.72) 
  param$Delta2_ocean <- 0 
  param$logDeltaSigma <- -0.412 # from Parken et al. 2006 where sig=0.662
  
  param$logNu1 <- 3
  param$logNu1_ocean <- 0
  param$logNu2 <- log(0.72)
  param$Nu2_ocean <- 0
  param$logNuSigma <- -0.412 #from Parken et al. 2006 where sig=0.66
  
  
  # 3. Estimate SR parameters from synoptic data set and SMSY and SREPs ----------
  mod <- "IWAM_Liermann" 
  # Compile model if changed:
  # dyn.unload(dynlib(paste("TMB_Files/", mod, sep="")))
  compile(here::here(paste("TMB_Files/", mod, ".cpp", sep="")))
    # Needs to be run to re-create the .dll and .o files from a new .cpp file
  dyn.load(dynlib(here::here(paste("TMB_Files/", mod, sep=""))))
  
  obj <- MakeADFun(data, param, DLL=mod, silent=TRUE, random = c("logA"))
  
  upper <- unlist(obj$par)
  upper[1:length(upper)]<- Inf
  
  lower <- unlist(obj$par)
  lower[1:length(lower)]<- -Inf
  
  
  #### RUN THE MODEL ---------------------------------------------------------
  # Required objects/inputs
    # obj created from MakeADFun function (TMB) that requires:
      # data,
      # parameters,
      # associated TMB file name (dll)
      # misc. information e.g. tracing, random effects parameters 
  
  opt <- nlminb(obj$par, 
                obj$fn, 
                obj$gr, 
                control = list(eval.max = 1e5, iter.max = 1e5), 
                lower=lower, 
                upper=upper)
  
  pl <- obj$env$parList(opt$par) # Gives the parameter estimates from the model
  #summary(sdreport(obj), p.value=TRUE)
  
  
  #### 4. Compile model outputs --------------------------------------------------
  # Create Table of outputs
  all_pars <- data.frame(summary(sdreport(obj)))
  all_pars$Param <- row.names(all_pars)
  # Rename parameter names
  all_pars$Param <- sapply(all_pars$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
  
  
  pars <- data.frame()
  pars <- all_pars %>% filter (Param %in% c("logA", 
                                            "logB", 
                                            "logSigma",  
                                            "SMSY", 
                                            "SREP"))
  
  stnum <- unique(srdat[, c("Stocknumber")])
  pars$Stocknumber <- rep(stnum)
  pars <- left_join(pars, unique(srdat[, c("Stocknumber", "Name")]))
  
  logDeltaSigma <- all_pars %>% filter (Param %in% c("logDeltaSigma")) 
  DeltaSigmaUCL <- exp(logDeltaSigma$Estimate + logDeltaSigma$Std..Error*1.96)
  DeltaSigmaLCL <- exp(logDeltaSigma$Estimate - logDeltaSigma$Std..Error*1.96) 
  DeltaSigma <- exp(logDeltaSigma$Estimate)
  
  # Combine again and rename
  # pars <- All_Est <- All_Ests_std
  pars$Param <- sapply(pars$Param, function(x) (unlist(strsplit(x, "[_]"))[[1]]))
  pars <- pars %>% left_join(stream, by="Stocknumber")
  
  
  all_Deltas <- data.frame()
  all_Deltas <- all_pars %>% filter (Param %in% c("logDelta1", 
                                                  "logDelta2",
                                                  "sigma_delta", 
                                                  "Delta2_bounded", 
                                                  "logDelta1_ocean", 
                                                  #"logDelta2ocean", # does not exist
                                                  "Delta2_ocean", 
                                                  "logNu1", 
                                                  "logNu2", 
                                                  "sigma_nu", 
                                                  "logNu1_ocean", 
                                                  "Nu2_ocean"))
  
  #### 5. Calculate diagnostics and plot SR curves, etc. -------------------------
  
  # Calculate AIC
    # No RE-SCALING
  nLL <- data.frame(nLL=obj$report()$nLL) %>% 
    add_column(Stocknumber=srdat$Stocknumber) %>% group_by(Stocknumber) %>% 
    summarize(CnLL=sum(nLL))
  AIC <- nLL %>% mutate(aic = 2 * 3 + 2*CnLL) 
  
  
  # Get predicted values and calculate r2
  pred_RS <- data.frame()
  pred_RS <- all_pars %>% filter (Param %in% c("logRS_pred"))
  
  # * all_pred previously Preds
  all_pred <- srdat %>% dplyr::select("Stocknumber",
                                   "yr_num", 
                                   "Sp", 
                                   "Rec", 
                                   "scale", 
                                   "Name") %>% 
              add_column(Pred=pred_RS$Estimate)
  
  # Check that length of vector of pred_RS is same as the length of spawner abundances in TMB data element
    # for spawners
  if (length(pred_RS$Estimate) == length(data$S)) {
    print("Lengths checked passed.")
  } else {
    print("WARNING: The output and inputs are not the same length.")
  }
  
  # mutate the predicted values with scale 
    # RE-SCALED VALUES
    # These Preds_stds are not used for plotting
  all_pred <- all_pred %>% mutate(ObslogRS = log ( (Rec / scale) / (Sp/scale) ) )
  r2 <- all_pred %>% group_by(Stocknumber) %>% summarize(r2=cor(ObslogRS,Pred)^2)
  
  
  # Get predicted values and their SEs to plot CIs
    # *These are not re-scaled*
    # They are used in the plotting functions and scaled within
    # pred_lnSMSY_S and pred_lnSMSY_O don't occur anywhere ************************************************************
  pred_lnSMSY <- data.frame() 
  pred_lnSMSY <- all_pars %>% filter (Param %in% c("pred_lnSMSY_stream", 
                                                  "pred_lnSMSY_ocean", 
                                                  "pred_lnSMSY_CI", 
                                                  "pred_lnSMSY_stream_CI", 
                                                  "pred_lnSMSY_ocean_CI"))
    # pred_lnSREP_S and pred_lnSREP_O don't occur elsewhere ***********************************************************
  pred_lnSREP <- data.frame() 
  pred_lnSREP <- all_pars %>% filter (Param %in% c("pred_lnSREP_stream", 
                                                  "pred_lnSREP_ocean", 
                                                  "pred_lnSREP_CI", 
                                                  "pred_lnSREP_stream_CI", 
                                                  "pred_lnSREP_ocean_CI"))
  
  # Calculate standardized residuals
    # These are RE-SCALED values
  SRes <- all_pred %>% mutate ( Res = ObslogRS- Pred) #%>% mutate (StdRes = Res/??)
  sigma <- pars %>% filter(Param=="logSigma") %>% dplyr::select(Stocknumber, Estimate, Name)
  SRes <- SRes %>% left_join(sigma) %>% rename(logSig = Estimate)
    # Slight problem naming with Estimate - can produce clang errors due to overlap
    # with function name.
  SRes <- SRes %>% mutate (StdRes = Res/exp(logSig))
  
  #### * Plot SR Curves ----------------------------------------------------------
  # Plot SR curves. linearized model, standardized residuals, autocorrleation plots for synoptic data set
  # if using a Liermann model, use srdat=srdat_std; otherwise srdat=srdat
  # plot <- TRUE
  # Plotted values are RE-SCALED either by plotting function or are already
    # scaled e.g., "SRes"
  
  # PlotSRCurve(srdat=srdat, pars=pars, r2=r2, removeSkagit = FALSE, mod=mod)
  
  if (plot==TRUE){
    png(paste("DataOut/SR_", mod, ".png", sep=""), width=7, height=7, units="in", res=500)
    PlotSRCurve(srdat=srdat, pars=pars, r2=r2, removeSkagit = FALSE, mod=mod)
    dev.off()
    png(paste("DataOut/SRLin_", mod, ".png", sep=""), width=7, height=7, units="in", res=1000)
    PlotSRLinear(srdat=srdat, pars=pars, r2=r2, removeSkagit = FALSE)
    dev.off()
    png(paste("DataOut/StdResid_", mod, ".png", sep=""), width=7, height=7, units="in", res=1000)
    PlotStdResid(SRes)
    dev.off()
    png(paste("DataOut/ACF_", mod, ".png", sep=""), width=7, height=7, units="in", res=1000)
    Plotacf(SRes)
    dev.off()
  }
  
  #### * Plot WA Regression ------------------------------------------------------
  # Plotted values are RE-SCALED within plot func()
  
  if(plot==TRUE){
    png(paste("DataOut/WAregSMSY_", mod, "_wBC.png", sep=""), width=7, height=7, units="in", res=500)
    par(mfrow=c(1,1), mar=c(4, 4, 4, 2) + 0.1)
    title_plot <- "Prior Ricker sigma and prior WA regression sigma"
    #title_plot <- "Separate life-histories: n=17\nFixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
    plotWAregressionSMSY (pars, all_Deltas, srdat, stream, WA, pred_lnSMSY, 
                          pred_lnWA = data$pred_lnWA, title1=title_plot, mod)
    dev.off()
    
    png(paste("DataOut/WAregSREP_", mod, "_wBC.png", sep=""), width=7, height=7, units="in", res=500)
    #png(paste("DataOut/WAreg_Liermann_SepRicA_UniformSigmaAPrior.png", sep=""), width=7, height=7, units="in", res=500)
    par(mfrow=c(1,1), mar=c(4, 4, 4, 2) + 0.1)
    title_plot <- "Prior Ricker sigmas and prior on WA regression sigma"
    #title_plot <- "Separate life-histories: n=17\nFixed-effect yi (logDelta1), \nFixed-effect slope (Delta2)"
    plotWAregressionSREP (pars, all_Deltas, srdat, stream, WA, pred_lnSREP, 
                          pred_lnWA = data$pred_lnWA, title1=title_plot, mod)
    dev.off()
    #plotWAregression (pars, all_Deltas, srdat, stream, WA, pred_lnSMSY, pred_lnWA = data$pred_lnWA, 
    # title1="Common, fixed yi (logDelta1), \nRandom slope (Delta2)")
  }
  
  #### 6. Calculate prediction intervals for SMSY and SREP for additional stocks ----
  
  # Get predicted values to estimate prediction intervals
    # These values are RE-SCALED to raw estimates during outputting in the TMB code
  pred_lnSMSY_pi <- data.frame()
  pred_lnSMSY_pi <- all_pars %>% filter (Param %in% c("pred_lnSMSY", "lnSMSY"))
  pred_lnSREP_pi <- data.frame()
  pred_lnSREP_pi <- all_pars %>% filter (Param %in% c("pred_lnSREP", "lnSREP"))
  
  pred_lnSMSY_pi$Stocknumber <- rep(stnum)
  pred_lnSREP_pi$Stocknumber <- rep(stnum)
  
  # To calculate prediction intervals, first get predicted and observed logSMSY 
  # and logSREP values for synoptic data set
  #   (actually only need observed logSMSY and logSREP values)
  
  #  First need to get the scale for each stock
  scale_pi <- srdat %>% dplyr::select(Stocknumber, scale) %>% distinct()
  pred_lnSMSY_pi <- pred_lnSMSY_pi %>% left_join(unique(srdat[, c("Stocknumber", "Name")])) %>% 
    left_join(scale_pi)
  pred_lnSREP_pi <- pred_lnSREP_pi %>% left_join(unique(srdat[, c("Stocknumber", "Name")])) %>% 
    left_join(scale_pi)
  
  # Then need to separate observed stream vs ocean type
  
  # pred_lSMSY_stream = predicted log SMSY for stream type
  # pred_lSMSY_ocean = predicted log SMSY for ocean type
  # obs_lSMSY_stream = observed log SMSY for stream type
  # obs_lSMSY_ocean = observed log SMSY for ocean type
  # and same for SREP
  
  # consider removing "_" between to match the names
  pred_lSMSY_stream <- pred_lnSMSY_pi %>% filter(Param=="pred_lnSMSY") %>% left_join(stream) %>% 
    filter(lh == 0) %>% pull(Estimate) #Predicted lnSMSY from WA regression- stream - PlSMSYs
  pred_lSMSY_ocean <- pred_lnSMSY_pi %>% filter(Param=="pred_lnSMSY") %>% left_join(stream) %>% 
    filter(lh == 1) %>% pull(Estimate) #Predicted lnSMSY from WA regression- ocean
  obs_lSMSY_stream <- pred_lnSMSY_pi %>% filter(Param=="lnSMSY") %>% left_join(stream) %>% 
    filter( lh== 0) %>% pull(Estimate) # "observed" lnSMSY data output from SR models- stream
  obs_lSMSY_ocean <- pred_lnSMSY_pi %>% filter(Param=="lnSMSY") %>% left_join(stream) %>% 
    filter(lh == 1) %>% pull(Estimate) # "observed" lnSMSY data output from SR models- ocean
  
  pred_SREP_stream <- pred_lnSREP_pi %>% filter(Param=="pred_lnSREP") %>% left_join(stream) %>% 
    filter(lh == 0) %>% pull(Estimate) #Predicted lnSMSY from WA regression- stream
  pred_SREP_ocean <- pred_lnSREP_pi %>% filter(Param=="pred_lnSREP") %>% left_join(stream) %>% 
    filter(lh == 1) %>% pull(Estimate) #Predicted lnSMSY from WA regression- ocean
  obs_SREP_stream <- pred_lnSREP_pi %>% filter(Param=="lnSREP") %>% left_join(stream) %>% 
    filter(lh == 0) %>% pull(Estimate) # "observed" lnSMSY data output from SR models- stream
  obs_SREP_ocean <- pred_lnSREP_pi %>% filter(Param=="lnSREP") %>% left_join(stream) %>% 
    filter(lh == 1) %>% pull(Estimate) # "observed" lnSMSY data output from SR models- ocean
  
  
  # Get watershed areas for synoptic data set to calculate PIs for stream and ocean
  wa_stream <- WA %>% left_join(stream) %>% filter(lh == 0) %>% pull(WA)
  wa_ocean <- WA %>% left_join(stream) %>% filter(lh == 1) %>% pull(WA)
  
  # Get names of WCVI stocks
    # using wcvistocks instead of sn now
    # * Requires Inlet aggregation information **************************************************************************
  # Create another triple if, elif, else statment
  if (all(sapply(c("Inlet","CU"), function(col) exists(col, where = wcvistocks)))) { # Complete aggregation
    stocknames <- c(as.vector(wcvistocks$Stock), as.vector(InletlnWA$Inlet), as.vector(CUlnWA$CU))
  } else if (all(sapply("CU", function(col) exists(col, where = wcvistocks)) &
             !sapply("Inlet", function(col) exists(col, where = wcvistocks, inherits = FALSE)))) { # Just CU's
    stocknames <- c(as.vector(wcvistocks$Stock), as.vector(CUlnWA$CU))
  } else { # Just overwriting at this point - this is overkill
    stocknames <- c(as.vector(wcvistocks$Stock))
  }
  # stocknames <- c(as.vector(wcvistocks$Stock), as.vector(InletlnWA$Inlet), as.vector(CUlnWA$CU)) # Original line
    # *******************************************************************************************************************
  
  # Get Predicted SMSY and SREP values for new "test/target" WCVI stocks and their Prediction Intervals
  # For single life-history events (stream OR ocean targets)
  # targetSMSY <- data.frame() 
  # targetSREP <- data.frame()
  
  # For instances of both life histories (stream AND ocean targets)
  target_SMSY_ocean <- data.frame()
  target_SREP_ocean <- data.frame()
  
  target_SMSY_stream <- data.frame()
  target_SREP_stream <- data.frame()
  
  condition_target <- cbind("ocean","stream") # choose if you want ocean, stream, or both
  
  if ("ocean" %in% condition_target) {
    # Step 1
    target_SMSY_ocean <- all_pars %>% filter (Param %in% c("target_lnSMSY_ocean")) %>% add_column(Stock = stocknames)
    target_SREP_ocean <- all_pars %>% filter (Param %in% c("target_lnSREP_ocean")) %>% add_column(Stock = stocknames)
    target_SMSY_pull_ocean <- target_SMSY_ocean %>% pull(Estimate)
    target_SREP_pull_ocean <- target_SREP_ocean %>% pull(Estimate)
    # Step 2
    target_SMSY_pi_ocean <- PredInt(x = log(wa_ocean), y = obs_lSMSY_ocean, Predy = target_SMSY_pull_ocean, Newx = data$target_lnWA_ocean)
    target_SREP_pi_ocean <- PredInt(x = log(wa_ocean), y = obs_SREP_ocean, Predy = target_SREP_pull_ocean, Newx = data$target_lnWA_ocean)
    # Step 3
    target_SMSY_ocean <- target_SMSY_ocean %>% add_column(LL = exp(target_SMSY_pi_ocean$lwr), UL = exp(target_SMSY_pi_ocean$upr))
    target_SREP_ocean <- target_SREP_ocean %>% add_column(LL = exp(target_SREP_pi_ocean$lwr), UL = exp(target_SREP_pi_ocean$upr))
    # Step 4
    target_SMSY_ocean <- target_SMSY_ocean %>% mutate (Estimate = exp(Estimate)) %>% dplyr::select(-Std..Error, - Param) %>% 
      add_column(Param = "SMSY")
    target_SREP_ocean <- target_SREP_ocean %>% mutate (Estimate = exp(Estimate)) %>% dplyr::select(-Std..Error, - Param) %>% 
      add_column(Param = "SREP")
    # Step 5
    target_estimates_SMSY_ocean <- target_SMSY_ocean %>% mutate(Estimate = round(Estimate, 0), LL = round(LL,0), UL = round(UL,0))
    target_estimates_SREP_ocean <- target_SREP_ocean %>% mutate(Estimate = round(Estimate, 0), LL = round(LL,0), UL = round(UL,0))
  }
  
  # WARNING
  length_check_params <- all_pars %>% filter (Param %in% c("target_lnSMSY_ocean"))
  if (length(stocknames) == length(length_check_params$Estimate)) {
    print("Lengths checked passed.")
  } else {
    print("WARNING: The output and inputs are not the same length.")
  }
  
  # Does not currently work * no data for stream
  # if ("stream" %in% condition_target){
  #   # Step 1
  #   target_SMSY_stream <- all_pars %>% filter (Param %in% c("target_lnSMSY_stream")) %>% add_column(Stock = stocknames)
  #   target_SREP_stream <- all_pars %>% filter (Param %in% c("target_lnSREP_stream")) %>% add_column(Stock = stocknames)
  #   target_SMSY_pull_stream <- target_SMSY_stream %>% pull(Estimate)
  #   target_SREP_pull_stream <- target_SREP_stream %>% pull(Estimate)
  #   # Step 2
  #   target_SMSY_pi_stream <- PredInt(x = log(wa_stream), y = obs_lSMSY_stream, Predy = target_SMSY_pull_stream, Newx = data$target_lnWA_stream)
  #   target_SREP_pi_stream <- PredInt(x = log(wa_stream), y = obs_SREP_stream, Predy = target_SREP_pull_stream, Newx = data$target_lnWA_stream)
  #   # Step 3
  #   target_SMSY_stream <- target_SMSY_stream %>% add_column(LL = exp(target_SMSY_pi_stream$lwr), UL = exp(target_SMSY_pi_stream$upr))
  #   target_SREP_stream <- target_SREP_stream %>% add_column(LL = exp(target_SREP_pi_stream$lwr), UL = exp(target_SREP_pi_stream$upr))
  #   # Step 4
  #   target_SMSY_stream <- target_SMSY_stream %>% mutate (Estimate = exp(Estimate)) %>% dplyr::select(-Std..Error, - Param) %>% 
  #     add_column(Param = "SMSY")
  #   target_SREP_stream <- target_SREP_stream %>% mutate (Estimate = exp(Estimate)) %>% dplyr::select(-Std..Error, - Param) %>% 
  #     add_column(Param = "SREP")
  #   # Step 5
  #   target_estimates_SMSY_stream <- target_SMSY_stream %>% mutate(Estimate = round(Estimate, 0), LL = round(LL,0), UL = round(UL,0))
  #   target_estimates_SREP_stream <- target_SREP_stream %>% mutate(Estimate = round(Estimate, 0), LL = round(LL,0), UL = round(UL,0))
  # }
  
      # STEPS
  # For inclusion into if loops above
  # Main goal is to allow for stream/ocean AND scenarios
  # Step 1: Filter out targets into new df and extract a single column estimates
  # targetSMSY <- all_pars %>% filter (Param %in% c("target_lnSMSY_ocean")) %>% add_column(Stock = stocknames) # OLD # target ocean example
  # targetSREP <- all_pars %>% filter (Param %in% c("target_lnSREP_ocean")) %>% add_column(Stock = stocknames) # OLD # target ocean example
  # targetSMSYpull <- targetSMSY %>% pull(Estimate) # OLD
  # targetSREPpull <- targetSREP %>% pull(Estimate) # OLD
  
  # Step 2: create new intervals object
  # Use custom function: PredInt() to estimate prediction intervals 
    # PredInt() defined in helperFunctions.R
  # targetSMSY_pi <- PredInt(x = log(wa_ocean), y = obs_lSMSY_ocean, Predy = targetSMSYpull, Newx = data$target_lnWA_ocean) # OLD
  # targetSREP_pi <- PredInt(x = log(wa_ocean), y = obs_SREP_ocean, Predy = targetSREPpull, Newx = data$target_lnWA_ocean) # OLD
  
  # Step 3: add a column with intervals to target dataframes (currently only ONE SCENARIO)
  # exp() bounds
  # targetSMSY <- targetSMSY %>% add_column(LL = exp(targetSMSY_pi$lwr), UL = exp(targetSMSY_pi$upr)) # OLD
  # targetSREP <- targetSREP %>% add_column(LL = exp(targetSREP_pi$lwr), UL = exp(targetSREP_pi$upr)) # OLD
  
  # Step 4: mutate exp() to estimate 
  # targetSMSY <- targetSMSY %>% mutate (Estimate = exp(Estimate)) %>% dplyr::select(-Std..Error, - Param) %>% 
  #   add_column(Param = "SMSY") # OLD
  # targetSREP <- targetSREP %>% mutate (Estimate = exp(Estimate)) %>% dplyr::select(-Std..Error, - Param) %>% 
  #   add_column(Param = "SREP") # OLD
  
  # Step 5: Mutate to round targets and bind together SMSY and SREP
  # target_estimates_SMSY <- targetSMSY %>% mutate(Estimate = round(Estimate, 0), LL = round(LL,0), UL = round(UL,0)) # OLD
  # target_estimates_SREP <- targetSREP %>% mutate(Estimate = round(Estimate, 0), LL = round(LL,0), UL = round(UL,0)) # OLD
  
  # Final combination of SMSY and SREP estimates into final df
  data_out_ocean <- target_estimates_SMSY_ocean %>% bind_rows(target_estimates_SREP_ocean)
  # data_out_stream <- target_estimates_SMSY_stream %>% bind_rows(target_estimates_SREP_stream)
  
  # data_out_combined <- data_out_ocean %>% bind_rows(data_out_stream)
  
  # data_out <- target_estimates_SMSY %>% bind_rows(target_estimates_SREP) # OLD
  
  # Step 6: Final writing step for outputting targets
  # Write SMSY and SREP with PIs to file
  if(remove.EnhStocks) write.csv(data_out_ocean, here::here("DataOut/FUNCTIONTEST_dataout_target_ocean_noEnh.csv"))
  if(remove.EnhStocks) datain <- c("DataOut/FUNCTIONTEST_dataout_target_ocean_noEnh.csv")
  if(!remove.EnhStocks) write.csv(data_out_ocean, here::here("DataOut/dataout_target_ocean_wEnh.csv"))
  if(!remove.EnhStocks) datain <- c("DataOut/dataout_target_ocean_wEnh.csv")
  
  # if(remove.EnhStocks) write.csv(data_out_stream, "DataOut/dataout_target_stream_noEnh.csv")
  # if(!remove.EnhStocks) write.csv(data_out_stream, "DataOut/dataout_target_stream_wEnh.csv")

  # if(remove.EnhStocks) write.csv(data_out, "DataOut/dataout_target_noEnh.csv")
  # if(!remove.EnhStocks) write.csv(data_out, "DataOut/dataout_target_wEnh.csv")
  
  #### End of IWAM Model ####
  
  #### run.bootstrap ####
    # Consider creating function objects for set.seed and for nBS
  # dfout <- NULL
  # dfout <- data.frame()
  
  if (run.bootstraps == TRUE){
    # set.seed(1) #10#12#13 (work for 1000), for 100, 200, 300, (for 5000trials), 1, 2, 3 (for 20000trials)
    set.seed(bs_seed)
    # nBS <- 10 # number trials for bootstrapping (original 20000), for testing use 10
    nBS <- bs_nBS
    outBench <- list()
    
    for (k in 1:nBS) {
      # datain must match the above Step 6 for writing output target estimates
      out <- Get.LRP.bs(datain = datain,# "DataOut/FUNCTIONTEST_dataout_target_ocean_noEnh.csv"
                        Bern_logistic=FALSE, 
                        prod="LifeStageModel", 
                        LOO = NA, 
                        run_logReg=FALSE) 
      outBench[[k]] <- out$bench
    }
    
    # # Is 200 enough trials? Yes
    # running.mean <- cumsum(LRP.bs$fit) / seq_along(LRP.bs$fit) 
    # plot(running.mean)
    
    
    # Compile bootstrapped estimates of Sgen, SMSY, and SREP, and identify 5th and 
    # 95th percentiles
    SGEN.bs <- select(as.data.frame(outBench), starts_with("SGEN"))
    # stockNames <- read.csv("DataOut/WCVI_SMSY_noEnh_wBC.csv") %>% 
    stockNames <- read.csv(here::here(datain)) %>% 
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
    
    write.csv(dfout, here::here("DataOut/1_FUNCTIONTEST_wcviCK-BootstrappedRPs.csv"))

  }
  # print(dfout)
  #### Table outputs ####
  
  # Can you store kables as objects?
    # Code taken from wcvi_workedexamples.Rmd
  # SREP_out <- NULL
  # SGEN_out <- NULL
  # SMSY_out <- NULL
  # 
  if (est.table == TRUE){
    # Output locations
    locations <- c('Barkley' , 'Clayoquot' , 'Kyuquot' , 'Quatsino' , 'Nootka/Esperanza')

    # SREP Store
    SREP_out <- data.frame()
    SREP_out <- dfout %>%
      filter(RP=='SREP') %>%
      filter(Stock %in% locations) %>%
      rename('Lower Quantile'=lwr, 'SREP'=Value, 'Upper Quantile'=upr) %>%
      mutate(RP = NULL) %>%
      relocate(Stock)

    # SGEN Store
    SGEN_out <- data.frame()
    SGEN_out <- dfout %>%
      filter(RP=='SGEN') %>%
      filter(Stock %in% locations) %>%
      rename('Lower Quantile'=lwr, 'SGEN'=Value, 'Upper Quantile'=upr) %>%
      mutate(RP = NULL) %>%
      relocate(Stock)

    # SMSY Store
    SMSY_out <-data.frame()
    SMSY_out <- dfout %>%
      filter(RP=='SMSY') %>%
      filter(Stock %in% locations) %>%
      rename('Lower Quantile'=lwr, 'SMSY'=Value, 'Upper Quantile'=upr) %>%
      mutate(RP = NULL) %>%
      relocate(Stock)

    # SMSY_store <- SMSY_out
    # storedtables <- list(SMSY_store, SREP_store, SGEN_store)
    # return(storedtables)
  }
  
  # Return function 
  return(list(dfout, SREP_out, SGEN_out, SMSY_out)) # WORKING
  #### End of IWAM_func ####
}

# # Check that the function runs: 
# store <- IWAM_func(WAin = "DataIn/WCVIStocks_NoAgg.csv", # insert Watershed areas file location within the base repository
#                    run.bootstraps = TRUE, # to turn on or off the bootstrap function added at the end
#                    bs_seed = 1, # seed for bootstrapping
#                    bs_nBS = 10, # trials for bootstrapping
#                    # mod = "IWAM_Liermann", # TMB model name for .cpp
#                    remove.EnhStocks = TRUE,
#                    plot = FALSE, # whether or not to create plots stored in DataOut/
#                    est.table = TRUE # store kable tables as per wcvi_workedexample.RMD
#                    )
# 
