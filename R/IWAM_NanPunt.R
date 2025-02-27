# IWAM Model Run for RPA on Nanaimo and Puntledge stocks ####

# Model steps available in the following files:
  # Model function: IWAM_model.R
  # TMB model: IWAM_Liermann.cpp
  # Bootstrapping method code: Get_LRP_bs.R

# Libraries and Sources ####
source(here::here("R/IWAM_model.R"))

# Function ####
    # 15-30 minutes of run time for 20,000 bootstrap iterations
iwamobj <- IWAM_func(WAinraw = c("DataIn/WCVIStocks_NanPunt.csv"),
                          targetname = "Target_name", # For user input in naming of outputs
                          bs_seed = 1, # Will default to 1
                          bs_nBS = 20000, # Number of bootstrapping iterations
                          bias.cor = TRUE, # Adding the sigma^2/2 bias correction term
                          # Remaining function parameters are default listed in IWAM_model.R
                          prod = "Parken")


# Understanding the outputs: ####

    # The final output file will be named with the following type:
      # "DataOut/Target_name_getLRP-BootstrappedRPs.csv"
    # This can also be accessed by iwamobj$dfout for processing.
    # The filename will always be saved within iwamobj$dataname.

    # The list of other available outputs is available above the return funcion
      # in the IWAM_model code.




