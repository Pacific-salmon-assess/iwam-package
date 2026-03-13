here::i_am("R/LambertWs.R") # For non-RStudio functionality
source(here::here("R/LambertWs.R")) # Lambert W function
source(here::here("R/helperFunctions.R")) # For bootstrapping
source(here::here("R/derived_post.R")) # For posterior extraction

Parkentable1 <- read.csv(here::here("DataIn/Parken_Table1n2.csv")) # Test stocks e.g. WCVI stocks
srdatwna <- read.csv(here::here("DataIn/SRinputfile.csv"))
WAbase <- read.csv(here::here("DataIn/WatershedArea.csv"))

# S <- srdat$Sp[srdat$Stocknumber == 1]
# logRS <- log(srdat$Rec[srdat$Stocknumber == 1] / srdat$Sp[srdat$Stocknumber == 1])
intercepts <- c()
slopes <- c()
sigma <- c()
stocks <- Parkentable1$Stock
Parkentable1$logalpha_new <- NA
Parkentable1$beta_new <- NA
Parkentable1$sigma_new <- NA

for (i in seq_along(stocks)){
	# Create lm(log(R/S) ~ S) for each population (25 in total)
	newdat <- srdat[srdat$Name == stocks[i],]
	# Save the model
	model <- lm(log(Rec/Sp) ~ Sp, data = newdat)
	# Get access to summary(model) and the parameters
	Parkentable1$logalpha_new[i] <- coef(model)[1] #model[[1]][[1]]
	Parkentable1$beta_new[i] <- -coef(model)[2] # model[[1]][[2]]
	Parkentable1$sigma_new[i] <- sigma(model)
}

# Write a code to fit habitat model and loop with simulated new values of alpha and beta

# Assume SMSY/SREP - 