library(tidyverse)
library(nimble)
library(coda)
library(dplyr)
library(MCMCvis)


## This seems like the best LambertsW function so far:
## https://github.com/cran/lamW/blob/master/R/lamW.R
nimLambertsW <- nimbleFunction(
  run = function(x = double()) {
    REXP <- 2.718281828459045090795598298427648842334747314453125 ## exp(1)
    REXPI <- 0.367879441171442334024277442949824035167694091796875 # exp(-1)
    EPS <- 2.2204460492503131e-16
      
    if (x == Inf) {
      return(Inf)
    } else if (x < -REXPI) {
      return(NaN)
    } else if (abs(x + REXPI) <= EPS) {
      return(-1.0)
    } else if (abs(x) <= 1e-16) {
      return(x)
    } else {
      if (abs(x) <= 6.4e-3) {
        ## When this close to 0 the Fritsch iteration may underflow. Instead,
        ## function will use degree-6 minimax polynomial approximation of Halley
        ## iteration-based values. Should be more accurate by three orders of
        ## magnitude than Fritsch's equation (5) in this range.
        ans <- (((((-1.0805085529250425e1 * x + 5.2100070265741278) * x -
               2.6666665063383532) * x + 1.4999999657268301) * x -
               1.0000000000016802) * x + 1.0000000000001752) * x +
               2.6020852139652106e-18       
        ## Minimax Approximation calculated using R package minimaxApprox 0.1.0
        return(ans);

      } else if (x <= REXP) {
        p = sqrt(2.0 * (REXP * x + 1.0));
        Numer = (0.2787037037037037 * p + 0.311111111111111) * p - 1.0;
        Denom = (0.0768518518518518 * p + 0.688888888888889) * p + 1.0;
        w = Numer / Denom;
      } else {
        w = log(x)
        L_2 = log(w);
        L_3 = L_2 / w;
        L_3_sq = L_3 * L_3;
        w <- w + -L_2 + L_3 + 0.5 * L_3_sq - L_3 / w + L_3 / (w * w) - 
             1.5 * L_3_sq / w + L_3_sq * L_3 / 3.0;
      }
      ## Fritsch Iteration for up to 5 iterations.
      MaxEval <- 5
      CONVERGED <- FALSE
      k <- 2.0 / 3.0;
      i <- 0;
      while (!CONVERGED & i < MaxEval){
        z <- log(x / w) - w
        w1 <- w + 1.0
        q <- 2.0 * w1 * (w1 + k * z)
        qmz <- q - z
        e <- z / w1 * qmz / (qmz - z)
        CONVERGED <- abs(e) <= EPS
        w <- w*(1.0 + e)
        i <- i + 1
      }
      return(w)
    }
    returnType(double())
  }
)

# Setup data:
srdatwna <- read.csv(here::here("DataIn/SRinputfile.csv"))
WAbase <- read.csv(here::here("DataIn/WatershedArea.csv"))

srdatwna <- srdatwna %>% 
  filter(!Name %in% c("Hoko","Hoh")) 
  
# Remove years with NAs and renumber.
srdat <- srdatwna %>% 
  filter(!Name %in% c("Hoko","Hoh")) %>% 
  filter(Rec != "NA") %>%
  filter( !(Name == "Cowichan" & (Yr < 1985 | Yr == 1986 | Yr == 1978))) %>%
  group_by(Name, Stocknumber, Stream) %>%
  arrange(Yr) %>%
  mutate(yr_num = 0:(n()-1)) %>%
  ungroup() %>%
  arrange(Stocknumber) %>%
  mutate(lh = factor(ifelse(Stream == 0, "stream", "ocean"), levels = c("stream", "ocean")))

names <- srdat %>% 
  dplyr::select (Stocknumber, Name, lh) %>% 
    distinct()

WAbase <- WAbase %>% 
  full_join(names, by="Name") %>% 
    arrange(Stocknumber) %>%
    mutate(logWA = log(WA)) 

## Shift log WA for the mean.
mean_logWA <- mean(WAbase$logWA)
WAbase$logWAshifted <- WAbase$logWA - mean_logWA

data <- list(
  logRS = log(srdat$Rec) - log(srdat$Sp)  
)

constants <- list()
constants$N_Stk = max(srdat$Stocknumber + 1)
constants$stk = srdat$Stocknumber + 1
constants$N_Obs = nrow(srdat)
constants$S = srdat$Sp
constants$betaPriorMean = c(10,0,0,0)

## Build a linear model:
constants$X <- model.matrix( ~ lh*logWAshifted, data = WAbase)
constants$nbeta <- ncol(constants$X)

inits <- function(){
  list(beta = c(10,0,0,0),
    logAlpha0 = 1.5,
    logESD = 1,
    logAlphaSD = 1)
}

stock_recruit_srep_biasCor <- nimbleCode({
  ## priors
  logAlpha0 ~ dnorm(mean=0.6, sd=0.45)
  # logAlpha0 ~ dnorm(mean=1.5, sd=2.5)
  logAlphaSD ~ dunif(0, 100)
  logESD ~ dunif(0, 100)
  ## Can we create regional mean alpha's?
  ## Can we use the distribution's used by Diane Dobson - see bootstrap

  ## Based on Parken et al. 2006 (log ( 1/beta ) has a linear relationship with log Watershed Area.
  for(i in 1:nbeta) beta[i] ~ dnorm(betaPriorMean[i], tau = 0.001)

  for( i in 1:N_Stk ) { # 25 iterations - number of base stocks
    logAlpha[i] ~ dnorm(logAlpha0, sd = logAlphaSD) ## Random effect for log alpha.
    logE0[i] ~ dnorm(mean = 0, sd = logESD) ## Stock level random effect 

    log(E[i]) <- inprod(beta[1:nbeta], X[i,1:nbeta]) + logE0[i] ## Stock level regression
    tauobs[i] ~ dgamma(0.001, 0.001) # stock specific precision
  }

  ## Model and residuals:  
  for(i in 1:N_Obs){ # 501 iterations - number of base observations
    logRS_pred[i] <- logAlpha[stk[i]]*(1 - S[i]/E[stk[i]])
    logRS[i] ~ dnorm( logRS_pred[i], tau = tauobs[stk[i]]) # stock specific precision
  }
})

model <- nimbleModel(stock_recruit_srep_biasCor, data = data, constants = constants,
  inits = inits(), buildDerivs=FALSE) # Build derivs if you want MLE
conf <- configureMCMC(model)
conf$setMonitors(c("logAlpha", "logAlpha0", "E", "logRS_pred",
                    "beta", "tauobs",
                    "logAlphaSD", "logESD", "logE0"))
mcmc <- buildMCMC(conf)
cmodel <- compileNimble(model)

cmcmc <- compileNimble(mcmc, project = model)
mcmc.out <- runMCMC(cmcmc, niter=50000, nburnin=5000, nchains=3, samplesAsCodaMCMC = TRUE)

MCMCtrace(object = mcmc.out, params = "beta")

## This is good to check how well your model works:
posteriorPredictiveFunction <- nimbleFunction(
  setup = function(model, mcmc){
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)
    cat("Stochastic parents of data are:", paste(parentNodes, collapse = ','), ".\n")
    simNodes <- model$getDependencies(parentNodes, self = FALSE)
    vars <- mcmc$mvSamples$getVarNames()  # need ordering of variables in mvSamples / samples matrix
    varsScalar <- model$expandNodeNames(vars, returnScalarComponents = TRUE)
    cat("Using posterior samples of:", paste(vars, collapse = ','), ".\n")
    n <- length(model$expandNodeNames(dataNodes, returnScalarComponents = TRUE))
    samples <- matrix(0, nrow = 2, ncol = 2)
    nMCMC <- 2
    Smsy <- matrix(0, nrow = 2, ncol = 2)
    Sgen <- matrix(0, nrow = 2, ncol = 2)
    Srep <- matrix(0, nrow = 2, ncol = 2)
  },
  run = function(){},
  methods = list(
    getVarNames = function(){
      returnType(character(1))
      return(varsScalar)
    },
    posteriorPredict = function(mcmcSamples = double(2)){
      samples <<- mcmcSamples
      nSamp <- dim(samples)[1]
      ppSamples <- matrix(nrow = nSamp, ncol = n)
      for(i in 1:nSamp) {
            values(model, vars) <<- samples[i, ]
            model$simulate(simNodes, includeData = TRUE)
            ppSamples[i, ] <- values(model, dataNodes)
      }
      returnType(double(2))
      return(ppSamples)
    },
    saveMCMC = function(mcmcSamples = double(2)){
      samples <<- mcmcSamples   
      nMCMC <<- dim(samples)[1]
    },
    calcRefPoints = function(X = double(2)){
      npred <- dim(X)[1]
      Smsy <<- matrix(0, nrow = nMCMC, ncol = npred)
      Sgen <<- matrix(0, nrow = nMCMC, ncol = npred)
      Srep <<- matrix(0, nrow = nMCMC, ncol = npred)

      for( i in 1:nMCMC ){
      values(model, vars) <<- samples[i, ]
        for( j in 1:npred){
          logalpha <- rnorm(1, mean = model$logAlpha0[1], sd = model$logAlphaSD[1])
          E <- exp(inprod(values(model, "beta"), X[j,]))
          beta <- logalpha/E
          Srep[i,j] <<- E
          Smsy[i,j] <<- (1-nimLambertsW(exp(1-logalpha)))/beta
          Sgen[i,j] <<- -1/beta*nimLambertsW(-beta*Smsy[i,j]/(exp(logalpha)))
        }
      }
    },
    getSmsy = function(){
      returnType(double(2))
      return(Smsy)
    },
    getSgen = function(){
      returnType(double(2))
      return(Sgen)      
    },
    getSrep = function(){
      returnType(double(2))
      return(Srep)
    }
  )
)

predict <- posteriorPredictiveFunction(model, mcmc)
cpredict <- compileNimble(predict)

## Prediction design matrix:
WAin <- read.csv(here::here("DataIn/WCVIStocks.csv"))
WAin$logWA <- log(WAin$WA)
WAin$logWAshifted <- WAin$logWA-mean_logWA
WAin <- WAin %>%
  mutate(lh = factor(ifelse(lh == 0, "stream", "ocean"), levels = c("stream", "ocean")))
Xpred <- model.matrix( ~ lh*logWAshifted, data = WAin)  ## This has to match model input.

## Make sure internal variable names match with output.
getOrder <- cpredict$getVarNames()

## Predict new sites:
Smsy <- Sgen <- Srep <- NULL
for( i in 1:length(mcmc.out)){
  samples <- mcmc.out[[i]][, getOrder] ## If this fails need to check your MCMC tracked vars.
  cpredict$saveMCMC(samples)
  cpredict$calcRefPoints(Xpred)
  Smsy <- rbind(Smsy, cpredict$getSmsy())
  Sgen <- rbind(Sgen, cpredict$getSgen())
  Srep <- rbind(Srep, cpredict$getSrep())
}

sum.smsy <- do.call("cbind", summary(mcmc(Smsy)))
colnames(sum.smsy) <- paste0("Smsy_",colnames(sum.smsy))
sum.sgen <- do.call("cbind", summary(mcmc(Sgen)))
colnames(sum.sgen) <- paste0("Sgen_",colnames(sum.sgen))
sum.srep <- do.call("cbind", summary(mcmc(Srep)))
colnames(sum.srep) <- paste0("Srep_",colnames(sum.srep))

## 
WAin <- WAin %>% 
  bind_cols(sum.smsy) %>% 
  bind_cols(sum.sgen) %>% 
  bind_cols(sum.srep)

## Make a prediction matrix for a line:
WAline <- expand.grid(lh = c("stream", "ocean"), logWA = 2:14)
WAline <- WAline %>%
  mutate(lh = factor(lh, levels = c("stream", "ocean")),
    logWAshifted = logWA - mean_logWA)
Xline <- model.matrix( ~ lh*logWAshifted, data = WAline) ## This has to match model input.
  
## Find the posterior predictions on the line:
Smsy_line <- Sgen_line <- Srep_line <- NULL
for( i in 1:length(mcmc.out)){
  samples <- mcmc.out[[i]]
  samples <- samples[,order(colnames(samples))]
  cpredict$saveMCMC(samples)
  cpredict$calcRefPoints(Xline) ## Design matrix for this prediction.
  Smsy_line <- rbind(Smsy_line, cpredict$getSmsy())
  Sgen_line <- rbind(Sgen_line, cpredict$getSgen())
  Srep_line <- rbind(Srep_line, cpredict$getSrep())
}

sum.smsy.line <- do.call("cbind", summary(mcmc(Smsy_line)))
colnames(sum.smsy.line) <- paste0("Smsy_",colnames(sum.smsy.line))
sum.sgen.line <- do.call("cbind", summary(mcmc(Sgen_line)))
colnames(sum.sgen.line) <- paste0("Sgen_",colnames(sum.sgen.line))
sum.srep.line <- do.call("cbind", summary(mcmc(Srep_line)))
colnames(sum.srep.line) <- paste0("Srep_",colnames(sum.srep.line))

## 
WAline <- WAline %>% 
  bind_cols(sum.smsy.line) %>% 
  bind_cols(sum.sgen.line) %>% 
  bind_cols(sum.srep.line)

## Posterior Predictive Checks:
library(DHARMa)
## Combine random amounts of the multiple chains.
samples <- NULL
for( i in 1:length(mcmc.out) ){
  samples <- rbind(samples, mcmc.out[[i]][sample(nrow(mcmc.out[[1]]), ceiling(0.25*nrow(mcmc.out[[1]]))), getOrder])
}
ppsim <- cpredict$posteriorPredict(samples)
logRS_predMean <- apply(samples[,grep("logRS_pred", colnames(samples)) ], 2, mean)
DHARMaRes <- createDHARMa(simulatedResponse = t(ppsim), observedResponse = data$logRS, 
             fittedPredictedResponse = logRS_predMean, integerResponse = T)
plot(DHARMaRes)

## Now can plot with whatever naming conventions you want.
##########################################################
ggplot(data = WAline, aes(x = logWA, y = log(Sgen_Mean), colour = lh, fill = lh)) +
  geom_line() + 
  geom_ribbon(aes(ymin = log(`Sgen_2.5%`), ymax = log(`Sgen_97.5%`)), alpha = .2) +
  geom_point(data = WAin, mapping = aes(logWA, log(Sgen_Mean), colour = lh)) +
  geom_errorbar(data = WAin, mapping = aes(ymin = log(`Sgen_2.5%`), ymax = log(`Sgen_97.5%`), colour = lh )) +
  theme_classic() +
  ylab("Log Sgen") +
  xlab("Log Watershed Area")
  
ggplot(data = WAline, aes(x = logWA, y = log(Smsy_Mean), colour = lh, fill = lh)) +
  geom_line() + 
  geom_ribbon(aes(ymin = log(`Smsy_2.5%`), ymax = log(`Smsy_97.5%`)), alpha = .2) +
  geom_point(data = WAin, mapping = aes(logWA, log(Smsy_Mean), colour = lh)) +
  geom_errorbar(data = WAin, mapping = aes(ymin = log(`Smsy_2.5%`), ymax = log(`Smsy_97.5%`), colour = lh )) +
  theme_classic() +
  ylab("Log Smsy") +
  xlab("Log Watershed Area")

ggplot(data = WAline, aes(x = logWA, y = log(Srep_Mean), colour = lh, fill = lh)) +
  geom_line() + 
  geom_ribbon(aes(ymin = log(`Srep_2.5%`), ymax = log(`Srep_97.5%`)), alpha = .2) +
  geom_point(data = WAin, mapping = aes(logWA, log(Srep_Mean), colour = lh)) +
  geom_errorbar(data = WAin, mapping = aes(ymin = log(`Srep_2.5%`), ymax = log(`Srep_97.5%`), colour = lh )) +
  theme_classic() +
  ylab("Log Srep") +
  xlab("Log Watershed Area")

