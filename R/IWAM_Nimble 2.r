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
  filter( !(Name == "Cowichan" & (Yr < 1985 | Yr == 1986 | Yr == 1987))) %>%
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
constants$betaPriorMean = c(10,10,0,0)
# constants$betaPriorMean = c(10, 0,0,0)

## Build a linear model:
# constants$X <- model.matrix( ~ -1+lh+lh:logWAshifted, data = WAbase)
constants$X <- model.matrix( ~ lh*logWAshifted, data = WAbase)
constants$nbeta <- ncol(constants$X)



inits <- function(){
  list(beta = c(10,0,0,0),
  # list(beta = c(10,10,0,0),
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
    logRS[i] ~ dnorm(logRS_pred[i], tau = tauobs[stk[i]]) # stock specific precision
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
mcmc.out <- runMCMC(cmcmc, niter=50000, nburnin=5000, nchains=3, samplesAsCodaMCMC=TRUE)

MCMCtrace(object = mcmc.out, params = "beta")
# Can I see the logAlpha traces
MCMCtrace(object = mcmc.out, params = "logAlpha")

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
    Smax <- matrix(0, nrow = 2, ncol = 2)
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
    ## Assumes X is sorted by stk when passed at obsLevel.
    calcRefPoints = function(X = double(2), obsLevel = integer(0, default = 0)){
      npred <- dim(X)[1]
      Smsy <<- matrix(0, nrow = nMCMC, ncol = npred)
      Sgen <<- matrix(0, nrow = nMCMC, ncol = npred)
      Srep <<- matrix(0, nrow = nMCMC, ncol = npred)
      Smax <<- matrix(0, nrow = nMCMC, ncol = npred)
      
      if(obsLevel == 1 & npred != dim(model$logAlpha)[1]){
        stop("Can only predict for the observed locations")
      }
      for( i in 1:nMCMC ){
      values(model, vars) <<- samples[i, ]
        for( j in 1:npred){
          if(obsLevel == 0 ) logalpha <- rnorm(1, mean = model$logAlpha0[1], sd = model$logAlphaSD[1])
          else logalpha <- model$logAlpha[j]
          E <- exp(inprod(values(model, "beta"), X[j,]) + obsLevel*model$logE0[j])
          beta <- logalpha/E
          Smax[i,j] <<- 1/beta
          Srep[i,j] <<- E
          Smsy[i,j] <<- (1-nimLambertsW(exp(1-logalpha)))/beta
          Sgen[i,j] <<- -1/beta*nimLambertsW(-beta*Smsy[i,j]/(exp(logalpha)))
        }
      }
    },
    getRefPostPred = function(type = character(0, default = "Smsy")){
      returnType(double(2))
      if(type == "Smsy") return(Smsy)  
      if(type == "Sgen") return(Sgen)
      if(type == "Smax") return(Smax)
      if(type == "Srep") return(Srep)
      print("Returning Srep")
      return(Srep)
    }
  )
)

predict <- posteriorPredictiveFunction(model, mcmc)
cpredict <- compileNimble(predict)

summarize_all_ref_pts <- function(predInfo){
  summarize <- function(mat){
    res <- data.frame()
    for( i in 1:ncol(mat) ){
      res <- rbind(res, data.frame(Mean = mean(mat[,i]), SD = sd(mat[,i]), 
        LowerCI = quantile(mat[,i], 0.025), Median = quantile(mat[,i], 0.5), UpperCI = quantile(mat[,i], 0.975), 
        row.names = i)
      )
    }
    return(res)
  }
  out <- data.frame()
  for( i in c("Srep", "Smax", "Smsy", "Sgen") ){
    tmp <- summarize(cpredict$getRefPostPred(i))
    tmp$type <- i
    tmp <- cbind(predInfo, tmp)
    out <- rbind(out, tmp)
  }  
  return(out)
}

## Prediction design matrix:
WAin <- read.csv(here::here("DataIn/WCVIStocks.csv"))
WAin$logWA <- log(WAin$WA)
WAin$logWAshifted <- WAin$logWA-mean_logWA
WAin <- WAin %>%
  mutate(lh = factor(ifelse(lh == 0, "stream", "ocean"), levels = c("stream", "ocean")))
Xpred <- model.matrix( ~ lh*logWAshifted, data = WAin)  ## This has to match model input.
# Xpred <- model.matrix( ~ -1+lh+lh:logWAshifted, data = WAin)
## Make sure internal variable names match with output.
getOrder <- cpredict$getVarNames()

## Predict new sites:
samples <- do.call("rbind", mcmc.out)[, getOrder] ## If this fails need to check your MCMC tracked vars.
cpredict$saveMCMC(samples)
cpredict$calcRefPoints(Xpred) ## 45,000 x 5 this is a big operation.
sum.pred <- summarize_all_ref_pts(WAin)

## Make a prediction matrix for a line:
WAline <- expand.grid(lh = c("stream", "ocean"), logWA = 2:14)
WAline <- WAline %>%
  mutate(lh = factor(lh, levels = c("stream", "ocean")),
    logWAshifted = logWA - mean_logWA)
Xline <- model.matrix( ~ lh*logWAshifted, data = WAline) ## This has to match model input.
# Xline <- model.matrix( ~ -1+lh+lh:logWAshifted, data = WAline)
cpredict$calcRefPoints(Xline) ## 45,000 x 5 this is a big operation.
sum.line <- summarize_all_ref_pts(WAline)
  
## Make a prediction matrix for observations
Xbase <- model.matrix( ~ lh*logWAshifted, data = WAbase) ## This has to match model input.
# Xbase <- model.matrix( ~ -1+lh+lh:logWAshifted, data = WAbase)
cpredict$calcRefPoints(Xbase, 1) ## 45,000 x 5 this is a big operation.
sum.base <- summarize_all_ref_pts(WAbase)



## Example plots:
ggplot(data = sum.line %>% filter(type == "Sgen"), 
    aes(x = logWA, y = log(Mean), colour = lh, fill = lh)) +
  geom_line() + 
  geom_ribbon(aes(ymin = log(LowerCI), ymax = log(UpperCI)), alpha = .2) +
  geom_point(data = sum.base %>% filter(type == "Sgen"), 
              mapping = aes(logWA, log(Mean), colour = lh)) +
  geom_errorbar(data = sum.base %>% filter(type == "Sgen"), 
                  mapping = aes(ymin = log(LowerCI), ymax = log(UpperCI), colour = lh )) + 
  theme_classic() +
  ylab("Log Sgen") +
  xlab("Log Watershed Area")

ggplot(data = sum.line %>% filter(type == "Srep"), 
    aes(x = logWA, y = log(Mean), colour = lh, fill = lh)) +
  geom_line() + 
  geom_ribbon(aes(ymin = log(LowerCI), ymax = log(UpperCI)), alpha = .2) +
  geom_point(data = sum.base %>% filter(type == "Srep"), 
              mapping = aes(logWA, log(Mean), colour = lh)) +
  geom_errorbar(data = sum.base %>% filter(type == "Srep"), 
                  mapping = aes(ymin = log(LowerCI), ymax = log(UpperCI), colour = lh )) + 
  theme_classic() +
  ylab("Log Srep") +
  xlab("Log Watershed Area")

ggplot(data = sum.line %>% filter(type == "Smsy"), 
    aes(x = logWA, y = log(Mean), colour = lh, fill = lh)) +
  geom_line() + 
  geom_ribbon(aes(ymin = log(LowerCI), ymax = log(UpperCI)), alpha = .2) +
  geom_point(data = sum.base %>% filter(type == "Smsy"), 
              mapping = aes(logWA, log(Mean), colour = lh)) +
  geom_errorbar(data = sum.base %>% filter(type == "Smsy"), 
                  mapping = aes(ymin = log(LowerCI), ymax = log(UpperCI), colour = lh )) + 
  theme_classic() +
  ylab("Log Smsy") +
  xlab("Log Watershed Area")


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
             fittedPredictedResponse = logRS_predMean, integerResponse = FALSE)
plot(DHARMaRes)
boxplot(DHARMaRes$scaledResiduals~srdat$Name)
abline(h = 0.5, col = 'red')
srdat$res <- DHARMaRes$scaledResiduals

## Which stocks do we not predict well?
ggplot(data = srdat, aes(x = Name, y = res)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.5, col = 'red') + 
  # add in tilt
  theme_classic()

## Look how alpha varies by stock?
logalpha.sum <- summary(mcmc.out[, grep('logAlpha\\[', colnames(mcmc.out[[1]]))])
srdat.alpha <- WAbase %>% bind_cols(data.frame(do.call("cbind", logalpha.sum)))
ggplot(data = srdat.alpha, aes(x = Name, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = X2.5., ymax = X97.5.)) +
  # add titl
  theme_classic()

## Are there bad years?
ggplot(data = srdat, aes(x = factor(Yr), y = res)) +
  geom_boxplot() +
  theme_classic()
## Probably some sort of overall annual impact on productivity to account for.




## Just Paul playing around....
### Parameterizing Lognormal state space models using moment matching  
### https://link.springer.com/article/10.1007/s10651-023-00570-x
stock_recruit_srep_momentmatching <- nimbleCode({
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
    tauobs[i] ~ dgamma(0.1, 0.1) # stock specific precision
  }

  ## Here we are going to do moment matching:
  for(i in 1:N_Obs){ # 501 iterations - number of base observations
    RS_pred[i] <- exp(logAlpha[stk[i]]*(1 - S[i]/E[stk[i]]))
    logRS_pred[i] <- log(RS_pred[i]^2/sqrt(RS_pred[i]^2 + 1/tauobs[stk[i]]))
    tau_obs[i] <- 1/log(1 + 1/(tauobs[stk[i]]*RS_pred[i]^2))
    logRS[i] ~ dnorm( logRS_pred[i], tau = tau_obs[i]) # stock specific precision
  }
})

inits <- function(){
  list(
    beta = c(10,0,0,0),
    logAlpha0 = 1.5,
    logESD = 1,
    logAlphaSD = 1,
    logE0 = rep(0, constants$N_Stk),
    logAlpha = rep(1.2, constants$N_Stk),
    tauobs = rep(0.001, constants$N_Stk)
    )
}

model <- nimbleModel(stock_recruit_srep_momentmatching, data = data, constants = constants,
  inits = inits(), buildDerivs=FALSE) # Build derivs if you want MLE

model$beta <- c(10,0,0,0)
model$calculate()
model$E

model$calculate("logRS[10]")
model$simulate()
model$logRS_pred[1]
model$tau_obs
hist(model$RS_pred)
sqrt(1/model$tauobs[1])
model$simulate("tauobs")

cmodel$RS_pred

conf <- configureMCMC(model)
conf$setMonitors(c("logAlpha", "logAlpha0", "E", "logRS_pred",
                    "beta", "tauobs",
                    "logAlphaSD", "logESD", "logE0"))
mcmc <- buildMCMC(conf)
cmodel <- compileNimble(model)

cmcmc <- compileNimble(mcmc, project = model)
mcmc.out.new <- runMCMC(cmcmc, niter=50000, nburnin=5000, nchains=3, samplesAsCodaMCMC = TRUE)

MCMCtrace(object = mcmc.out, params = "beta")

summary(mcmc.out.new[, c("beta[1]", "beta[2]", "beta[3]",  "beta[4]")])
summary(mcmc.out[, c("beta[1]", "beta[2]", "beta[3]",  "beta[4]")])

summary(mcmc.out.new[, "logAlpha0"])
summary(mcmc.out[, "logAlpha0"])
