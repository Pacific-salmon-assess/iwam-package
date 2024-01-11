# Inverse logit and logit funcitons can come in handy =====================================================================
inv_logit <- function(x){
  exp(x)/(1+exp(x))
}

logit <- function(x){
  log(x/(1-x))
}

# Scaling by digits
# Calculate scale for each stock as a tibble (tidyverse df)
digit_scaling <- function(input_data){ # THIS input_data IS SPECIFICALLY FOR SRDAT - with many years of data that need to be grouped
  digits <- input_data %>% group_by(Stocknumber) %>% 
    summarize(maxDigits = count.dig(max(Sp)))
  # count.dig() Creates a count [numeric] of the max number of digits 
  # of spawners as digits per stock
  # the function count.dig() can be found in the script: helperFunctions.R
  input_data <- left_join(input_data, digits)
  input_data <- input_data %>% mutate(scale = 10^(maxDigits-1))
  # data <- data %>% mutate(scale = 10^4)) # Alternate scale
  # using mutate; creates a per stock scale by taking the number of digits - 1,
  # as the exponent on a base 10 log scale
}
# 1) group dataset by stock id
# 2) summarize(maxDigits = count.dig(max()))


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


count.dig <- function(x) {floor(log10(x)) + 1}

'%not in%' <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

# PredInt = a function to calculate prediction intervals
# x = independent variable
# y = dependenet variable
# Newx = x variables for which you'd like the prediction intervals
# Predy = Predicted y variable at those Newx's
PredInt <- function(x,y,Newx=x, Predy){
  sumXErr <- sum( (x-mean(x))^2 )
  sumYErr <- sum( (y-mean(y))^2 )
  sumXYErr <- sum( (y - mean(y)) * (x - mean(x)) ) ^2
  STEYX <- sqrt( (1/(length(x)-2)) * (sumYErr - sumXYErr/sumXErr) )
  
  # SE of the prediction from http://www.real-statistics.com/regression/confidence-and-prediction-intervals/
  SE.pred <- STEYX * sqrt( (1 + 1/length(x) + ((Newx - mean(x))^2)/sumXErr) ) 
  t.crit <- qt(0.975,df=length(x)-2) #95% intervals
  
  upr <- Predy + SE.pred*t.crit
  lwr <- Predy - SE.pred*t.crit
  PI <- list()
  PI$upr <- upr
  PI$lwr <- lwr
  return(PI)
  
}



sGenOptimum <- function ( S, theta ) {
  # Function called from sGenSolver 
  loga <- theta[1]
  b <- theta[2]
  prt <- S * exp( loga - b * S)
  sMSY <- ( 1 - gsl::lambert_W0 (exp ( 1 - loga) ) ) / b
  epsilon <- log(sMSY) - log(prt)
  nlogLike <- - sum( dnorm ( epsilon, 0, 1, log = T))
  return( nlogLike )
}


sGenSolver <- function (loga, b) {
  # Function to estimate Sgen from loga and b Ricker parameters
  theta <- c(loga, b)
  sMSY <- (1 - gsl::lambert_W0(exp(1 - loga))) / b
  fit <- optimize(f = sGenOptimum, interval = c(0, sMSY),
                  theta = theta)
  return(fit$minimum)
}


bEst <- function(b, SMSY, SREP){
  # Function to estimate Ricker b paramter from SMSY and SREP
  # (1 − b * SMSY) * exp( 1 − b * SMSY) = exp(1 − loga) (Explicit formula for SMSY; Eqn 10 Scheueurell 2014)
  # Rearranging:
  # log (1 − b * SMSY) + (1 − b * SMSY) = 1 - loga
  # log (1 − b * SMSY) - b * SMSY = - loga
  # loga <-  b*SMSY - log(1-b*SMSY)  
  
  # SREP <- loga / b # from Ricker 1975
  # loga <- SREP * b # Rearrange and substitue this into equation above
  # SREP * b <- b*SMSY - log(1-b*SMSY)  
  # SREP <- (b * SMSY - log( 1 - b * SMSY) ) / b
  # Use this equation to solve for b
  
  pred <- (b * SMSY - log( 1 - b * SMSY) ) / b
  epsilon <- log(SREP) - log(pred)
  nlogLike <- - sum( dnorm ( epsilon, 0, 1, log = T))
  return(nlogLike)
}
bSolver <- function (SMSY, SREP) {
  # Function to estimate b from SMSY and SREP, where SMAX is bounded between 1/3 of SREP and SREP
  fit <- optimize(f = bEst, interval = c(1/SREP, 3/SREP), SMSY=SMSY, SREP=SREP)
  return(fit$minimum)
}

Sgen.fn <- function ( SMSY, SREP, half.a = FALSE, const.SMAX =FALSE, explicit = TRUE , plot=FALSE) {
  # Function to convert SMSY and SREP from watershed-area model into Sgen
  # half.a = should we assume productivity has dropped by half?
  # explicit = should we use the explicit relationship between SMSY and Ricker parameters as in Scheuerell 2014?
  
  
  if (explicit){
    a.par <- exp( bSolver(SMSY, SREP) * SREP )
    b.par <- log( a.par ) / SREP
    if (half.a) a.par <- 0.5 * exp( bSolver(SMSY, SREP) * SREP )
    if (half.a) b.par <- log( a.par ) / SREP # b parameter adjusted based on half alpha, assuming SREP is unchanged
    # Alternaively, could assume SMAX (= 1/b) remains constant, but SREP changes
    if (half.a) if(const.SMAX) b.par <- bSolver(SMSY, SREP) #(same b as for constant prod case)
    if (half.a) if(const.SMAX) SREP <- log(a.par) / b.par
    if (half.a) SMSY <- (1 - gsl::lambert_W0(exp(1 - log(a.par) ))) / b.par
  }
  
  if( !explicit ){
    a.par <- exp( (0.5 - SMSY / SREP ) / 0.07 )                
    if (half.a) a.par <- (1/2) * exp( ( 0.5 - SMSY / SREP) / 0.07 )  ## Half alpha calculated here
    b.par <- log (a.par) / SREP  ## b parameter adjusted based on half alpha, assuming SREP is unchanged
    # Alternaively, could assume SMAX (= 1/b) remains constant, but SREP changes
    if (half.a) if(const.SMAX) b.par <- exp( (0.5 - SMSY / SREP) / 0.07) / SREP #(same b as for constant prod case)
    if (half.a) if(const.SMAX) SREP <- log(a.par) / b.par
    if (half.a) sMSY <-(1 - gsl::lambert_W0(exp(1 - log(a.par) ))) / b.par
  }
  
  sgen.out <- sGenSolver( log(a.par), b.par )
  
  if(plot){
    Rpred <- NA
    for (i in 1:1000){ Rpred[i]<- a.par * i * exp (- b.par * i)}
    if (const.SMAX) xlab <- "Spawners" else xlab <- ""
    plot(1:1000, Rpred, type="l", ylim = c (0, 1400), xlab = xlab,  ylab = "Recruits", lwd=2 )
    abline(a=0, b=1)
    abline(v=sgen.out, lty="dotted")
    abline(v=SMSY, lty="dashed")
    abline(v=(1/b.par), lty="dotdash")
    abline(v=SREP)
    legend( x = 700, y = 600, legend = c("Sgen", "SMSY", "SMAX", "SREP" ), lty=c("dotted","dashed", "dotdash", "solid"), bty="n" )
    if (!half.a) title("Constant productivity")
    if (half.a) if(!const.SMAX) title("Half productivity; constant SREP")
    if (half.a) if(const.SMAX) title("Half productivity; constant SMAX")
  }
  
  return( list( SGEN = sgen.out , SMSY = SMSY, SREP = SREP, apar = a.par, bpar = b.par) )
  
}

# # Example: Artlish
# SMSY <- 345 
# SREP <- 971
# 
# png(paste("DataOut/Artlish_WCVI_SRcurve.png", sep=""), width=4, height=7, units="in", res=500)
#   par(mfcol=c(3,1),  mar = c(4, 4, 2.5, 2) + 0.1)
#   Sgen.fn(SMSY, SREP, half.a = FALSE, const.SMAX = FALSE, plot=TRUE)
#   Sgen.fn(SMSY, SREP, half.a = TRUE, const.SMAX = FALSE, plot=TRUE)
#   Sgen.fn(SMSY, SREP, half.a = TRUE, const.SMAX = TRUE, plot=TRUE)
# dev.off()

Sgen.fn2 <- function ( a.par, SREP,  explicit = TRUE , plot=FALSE) {
  # Function to convert SREP from watershed-area model with independent alpha into Sgen
  # explicit = should we use the explicit relationship between SMSY and Ricker parameters as in Scheuerell 2014?
  # b = log(a.par)/SREP
  b.par <- log(a.par)/SREP
  if (explicit){
    SMSY <- (1 - gsl::lambert_W0(exp(1 - log(a.par) ))) / (b.par)
    
    }
  
  if( !explicit ){
    SMSY <- log(a.bar)/ b.par * (0.5 - 0.07 * log(a.par))
  }
  
  sgen.out <- sGenSolver( log(a.par), b.par )
  
  if(plot){
    Rpred <- NA
    for (i in 1:1000){ Rpred[i]<- a.par * i * exp (- b.par * i)}
    if (const.SMAX) xlab <- "Spawners" else xlab <- ""
    plot(1:1000, Rpred, type="l", ylim = c (0, 1400), xlab = xlab,  ylab = "Recruits", lwd=2 )
    abline(a=0, b=1)
    abline(v=sgen.out, lty="dotted")
    abline(v=SMSY, lty="dashed")
    abline(v=(1/b.par), lty="dotdash")
    abline(v=SREP)
    legend( x = 700, y = 600, legend = c("Sgen", "SMSY", "SMAX", "SREP" ), lty=c("dotted","dashed", "dotdash", "solid"), bty="n" )
    if (!half.a) title("Constant productivity")
    if (half.a) if(!const.SMAX) title("Half productivity; constant SREP")
    if (half.a) if(const.SMAX) title("Half productivity; constant SMAX")
  }
  
  return( list( SGEN = sgen.out , SMSY = SMSY, SREP = SREP, apar = a.par, bpar = b.par) )
  
}

ggplot.corr <- function(data, lag.max = 10, ci = 0.95, title="") {
  # Adapted from https://rh8liuqy.github.io/ACF_PACF_by_ggplot2.html

  list.acf <- acf(data, lag.max = lag.max, type = "correlation", plot = FALSE)
  N <- as.numeric(list.acf$n.used)
  df1 <- data.frame(lag = list.acf$lag, acf = list.acf$acf)
  df1$lag.acf <- dplyr::lag(df1$acf, default = 0)
  df1$lag.acf[2] <- 0
  df1$lag.acf.cumsum <- cumsum((df1$lag.acf)^2)
  df1$acfstd <- sqrt(1/N * (1 + 2 * df1$lag.acf.cumsum))
  df1$acfstd[1] <- 0
  df1 <- df1 %>% dplyr::select(lag, acf, acfstd)
  
  plot.acf <- ggplot(data = df1, aes( x = lag, y = acf)) +
    geom_col(fill = "#4373B6", width = 0.7) +
    geom_hline(yintercept = qnorm((1+ci)/2)/sqrt(N), 
               colour = "sandybrown",
               linetype = "dashed") + 
    geom_hline(yintercept = - qnorm((1+ci)/2)/sqrt(N), 
               colour = "sandybrown",
               linetype = "dashed") + 
    scale_x_continuous(breaks = seq(0,max(df1$lag),6)) +
    scale_y_continuous(name = element_blank(), 
                       limits = c(min(df1$acf, - qnorm((1+ci)/2)/sqrt(N)) ,1)) +
    ggtitle(paste("ACF for", title)) +
    theme_bw() + 
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(size = 14))
          
  
  return(plot.acf)
 
}
