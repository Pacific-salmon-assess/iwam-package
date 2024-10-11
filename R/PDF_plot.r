# Code to plot prior distributions for sigma

library(invgamma)
library(viridis)
library(latex2exp)

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


plotPriors <- function (plot_inv_gamma_only, Delta, modelobject){
  #plot_inv_gamma_only TRUE = then only plot inverse gamma distributions and not cauchy, normal or uniform
  #Delta TRUE = plot priors for the SD of the Watershed-area model
  #Delta FALSE = plot priors for the SD of the Ricker model 
  
  test <- seq(0.00001,9,len=10000)
  
  #plot(x=test, y=abs(dcauchy(test,0,1)), type="n", ylab="Probability density", xlab="Ricker Sigma (or LogA sigma)" , ylim=c(0,0.8), xlim=c(0,2.5))
  plot(x=test, y=abs(dcauchy(test,0,1)), type="n", ylab="Probability density", xlab="Sigma for watershed-area regression" , ylim=c(0,0.8), xlim=c(0,2.5))
  cols<-viridis(4, alpha=0.9)
  
  if(!Delta){ # FLAG ***********************************************************
    # Add histogram of Ricker sigmas
    
    # Ricker sigma from Parken data, excluding Upper Columbia and Siluetz, where they included AR(1) term in model
    RicSigPark <- read.csv(here::here("DataIn/ParkenRicSig.csv"))$RicSig 
    
    # Ricker sigma from PSE Chinook SR data extracted 15 Oct 2020. See sigR_metaanalysis.R
    RicSig <- exp(read.table(here::here("DataOut/PSE_sigma.csv"))$Estimate) 
    
    hist_out <- hist(RicSig, plot=FALSE)
    hist_outPark <- hist(RicSigPark, plot=FALSE)
    lowlimit <- hist_out$breaks[1] # Indicates where bars should start on the plot, assuming units of 0.1
    
    barplot(height =  c( rep(0,lowlimit*10),(hist_out$density)*0.2), width=0.1, col="light blue", border=grey(0.7), space=0, add=TRUE)
    barplot(height= hist_outPark$density*0.2, width=0.1, col=grey(0.8, alpha=0.5), border=grey(0.7, alpha=0.5), space=0, add=TRUE)
    
    RicSig_Thorson <- 0.69 #Thorson et al. 2014 marginal sigma from predictive distribution from hierarchical Ricker model of 20 salmonid stocks from Myers et al 1995
    RicSigSD_Thorson <- 0.27
    abline(v=RicSig_Thorson, lwd=2, lty="dotted", col=grey(0.8)) # RicSig_Thorson
    #polgyon(x=c(RicSig_Thorson-RigSig_Thorson))
  }
  
  # Better to use external data to form priors, 
  # E.g., NCEAS State of Alaska Salmon and People project:	
  # Brendan Connors: Chinook (n = 75) most from US West Coast and of questionable quality, remaining ~20 stocks are from AK and of higher quality
  
  
  
  #Inverse gamma on sqrt(variance)=sigma
  if(!Delta){
    shape<-rate<-0.001
    lines(sqrt(test), sqrt(dinvgamma(test,shape=shape, rate=rate)), lwd=4, col=t_col(color=cols[1], percent=90))
  }
  
  # shape<-rate<-0.01
  # lines(sqrt(test), sqrt(dinvgamma(test,shape=shape, rate=rate)), lwd=4, col=t_col(color=cols[1], percent=70))
  # shape<-rate<-0.1
  # lines(sqrt(test), sqrt(dinvgamma(test,shape=shape, rate=rate)), lwd=4, col=t_col(color=cols[1], percent=50))
  # shape<-rate<-1
  # lines(sqrt(test), sqrt(dinvgamma(test,shape=shape, rate=rate)), lwd=4, col=t_col(color=cols[1], percent=30))
  
  if (!plot_inv_gamma_only){
    
    if(!Delta){
      #Half-Normal on sigma
      lines(x=test, y=abs(dnorm(test,0,1)), col=cols[2], lwd=4)
    
    
      # Half-cauchy on sigma, as implemented in TMB
      # student T, reduced to cauchy with scale=1. Adapted from Kasper Kristensen
      # https://www.stockassessment.org/datadisk/stockassessment/userdirs/user80/WBSS_HAWG_2018/TMB/include/distributions_R.hpp
      dt<- function(x, df=1, give_log=FALSE)
      {
        logres <- lgamma((df+1)/2) - 1/2*log(df*pi) -lgamma(df/2) - (df+1)/2*log(1+x*x/df)
        #logres <- log(1/(pi*(1+x^2)))  #Equivalent ( see http://www.math.wm.edu/~leemis/chart/UDR/PDFs/TStandardcauchy.pdf)
        if(!give_log) return (exp(logres))
        else return (logres)  
      }
    
      lines(x=test, y=dt(test), col=cols[3], lwd=4)
    
      #Altrenative, equivalent half cauchy using R function
      #lines(x=test, y=abs(dcauchy(test,0,1)), lwd=2, col="red")
    
      #Uniform 0-1
      lines(x=c(0,2,2,3), y=c(0.5,0.5,0,0), col=cols[4], lwd=2)
    }
    
    if(Delta){
      #Uniform 0-1
      lower <- 0.21 #See KFrun.R
      upper <- sd(log(read.csv("DataIn/ParkenSMSY.csv")$SMSY))
      lines(x=c(0,lower,lower, upper,upper,3), y=c(0,0,0.8,0.8,0,0), col=cols[4], lwd=4)
      #Normal bounded
      norm_scalar <- 1.8
      testb <- seq(lower, upper,len=100)
      lines(x=testb, y=abs(dnorm(testb,0.8,0.28))/norm_scalar, col=cols[2], lwd=4)#See KRrun.R for N(0.8,0.28)
      lines(x=c(0,lower, lower), y=c(0, 0, dnorm(lower,0.8,0.28))/norm_scalar, col=cols[2], lwd=4)
      lines(x=c(upper,upper, 3), y=c(dnorm(upper,0.8,0.28)/norm_scalar, 0, 0), col=cols[2], lwd=4)
      
    }
    
    if(Delta){
      # For Delta sigma, vertical line SD of SMSY among 25 stocks
      abline(v= sd(log(read.csv("DataIn/ParkenSMSY.csv")$SMSY)), col=grey(0.5))
      abline(v= sqrt((0.293+0.146)/2), col=grey(0.5), lty="dashed")#See Parken et al. (2006)
      abline(v= 0.21, col=grey(0.5), lty="dotted")#See KFrun.R
    }
    
    if(!Delta){
      legend(x=1.1, y=0.78, legend=c("Inverse gamma(1,1)", "Inverse gamma(0.1,0.1)", "Inverse gamma(0.01,0.01)", 
                                     "Inverse gamma(0.001,0.001)", "Half Normal (0,1)", "Half Cauchy (0,1)", "Uniform (0,2)"),
             col=c(t_col(cols[1], 30), t_col(cols[1], 50), t_col(cols[1], 70), t_col(cols[1], 90), cols[2:4]), bty="n", lwd=2) 
    }
    if(Delta){
      # legend(x=1.4, y=0.78, legend=c("Inverse gamma(1,1)", "Inverse gamma(0.1,0.1)", "Inverse gamma(0.01,0.01)",
      #                               "Half Normal (0,1)", "Half Cauchy (0,1)", "SD of ln(SMSY) Parken et al.",
      #                               "sigma WA regression Parken et al.", "Med. SD of time-varying SMSYs"),
      #        col=c(t_col(cols[1], 30), t_col(cols[1], 50), t_col(cols[1], 70), cols[2:3], rep(grey(0.5),3)),
      #        lty=c(rep("solid", 6), "dashed", "dotted"), bty="n", lwd=2, cex=0.8)
      legend(x=1.4, y=0.78, legend=c("Inverse gamma(1,1)", "Inverse gamma(0.1,0.1)", "Inverse gamma(0.01,0.01)",
                                     "Normal bounded", "Uniform bounded", "SD of ln(SMSY) Parken et al.",
                                     "sigma WA regression Parken et al.", "Med. SD of time-varying SMSYs"),
             col=c(t_col(cols[1], 30), t_col(cols[1], 50), t_col(cols[1], 70), cols[2], cols[4], rep(grey(0.5),3)),
             lty=c(rep("solid", 6), "dashed", "dotted"), bty="n", lwd=2, cex=0.8)
      
    }
  }
  
  if (plot_inv_gamma_only){
    if(!Delta){
      legend(x=1.05, y=0.8, legend=c("Inverse gamma(1,1)", # 1
                                   "Inverse gamma(0.1,0.1)", # 2 
                                   "Inverse gamma(0.01,0.01)", # 3
                                   "Inverse gamma(0.001,0.001)", # 4
                                   "Thorson et al. 2014 marginal sigma", # 5
                                   "Distribution of Ricker Sigma Parken", # 6
                                   "Distribution of Ricker Sigma PSE"), # 7
             col=c(t_col(cols[1], 30), 
                   t_col(cols[1], 50), 
                   t_col(cols[1], 70), 
                   t_col(cols[1], 90), 
                   grey(0.7, alpha = 0.5),
                   grey(0.7, alpha = 0.5), 
                   "light blue"), 
             bty="n", 
             pch=c(NA,NA,NA,NA,NA, 22,22), 
             pt.bg= c(NA, NA, NA, NA, NA,
                      grey(0.7, alpha=0.5), 
                     "light blue"), 
             lty = c("solid", "solid", "solid", "solid", "dashed", NA, NA),
             lwd = c(rep(2,5),NA,NA) )
    }
    if(Delta){
      # NEW LINES
      mleest_smsy <- exp(modelobject[5,1])
      # mleest_srep <- modelobject[10,1]
      smsy95 <- mleest_smsy + 1.96 * modelobject[5,2]
      smsy5 <- mleest_smsy - 1.96 * modelobject[5,2]
        # MLE line for estimate of WA regression sigma
        # Boundaries for 95% confidence intervals of sigma
      # abline(v = mleest_smsy, col = grey(0.5)) # solid
      # abline(v = smsy5, col = grey(0.5), lty = "dotted")
      # abline(v = smsy95, col = grey(0.5), lty = "dotted")
      
      # as a Rectangle
      y_minplot <- par("usr")[3]  # Minimum of the y-axis (automatically determined by R)
      y_maxplot <- par("usr")[4]  # Maximum of the y-axis (automatically determined by R)
      
      rect(smsy5, y_minplot, smsy95, y_maxplot - 0.001, col = grey(0.9), border = NA)
      abline(v = mleest_smsy, col = grey(0.5), lwd = 2) # solid
      
      # For Delta sigma, vertical line SD of SMSY among 25 stocks
      abline(v= sd(log(read.csv("DataIn/ParkenSMSY.csv")$SMSY)), col=t_col(cols[2], 50), lwd = 2)
      abline(v= sqrt((0.293+0.146)/2), col=t_col(cols[2], 50), lty="dashed", lwd = 2) #See Parken et al. (2006)
      abline(v= 0.21, col=t_col(cols[2], 50), lty="dotted", lwd = 2) #See KFrun.R
      
      # Adjust legend for new lines - make sure they are short enough titles
      legend(x=1.4, y=0.78, legend=TeX(c("Inverse gamma(1,1)", 
                                     "Inverse gamma(0.1,0.1)", 
                                     "Inverse gamma(0.01,0.01)", 
                                     r"($\sigma$ WA-SMSY reg. MLE)", 
                                     "95% CI",
                                     "SD of ln(SMSY) Parken et al.",
                                     r"($\sigma$ WA reg. Parken et al.)",
                                     "Med. SD of time-varying SMSYs")),
             # col=c(t_col(cols[1], 30), t_col(cols[1], 50), t_col(cols[1], 70), rep(grey(0.5),3)), 
             col = c(t_col(cols[1], 30), t_col(cols[1], 50), t_col(cols[1], 70), 
                     rep(grey(0.5),2),
                     rep(t_col(cols[2], 50), 3)), 
             # lty=c(rep("solid", 4), "dashed", "dotted"),
             lty = c(rep("solid", 4), NA, "solid", "dashed", "dotted"),
             pch = c(rep(NA, 4), 15, rep(NA, 3)),
             pt.bg = c(rep(NA,4), grey(0.5), rep(NA,3)),
             pt.cex = 2,
             bty = "n", 
             lwd=2, 
             cex=0.8)
    }
  }
  
  shape<-rate<-0.01
  lines(sqrt(test), sqrt(dinvgamma(test,shape=shape, rate=rate)), lwd=4, col=t_col(color=cols[1], percent=70))
  shape<-rate<-0.1
  lines(sqrt(test), sqrt(dinvgamma(test,shape=shape, rate=rate)), lwd=4, col=t_col(color=cols[1], percent=50))
  shape<-rate<-1
  lines(sqrt(test), sqrt(dinvgamma(test,shape=shape, rate=rate)), lwd=4, col=t_col(color=cols[1], percent=30))
  
}

# cols <- c("#440154E6" "#31688EE6" "#35B779E6" "#FDE725E6")
# Cauchy = normal/sqrt(chi^2), Gelman et al. 2006. Doesnt work
# prior.scale <- 1
# xi <- dnorm (test, 0, prior.scale)
# tau.eta <-  dgamma (test,0.5,0.5) # chi^2 with 1 d.f.
# sigma.theta <- abs(xi)/sqrt(tau.eta)
# plot(x=test, y=sigma.theta, type="l")

# Plotting for: Sigma for Watershed Area Regression ####
# TK: Shows the probability density function for each of the prior distributions
# TK: For examples see TWG LRPs Mtg 7 ppt
# plotPriors only has two inputs - so 4 possible plots
# Plot should be the same between prior sets since it should have no dependence
# on the data itself, just the density functions?
  # E.g. one only has to change the labelling between results sections for
  # the sensitivity analysis

# Copy the below code whenever you want to print the actual plot
    # par(mfrow=c(1,1))  
    # 
    #   # Ricker Priors - InvGamma
    # png(paste("DataOut/RicPriors_InvGamma.png", sep=""), width=7, height=7, units="in", res=500)
    # plotPriors(plot_inv_gamma_only=TRUE, Delta=FALSE)
        # Dashed line?
        # X-title is wrong
    # dev.off()
    # 
    #   # WA sigma priors - InvGamma
    # #png(paste("DataOut/DeltaPriors_InvGamma.png", sep=""), width=7, height=7, units="in", res=500)
    # plotPriors(plot_inv_gamma_only=TRUE, Delta=TRUE, modelobject = iwam_default$all_Deltas)
    # dev.off()
    # 
    #   # Random.
    # plotPriors(plot_inv_gamma_only=FALSE, Delta=TRUE)
    # 
    #   # Ricker Priors 
    # #png(paste("DataOut/RicPriors_sm.png", sep=""), width=7, height=7, units="in", res=500)
    # plotPriors(plot_inv_gamma_only=FALSE, Delta=FALSE)
    # #dev.off()


#png(paste("DataOut/DeltaPriors.png", sep=""), width=7, height=7, units="in", res=500)
# plotPriors(plot_inv_gamma_only=FALSE, Delta=TRUE)
#dev.off()
