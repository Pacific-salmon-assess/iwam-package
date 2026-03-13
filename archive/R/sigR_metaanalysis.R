#Code to estimate sigR from PSE data

#---------------------------------------------------------
# Libaries

library(rsample)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(reshape2)
library(TMB)
library(zoo)
library(viridis)
library(hrbrthemes)

# Functions
count.dig <- function(x) {floor(log10(x)) + 1}
'%not in%' <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))


# Read in data
SRDat <- read.csv("DataIn/PSE_ChinookData.csv")

# Summary
summary_df <- SRDat %>%
     group_by(Name) %>%
     summarize(
           num_years = n(),
           year_range = paste0(min(Yr), "-", max(Yr)),
          .groups = "drop"
    )

# Calculate scale for each stock
digits <- SRDat %>% group_by(Stocknumber) %>% summarize(maxDigits = count.dig(max(Sp)))
SRDat <- left_join(SRDat, digits)
SRDat <- SRDat %>% mutate(Scale = 10^(maxDigits-1))
Scale.stock <- (SRDat %>% group_by(Stocknumber) %>% summarize(Scale.stock = max(Scale)))$Scale.stock

data <- list()
Scale <- SRDat$Scale 
data$S <- SRDat$Sp/Scale 
data$logR <- log( (SRDat$Rec/Scale) )
data$stk <- as.numeric(SRDat$Stocknumber)
N_Stocks <- length(unique(SRDat$Name))
# data$yr <- SRDat$yr_num

param <- list()
param$logA <- ( SRDat %>% group_by (Stocknumber) %>% summarise(yi = lm(log( Rec / Sp) ~ Sp )$coefficients[1] ) )$yi
B <- SRDat %>% group_by(Stocknumber) %>% summarise( m = - lm(log( Rec / Sp) ~ Sp )$coefficients[2] )
param$logB <- log ( 1/ ( (1/B$m)/Scale.stock ))
param$logSigma <- rep(-2, N_Stocks)

# dyn.unload(dynlib("TMB_Files/Ricker"))
# compile("TMB_Files/Ricker.cpp")
dyn.load(dynlib("TMB_Files/Ricker"))

obj <- MakeADFun(data, param, DLL="Ricker", silent=TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e5, iter.max = 1e5))
pl <- obj$env$parList(opt$par)

All_Ests <- data.frame(summary(sdreport(obj)))
All_Ests$Param <- row.names(All_Ests)


# Put together readable data frame of values
All_Ests$Param <- sapply(All_Ests$Param, function(x) (unlist(strsplit(x, "[.]"))[[1]]))
All_Ests <- All_Ests %>% filter (Param %in% c("logA", "logB", "logSigma"))

SN <- unique(SRDat[, c("Stocknumber")])
All_Ests$Stocknumber <- rep(SN)
All_Ests <- left_join(All_Ests, unique(SRDat[, c("Stocknumber", "Name")]))

All_Est <- All_Ests

LogSigma <- All_Ests %>% filter (Param %in% c("logSigma"))
PSE_sigma<-exp(LogSigma$Estimate)
write.table(LogSigma,"DataOut/PSE_sigma.csv")

# read.table("DataOut/PSE_sigma.csv")

# Plot SR curves

Stks <- unique(SRDat$Stocknumber)
NStks <- length(Stks)
par(mfrow=c(5,5), mar=c(2, 2, 1, 0.1) + 0.1)
  
for (i in Stks){
  names <- All_Est %>% select ("Name", "Stocknumber") %>% distinct()
  name <- All_Est %>% filter (Stocknumber==i) %>% select ("Name") %>% distinct()
    
  R <- SRDat %>% filter (Stocknumber==i) %>% select(Rec) 
  S <- SRDat %>% filter (Stocknumber==i) %>% select(Sp) 
  # what is the scale of Ricker b estimate?
  Sc <- SRDat %>% filter (Stocknumber==i) %>% select(Scale) %>% distinct() %>% as.numeric()
    
  plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)), ylim=c(0,max(R$Rec) ) )
  a <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="logA") %>% 
    summarise(A=exp(Estimate)) %>% as.numeric()
  # Divide b by scale
  b <- All_Est %>% filter (Stocknumber==i) %>% filter(Param=="logB") %>% 
    summarise(B=exp(Estimate)/Sc) %>% as.numeric()
    
    
  for (j in 1:100){
    SS[j] <- j*(max(S$Sp)/100)
    RR[j] <- a * SS[j] * exp(-b * SS[j])
    }
     
  col.use <- "black"
  lines(x=SS, y=RR, col=col.use) 
    
  mtext(name$Name, side=3, cex=0.8)
    
  }

