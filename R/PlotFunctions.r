# Plot SR models: curves and linearized mddel

# Functions ####
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
  t.col
  #invisible(t.col)
}

# Plot SR curves ####

PlotSRCurve <- function(srdat, pars, SMSY_std=NULL, stksNum_ar=NULL, stksNum_surv=NULL, stks_surv, r2, removeSkagit, mod) {
  Stks <- unique(srdat$Stocknumber)
  NStks <- length(Stks)
  par(mfrow=c(5,5), mar=c(2, 2, 1, 0.1) + 0.1)
  
  for (i in Stks){
    names <- pars %>% dplyr::select ("Name", "Stocknumber") %>% distinct()
    name <- pars %>% filter (Stocknumber==i) %>% dplyr::select ("Name") %>% distinct()
    
    R <- srdat %>% filter (Stocknumber==i) %>% dplyr::select(Rec) 
    S <- srdat %>% filter (Stocknumber==i) %>% dplyr::select(Sp) 
    # what is the scale of Ricker b estimate?
    Sc <- srdat %>% filter (Stocknumber==i) %>% dplyr::select(scale) %>% distinct() %>% as.numeric()
    if(name$Name != "Skagit" & name$Name != "KSR") plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)), ylim=c(0,max(R$Rec) ) )
    if(name$Name == "Skagit") plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)*3), ylim=c(0,max(R$Rec) ) )
    if(name$Name == "KSR") plot(x=S$Sp, y=R$Rec, xlab="", ylab="", pch=20, xlim=c(0,500), ylim=c(0,max(R$Rec) ) )
    
    a <- pars %>% filter (Stocknumber==i) %>% filter(Param=="logA") %>% 
      summarise(A=exp(Estimate)) %>% as.numeric()
    # Divide b by scale
    b <- pars %>% filter (Stocknumber==i) %>% filter(Param=="logB") %>% 
      summarise(B=exp(Estimate)/Sc) %>% as.numeric()
    
# Test removal 
  #   if(mod!="IWAM_FixedSep_RicStd" & mod!="Liermann" & mod!= 'IWAM_Liermann' & mod!="Liermann_PriorRicSig_PriorDeltaSig" & mod!="Liermann_HalfNormRicVar_FixedDelta"){
  #     if (i %in% stksNum_surv) {  #stocknumber 0 and either 22 or 23, depending on if Skagit is removed
  #   surv.dat <- as.data.frame(read.csv(here::here("DataIn/Surv.csv"))) %>% filter(Name==name$Name)
  #   #if(i==0) surv.dat <- as.data.frame(read.csv("DataIn/Surv.csv")) %>% filter(Name==name$Name)
  #   #if(i==22|i==23)surv.dat <- as.data.frame(read.csv("DataIn/Surv.csv")) %>% filter(Name=="Cowichan")
  #   #if(i==22|i==23) { surv.dat <- surv.dat %>% filter(Yr >= 1985 & Yr !=1986 & Yr != 1987) }
  #   if(name$Name == "Cowichan") { surv.dat <- surv.dat %>% filter(Yr >= 1985 & Yr !=1986 & Yr != 1987) }
  #   mean.log.surv <- surv.dat %>% summarize(mean = mean(log(Surv)))
  #   gamma <- pars %>% filter (Stocknumber==i) %>% filter(Param=="gamma") %>% dplyr::select(Estimate)
  #   a <- pars %>% filter (Stocknumber==i) %>% filter(Param=="logA") %>%
  #     summarise(A=exp(Estimate + gamma$Estimate*mean.log.surv$mean)) %>% as.numeric()
  # }
  #     
  #   }
    #Parken values for skagit
    skagit_alpha <- 7.74
    skagit_beta <- 0.0000657
    RR_skagit <- NA
    SS <- RR<- NA
    #RR_std <- NA
    
    if(removeSkagit==FALSE){
      for (j in 1:100){
        if (i!=22 & i!=7) SS[j] <- j*(max(S$Sp)/100)
        if (i==22) SS[j] <- j*(max(S$Sp*3)/100)
        if (i==7) SS[j] <- j*(500/100)
        RR[j] <- a * SS[j] * exp(-b * SS[j])
        if(i==22) {RR_skagit[j] <- skagit_alpha * SS[j] * exp(-skagit_beta * SS[j])}
        #if (i %in% stks_ar) {RR_std[j] <- A_std$A[which(A_std$Stocknumber==i)] * SS[j] *  exp(-B_std$B[which(B_std$Stocknumber==i)] * SS[j])}
      }
    } 
    if(removeSkagit==TRUE){
      for (j in 1:100){
        SS[j] <- j*(max(S$Sp)/100)
        RR[j] <- a * SS[j] * exp(-b * SS[j])
      }
    } 
    

    # if (i %not in% c(stksNum_ar, stksNum_surv)) col.use <- "black"
    # if (i %in% stksNum_ar) col.use <- "red"
    # if (i %in% stksNum_surv) col.use <- "blue"
    # if(mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=="Liermann_PriorRicSig_PriorDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta") 
    col.use <- "black"
      # removed above give only one mod usage and no more _ar or _surv model order usage
    lines(x=SS, y=RR, col=col.use) 
    
    #For Skagit, add Parken et al. 2006 model curve
    if(removeSkagit==FALSE){if(i==22) lines(x=SS, y=RR_skagit, lty="dashed")}
    
    mtext(name$Name, side=3, cex=0.8)
    
    # Plot SMSY_stream (black for std, red for AR(1), and dashed for Parken et al. 2006)
    SMSY <- pars %>% filter (Stocknumber==i) %>% filter(Param=="SMSY") %>% 
      summarise(SMSY = Estimate * Sc) %>% as.numeric()
    SMSY_ul <- pars %>% filter (Stocknumber==i) %>% filter(Param=="SMSY") %>% 
      summarise(SMSY_ul = Estimate * Sc + 1.96 * Std..Error * Sc ) %>% as.numeric()
    SMSY_ll <- pars %>% filter (Stocknumber==i) %>% filter(Param=="SMSY") %>% 
      summarise(SMSY_ul = Estimate * Sc - 1.96 * Std..Error * Sc ) %>% as.numeric()
    
    
    abline(v = SMSY, col=col.use)

    # if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedCombined"|mod=="Ricker_AllMod"){
    #   if (i %in% stksNum_ar) polygon(x=c(SMSY_ul, SMSY_ll, SMSY_ll, SMSY_ul), y=c(-10000,-10000,max(R$Rec),max(R$Rec)), col=rgb(1,0,0, alpha=0.1), border=NA ) 
    #   if (i %in% stksNum_surv) polygon(x=c(SMSY_ul, SMSY_ll, SMSY_ll, SMSY_ul), y=c(-10000,-10000,max(R$Rec),max(R$Rec)), col=rgb(0,0,1, alpha=0.1), border=NA ) 
    #   if (i %not in% c(stksNum_ar, stksNum_surv))  polygon(x=c(SMSY_ul, SMSY_ll, SMSY_ll, SMSY_ul), y=c(-10000,-10000,max(R$Rec),max(R$Rec)), col=grey(0.8, alpha=0.4), border=NA )
    # }
    
    if(mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=='IWAM_Liermann'|mod=="Liermann_PriorRicSig_PriorDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta")  polygon(x=c(SMSY_ul, SMSY_ll, SMSY_ll, SMSY_ul), y=c(-10000,-10000,max(R$Rec),max(R$Rec)), col=grey(0.8, alpha=0.4), border=NA )
    #else polygon(x=c(SMSY_ul, SMSY_ll, SMSY_ll, SMSY_ul), y=c(0,0,max(R$Rec),max(R$Rec)), col=grey(0.8, alpha=0.4), border=NA )
  
    if(!is.null(SMSY_std)) {
      SMSY_std <- SMSY_std %>% right_join(names) %>% filter(Name==name$Name)#filter(Stocknumber != 22) 
      if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedCombined"|mod=="Ricker_AllMod"){
        if(i %in% stksNum_ar) abline(v=SMSY_std$Estimate[which(SMSY_std$Stocknumber==i)]*scale.stock[i+1] , col="black")
        if(i %in% stksNum_surv) abline(v=SMSY_std$Estimate[which(SMSY_std$Stocknumber==i)]*scale.stock[i+1] , col="black")
      }
    }
    
    # PARKEN SMSY REMOVAL *******************************************************************************************************
    # ParkenSMSY <- read.csv(here::here("DataIn/ParkenSMSY.csv"))
    #if (removeSkagit==TRUE) ParkenSMSY <- ParkenSMSY %>% filter(Name != "Skagit")
    # ParkenSMSY <- ParkenSMSY %>% filter(Name==as.character(name$Name)) %>% dplyr::select (SMSY) %>% as.numeric()
    # abline(v=ParkenSMSY, lty="dashed")
    
    if(is.data.frame(r2)==TRUE) {
      lab <-  r2 %>% filter(Stocknumber==i) %>% dplyr::select(r2) %>% as.numeric() %>% round(2)
      legend("topright", legend = "", title= paste0("r2=",lab), bty="n")
    }
  }
  
}

# ------------------------------------------------------------------------------
# Plot SR linearized model
PlotSRLinear <- function(srdat, pars, SMSY_std=NULL, stksNum_ar=NULL, stksNum_surv=NULL, r2, removeSkagit) {
  Stks <- unique(srdat$Stocknumber)
  NStks <- length(Stks)
  par(mfrow=c(5,5), mar=c(3, 2, 2, 1) + 0.1)
  
  for (i in Stks){
    R <- srdat %>% filter (Stocknumber==i) %>% dplyr::select(Rec) 
    S <- srdat %>% filter (Stocknumber==i) %>% dplyr::select(Sp) 
    logRS <- log(R$Rec/S$Sp)
    
    # what is the scale of Ricker b estimate?
    Sc <- srdat %>% filter (Stocknumber==i) %>% dplyr::select(scale) %>% distinct() %>% as.numeric()
    plot(x=S$Sp, y=logRS, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)), ylim=c(0,max(logRS) ) )
    #if(i !=22 & i!=0) plot(x=S$Sp, y=logRS, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)), ylim=c(0,max(logRS) ) )
    #if(i ==22) plot(x=S$Sp, y=logRS, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)*3), ylim=c(0,max(logRS) ) )
    #if(i ==0) plot(x=S$Sp, y=logRS, xlab="", ylab="", pch=20, xlim=c(0,max(S$Sp)*3), ylim=c(min(logRS),max(logRS) ) )
    
    name <- pars %>% filter (Stocknumber==i) %>% dplyr::select ("Name") %>% distinct()
    mtext(name$Name, side=3)
    
     
    logA <- pars %>% filter (Stocknumber==i) %>% filter(Param=="logA") %>% 
      summarise(A=Estimate) %>% as.numeric()
    # Divide b by scale
    B <- pars %>% filter (Stocknumber==i) %>% filter(Param=="logB") %>% 
      summarise(B=exp(Estimate)/Sc) %>% as.numeric()

    # if(mod!="IWAM_FixedSep_RicStd" & mod!="Liermann" & mod!="IWAM_Liermann" & mod!="Liermann_PriorRicSig_PriorDeltaSig" & mod!="Liermann_HalfNormRicVar_FixedDelta"){
    #   if (i %in% stksNum_surv) {  #stocknumber 0 and either 22 or 23, depending on if Skagit is removed
    #     
    #     if(i==0) surv.dat <- as.data.frame(read.csv(here::here("DataIn/Surv.csv"))) %>% filter(Name=="Harrison") 
    #     if(i==22|i==23)surv.dat <- as.data.frame(read.csv(here::here("DataIn/Surv.csv"))) %>% filter(Name=="Cowichan") 
    #     if(i==22|i==23) { surv.dat <- surv.dat %>% filter(Yr >= 1985 & Yr !=1986 & Yr != 1987) }
    #     mean.log.surv <- surv.dat %>% summarize(mean = mean(log(Surv)))
    #     gamma <- pars %>% filter (Stocknumber==i) %>% filter(Param=="gamma") %>% dplyr::select(Estimate)
    #     logA <- pars %>% filter (Stocknumber==i) %>% filter(Param=="logA") %>% 
    #       summarise(logA.adj=(Estimate + gamma$Estimate*mean.log.surv$mean)) %>% as.numeric()
    #   }
    # }
    

    # if (i %in% stksNum_ar) col.use <- "red"
    # if (i %in% stksNum_surv) col.use <- "blue"
    # if (i %not in% c(stksNum_ar, stksNum_surv)) col.use <- "black"
    # if (mod=="IWAM_FixedSep_RicStd"|mod=="Liermann"|mod=="Liermann_PriorRicSig_PriorDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta") 
    col.use <- "black"
      # removed above give only one mod usage and no more _ar or _surv model order usage
    
    abline(a=logA, b=-B, col=col.use)
    
    # For Skagit, add the best fit based on Parken et al. 2006 Ricker estimates
    skagit_alpha <- 7.74
    skagit_beta <- 0.0000657
    if(removeSkagit==FALSE) {if (i==22)     abline(a=log(skagit_alpha), b=-skagit_beta, col="black", lty="dashed")}

    lab <-  r2 %>% filter(Stocknumber==i) %>% dplyr::select(r2) %>% as.numeric() %>% round(2)
    legend("topright", legend = "", title= paste0("r2=",lab), bty="n")

  }
  
  
}

#-------------------------------------------------------------------------------
#Plot standardized residuals ####

PlotStdResid <- function(SRes){
  
  Stks <- unique(SRes$Stocknumber)
  NStks <- length(Stks)
  par(mfrow=c(5,5), mar=c(3, 2, 2, 1) + 0.1)
  
  for (i in Stks){
    SR <- SRes %>% filter (Stocknumber==i) %>% dplyr::select(Res) 
    plot(x=1:length(SR$Res), y=SR$Res, xlab="", ylab="", pch=20, ylim=c(min(-4, min(SR$Res)),max(4, max(SR$Res)) ) ) 
    abline(h=c(3,-3), col="red")
    abline(h=0, col="black")
    name <- SRes %>% filter (Stocknumber==i) %>% dplyr::select ("Name") %>% distinct()
    mtext(name$Name, side=3)
  }
}

Plotacf <- function(Preds){
  Stks <- unique(Preds$Stocknumber)
  NStks <- length(Stks)
  par(mfrow=c(5,5), mar=c(2, 2, 1.5, 1) + 0.1)
  for (i in Stks){
    Res <- Preds %>% filter (Stocknumber==i) %>% dplyr::select(Res)
    acf(Res$Res, plot=T)
    acf1 <- acf(Res$Res, plot=F)$acf[2]
    len <- length(Res$Res)
    acf1.ci <- qnorm((1 + 0.95)/2)/sqrt(len)
    if (abs(acf1)>acf1.ci) col="red" else col="black"
    name <- Preds %>% filter (Stocknumber==i) %>% dplyr::select ("Name") %>% distinct()
    mtext(name$Name, side=3, col=col)
  }
  
}
# What are the upper and lower bounds of KSR, Stikine and Cowichan?

# KSR.SMSY <- pars %>% filter(Name=="KSR") %>% filter(Param=="SMSY") %>% select(Estimate) %>% as.numeric()
# KSR.SMSY.ul <- KSR.SMSY + 1.96 * (pars %>% filter(Name=="KSR") %>% filter(Param=="SMSY") %>% select(Std..Error) %>% as.numeric())
# KSR.SMSY.ul <- KSR.SMSY.ul * srdat %>% filter (Name=="KSR") %>% select(scale) %>% distinct() %>% as.numeric()
# 
# Stikine.SMSY <- pars %>% filter(Name=="Stikine") %>% filter(Param=="SMSY") %>% select(Estimate) %>% as.numeric()
# Stikine.SMSY.ul <- Stikine.SMSY + 1.96 * (pars %>% filter(Name=="Stikine") %>% filter(Param=="SMSY") %>% select(Std..Error) %>% as.numeric())
# Stikine.SMSY.ul <- Stikine.SMSY.ul * srdat %>% filter (Name=="Stikine") %>% select(scale) %>% distinct() %>% as.numeric()
# 
# Cow.SMSY <- pars %>% filter(Name=="Cowichan") %>% filter(Param=="SMSY") %>% select(Estimate) %>% as.numeric()
# Cow.SMSY.ul <- Cow.SMSY + 1.96 * (pars %>% filter(Name=="Cowichan") %>% filter(Param=="SMSY") %>% select(Std..Error) %>% as.numeric())
# Cow.SMSY.ul <- Cow.SMSY.ul * srdat %>% filter (Name=="Cowichan") %>% select(scale) %>% distinct() %>% as.numeric()
# Cow.SMSY.ll <- Cow.SMSY - 1.96 * (pars %>% filter(Name=="Cowichan") %>% filter(Param=="SMSY") %>% select(Std..Error) %>% as.numeric())
# Cow.SMSY.ll <- Cow.SMSY.ll * srdat %>% filter (Name=="Cowichan") %>% select(scale) %>% distinct() %>% as.numeric()


# ------------------------------------------------------------------------------
# Plot  WA regression ####
  # SREP is currently modded for loglog plotting
  # 11/12/2024 - Working on SMSY

plotWAregressionSMSY <- function (pars, all_Deltas, srdat, lifehist, WAbase,  
                                  pred_lnSMSY=NA, pred_lnWA, title1, mod) {

  options(scipen = 5) # for ticks without sci. notation
  
  SMSY <- pars %>% filter(Param=="SMSY") %>% 
    mutate(ModelOrder=0:(length(unique(pars$Stocknumber))-1))
  # what is scale of SMSY?
  Sc <- srdat %>% dplyr::select(Stocknumber, scale) %>% 
    distinct()
  SMSY <- SMSY %>% left_join(Sc, by="Stocknumber") %>% 
    mutate(rawSMSY=Estimate*scale)
  lnSMSY <- log(SMSY$rawSMSY)
  lnWA <- log(WAbase$WA)
  
  par(cex=1.5, mar = c(5, 4, 2, 2) + 0.1)
  # If you want a title - restore to default margins or just take "mar" out
  col.use <- NA
  for(i in 1:length(SMSY$lh)) {
    if (SMSY$lh[i]==0) col.use[i] <- "forestgreen" 
    else col.use[i] <- "dodgerblue3"
    }
  
  # plot(y=lnSMSY, x=lnWA, pch=20, col=col.use, xlab="log(Watershed Area, km2)", ylab="log(SMSY)")
  # points(y=lnSMSY, x=lnWA, pch=20, col=col.use, cex=1.5)
  #points(y=lnSMSY[18:length(SMSY$lh)], x=lnWA[18:length(SMSY$lh)], pch=3, col=col.use[18:length(SMSY$lh)], cex=1.5)
  
  # LOG LOG 
  plot(y=exp(lnSMSY), x=exp(lnWA), pch=20, col=col.use, 
       xlab= expression("Accessible Watershed Area, km"^2), 
       ylab= expression(S[MSY]), log = 'xy',
       # xlim = c(50,500000) , ylim=c(100,5000000))
       xlim = c(50,200000) , ylim=c(200,2000000))
  # 114434 is the max watershed area - shouldn't go further over that on the x axis
  points(y=exp(lnSMSY), x=exp(lnWA), pch=20, col=col.use, cex=1.5)
  
  logD1 <- all_Deltas %>% 
    filter(Param=="logDelta1") %>% 
    dplyr::select(Estimate) %>% 
    pull()
  logD2 <- all_Deltas %>% 
    filter(Param=="logDelta2") %>% 
    dplyr::select(Estimate) %>% 
    pull()
  
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|
     mod=="IWAM_FixedSep_Constm"|mod=="Liermann"|
     mod=="IWAM_Liermann"|mod=="Liermann_PriorRicSig_PriorDeltaSig"|
     mod=="Liermann_HalfNormRicVar_FixedDelta") {
    logD1o <- all_Deltas %>% 
      filter(Param=="logDelta1_ocean") %>% 
      dplyr::select(Estimate) %>% 
      pull() + logD1}
  
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|
     mod=="IWAM_FixedSep_Constyi"|mod=="Liermann"|
     mod=="IWAM_Liermann"|mod=="Liermann_PriorRicSig_PriorDeltaSig"|
     mod=="Liermann_HalfNormRicVar_FixedDelta") {
    D2o <- exp(all_Deltas %>% filter(Param=="logDelta2_ocean") %>% 
                 dplyr::select(Estimate) %>% 
                 pull() ) + exp(logD2)
    if(nrow(all_Deltas %>% 
            filter(Param=="Delta2_ocean"))>=1)  D2o <- (all_Deltas %>% 
                                                          filter(Param=="Delta2_ocean") %>% 
                                                          dplyr::select(Estimate) %>% 
                                                          pull() ) + exp(logD2)
  }
  
  #if (mod=="IWAM_FixedCombined") D2 <- all_Deltas %>% filter(Param=="Delta2_bounded") %>% dplyr::select(Estimate) %>% pull()
  if(length(logD1)==1&length(logD2)==2){
    abline(a=logD1[1], b=exp(logD2[1]), col="forestgreen", lwd=2)
    abline(a=logD1[1], b=exp(logD2[2]), col="dodgerblue3", lwd=2)
  }
  
  if(length(logD1)==2&length(logD2)==2){
    abline(a=logD1[1], b=exp(logD2[1]), col="forestgreen", lwd=2)
    abline(a=logD1[2], b=exp(logD2[2]), col="dodgerblue3", lwd=2)
  }
  
  if(mod=="IWAM_FixedCombined") abline(a=logD1, b=exp(logD2), col="maroon", lwd=2) # Actually pulls Delta2_bounded, so no need to log
  
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|
     mod=="Liermann"|mod=="IWAM_Liermann"|
     mod=="Liermann_PriorRicSig_PriorDeltaSig"|
     mod=="Liermann_HalfNormRicVar_FixedDelta") {
    # abline(a=logD1, b=exp(logD2), col="forestgreen", lwd=2)
    # abline(a=logD1o, b=D2o, col="dodgerblue3", lwd=2)
    
    simWA <- seq(3, 13, 0.5)
    Preds <-  logD1 + simWA*exp(logD2) 
    Predso <- logD1o + simWA*D2o
    lines(x=exp(simWA), y=exp(Preds), col="forestgreen", lwd=2)
    lines(x=exp(simWA), y=exp(Predso), col="dodgerblue3", lwd=2)
  }
  
  if(mod=="IWAM_FixedSep_Constm") {
    abline(a=logD1, b=exp(logD2), col="forestgreen", lwd=2)
    abline(a=logD1o, b=exp(logD2), col="dodgerblue3", lwd=2)
  }
  
  if(mod=="IWAM_FixedSep_Constyi") {
    abline(a=logD1, b=exp(logD2), col="forestgreen", lwd=2)
    abline(a=logD1, b=D2o, col="dodgerblue3", lwd=2)
  }
  
  if(exists("pred_lnSMSY")){
    pred_lnSMSY <- pred_lnSMSY %>% 
      mutate (up = Estimate + 1.96 * Std..Error, lo=Estimate - 1.96*Std..Error) 
    #up_S <- pred_lnSMSY %>% filter(Param== "pred_lnSMSY_S") %>% 
    # dplyr::select(up) %>% pull()
    #lo_S <- pred_lnSMSY %>% filter(Param== "pred_lnSMSY_S") %>% 
    # dplyr::select(lo) %>% pull()
    #up_O <- pred_lnSMSY %>% filter(Param== "pred_lnSMSY_O") %>% 
    # dplyr::select(up) %>% pull()
    #lo_O <- pred_lnSMSY %>% filter(Param== "pred_lnSMSY_O") %>% 
    # dplyr::select(lo) %>% pull()
    up_S <- pred_lnSMSY %>% 
      filter(Param== "pred_lnSMSY_stream_CI") %>% 
      dplyr::select(up) %>% 
      pull()
    lo_S <- pred_lnSMSY %>% 
      filter(Param== "pred_lnSMSY_stream_CI") %>% 
      dplyr::select(lo) %>% 
      pull()
    up_O <- pred_lnSMSY %>% 
      filter(Param== "pred_lnSMSY_ocean_CI") %>% 
      dplyr::select(up) %>% 
      pull()
    lo_O <- pred_lnSMSY %>% 
      filter(Param== "pred_lnSMSY_ocean_CI") %>% 
      dplyr::select(lo) %>% 
      pull()
    up <- pred_lnSMSY %>% 
      filter(Param== "pred_lnSMSY_CI") %>% 
      dplyr::select(up) %>% 
      pull()
    lo <- pred_lnSMSY %>% 
      filter(Param== "pred_lnSMSY_CI") %>% 
      dplyr::select(lo) %>% 
      pull()
    if(is.na(up_S[1])==FALSE) polygon(x=c(exp(pred_lnWA), exp(rev(pred_lnWA))), 
                                      y=c(exp(up_S), exp(rev(lo_S))), 
                                      col=rgb(0,0.4,0, alpha=0.2), border=NA) 
    if(is.na(up_O[1])==FALSE) polygon(x=c(exp(pred_lnWA), exp(rev(pred_lnWA))), 
                                      y=c(exp(up_O), exp(rev(lo_O))), 
                                      col=rgb(0,0.2,0.4, alpha=0.2), border=NA) 
    if(is.na(up[1])==FALSE) polygon(x=c(exp(pred_lnWA), exp(rev(pred_lnWA))), 
                                    y=c(exp(up), exp(rev(lo))), 
                                    col=rgb(0.6,0.2,0.4, alpha=0.2), border=NA)
  }
  
  if(length(logD1)==2&length(logD2)==2){
    text(x=6, y=10.5,labels= paste0("Ocean-type\nlog(Delta1)=",
                                    round(logD1[1],2), ", \nDelta2=", 
                                    round(exp(logD2[1]),2)), col="dodgerblue3", cex=0.8)
    text(x=9, y=6.5, labels= paste0("Stream-type\nlog(Delta1)=",
                                    round(logD1[3],2), ", \nDelta2=", 
                                    round(exp(logD2[2]),2)), col="forestgreen", cex=0.8)
  }
  
  if(length(logD1)==1&length(logD2)==2){
    text(x=6, y=10.5,labels= paste0("log(Delta1)=",
                                    round(logD1[1],2), ", \nDelta2=", 
                                    round(exp(logD2[1]),2)), col="dodgerblue3", cex=0.8)
    text(x=9, y=6.5, labels= paste0("log(Delta1)=",
                                    round(logD1[1],2), ", \nDelta2=", 
                                    round(exp(logD2[2]),2)), col="forestgreen", cex=0.8)
  }
  
  if(mod=="IWAM_FixedCombined"){
    text(x=6, y=10.5,labels= paste0("log(Delta1)=",
                                    round(logD1[1],2), ", \nDelta2=", 
                                    round(exp(logD2[1]),2)), col="maroon", cex=0.8)
  
  }
  
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|
     mod=="Liermann"|mod=="IWAM_Liermann"|
     mod=="Liermann_PriorRicSig_PriorDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta"){
    text(x= 50000, y= 500, labels= paste0("log(Delta1)=",
                                 round(logD1[1],2), ", \nDelta2=",
                                 round(exp(logD2[1]),2)), col="forestgreen", cex=0.8)
    text(x= 500, y= 50000, labels= paste0("log(Delta1)=",
                                   round(logD1o[1],2), ", \nDelta2=",
                                   round(D2o[1],2)), col="dodgerblue3", cex=0.8)
  }
  
  if(mod=="IWAM_FixedSep_Constm"){
    text(x=9, y=7,labels= paste0("log(Delta1)=",
                                 round(logD1[1],2), ", \nDelta2=", 
                                 round(exp(logD2[1]),2)), col="forestgreen", cex=0.8)
    text(x=6, y=9.5,labels= paste0("log(Delta1)=",
                                   round(logD1o[1],2), ", \nDelta2=", 
                                   round(exp(logD2[1]),2)), col="dodgerblue3", cex=0.8)
  }
  
  if(mod=="IWAM_FixedSep_Constyi"){
    text(x=9, y=7,labels= paste0("log(Delta1)=",
                                 round(logD1[1],2), ", \nDelta2=", 
                                 round(exp(logD2[1]),2)), col="forestgreen", cex=0.8)
    text(x=6, y=9.5,labels= paste0("log(Delta1)=",
                                   round(logD1[1],2), ", \nDelta2=", 
                                   round(D2o[1],2)), col="dodgerblue3", cex=0.8)
  }
  
  # Title if you want it
  # title(title1, cex.main=0.9)
  
}

# ------------------------------------------------------------------------------
plotWAregressionSREP <- function (pars, all_Deltas, srdat, lifehist, WAbase,  
                                  pred_lnSREP=NA, pred_lnWA, title1, mod) {
  
  options(scipen = 5) # for ticks without sci. notation
  
  SREP <- pars %>% filter(Param=="SREP") %>% 
    mutate(ModelOrder=0:(length(unique(pars$Stocknumber))-1))
  # what is scale of SREP?
  Sc <- srdat %>% 
    select(Stocknumber, scale) %>% 
    distinct()
  SREP <- SREP %>% 
    left_join(Sc) %>% 
    mutate(rawSREP=Estimate*scale)
  lnSREP <- log(SREP$rawSREP)
  lnWA <- log(WAbase$WA)
  
  par(cex=1.5, mar = c(5, 4, 2, 2) + 0.1)
    # If you want a title - restore to default margins or just take "mar" out
  col.use <- NA
  for(i in 1:length(SREP$lh)) {
    if (SREP$lh[i]==0) col.use[i] <- "forestgreen" 
    else col.use[i] <- "dodgerblue3"
    }
  
  # plot(y=lnSREP, x=lnWA, pch=20, col=col.use, xlab="log(Watershed Area, km2)", ylab="log(SREP)") # ORIGINAL
  # points(y=lnSREP, x=lnWA, pch=20, col=col.use, cex=1.5) # ORIGINAL
  # points(y=lnSREP[18:length(SREP$lh)], x=lnWA[18:length(SREP$lh)], pch=3, col=col.use[18:length(SREP$lh)], cex=1.5)
  
  # LOG LOG
  plot(y=exp(lnSREP), x=exp(lnWA), pch=20, col=col.use, 
       xlab= expression("Accessible Watershed Area, km"^2), 
       ylab=expression(S[REP]), log = 'xy',
       # xlim =c(50,500000) , ylim=c(500,1000000))
       xlim = c(50,200000) , ylim=c(200,2000000))
  points(y=exp(lnSREP), x=exp(lnWA), pch=20, col=col.use, cex=1.5)
    # Works for getting points on there - but lines are not showing properly
    # Axes numbers are in sci-notif - can that be changed?

  logN1 <- all_Deltas %>% 
    filter(Param=="logNu1") %>% 
    dplyr::select(Estimate) %>% 
    pull()
  logN2 <- all_Deltas %>% 
    filter(Param=="logNu2") %>% 
    dplyr::select(Estimate) %>% 
    pull()
  
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|
     mod=="IWAM_FixedSep_Constm"|mod=="Liermann"|
     mod=="IWAM_Liermann"|mod=="Liermann_PriorRicSig_PriorDeltaSig"|
     mod=="Liermann_HalfNormRicVar_FixedDelta") logN1o <- all_Deltas %>% 
    filter(Param=="logNu1_ocean") %>% 
    dplyr::select(Estimate) %>% 
    pull() + logN1
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|
     mod=="IWAM_FixedSep_Constyi"|mod=="Liermann"|
     mod=="IWAM_Liermann"|mod=="Liermann_PriorRicSig_PriorDeltaSig"|
     mod=="Liermann_HalfNormRicVar_FixedDelta") {
    N2o <- exp(all_Deltas %>% 
                 filter(Param=="logNu2_ocean") %>% 
                 dplyr::select(Estimate) %>% 
                 pull() ) + exp(logN2)
    if(nrow(all_Deltas %>% 
            filter(Param=="Nu2_ocean"))>=1)  N2o <- (all_Deltas %>% filter(Param=="Nu2_ocean") %>% 
                                                       dplyr::select(Estimate) %>% 
                                                       pull() ) + exp(logN2)
  }
  
  #if (mod=="IWAM_FixedCombined") D2 <- all_Deltas %>% filter(Param=="Delta2_bounded") %>% dplyr::select(Estimate) %>% pull()
  if(length(logN1)==1&length(logN2)==2){ 
    abline(a=logN1[1], b=exp(logN2[1]), col="forestgreen", lwd=2)
    abline(a=logN1[1], b=exp(logN2[2]), col="dodgerblue3", lwd=2)
  }
  if(length(logN1)==2&length(logN2)==2){
    abline(a=logN1[1], b=exp(logN2[1]), col="forestgreen", lwd=2)
    abline(a=logN1[2], b=exp(logN2[2]), col="dodgerblue3", lwd=2)
  }
  
  if(mod=="IWAM_FixedCombined") abline(a=logN1, b=exp(logN2), col="maroon", lwd=2) #Actually pulls Delta2_bounded, so no need to log
  
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|
     mod=="Liermann"|mod=="IWAM_Liermann"|
     mod=="Liermann_PriorRicSig_PriorDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta") {
    # abline(a=logN1, b=exp(logN2), col="forestgreen", lwd=2)
    # abline(a=logN1o, b=N2o, col="dodgerblue3", lwd=2)
    
    simWA <-  seq(3, 13, 0.5)
    Preds <-  logN1 + simWA*exp(logN2) 
    Predso <- logN1o + simWA*N2o
    lines(x=exp(simWA), y=exp(Preds), col="forestgreen", lwd=2)
    lines(x=exp(simWA), y=exp(Predso), col="dodgerblue3", lwd=2)
  }
  
  if(mod=="IWAM_FixedSep_Constm") {
    abline(a=logN1, b=exp(logN2), col="forestgreen", lwd=2)
    abline(a=logN1o, b=exp(logN2), col="dodgerblue3", lwd=2)
  }
  if(mod=="IWAM_FixedSep_Constyi") {
    abline(a=logN1, b=exp(logN2), col="forestgreen", lwd=2)
    abline(a=logN1, b=N2o, col="dodgerblue3", lwd=2)
  }
  
  if(exists("pred_lnSREP")){
    pred_lnSREP <- pred_lnSREP %>% mutate (up = Estimate + 1.96 * Std..Error, lo=Estimate - 1.96*Std..Error) 
    #up_S <- pred_lnSREP %>% filter(Param== "pred_lnSREP_S") %>% 
    # dplyr::select(up) %>% pull()
    #lo_S <- pred_lnSREP %>% filter(Param== "pred_lnSREP_S") %>% 
    # dplyr::select(lo) %>% pull()
    #up_O <- pred_lnSREP %>% filter(Param== "pred_lnSREP_O") %>% 
    # dplyr::select(up) %>% pull()
    #lo_O <- pred_lnSREP %>% filter(Param== "pred_lnSREP_O") %>% 
    # dplyr::select(lo) %>% pull()
    up_S <- pred_lnSREP %>% 
      filter(Param== "pred_lnSREP_stream_CI") %>% 
      dplyr::select(up) %>% 
      pull()
    lo_S <- pred_lnSREP %>% 
      filter(Param== "pred_lnSREP_stream_CI") %>% 
      dplyr::select(lo) %>% 
      pull()
    up_O <- pred_lnSREP %>% 
      filter(Param== "pred_lnSREP_ocean_CI") %>% 
      dplyr::select(up) %>% 
      pull()
    lo_O <- pred_lnSREP %>% 
      filter(Param== "pred_lnSREP_ocean_CI") %>% 
      dplyr::select(lo) %>%  
      pull()
    up <- pred_lnSREP %>% 
      filter(Param== "pred_lnSREP_CI") %>% 
      dplyr::select(up) %>% 
      pull()
    lo <- pred_lnSREP %>% 
      filter(Param== "pred_lnSREP_CI") %>% 
      dplyr::select(lo) %>% 
      pull()
    
    if(is.na(up_S[1])==FALSE) polygon(x=c(exp(pred_lnWA), exp(rev(pred_lnWA))), 
                                      y=c(exp(up_S), exp(rev(lo_S))), 
                                      col=rgb(0,0.4,0, alpha=0.2), border=NA)
    if(is.na(up_O[1])==FALSE) polygon(x=c(exp(pred_lnWA), exp(rev(pred_lnWA))), 
                                      y=c(exp(up_O), exp(rev(lo_O))), 
                                      col=rgb(0,0.2,0.4, alpha=0.2), border=NA)
    if(is.na(up[1])==FALSE) polygon(x=c(exp(pred_lnWA), exp(rev(pred_lnWA))), 
                                    y=c(exp(up), exp(rev(lo))), 
                                    col=rgb(0.6,0.2,0.4, alpha=0.2), border=NA)
  }
  
  if(length(logN1)==2&length(logN2)==2){
    text(x=6, y=10.5,labels= paste0("Ocean-type\nlog(Nu1)=",
                                    round(logN1[1],2), ", \nNu2=", 
                                    round(exp(logN2[1]),2)), col="dodgerblue3", cex=0.8)
    text(x=9, y=6.5, labels= paste0("Stream-type\nlog(Nu1)=",
                                    round(logN1[3],2), ", \nNu2=", 
                                    round(exp(logN2[2]),2)), col="forestgreen", cex=0.8)
  }
  
  if(length(logN1)==1&length(logN2)==2){
    text(x=6, y=10.5,labels= paste0("log(Nu1)=",
                                    round(logN1[1],2), ", \nNu2=", 
                                    round(exp(logN2[1]),2)), col="dodgerblue3", cex=0.8)
    text(x=9, y=6.5, labels= paste0("log(Nu1)=",
                                    round(logN1[1],2), ", \nNu2=", 
                                    round(exp(logN2[2]),2)), col="forestgreen", cex=0.8)
  }
  
  if(mod=="IWAM_FixedCombined"){
    text(x=6, y=10.5,labels= paste0("log(Nu1)=",
                                    round(logN1[1],2), ", \nNu2=", 
                                    round(exp(logN2[1]),2)), col="maroon", cex=0.8)
  }
  
  if(mod=="IWAM_FixedSep"|mod=="IWAM_FixedSep_RicStd"|
     mod=="Liermann"|mod=="IWAM_Liermann"|
     mod=="Liermann_PriorRicSig_PriorDeltaSig"|mod=="Liermann_HalfNormRicVar_FixedDelta"){
    text(x = 50000, y = 500, labels = paste0("log(Nu1)=",
                                 round(logN1[1],2), ", \nNu2=", 
                                 round(exp(logN2[1]),2)), col="forestgreen", cex=0.8)
    text(x = 500, y = 50000, labels = paste0("log(Nu1)=",
                                    round(logN1o[1],2), ", \nNu2=", 
                                    round(N2o[1],2)), col="dodgerblue3", cex=0.8)
    # text(x= 50000, y= 500)
    # text(x= 500, y= 50000)
    # x = exp(9), y = exp(8)
    # x = exp(6), y = exp(10.5)
  }
  
  if(mod=="IWAM_FixedSep_Constm"){
    text(x=9, y=7,labels= paste0("log(Nu1)=",
                                 round(logN1[1],2), ", \nNu2=", 
                                 round(exp(logN2[1]),2)), col="forestgreen", cex=0.8)
    text(x=6, y=9.5,labels= paste0("log(Nu1)=",
                                   round(logN1o[1],2), ", \nNu2=", 
                                   round(exp(logN2[1]),2)), col="dodgerblue3", cex=0.8)
  }
  
  if(mod=="IWAM_FixedSep_Constyi"){
    text(x = 9, y = 7, labels = paste0("log(Nu1)=",
                                 round(logN1[1],2), ", \nNu2=", 
                                 round(exp(logN2[1]),2)), col="forestgreen", cex=0.8)
    text(x = 6, y = 9.5, labels = paste0("log(Nu1)=",
                                   round(logN1[1],2), ", \nNu2=", 
                                   round(N2o[1],2)), col="dodgerblue3", cex=0.8)
  }
  
  # If you want titles
  # title(title1, cex.main=0.9)
  
}

#-------------------------------------------------------------------------------
# Plot WA regression for SREP with WCVI Chinook stock superimposed
plotWAregressionSREP_withWCVI <- function (All_Est, All_Deltas, SRDat, Stream, 
                                           WA,  PredlnSREP=NA, PredlnWA, 
                                           WCVISREP, title1, mod,
                                           WCVIlnWA_file = 
                                             "DataIn/WCVIStocks.csv") {
  
  options(scipen=5) # This removes the tendancy for ticks with sci. notation.
  SREP <- All_Est %>% filter(Param=="SREP") %>% 
    mutate(ModelOrder=0:(length(unique(All_Est$Stocknumber))-1))
  # Only ocean-type WAs
  SREP.x <- SREP %>% select(Stocknumber, lh)
  WA.x <- left_join(SREP.x,WA, by=c("Stocknumber"))
  WA.x <- WA.x %>% filter(lh==1)
  lnWA <- log(WA.x$WA)
  
  # Only ocean-type SREP:
  
  SREP <- SREP %>% filter(lh==1)
  # what is scale of SMSY?
  Sc <- SRDat %>% dplyr::select(Stocknumber, Scale) %>% distinct()
  SREP <- SREP %>% left_join(Sc, by="Stocknumber") %>% 
    mutate(rawSREP=Estimate*Scale)
  lnSREP <- log(SREP$rawSREP)
  
  par(cex=1.5)
  col.use <- NA
  for(i in 1:length(SREP$lh)) {
    if (SREP$lh[i]==0) col.use[i] <- "forestgreen" 
        else col.use[i] <- "dodgerblue3"
  }
  plot(y = exp(lnSREP), x = exp(lnWA), log="xy", pch = 20, col = col.use, 
       xlab = expression("Accessible watershed area, km"^2), 
       ylab=expression(S[REP]), 
       xlim =c(10,100000) , ylim=c(50,100000))
  
  logD1 <- All_Deltas %>% filter(Param=="logNu1") %>% 
    dplyr::select(Estimate) %>% pull()
  logD2 <- All_Deltas %>% filter(Param=="logNu2") %>% 
    dplyr::select(Estimate) %>% pull()
  if(mod=="Liermann_PriorRicSig_PriorDeltaSig") {
    logD1o <- All_Deltas %>% 
      filter(Param=="logNu1ocean") %>% 
      dplyr::select(Estimate) %>% pull() + logD1
  }
  if(mod=="Liermann_PriorRicSig_PriorDeltaSig") {
    D2o <- exp(All_Deltas %>% filter(Param=="logNu2ocean") %>% 
                 dplyr::select(Estimate) %>% pull() ) + exp(logD2)
    
    if(nrow(All_Deltas %>% filter(Param=="Nu2ocean"))>=1)  {
      D2o <- (All_Deltas %>% filter(Param=="Nu2ocean") %>% 
                dplyr::select(Estimate) %>% pull() ) + exp(logD2) 
    }
  }
  if(mod=="Liermann_PriorRicSig_PriorDeltaSig") {
    simWA <-  seq(2,12,0.5)
    Preds <-  logD1 + simWA*exp(logD2) 
    Predso <- logD1o + simWA*D2o
    # lines(x=exp(simWA), y=exp(Preds), col="forestgreen", lwd=2)
    lines(x=exp(simWA), y=exp(Predso), col="dodgerblue3", lwd=2)
  }
  if(exists("PredlnSREP")){
    PredlnSREP <- PredlnSREP %>% 
      mutate (up = Estimate + 1.96 * Std..Error, lo=Estimate - 1.96*Std..Error) 
    up_S <- PredlnSREP %>% filter(Param== "PredlnSREPs_CI") %>% 
      dplyr::select(up) %>% pull()
    lo_S <- PredlnSREP %>% filter(Param== "PredlnSREPs_CI") %>% 
      dplyr::select(lo) %>% pull()
    up_O <- PredlnSREP %>% filter(Param== "PredlnSREPo_CI") %>% 
      dplyr::select(up) %>% pull()
    lo_O <- PredlnSREP %>% filter(Param== "PredlnSREPo_CI") %>% 
      dplyr::select(lo) %>% pull()
    up <- PredlnSREP %>% filter(Param== "PredlnSREP_CI") %>% 
      dplyr::select(up) %>% pull()
    lo <- PredlnSREP %>% filter(Param== "PredlnSREP_CI") %>% 
      dplyr::select(lo) %>% pull()
    # if(is.na(up_S[1])==FALSE) polygon(x = c(exp(PredlnWA), exp(rev(PredlnWA))), 
    #                                   y = c(exp(up_S), exp(rev(lo_S))), 
    #                                   col = rgb(0,0.4,0, alpha=0.2), border=NA)
    if(is.na(up_O[1])==FALSE) polygon(x = c(exp(PredlnWA), exp(rev(PredlnWA))), 
                                      y = c(exp(up_O), exp(rev(lo_O))), 
                                      col = rgb(0,0.2,0.4, alpha=0.2), 
                                      border=NA)
    if(is.na(up[1])==FALSE) polygon(x = c(PredlnWA, rev(PredlnWA)), 
                                    y = c(up, rev(lo)), 
                                    col = rgb(0.6,0.2,0.4, alpha=0.2), 
                                    border=NA)
  }
  
  if(mod=="Liermann_PriorRicSig_PriorDeltaSig"){
    text(x=9, y=7,labels= 
           paste0( "log(Delta1)=", round(logD1[1],2), ", \nDelta2=", 
                   round(exp(logD2[1]),2)), col="forestgreen", cex=0.8)
    text(x=6, y=9.5,labels= 
           paste0("log(Delta1)=", round(logD1o[1],2), ", \nDelta2=", 
                  round(D2o[1],2)), col="dodgerblue3", cex=0.8)
  }
  
  WCVIlnWA <- read.csv("DataIn/WCVIStocks.csv") %>% mutate (lnWA=log(WA)) %>% 
    filter(lh==1) %>% select(Stock, lnWA) %>% filter(Stock!="Cypre")
  WCVISREP.x <- WCVISREP %>% filter(Param=="SREP") %>% 
    select(Stock, Estimate, LL, UL) %>% filter(Stock!="Cypre")
  WCVIplot.df <- left_join(WCVIlnWA, WCVISREP.x)
  WCVIplot.df <- WCVIplot.df[complete.cases(WCVIplot.df),]
  WCVIplot.df <- arrange(WCVIplot.df, lnWA)
  points(x=exp(WCVIplot.df$lnWA), y=WCVIplot.df$Estimate, pch = 3, cex=0.5)
  # polygon(x = c(exp(WCVIplot.df$lnWA), exp(rev(WCVIplot.df$lnWA))), 
  #         y = c(WCVIplot.df$UL, rev(WCVIplot.df$LL)), 
  #         col = grey(0.3, alpha=0.2), border=NA)
  
  title(title1, cex.main=0.9)
  
}

# -----------------------------------------------------------------------------
# plotWAregression_Parken <- function(data, all_Deltas){
#   par(cex=1.5)
#   col.use <- NA
#   Parken_data <- as.data.frame(read.csv(here::here("DataIn/WA_Parken.csv")))
#   for(i in 1:length(Parken_data$lh)) {if (Parken_data$lh[i]==0) col.use[i] <- "forestgreen" else col.use[i] <- "dodgerblue3"}
#   plot(x=log(data$WAbase), y=log(data$SMSY*data$scale),pch=20, col=col.use, xlab="log(Watershed Area, km2)", ylab="log(SMSY)")
#   points(x=log(data$WAbase), y=log(data$SMSY*data$scale), pch=20, col=col.use, cex=1.5)
#   logD1 <- all_Deltas %>% filter(Param=="logDelta1") %>% dplyr::select(Estimate) %>% pull()
#   D2 <- all_Deltas %>% filter(Param=="Delta2_bounded") %>% dplyr::select(Estimate) %>% pull()
#   abline(a=logD1, b=D2, col="maroon", lwd=2)
#   
#   #abline(lm(log(data$SMSY*data$scale) ~ log(data$WAbase)))
#   #text(x=9, y=8,labels= paste0("log(Delta1)=",round(logD1[1],2), ", \nDelta2=", round(exp(logD2[1]),2)), col="forestgreen", cex=0.8)
#   text(x=6, y=10.5,labels= paste0("log(Delta1)=",round(logD1[1],2), ", \nDelta2=", round(D2[1],2)), col="maroon", cex=0.8)
#   
#   title("Parken Data: Combined")
#   
# }

# plotWAregression_ParkenSep <- function(data, all_Deltas){
#   par(cex=1.5)
#   col.use <- NA
#   Parken_data <- as.data.frame(read.csv(here::here("DataIn/WA_Parken.csv")))
#   for(i in 1:length(Parken_data$lh)) {if (Parken_data$lh[i]==0) col.use[i] <- "forestgreen" else col.use[i] <- "dodgerblue3"}
#   plot(x=log(data$WAbase), y=log(data$SMSY*data$scale),pch=20, col=col.use, xlab="log(Watershed Area, km2)", ylab="log(SMSY)")
#   points(x=log(data$WAbase), y=log(data$SMSY*data$scale), pch=20, col=col.use, cex=1.5)
#   logD1 <- all_Deltas %>% filter(Param=="logDelta1") %>% dplyr::select(Estimate) %>% pull()
#   logD1o <- all_Deltas %>% filter(Param=="logDelta1_ocean") %>% dplyr::select(Estimate) %>% pull() + logD1
#   logD2 <- all_Deltas %>% filter(Param=="logDelta2") %>% dplyr::select(Estimate) %>% pull()
#   D2o <- exp(all_Deltas %>% filter(Param=="logDelta2_ocean") %>% dplyr::select(Estimate) %>% pull() ) + exp(logD2)
#   abline(a=logD1, b=exp(logD2), col="forestgreen", lwd=2)
#   abline(a=logD1o, b=D2o, col="dodgerblue3", lwd=2)
#   
#   #abline(lm(log(data$SMSY*data$scale) ~ log(data$WAbase)))
#   text(x=9, y=7,labels= paste0("log(Delta1)=",round(logD1[1],2), ", \nDelta2=", round(exp(logD2[1]),2)), col="forestgreen", cex=0.8)
#   text(x=6, y=10.5,labels= paste0("log(Delta1)=",round(logD1o[1],2), ", \nDelta2=", round(D2o[1],2)), col="dodgerblue3", cex=0.8)
#   
#   title("Parken Data: Separated")
#   
# }


# Not currently adjusted for IWAM_Liermann * MOST UCRRENT TMB MODEL
plotRicA <- function (){#pars_Liermann, pars_Ricker_AllMod, pars_Liermann_SepRicA=NA){
  #par(mfrow=c(1,1), mar=c(4, 4, 4, 2) + 0.1)
  pars_Liermann <- readRDS("DataOut/pars_Liermann.RDS")
  pars_Ricker_AllMod <- readRDS("DataOut/pars_Ricker_AllMod.RDS")
  pars_Ricker_std <- readRDS("DataOut/pars_Ricker_std_wBC.RDS")
  pars_Ricker_std_noWAreg <- readRDS("DataOut/pars_Ricker_std_noWAreg_wBC.RDS")#From CheckAR1.r
  pars_Liermann_SepRicA <- readRDS("DataOut/pars_Liermann_PriorRicSig_PriorDeltaSig_invGamma0.01_wBC.RDS")#readRDS("DataOut/pars_Liermann_SepRicA.RDS")#r
  pars_Liermann_HalfNormRicVar <-readRDS( "DataOut/pars_Liermann_PriorRicSig_PriorDeltaSig_HalfNormRicVar_wBC.RDS")#readRDS("DataOut/pars_Liermann_HalfNormRicVar.RDS")
  pars_Liermann_HalfCauchyRicVar <- readRDS( "DataOut/pars_Liermann_PriorRicSig_PriorDeltaSig_HalfCauchyRicVar_wBC.RDS")#readRDS("DataOut/pars_Liermann_HalfCauchyRicVar.RDS")
  pars_Liermann_SepRicA_uniformSigmaPrior <- readRDS( "DataOut/pars_Liermann_PriorRicSig_PriorDeltaSig_uniformSigmaPrior_wBC.RDS")#readRDS("DataOut/pars_Liermann_SepRicA_uniformSigmaPrior.RDS")
  pars_Liermann_SepRicA_noSigmaPrior <- readRDS( "DataOut/pars_Liermann_PriorRicSig_PriorDeltaSig_noPriorRicVar_wBC.RDS")# readRDS("DataOut/pars_Liermann_SepRicA_noSigmaPrior.RDS")
  pars_Liermann_SepRicA_invGamma0.001 <- readRDS("DataOut/pars_Liermann_PriorRicSig_PriorDeltaSig_invGamma0.001_wBC.RDS")#readRDS("DataOut/pars_Liermann_SepRicA_invGamma0.001.RDS")
  #pars_Liermann_SepRicA_invGamma0.1 <- readRDS("DataOut/pars_Liermann_SepRicA_invGamma0.1.RDS")
  pars_Liermann_SepRicA_invGamma0.1 <- readRDS("DataOut/pars_Liermann_PriorRicSig_PriorDeltaSig_invGamma0.1_wBC.RDS")
  pars_Liermann_SepRicA_invGamma0.01_invGammaA0.001 <- readRDS("DataOut/pars_Liermann_SepRicA_invGamma0.01_invGammaA0.001.RDS")
  pars_Liermann_SepRicA_invGamma0.001_invGammaA0.01 <- readRDS("DataOut/pars_Liermann_SepRicA_invGamma0.001_invGammaA0.01.RDS")
  
  RicAL <- pars_Liermann %>% filter (Param=="logA")
  RicALsep <- pars_Liermann_SepRicA %>% filter (Param=="logA")
  RicALhNRV <- pars_Liermann_HalfNormRicVar %>% filter (Param=="logA")
  RicALhCRV <- pars_Liermann_HalfCauchyRicVar %>% filter (Param=="logA")
  RicALsep_uni <- pars_Liermann_SepRicA_uniformSigmaPrior %>% filter (Param=="logA")
  RicALsep_0.1 <- pars_Liermann_SepRicA_invGamma0.1 %>% filter (Param=="logA")
  RicALsep_0.001 <- pars_Liermann_SepRicA_invGamma0.001 %>% filter (Param=="logA")
  RicALsep_0.001_0.01 <- pars_Liermann_SepRicA_invGamma0.001_invGammaA0.01 %>% filter (Param=="logA")
  RicALsep_0.01_0.001 <- pars_Liermann_SepRicA_invGamma0.01_invGammaA0.001 %>% filter (Param=="logA")
  RicALsep_none <- pars_Liermann_SepRicA_noSigmaPrior %>% filter (Param=="logA")
  RicAR <- pars_Ricker_AllMod %>% filter (Param=="logA")
  RicARstd <- pars_Ricker_std %>% filter (Param=="logA")#(Param=="logA_std")
  RicARstd_noWAreg <- pars_Ricker_std_noWAreg %>% filter (Param=="logA")#(Param=="logA_std")
  nRicAL <- length(RicALsep$Estimate)
  nRicAstd <- length(RicARstd$Estimate)

  
  #box.data <- data.frame(RiclogA=c(RicAL$Estimate,RicALsep$Estimate, RicAR$Estimate), Model=c(rep("Hierarchical_logA_combined",nRicA),rep("Hierarchical_logA_sep",nRicA), rep("Fixed_logA",nRicA)))
  box.data <- data.frame(RiclogA = c(RicARstd_noWAreg$Estimate, RicARstd$Estimate, RicALsep_0.1$Estimate, RicALsep$Estimate, RicALsep_0.001$Estimate, 
                                   RicALhNRV$Estimate, RicALhCRV$Estimate, RicALsep_uni$Estimate, RicALsep_none$Estimate), 
                         Model = c(rep("Fixed\neffects\nNoWAreg",nRicAstd), rep("Fixed\neffects",nRicAstd), rep("InvGamma\n(0.1, 0.1)", nRicAL), 
                                 rep("InvGamma\n(0.01,0.01)",nRicAL), rep("Inv\nGamma\n(0.001, 0.001)", nRicAL),
                                 rep("Half\nNormal\n(0,1)",nRicAL), 
                                 rep("Half\nCauchy\n(0,1)",nRicAL),  
                                 rep("Uniform\n(0,2)",nRicAL), rep("No priors\nSigmas",nRicAL)), 
                         lh = c(RicARstd$lh, RicARstd$lh, RicALsep_0.1$lh, RicALsep$lh, RicALsep_0.001$lh, 
                                      RicALhNRV$lh, RicALhCRV$lh, RicALsep_uni$lh, RicALsep_none$lh) )
  
  order <- data.frame(Model=unique(box.data$Model), order=1:9)
  box.data <- box.data %>% left_join(order, by="Model")
  #Add reorder to as.factor(Model), I think, but need to specify order as above 1:8. Add column to dataframe to do this, of length 200.
  # see my old dplyr code on adding new columns with string of numbers aligned with a factor (IWAM.r)
  cols<-viridis(4, alpha=0.9)
  cols.order <- c(grey(0.2), grey(0.5), cols[3], cols[2], t_col(color=cols[1], percent=70), t_col(cols[1], percent=50), t_col(cols[1], percent=30), "white", cols[4])
  #print cols.order and then cut and past into line: scale_fill_manual below. These are the colours I used in the plot of priors. can add a legend to quickly see how ggplot reorders this befor plotting
  
  ggplot(box.data, aes(x=reorder(as.factor(Model),order), y=RiclogA, fill=as.factor(Model))) +  
    geom_boxplot(outlier.size=-1) + 
    scale_fill_manual(values=c("#808080", "#808080", "#35B779E6", "#31688EE6", "#4401544C", "#4401547F", "#440154B2", "white", "#FDE725E6" )) +
    #scale_fill_viridis(discrete = TRUE, alpha=0.6) + 
    geom_jitter(color="black", shape=as.factor(box.data$lh), size=2, alpha=0.9) + #,aes(shape=lh)?. I like shape 16,17
    scale_shape_manual(values=c(16,17)) + #Does not work 
    ggtitle("Distribution of Ricker logA") +
    theme(legend.position="none") +
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=20,face="bold"),
          plot.title = element_text(size = 20)) + 
    xlab("Model")
  #boxplot(list(Hierarchical=RicAL$Estimate, Fixed=RicAR$Estimate), ylab="Ricker log(a)", outline=FALSE, col="forestgreen")
  #points(x=jitter(rep(1,nRicA)), y=RicAL$Estimate, pch=20)
  #points(x=jitter(rep(2,nRicA)), y=RicAR$Estimate, pch=20)
}

# png(paste("DataOut/RicADist_ComparePriors_wBC.png", sep=""), width=7, height=7, units="in", res=500)
# plotRicA()
# dev.off()

plotRicA_reduc <- function (){
  # New read-ins
    # TK: These are no broken due to the changes in naming conventions by IWAM_model.R
  All_Est_Fixed <- readRDS("DataOut/pars_techreport_ricgamma_0.1_wagamma_1_IWAM_FixedEffects")
  # pars_techreport_ricgamma_0.1_wagamma_1_IWAM_FixedEffects - new
  # pars_ricgamma_0.1_wagamma_1_IWAM_FixedEffects - original
  RicFixed <- All_Est_Fixed %>% 
    filter(Param=="logA")
  
  # Inv Gamma 0.1
  All_Est_Default <- readRDS("DataOut/pars_ricgamma_0.1_wagamma_1_wagamma_IWAM_Liermann")
  RicDefault <-  All_Est_Default %>%  
    filter(Param=="logA")
  
  # Halfcauchy
  All_Est_riccauchy <- readRDS("DataOut/pars_riccauchy_0.1_wagamma_1_wagamma_IWAM_Liermann")
  Riccauchy <-  All_Est_riccauchy %>%  
    filter(Param=="logA")
  # HalfNormal
  All_Est_richalfnorm <- readRDS("DataOut/pars_richalfnorm_0.1_wagamma_1_wagamma_IWAM_Liermann")
  Richalfnorm <-  All_Est_richalfnorm %>%  
    filter(Param=="logA")
  # Inv Gamma 0.001
  All_Est_ricg0.001 <- readRDS("DataOut/pars_ricgamma_0.001_wagamma_1_wagamma_IWAM_Liermann")
  Ric0.001 <-  All_Est_ricg0.001 %>%  
    filter(Param=="logA")
  # Inv Gamma 0.01
  All_Est_ricg0.01 <- readRDS("DataOut/pars_ricgamma_0.01_wagamma_1_wagamma_IWAM_Liermann")
  Ric0.01 <-  All_Est_ricg0.01 %>%  
    filter(Param=="logA")
  
  # Unused - specific to WA regression priors
  # All_Est_wag0.01 <- readRDS("DataOut/pars_ricgamma_0.1_wagamma_0.01_wagamma_IWAM_Liermann")
  # WA0.01 <-  All_Est_wag0.01 %>%  filter(Param=="logA")
  # All_Est_wag0.1 <- readRDS("DataOut/pars_ricgamma_0.1_wagamma_0.1_wagamma_IWAM_Liermann")
  # WA0.01 <-  All_Est_wag0.1 %>%  filter(Param=="logA")
  # 
  # All_Est_wacauchy <- readRDS("DataOut/pars_ricgamma_wacauchy_IWAM_Liermann")
  # WAcauchy <-  All_Est_wacauchy %>%  filter(Param=="logA")
  # All_Est_wahalfnorm <- readRDS("DataOut/pars_ricgamma_wahalfnorm_IWAM_Liermann")
  # WAhalfnorm <- All_Est_wahalfnorm %>%  filter(Param=="logA")
  
  # Calls and transforms
  # RicAL <- All_Est_Liermann %>% filter (Param=="logA")
  # nRicAL <- length(RicALsep$Estimate)
  nRicAL <- length(RicFixed$Estimate) # 25 stocks - from synoptic set?
  # nRicAstd <- length(RicARstd$Estimate)
  
  # Labelling:
      # Fixed effects: nRicAstd
      # InvGamma: nRicAL (all 3 variants)
      # Half normal + Half cauchy: nRicAL (why these also?)
  box.data <- data.frame(RicLogA = c(RicFixed$Estimate, 
                                     RicDefault$Estimate, 
                                     Ric0.01$Estimate, 
                                     Ric0.001$Estimate, 
                                     Richalfnorm$Estimate, 
                                     Riccauchy$Estimate), 
                         Model = c(rep("Single\nPopulation\nModel",nRicAL), # Originally Fixed\neffects
                                   rep("InvGamma\n(0.1, 0.1)", nRicAL), 
                                   rep("InvGamma\n(0.01,0.01)",nRicAL), 
                                   rep("InvGamma\n(0.001, 0.001)", nRicAL),
                                   rep("Half\nNormal\n(0,1)",nRicAL), 
                                   rep("Half\nCauchy\n(0,1)",nRicAL)), 
                         lh = c(RicFixed$lh,
                                RicFixed$lh,
                                RicFixed$lh,
                                RicFixed$lh,
                                RicFixed$lh,
                                RicFixed$lh)
                         )
  
  # order <- data.frame(Model=unique(box.data$Model), order=1:5)
  # box.data <- box.data %>% left_join(order, by="Model")
  
  #Add reorder to as.factor(Model), I think, but need to specify order as above 1:8. Add column to dataframe to do this, of length 200.
  # see my old dplyr code on adding new columns with string of numbers aligned with a factor (IWAM.r)
  cols<-viridis(4, alpha=0.9)
  cols.order <- c(grey(0.2), grey(0.5), cols[3], cols[2], 
                  t_col(color=cols[1], percent=70), 
                  t_col(cols[1], percent=50), 
                  t_col(cols[1], percent=30), "white", cols[4])
  #print cols.order and then cut and past into line: scale_fill_manual below. These are the colours I used in the plot of priors. can add a legend to quickly see how ggplot reorders this befor plotting
  
  # fmodel <- function(x) { # used to quickly ifll model - could also do it internally
  #   factor(x, levels = unique(x))
  # }
  
  ggplot(box.data, aes(x=factor(Model, levels = unique(Model)), y=RicLogA, fill=factor(Model, levels = unique(Model)))) +  
    geom_boxplot(outlier.size=-1) + 
    scale_fill_manual(values=c("#808080", "#35B779E6", "#31688EE6", 
                               "#4401544C", "#4401547F", "#440154B2")) +
    geom_jitter(color="black", shape=as.factor(box.data$lh), size=2, alpha=0.9) + #,aes(shape=lh)?. I like shape 16,17
    scale_shape_manual(values=c(16,17)) +
    # ggtitle("Distribution of Ricker LogA") + # Remove title
    theme(legend.position="none") +
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=20,face="bold"),
          plot.title = element_text(size = 20)) + 
    xlab("") + # blank 
    ylab("Productivity") 
}

# png(paste("DataOut/RicADist_ComparePriors_wBC.png", sep=""), width=7, height=7, units="in", res=500)
# plotRicA()
# dev.off()

plotTestStocks <- function(data = ParkenTestSMSY){
  
  pd <- position_dodge(0.4)
  #data$Stock <- as.factor(data$Stock)
  ggplot(data, aes(x=Stock, y=SMSY, group=Source, colour=Source)) + 
    geom_point(position=pd, size=1.5) +
    geom_errorbar(aes(ymin=LL, ymax=UL), width=0.1, size=1, position=pd) + 
    coord_cartesian(ylim = c(0, 25000)) +
    #ylim(0,51000) + 
    theme(axis.text = element_text(angle = 90, vjust = 0.5, hjust=1, size=12) ) + 
    xlab("") + ylab("SMSY") + 
    theme(
      #plot.title = element_text(color="black", size=14, face="bold.italic"),
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14))
  

}
 
#South Thompson values are > y-axis upper limit
#TestSMSY %>% filter(Stock=="S.Thompson" ) 
# For Inv Gammma (0.1, 0.1) Ric Sig and InvGamma(1,1) deltaSig
#Stock     LL    UL    SMSY       Source
# S.Thompson 31092.91 126981.7 62834.95    IWAM
# S.Thompson 37300.00  87000.0 57600.00 Parken

# With bias correction invGamm0.1 for Rsig and inGamm1 deltaSig
#Stock       LL       UL     SMSY Source
#TestlnSMSY_ocean.9 S.Thompson 42839.86 150375.6 80262.51   IWAM

# png(paste("DataOut/TestStock_PI_wBC.png", sep=""), width=9, height=7, units="in", res=500)
# plotTestStocks()
# dev.off()

plotSMSY <- function(data = WCVISMSY){
  
  data <- data %>% filter(Param == "SMSY") %>% rename(SMSY = Estimate)
  #model order: 2  6 11 12 13 14 15 16 17 20  1  3  5  7  8  9 18 19 21  4 10
  order <- c(2,  6, 11, 12, 13, 14, 15, 16, 17, 20 , 1 , 3 , 5,  7,  8 , 9, 18 ,19, 21,  4, 10)
  #order <- c("Bedwell/Ursus","Cypre", "Megin","Moyeha", "Nahmint", "Nitinat", "San Juan", "Sarita", "Somass","Tranquil", "Artlish",
  #           "Burman","Conuma","Gold","Kaouk", "Leiner","Tahsis", "Tahsish","Zeballos","Cayeghle","Marble")
  #ggplot(data, aes(x=reorder(as.factor(Stock),order), y=SMSY, colour=CU)) + 
  ggplot(data, aes(x=Stock, y=SMSY, colour=CU)) + 
    geom_point(size=2) +
    geom_errorbar(aes(ymin=LL, ymax=UL), width=0.1, size=1) + 
    #coord_cartesian(ylim = c(0, 5000)) +
    #ylim(0,51000) + 
    theme(axis.text = element_text(angle = 90, vjust = 0.5, hjust=1, size=12) ) + 
    xlab("") + ylab("SMSY") + 
    theme(
      #plot.title = element_text(color="black", size=14, face="bold.italic"),
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14))
}

# png(paste("DataOut/WCVI_MinFixedDeltaPI.png", sep=""), width=9, height=7, units="in", res=500)
# plotSMSY()
# dev.off()

#========================================================== #
# Get Quantiles of escapement to plot against benchmarks ####
# 
# WCVIEsc <- as.data.frame(read.csv("DataIn/WCVIEsc.csv", row.names="Yr"))
# WCVI_RPs <- as.data.frame(read.csv("DataOut/WCVI_SMSY.csv")) %>% select (-X)
# 
# WCVIEscQuant <- data.frame(t(data.frame (apply(WCVIEsc, 2, quantile, probs = c(0.025,0.25, 0.5,0.75, 0.975), na.rm=T))))
# WCVIEscQuant$Stock <- row.names(WCVIEscQuant)
# WCVIEscQuant$Stock <- sapply(WCVIEscQuant$Stock, function(x) (gsub(".", " ", x, fixed=TRUE) ) )
# WCVIEscQuant$Stock <- sapply(WCVIEscQuant$Stock, function(x) (gsub("Bedwell Ursus", "Bedwell/Ursus", x, fixed=TRUE) ) )
# WCVIEscQuant$Stock <- as.character(WCVIEscQuant$Stock)
# stockCU <- WCVI_RPs %>% filter (Param == "SMSY") %>% select (Stock)#, CU)
# stockCU$Stock <- as.character(stockCU$Stock)
# 
# WCVIEscQuant <- WCVIEscQuant %>%  left_join(stockCU, by = "Stock") 
# WCVIEscQuant <- WCVIEscQuant %>% rename(Estimate=X50., LL=X25., UL=X75.) %>% select (-X2.5., -X97.5.) %>% add_column(Param="Quant")
# 
# WCVIEscQuant <- WCVIEscQuant %>% bind_rows(WCVI_RPs)
# WCVIEscQuant <- WCVIEscQuant %>% mutate(Estimate = round(Estimate, 0), LL = round(LL, 0), UL = round(UL, 0))
#write.csv(WCVIEscQuant, "DataOut/WCVIbenchmarks.csv")
#read.csv("DataOut/WCVIbenchmarks.csv")

# WCVI Specific ================================================================

plotWCVIBenchmarks <- function(data = WCVIEscQuant){
  
  pd <- position_dodge(0.4)
  
  ggplot(data, aes(x=Stock, y=Estimate, group=Param, colour=Param)) + 
    geom_point(position=pd, size=2) +
    geom_errorbar(aes(ymin=LL, ymax=UL), position=pd, width=0.1, size=1) + 
    #coord_cartesian(ylim = c(0, 5000)) +
    #ylim(0,51000) + 
    theme(axis.text = element_text(angle = 90, vjust = 0.5, hjust=1, size=12) ) + 
    xlab("") + ylab("Benchmarks") + 
    theme(
      #plot.title = element_text(color="black", size=14, face="bold.italic"),
      axis.title.x = element_text(size=14),
      axis.title.y = element_text(size=14))
}

#png(paste("DataOut/WCVIbenchmarks.png", sep=""), width=9, height=7, units="in", res=500)
#plotWCVIBenchmarks()
#dev.off()

# If you desire timeseries - make sure to update the below code with the CORRECT bootstrapped estimate csv
plotWCVI_timeseries <- function(WCVIEsc=WCVIEsc, remove.EnhStocks=TRUE, 
                                prod="LifeStageModel", boot=NULL){
  
 
  Years <- rownames(WCVIEsc)
  
  #EnhStocks <- c("Burman",  "Conuma", "Leiner", "Nitinat", "Sarita",  "Somass",  
  # "Zeballos", "San Juan")
  EnhStocks <- data.frame(read.csv(here::here("DataIn/WCVIStocks.csv"))) %>% 
    filter (Enh==1) %>% pull(Stock)
  EnhStocks <- as.character(EnhStocks)
  
  # remove Artlish: 23 Dec 2020
  
  if(prod=="LifeStageModel"){
    if(remove.EnhStocks) WCVI_RPs <- 
        as.data.frame(read.csv(here::here("DataOut/wcviRPs_noEnh.csv"))) %>% 
        rename(StockNames=Stock)  %>% filter (StockNames != "Little Zeballos")
    
    if(!remove.EnhStocks) WCVI_RPs <- 
        as.data.frame(read.csv(here::here("DataOut/wcviRPs_wEnh.csv"))) %>% 
        rename(StockNames=Stock)  %>% filter (StockNames != "Little Zeballos")
  }#End of if(prod=="LifeStageModel")
  
  if(prod=="RunReconstruction"){
    if(remove.EnhStocks) WCVI_RPs <- 
        as.data.frame(read.csv(here::here("DataOut/wcviRPs_noEnh_lowSgen.csv"))) %>% 
        rename(StockNames=Stock)  %>% filter (StockNames != "Little Zeballos")
    if(!remove.EnhStocks) print("Error. wcviRPs_wEnh_lowSgen.csv does not exist. 
                                Need to generate it for with Enhancement 
                                scenario")
  }# End of if(prod=="RunReconstruction")
 
 if (remove.EnhStocks) WCVI_RPs <- WCVI_RPs %>% filter(StockNames %not in% c(EnhStocks)) 
  
  StockNames <- row.names(t(WCVIEsc))
  StockNames <- sapply(StockNames, function(x) (gsub(".", " ", x, fixed=TRUE) ) )
  StockNames <- sapply(StockNames, function(x) (gsub("Bedwell Ursus", "Bedwell/Ursus", x, fixed=TRUE) ) )
  StockNames <- sapply(StockNames, function(x) (gsub("Nootka Esperanza", "Nootka/Esperanza", x, fixed=TRUE) ) )
  StockNames <- data.frame(StockNames) %>% filter(StockNames !="Little Zeballos")
  
  # Reorder in WCVI_RPs to be aligned with stock order in escapement data (=order in StockNames df)
  WCVI_rps <- suppressWarnings(StockNames %>% full_join( WCVI_RPs, by = "StockNames") %>% dplyr::select (-X))
    
  col.pal <- viridis(4)
  col.pal.light <- viridis(4, alpha=0.2)
  ymax <- NA
  N_stk <- length(StockNames$StockNames) - 5 #omit the 5 inlets
  
  # par(mfrow=c(7,3), mar = c(2, 3, 2, 1) + 0.1)
  # for (i in 1:N_stk){
  #   ymax[i] <- max(WCVI_SREP$Estimate[i], WCVIEsc[,i], na.rm=T)#WCVI_sgen$SGEN[i]*4
  #   plot(x= Years, y= WCVIEsc[,i], type="l", lwd=3, col=col.pal[1], ylab="", xlab="", las=0, xaxt="n", bty="n", ylim=c(0,ymax[i]))
  #   axis(side=1, tick=FALSE, padj=-1.8)
  #   abline(h=0)
  #   axis(side=2)
  #   mtext(StockNames$StockNames[i], side=3, line=0.5, at="1960")
  #   legend(x="topleft", legend=NA, title= paste( "   Sgen=", WCVI_sgen$SGEN[i] ), bty="n" )
  #   abline(h=WCVI_SMSY$Estimate[i], col=col.pal[4], lwd=2)
  #   abline(h=WCVI_SREP$Estimate[i], col=col.pal[3], lwd=2)
  #   abline(h=WCVI_sgen$SGEN[i], col=col.pal[2], lwd=2)
  #   polygon(x=c(range(Years), rev(range(Years))), y=c(rep(WCVI_SMSY$LL[i],2), rep(WCVI_SMSY$UL[i],2)), col=col.pal.light[4], border=NA)#col.pal[4])
  #   polygon(x=c(range(Years), rev(range(Years))), y=c(rep(WCVI_SREP$LL[i],2), rep(WCVI_SREP$UL[i],2)), col=col.pal.light[3], border=NA)#col.pal[3])
  # }
  
  N_inlets <- 5 + N_stk

  if (remove.EnhStocks) par(mfrow=c(3,2), mar = c(2, 3, 2, 1) + 0.1)
  if (!remove.EnhStocks) par(mfrow=c(4,2), mar = c(2, 3, 2, 1) + 0.1)

  # For PPT presentation
  ppt <- FALSE
  if (ppt) par(mfrow=c(5,1), mar = c(2, 3, 2, 1) + 0.1)
  
  # Plot in order of inlets from south to north
  if(!remove.EnhStocks) inlet.list <- c(which(StockNames == "San Juan"), which(StockNames == "Nitinat"),
                                        which(StockNames == "Barkley"), which(StockNames == "Clayoquot"), 
                                        which(StockNames == "Nootka/Esperanza"), which(StockNames == "Kyuquot"),
                                        which(StockNames == "Quatsino") )
  if(remove.EnhStocks) inlet.list <- c(which(StockNames == "Barkley"), which(StockNames == "Clayoquot"), 
                                       which(StockNames == "Nootka/Esperanza"), which(StockNames == "Kyuquot"),
                                       which(StockNames == "Quatsino") ) 
  if(ppt) inlet.list <- rev(inlet.list)

    for (i in c(inlet.list)){
    ymax[i] <- max(WCVI_rps$SREP[i], WCVIEsc[,i], na.rm=T)
    plot(x= Years, y= WCVIEsc[,i], type="l", lwd=3, col=col.pal[1], ylab="", xlab="", las=0, xaxt="n", bty="n", ylim=c(0,ymax[i]))
    axis(side=1, tick=FALSE, padj=-1.8)
    abline(h=0)
    axis(side=2)
    mtext(WCVI_rps$StockNames[i], side=3, line=0.5, at="1960")
    if(is.null(boot)) legend(x="topleft", legend=NA, title= paste( "   Sgen=", WCVI_rps$SGEN[i] ), bty="n" )
    if(!is.null(boot)) {if(is.list(boot)) abline(h=WCVI_rps$SMSY[i], col=col.pal[4], lwd=2) }
    if(is.null(boot)) {abline(h=WCVI_rps$SMSY[i], col=col.pal[4], lwd=2) }
    abline(h=WCVI_rps$SREP[i], col=col.pal[3], lwd=2)
    abline(h=WCVI_rps$SGEN[i], col=col.pal[2], lwd=2)
    # CIs for SMSY from WA-model with previous assumption about productivity
    #polygon(x=c(range(Years), rev(range(Years))), y=c(rep(WCVI_rps$SMSYLL[i],2), rep(WCVI_rps$SMSYUL[i],2)), col=col.pal.light[4], border=NA)#col.pal[4])
    polygon(x=c(range(Years), rev(range(Years))), y=c(rep(WCVI_rps$SREPLL[i],2), rep(WCVI_rps$SREPUL[i],2)), col=col.pal.light[3], border=NA)#col.pal[3])
   
      if(!is.null(boot)){
        if(is.list(boot)){
          SMSY.U <- boot$SMSY.boot %>% add_column(stock=row.names(boot$SMSY.boot)) %>% filter(stock==    WCVI_rps$StockNames[i]) %>% pull(upr)
          SMSY.L <- boot$SMSY.boot %>% add_column(stock=row.names(boot$SMSY.boot)) %>% filter(stock==    WCVI_rps$StockNames[i]) %>% pull(lwr)
          polygon(x=c(range(Years), rev(range(Years))), y=c(rep(SMSY.L,2), rep(SMSY.U,2)), col=col.pal.light[4], border=NA)#col.pal[3])
          sgen.U <- boot$SGEN.boot %>% add_column(stock=row.names(boot$SGEN.boot)) %>% filter(stock==    WCVI_rps$StockNames[i]) %>% pull(upr)
          sgen.L <- boot$SGEN.boot %>% add_column(stock=row.names(boot$SGEN.boot)) %>% filter(stock==    WCVI_rps$StockNames[i]) %>% pull(lwr)
          polygon(x=c(range(Years), rev(range(Years))), y=c(rep(sgen.L,2), rep(sgen.U,2)), col=col.pal.light[2], border=NA)#col.pal[3])
          SREP.U <- boot$SREP.boot %>% add_column(stock=row.names(boot$SREP.boot)) %>% filter(stock==    WCVI_rps$StockNames[i]) %>% pull(upr)
          SREP.L <- boot$SREP.boot %>% add_column(stock=row.names(boot$SREP.boot)) %>% filter(stock==    WCVI_rps$StockNames[i]) %>% pull(lwr)
          polygon(x=c(range(Years), rev(range(Years))), y=c(rep(SREP.L,2), rep(SREP.U,2)), col=col.pal.light[3], border=NA)#col.pal[3])
        }
        if(!is.list(boot)){
          # SMSY.U <- read.csv("DataOut/wcviCK-BootstrappedRPs.csv") %>% filter (RP=="SMSY") %>% filter(stock==    WCVI_rps$StockNames[i]) %>% pull(upr)
          # SMSY.L <- boot$SMSY.boot %>% add_column(stock=row.names(boot$SMSY.boot)) %>% filter(stock==    WCVI_rps$StockNames[i]) %>% pull(lwr)
          # polygon(x=c(range(Years), rev(range(Years))), y=c(rep(SMSY.L,2), rep(SMSY.U,2)), col=col.pal.light[4], border=NA)#col.pal[3])
          sgen.U <- read.csv("DataOut/wcviCK-BootstrappedRPs.csv") %>% filter (RP=="SGEN") %>% filter(Stock== WCVI_rps$StockNames[i]) %>% pull(upr)
          sgen.L <- read.csv("DataOut/wcviCK-BootstrappedRPs.csv") %>% filter (RP=="SGEN") %>% filter(Stock== WCVI_rps$StockNames[i]) %>% pull(lwr)
          polygon(x=c(range(Years), rev(range(Years))), y=c(rep(sgen.L,2), rep(sgen.U,2)), col=col.pal.light[2], border=NA)#col.pal[3])
          SREP.U <- read.csv("DataOut/wcviCK-BootstrappedRPs.csv") %>% filter (RP=="SREP") %>% filter(Stock== WCVI_rps$StockNames[i]) %>% pull(upr)
          SREP.L <- read.csv("DataOut/wcviCK-BootstrappedRPs.csv") %>% filter (RP=="SREP") %>% filter(Stock== WCVI_rps$StockNames[i]) %>%  pull(lwr)
          polygon(x=c(range(Years), rev(range(Years))), y=c(rep(SREP.L,2), rep(SREP.U,2)), col=col.pal.light[3], border=NA)#col.pal[3])
          
        }
       
       #WHY IS BOOTSTRPPED interval for SMSY so much higher than the point estimate? Am I using the correct mean/sig to bootstrap SMSY (wBC), ie.., comparing apples to apples?
       
     }
  }
  
  
}

# xx <- Get.LRP(remove.EnhStocks=TRUE, Bern_logistic=TRUE,
#               prod="LifeStageModel")
# png(paste("DataOut/WCVI_inlet_timeseries_nEnh_bs.png", sep=""), width=9, height=7, units="in", res=500)
# plotWCVI_timeseries(WCVIEsc=xx$WCVIEsc, remove.EnhStocks = TRUE,
#                     prod="LifeStageModel", boot=TRUE)
# dev.off()

## yy <- Get.LRP(remove.EnhStocks=FALSE)
# png(paste("DataOut/WCVI_inlet_timeseries_wEnh.png", sep=""), width=9, height=7, units="in", res=500)
# plotWCVI_timeseries(WCVIEsc=yy$WCVIEsc, remove.EnhStocks = FALSE) 
# dev.off()

# xx <- Get.LRP(remove.EnhStocks=TRUE)
# png(paste("DataOut/WCVI_inlet_timeseries_noEnh.png", sep=""), width=6, height=9, units="in", res=500)
# plotWCVI_timeseries(WCVIEsc=xx$WCVIEsc, remove.EnhStocks = TRUE)
# dev.off()

#png(paste("DataOut/WCVItimeseries%02d.png", sep=""), width=9, height=7, units="in", res=500)
# png(paste("DataOut/WCVItimeseries.png", sep=""), width=9, height=7, units="in", res=500)
# plotWCVI_timeseries(WCVIEsc) 
# dev.off()





plotWCVI_SMUtimeseries <- function(SMU_Esc=SMU_Esc, out=out$LRP, WCVI_Esc=WCVIEsc){
  
  Years <- rownames(WCVI_Esc)


  col.pal <- viridis(4)
  col.pal.light <- viridis(4, alpha=0.2)
  ymax <- NA
  #N_stk <- length(StockNames$StockNames) - 5 #omit the 5 inlets
  #N_stk <- length(StockNames) - 5 #omit the 5 inlets
  

  par(mfrow=c(1,1), mar = c(2, 4, 2, 1) + 0.1, cex.lab=1.3)

  ymax <- 20000#max(out$upr, out$fit, SMU_Esc, na.rm=T)#55000
  plot(x= Years, y= SMU_Esc, type="l", lwd=3, col=col.pal[1], ylab="Aggregate Spawner abundances", xlab="", las=0, xaxt="n", bty="n", ylim=c(0,ymax))
  points(x= Years, y= SMU_Esc, col=col.pal[1], pch=19)
  axis(side=1, tick=TRUE, pos=0, padj=-0.9)#-1.8
  abline(h=0)
  axis(side=2)
  mtext("WCVI SMU", side=3, line=0.5, at="1960", cex=1.5)
  #legend(x="topleft", legend=NA, title= paste( "   LRP=", signif(out$fit,4) ), bty="n" )
  abline(h=out$fit, col=col.pal[2], lwd=2)
  # #Projected LRP from "ProjectedLRPs.csv"
  # abline(h=projLRPa, col=col.pal[4], lwd=2)#abline(h=projLRP$fit, col=col.pal[4], lwd=2)
  # CIs for LRP from TMB logistic regression
  #polygon(x=as.numeric(c(range(Years), rev(range(Years)))), y=c(rep(out$lwr,2), rep(out$upr,2)), col=col.pal.light[3], border=NA)
}

 # xx <- Get.LRP(remove.EnhStocks=TRUE)
# png(paste("DataOut/WCVI_SMUtimeseries_noEnh.png", sep=""), width=9, height=4, units="in", res=500)
# plotWCVI_SMUtimeseries(SMU_Esc=xx$SMU_Esc, out=xx$out, WCVI_Esc=xx$WCVIEsc)
# dev.off()

# dfOut <- data.frame(SiteEsc=xx$SMU_Esc, Year=as.numeric(rownames(xx$WCVIEsc)))
# write.csv(dfOut, "DataOut/SMU_Sum.csv")

## yy <- Get.LRP(remove.EnhStocks=FALSE)
# png(paste("DataOut/WCVI_SMUtimeseries_wEnh.png", sep=""), width=9, height=4, units="in", res=500)
# plotWCVI_SMUtimeseries(SMU_Esc=yy$SMU_Esc, out=yy$out, WCVI_Esc=yy$WCVIEsc)
# dev.off()

# # Code to plot projection based LRPs
# # Need to run the function lines indepently to plot both empirical and projection LRPs
#   projLRPa <- data.frame(read.csv("c:/github/SalmonLRP_RetroEval/WCVIChinookStudy/DataOut/ProjectedLRPs/ProjectedLRPscvER0.21baseERn10000_ALLp.csv")) %>%
#     filter(ProbThresh=="0.5") %>% pull(LRP)
#   projLRPb <- data.frame(read.csv("c:/github/SalmonLRP_RetroEval/WCVIChinookStudy/DataOut/ProjectedLRPs/ProjectedLRPscvER0.42_ALLp.csv")) %>%
#     filter(ProbThresh=="0.5") %>% pull(LRP)
#   projLRPc <- data.frame(read.csv("c:/github/SalmonLRP_RetroEval/WCVIChinookStudy/DataOut/ProjectedLRPs/ProjectedLRPscvER0_Allp.csv")) %>%
#     filter(ProbThresh=="0.5") %>% pull(LRP)
# 
# projLRPAllp <- data.frame(read.csv("c:/github/SalmonLRP_RetroEval/WCVIChinookStudy/DataOut/ProjectedLRPs/ProjectedLRPsbaseER_ALLp.csv"))
# projLRPg <- projLRPAllp %>% filter(ProbThresh=="0.5") %>% pull(LRP)
# projLRPh <- projLRPAllp %>% filter(ProbThresh=="0.66") %>% pull(LRP)
# projLRPi <- projLRPAllp %>% filter(ProbThresh=="0.9") %>% pull(LRP)
# projLRPj <- projLRPAllp %>% filter(ProbThresh=="0.99") %>% pull(LRP)
# # # # # #
# # # # # # projLRP <- data.frame(fit=projLRPa)
# png(paste("DataOut/WCVI_SMUtimeseries_projLRPbaseERpValues_noEnh.png", sep=""), width=9, height=4, units="in", res=500)
# plotWCVI_SMUtimeseries(SMU_Esc=xx$SMU_Esc[38:length(xx$SMU_Esc)], out=xx$out$LRP, WCVI_Esc=xx$WCVIEsc[38:length(xx$SMU_Esc),])
# # abline(h=projLRPa+100, col="salmon", lwd=2)
# # abline(h=projLRPb, col="aquamarine3", lwd=2)
# # abline(h=projLRPc, col="black", lwd=2)
# # abline(h=projLRPd, col=viridis(3)[1], lwd=2)
# # abline(h=projLRPe, col=viridis(3)[2], lwd=2)
# # abline(h=projLRPf, col=viridis(3)[3], lwd=2)
# abline(h=projLRPg, col="orange", lwd=2)
# abline(h=projLRPh, col=viridis(4, alpha=0.2)[3], lwd=2)
# # abline(h=projLRPi, col=viridis(4, alpha=0.)[2], lwd=2)
# # abline(h=projLRPj, col=viridis(4, alpha=0.2)[1], lwd=2)
# dev.off()

# Current geometric mean?
# dum <- xx$WCVIEsc[61:64,c("Kyuquot", "Clayoquot", "Quatsino", "Barkley", "Nootka/Esperanza", "WCVI Nootka & Kyuquot", "WCVI South", "WCVI North")] 
# geoMean <- function(x){prod(x)^(1/length(x))}
# apply(dum,2,geoMean)
#Clayoquot: (prod(dum[c(1,2,4),2]))^(1/3)

#================================================================== #
# Plot time-series of CU status ----
#================================================================== #

## xx <- Get.LRP(remove.EnhStocks=TRUE)
# png("DataOut/WCVI_CUstatus_noEnh.png", width=4.5, height=6.5, units="in", res=200)
# par(mfrow=c(3,1), cex=1.1,  mar = c(2, 5, 2, 1) + 0.1)
# Years <- rownames(xx$WCVIEsc)
# 
# for (i in c(3,1,2)){
#   plot(x=Years, y=xx$CU_Status[,i], pch=20, xlab="", ylab="", ylim=c(0,1), las=1 , yaxp=c(0,1,1))
#   mtext(colnames(xx$CU_Status)[i], side=3, line=0.5, at="1952", adj=0, cex=1.4)
#   if(i==1) mtext("CU Status = Red(0) or Amber/Green(1)", side=2, line=2, cex=1.8)
# 
# }
# dev.off()



#================================================================== #
# Plot time-series of SMU ppns without enh ----
#================================================================== #
## xx <- Get.LRP(remove.EnhStocks=TRUE)

# png("DataOut/WCVI_SMUstatus_noEnh.png", width=4.5, height=3.5, units="in", res=200)
# par(mfrow=c(1,1), cex=1.1,  mar = c(2, 3, 2, 1) + 0.1)
# 
# plot(x=Years, y=xx$SMU_ppn, pch=20, xlab="", ylab="", ylim=c(0,1), las=1)
# mtext("WCVI SMU", side=3, line=0.5, at="1952", adj=0, cex=1.4)
# 
# dev.off()

# Plot data for log regression: Need df: SMU_Esc from "WCVILRPs.R"
# plot(x=SMU_Esc, y=SMU_ppn, xlab="Aggregate Escapement", ylab="Proportion CUs above red zone", las=1, ylim=c(0,1), xlim=c(0, max(SMUEsc, na.rm=T)))

#================================================================== #
# Plot logistic model from ####
#============================================================== #

plotLogistic <- function(Data, Preds, LRP, useGenMean = F, plotName, outDir, p=0.95, useBern_Logistic = F){
  
  # y label different depending on model
  if(useBern_Logistic){
    Ylab = "Pr(All CUs > Lower Benchmark)"
  } else {
    Ylab = "Prop. CUs > Lower Benchmark"
  }
  
  if(useGenMean){
    Xlab = "Gen. Mean Agg. Escapement"
  }else {
    Xlab = "Aggregate Escapement"
  }
  
  # If LRP$lwr is <0, assume can't fit LRP
  
  # Create plot
  annual_LRP_plot <- ggplot(data=Preds, mapping=aes(x=xx,y=fit)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr, x=xx), fill = "grey85") +
    geom_line(mapping=aes(x=xx, y=fit), col="red", size=1) +
    geom_line(mapping=aes(x=xx, y=upr), col="grey85") +
    geom_line(mapping=aes(x=xx, y=lwr), col="grey85") +
    geom_point(data=Data, aes(x = xx, y = yy)) +
    xlab(Xlab) + ylab(Ylab) +
    #  ggtitle(paste("Retrospective Year", max(Dat.LRP$yr), sep=" ")) +
    theme_classic()
  
  if(LRP$lwr > 0) {
    annual_LRP_plot <- annual_LRP_plot + 
      geom_vline(xintercept=LRP$fit, linetype="dashed", color="orange", size = 1) +
      geom_hline(yintercept= p, linetype="dashed", color="orange", size = 1) +
      geom_vline(xintercept = LRP$lwr, linetype = "dotted", color = "orange", size = 1) +
      geom_vline(xintercept = LRP$upr, linetype = "dotted", color = "orange", size = 1) 
  }
  
  
  annual_LRP_plot 
  # # Save plot
  # ggsave(paste(outDir, "/",plotName,".png",sep=""), plot = annual_LRP_plot,
  #        width = 4, height = 3, units = "in")      
  
}  

## xx <- Get.LRP(remove.EnhStocks=TRUE)
#plotLogistic(Data=xx$out$Logistic_Data, Preds=xx$out$Preds, LRP=xx$out$LRP, useGenMean = F, plotName="WCVI_logReg_noEnh", outDir="DataOut", p=0.95)

## yy <- Get.LRP(remove.EnhStocks=FALSE)
#plotLogistic(Data=yy$out$Logistic_Data, Preds=yy$out$Preds, LRP=yy$out$LRP, useGenMean = F, plotName="WCVI_logReg_wEnh", outDir="DataOut", p=0.95)

# Data <- out$Logistic_Data
# Preds <- out$Preds 
# LRP <- out$LRP
# useGenMean <- F
# p <- 0.95
# useBern_Logistic <- F

#if(xx$out$LRP$lwr<=0) xx$out$LRP$lwr <-1
#if(yy$out$LRP$lwr<=0) yy$out$LRP$lwr <-1


#=============================================================================== #
# Plot annual indicator abundance plots (WCVIO) ####
#=============================================================================== #

PlotAnnualIndicator <- FALSE
if(PlotAnnualIndicator){
  WCVIEsc <- data.frame(read.csv("DataIn/WCVIEsc.csv", row.names="Yr")) %>% 
    dplyr::select (-"Little.Zeballos")
  
  WCVIEsc <- WCVIEsc %>% add_column(Year=as.numeric( rownames(WCVIEsc) ) )
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
  
  WCVIEsc_long <- pivot_longer(WCVIEsc, cols = WCVIStocks$Stock, 
                               names_to = "Stock" , values_to = "Spawners")
  Enh_ <- WCVIStocks %>% select(c(Stock, Enh))
  WCVIEsc_long <- left_join(WCVIEsc_long, Enh_)
  
  WCVIEsc_long <- WCVIEsc_long %>% mutate(Spawners=Spawners/1000)
  
  
  PNI <- as.data.frame(read.csv("DataIn/PNI_indStocks.csv"))
  PNI <- PNI %>% select(c(-Stock.noHyphens,-startYr, -endYr))
  #WCVIEsc_long %>% left_join(PNI, by="Stock")
  
  #Add PNI to facets
  # e.g., for code to do this see 3rd and 4th plots here: https://stackoverflow.com/questions/11889625/annotating-text-on-individual-facet-in-ggplot2
  # 
  dat.text <- data.frame(label=PNI$PNI, Stock=PNI$Stock, Enh=PNI$Enh)
  
  IndicatorTimeSeries <- ggplot (WCVIEsc_long, aes (x=Year, y=Spawners, group=Enh)) + 
    geom_line( aes(colour = Enh)) + 
    facet_wrap(~Stock, scales="free_y") + 
    theme(legend.position = "none") +
    # xlab("Année") + 
    # ylab("Géniteurs") +
    geom_text(data=dat.text,
    mapping = aes(x = -Inf, y = -Inf, label = label),
    hjust   = -5.8,#0.1,
    vjust   = -11.5,
    size    =2,
    colour = grey(0.5)
  )

   
  
  ggsave("DataOut/IndicatorTimeSeriesFR.png", plot = IndicatorTimeSeries,
         height = 5, width = 8, units = "in")
  # ggsave("DataOut/IndicatorTimeSeriesFR.png", plot = IndicatorTimeSeries, 
  #        height = 5, width = 8, units = "in") 
  
  
}

#================================================================== #
# Plot annual inlet status plot (WCVI) ####
#============================================================== #

library(patchwork)

# ***must run code within WCVILRPs.R code to get variable, Inlet_Status

# Inlet_StatusNumeric <- t(apply(X= Inlet_Status, MARGIN = 1, FUN=as.numeric))
# 
# colnames(Inlet_StatusNumeric)  <- colnames(Inlet_Status)
# 
# Inlet_StatusNumeric <- as.data.frame(Inlet_StatusNumeric) %>% 
#   add_column(Years=as.numeric(Years))
# 
# colnames(Inlet_StatusNumeric)[5] <- "NootkaEsperanza"
# 
# Inlet_StatusNumeric
# g1 <- ggplot(Inlet_StatusNumeric) + geom_point(aes(x=Years, y=Barkley)) +
#   scale_x_continuous(name ="") + scale_y_continuous(breaks=c(0,1), limits=c(0,1)) + 
#   theme(text = element_text(size=17))
# 
# g2 <- ggplot(Inlet_StatusNumeric) + geom_point(aes(x=Years, y=Clayoquot))+
#   scale_x_continuous(name ="") + scale_y_continuous(breaks=c(0,1), limits=c(0,1)) + 
#   theme(text = element_text(size=17))
# g3 <- ggplot(Inlet_StatusNumeric) + geom_point(aes(x=Years, y=NootkaEsperanza))+
#   scale_x_continuous(name ="") + scale_y_continuous(breaks=c(0,1), limits=c(0,1))+ 
#   theme(text = element_text(size=17))
# g4 <- ggplot(Inlet_StatusNumeric) + geom_point(aes(x=Years, y=Kyuquot))+
#   scale_x_continuous(name ="") + scale_y_continuous(breaks=c(0,1), limits=c(0,1))+ 
#   theme(text = element_text(size=17))
# g5 <- ggplot(Inlet_StatusNumeric) + geom_point(aes(x=Years, y=Quatsino))+
#   scale_x_continuous(name ="") + scale_y_continuous(breaks=c(0,1), limits=c(0,1))+ 
#   theme(text = element_text(size=17))
# 
# Inlet_Status_plot <- g5 + g4 + g3 + g2 + g1 + plot_layout(ncol = 1)
# plotName <- "InletStatus"
# ggsave(paste("DataOut/",plotName,".png",sep=""), plot = Inlet_Status_plot,
#        width = 7, height = 10, units = "in")  

