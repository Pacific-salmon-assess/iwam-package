#include <TMB.hpp >

// Set up Lambert's W function to use to calculate SMSY
// Code taken from https://kaskr.github.io/adcomp/lambert_8cpp_source.html
// Step 1: Code up a plain C version
// Double version of Lambert W function
double LambertW(double x) {
  double logx = log(x);
  double y = (logx > 0 ? logx : 0);
  int niter = 100, i=0;
  for (; i < niter; i++) {
    if ( fabs( logx - log(y) - y) < 1e-9) break;
    y -= (y - exp(logx - y)) / (1 + y);
  }
  if (i == niter) Rf_warning("W: failed convergence");
  return y;
}

TMB_ATOMIC_VECTOR_FUNCTION(
  // ATOMIC_NAME
  LambertW
  ,
  // OUTPUT_DIM
  1,
  // ATOMIC_DOUBLE
  ty[0] = LambertW(tx[0]); // Call the 'double' version
,
// ATOMIC_REVERSE
Type W  = ty[0];                    // Function value from forward pass
Type DW = 1. / (exp(W) * (1. + W)); // Derivative
px[0] = DW * py[0];                 // Reverse mode chain rule
)
  
  // Scalar version
  template<class Type>
  Type LambertW(Type x){
    CppAD::vector<Type> tx(1);
    tx[0] = x;
    return LambertW(tx)[0];
  }
  
// Opening calls ---------------------------------------------------------------
template<class Type>
Type objective_function<Type>:: operator() ()
{
  DATA_VECTOR(S); 
  DATA_VECTOR(logRS); 
  
  DATA_IVECTOR(stk); // stock number 
  DATA_IVECTOR(yr); // only occurs once
  
  DATA_SCALAR(logMuA_stream_mean);
  DATA_SCALAR(logMuA_stream_sig);
  DATA_SCALAR(logMuA_ocean_mean);
  DATA_SCALAR(logMuA_ocean_sig);
  DATA_SCALAR(HalfNormMean);
  DATA_SCALAR(HalfNormSig);
  DATA_SCALAR(HalfNormMeanA);
  DATA_SCALAR(HalfNormSigA);
  
  DATA_VECTOR(WAbase); 
  DATA_VECTOR(scale);
  DATA_IVECTOR(lifehist);

  DATA_INTEGER(biasCor);
    
  DATA_INTEGER(SigRicPriorNorm); // on/off
  DATA_INTEGER(SigRicPriorGamma); // on/off
  DATA_INTEGER(SigRicPriorCauchy); // on/off
  
  DATA_INTEGER(SigDeltaPriorNorm); // on/off
  DATA_INTEGER(SigDeltaPriorGamma); // on/off
  DATA_INTEGER(SigDeltaPriorCauchy); // on/off
  
  // Alternative's
  DATA_INTEGER(SigRicPenal); // Gamma Ricker penalty on sigma
  DATA_INTEGER(SigDeltaPenal); // Gamma WA penalty on precision
  DATA_INTEGER(SigDeltaPenal_Jac); // Invgamma WA penalty on variance
  
  // Depreciated controls for shape and rate
  DATA_SCALAR(Tau_dist); // Turn off to manually change penalties
  DATA_SCALAR(Tau_D_dist); // Turn off to manually change penalties
  
  // New penalty term controls for shape and rate (scale = 1/rate)
  // DATA_SCALAR(Ric_sigshape);
  // DATA_SCALAR(Ric_sigrate);
  // DATA_SCALAR(WA_sigshape);
  // DATA_SCALAR(WA_sigscale);
  // DATA_SCALAR(WA_sigshapeJac); // Specific shape for Jacobian alt.
  
  //DATA_SCALAR(logDeltaSigma);
  //DATA_SCALAR(logNuSigma);
  
  DATA_SCALAR(SigDelta_mean);
  DATA_SCALAR(SigDelta_sig);
  DATA_SCALAR(SigNu_mean);
  DATA_SCALAR(SigNu_sig);
  
  DATA_VECTOR(pred_lnWA);
  DATA_VECTOR(target_lnWA_ocean);
  DATA_VECTOR(target_lnWA_stream);
  
  PARAMETER_VECTOR(logA); 
  PARAMETER_VECTOR(logB); 
  PARAMETER_VECTOR(logSigma); 
  
  PARAMETER(logMuA_stream);
  PARAMETER(logSigmaA);
  PARAMETER(logMuA_ocean);
  PARAMETER(logDelta1);
  PARAMETER(logDelta1_ocean);
  PARAMETER(logDelta2);
  PARAMETER(Delta2_ocean);
  PARAMETER(logDeltaSigma);
  
  PARAMETER(logNu1);
  PARAMETER(logNu1_ocean);
  PARAMETER(logNu2);
  PARAMETER(Nu2_ocean);
  PARAMETER(logNuSigma);
  
  
  
  Type ans=0.0; // ans is the log-likelihood - is then additive for each of the distributions
  int N_Obs = S.size(); //size() gives the size of the vector - TMB function
  int N_stks = scale.size(); 
  
  vector <Type> logRS_pred(N_Obs); 
  vector <Type> sigma = exp(logSigma); 
  Type sigmaA = exp(logSigmaA);
  vector <Type> nLL(N_Obs); // negative log-likelihood to calculate AIC - not need for est.

  // Standard Ricker model: ----------------------------------------------------
  for (int i = 0; i<N_Obs; i++){
    if(biasCor == 0) {
      logRS_pred(i) = logA(stk(i)) - exp(logB(stk(i))) * S(i); // S is: Sp/scale
    }
    if(biasCor == 1) { // correcting for the back-calculation bias - from log transform to raw
      logRS_pred(i) = logA(stk(i)) - exp(logB(stk(i))) * S(i) - pow(sigma(stk(i)),2) / Type(2);
    } // power function - squared sigma / 2
    ans += -dnorm(logRS_pred(i), logRS(i),  sigma(stk(i)), true);    
    nLL(i) = -dnorm(logRS_pred(i), logRS(i),  sigma(stk(i)), true);
  }
  
  // Add hierarchical structure to A: ------------------------------------------
  for(int i=0; i<N_stks; i++){
    // add prior on logA, 
    ans += -dnorm(logA(i), logMuA_stream + logMuA_ocean * lifehist(i), sigmaA, true );
    
     // add prior on sigma, 
    if (SigRicPriorGamma == 1) {
      ans += -dgamma(pow(sigma(i),-2), Tau_dist, 1/Tau_dist, true);
        // Change prior to 0.1, 0.1 - because rate should be 1/term
        // (for gamma --> invgamma)
        // which means scale is 1/(1/term) or in this case just term
        // dgamma on precision
        // invgamma on tau
    }
    
    // Alternative penalty directly on sigma        
    if (SigRicPenal == 1) {
    ans += -dgamma(sigma(i), Type(7.5), Type(0.1), true);
    // Type(7.5), Type(0.1)
    }
    
    // Half normal
    if (SigRicPriorNorm == 1) {
      //ans += -abs( dnorm( sigma(i), HalfNormMean, HalfNormSig, true) );
      //3 June 2021. abs() function no longer works with TMB, so have removed
      ans += -(dnorm(sigma(i), HalfNormMean, HalfNormSig, true) );
    }
    
    // // Half cauchy
    if (SigRicPriorCauchy == 1) {
      //ans += - abs( dt( sigma(i), Type(1), true ));
      //3 June 2021. abs() function no longer works with TMB, so have removed
      ans += - (dt(sigma(i), Type(1), true ));
    }
  }
  
  
  
  // Add priors for hyperpars: -------------------------------------------------
  // MuA prior for stream type
  ans += -dnorm(logMuA_stream, logMuA_stream_mean, logMuA_stream_sig, true);
  // MuA prior for ocean type
  ans += -dnorm(logMuA_ocean, logMuA_ocean_mean, logMuA_ocean_sig, true);
  
  // sigmaA prior
  if (SigRicPriorGamma == 1) {
    ans += -dgamma(pow(sigmaA,-2), Tau_dist, 1/Tau_dist, true);
  }
  
  // Alternative penalty directly on sigma
  if (SigRicPenal == 1) {
    ans += -dgamma(sigmaA, Type(7.5), Type(0.1), true);
  }
  
  // // Half Normal
  if (SigRicPriorNorm == 1) {
    //3June 2021. abs() functin no longer works in TMB
    //ans += -abs( dnorm( sigmaA, HalfNormMeanA, HalfNormSigA, true) );
    ans += -(dnorm(sigmaA, HalfNormMeanA, HalfNormSigA, true) );
  }
  
  // Half cauchy
  if (SigRicPriorCauchy == 1) {
    //ans += - abs(dt( sigmaA, Type(1), true));
    ans += -(dt(sigmaA, Type(1), true));
  }
  
  
  
  //Calculate SMSY and SREP ----------------------------------------------------
  vector <Type> SMSY(N_stks); // Removed _std 
  vector <Type> SREP(N_stks); // Removed _std
  
  //For SMSY calculation, 
  for(int i=0; i<N_stks; i++){
    SMSY(i) =  (1 - LambertW( exp (1- logA(i)) ) ) / exp(logB(i)) ; // scaled SMSY
  }
  SREP = logA / exp(logB); // scaled SREP
  // SMSY(i) =  (1 - LambertW( exp (1- logA(i) - sigma^2/2) ) ) / exp(logB(i)) ; 
    // non-bias correction version
   
  
  
  // Liermann's model with both stream and ocean type --------------------------
  vector <Type> pred_lnSMSY(N_stks);
  vector <Type> pred_lnSREP(N_stks);
  Type sigma_delta = exp(logDeltaSigma);
  Type sigma_nu = exp(logNuSigma);
  
  // for (int i=0; i<N_stks; i++){ // THE ACTUAL WATERSHED MODEL
  //   pred_lnSMSY(i) = logDelta1 + logDelta1_ocean * lifehist(i) + ( exp(logDelta2) + Delta2_ocean * lifehist(i) ) * log(WAbase(i)) ;
  //     // Confusion about log-space vs non log-space
  //     // From Parken model (allometric equation)
  //   ans += -dnorm( pred_lnSMSY(i), log(SMSY(i) * scale(i) ),  sigma_delta, true);
  //   pred_lnSREP(i) = logNu1 + logNu1_ocean * lifehist(i) + ( exp(logNu2) + Nu2_ocean * lifehist(i) ) * log(WAbase(i)) ;
  //   ans += -dnorm( pred_lnSREP(i), log(SREP(i) * scale(i) ),  sigma_nu, true);
  // }
  // Stream-type is the base and deviation for the ocean
  // How is process error shown in this model?
  
  for (int i=0; i<N_stks; i++){
    if(biasCor == 0) {
      pred_lnSMSY(i) = logDelta1 + logDelta1_ocean * lifehist(i) + ( exp(logDelta2) + Delta2_ocean * lifehist(i) ) * log(WAbase(i)) ;
    }
    if(biasCor == 1) { // Bias corrected
      pred_lnSMSY(i) = logDelta1 + logDelta1_ocean * lifehist(i) + ( exp(logDelta2) + Delta2_ocean * lifehist(i) ) * log(WAbase(i)) - pow(sigma_delta,2) / Type(2);
    }
    ans += -dnorm( pred_lnSMSY(i), log(SMSY(i) * scale(i) ),  sigma_delta, true);
    
    if(biasCor == 0) {
      pred_lnSREP(i) = logNu1 + logNu1_ocean * lifehist(i) + ( exp(logNu2) + Nu2_ocean * lifehist(i) ) * log(WAbase(i)) ;
    }
    if(biasCor == 1) { // Bias corrected
      pred_lnSREP(i) = logNu1 + logNu1_ocean * lifehist(i) + ( exp(logNu2) + Nu2_ocean * lifehist(i) ) * log(WAbase(i))  - pow(sigma_nu,2) / Type(2);
    }
    ans += -dnorm( pred_lnSREP(i), log(SREP(i) * scale(i) ),  sigma_nu, true);
  }
  
  
  
  // Normal prior on sigma_delta and sigma_nu
  if (SigDeltaPriorNorm == 1) {
    ans += -dnorm(sigma_delta, SigDelta_mean, SigDelta_sig, true);
    ans += -dnorm(sigma_nu, SigNu_mean, SigNu_sig, true);
  }
  
  // Inverse gamma prior on sigma_delta and sigma_nu
  if (SigDeltaPriorGamma == 1) {
    ans += -dgamma(pow(sigma_delta,-2), Tau_D_dist, 1/Tau_D_dist, true);
    ans += -dgamma(pow(sigma_nu,-2), Tau_D_dist, 1/Tau_D_dist, true);
  }
  
  // Alternative: Gamma penalty on precision
    // Shape = 3, scale = rate = 1
  if (SigDeltaPenal == 1) {
    ans += -dgamma(pow(sigma_delta, -2), Type(3), Type(1), true);
    ans += -dgamma(pow(sigma_nu, -2), Type(3), Type(1), true);
  }
  
  // Alternative: Invgamma on variance w/ Jacobian: 
    // Shape 0.75, scale = rate = 1
  if (SigDeltaPenal_Jac == 1) {
    ans += -dgamma(pow(sigma_delta,-2), Type(0.75), Type(1), true);
    ans += Type(2)*log(pow(sigma_delta,2)); //Jacobian adjustment
    ans += -dgamma(pow(sigma_nu,-2), Type(0.75), Type(1), true);
    ans += Type(2)*log(pow(sigma_nu,2)); //Jacobian adjustment
  }
  
  // Alternative: Another Normal
    // ans += -dnorm(sigma_delta, Type(1), Type(0.1), true);
    // ans += -dnorm(sigma_nu, Type(1), Type(0.1), true);
    
  // Alternative: STATIC VALUE
    // Comment out all of the above and set logNuSigma and logDeltaSigma 
    // within IWAM_model.R
  
  // Half cauchy prior on sigma_delta and sigma_nu
  if (SigDeltaPriorCauchy == 1) {
    //3 June 2021. abs() no longer works in TMB
    //ans += -abs( dt( sigma_delta, Type(1), true));
    //ans += - abs( dt( sigma_nu, Type(1), true ));
    ans += -( dt( sigma_delta, Type(1), true));
    ans += -( dt( sigma_nu, Type(1), true ));
  }
  
  
  
  // Get predicted values for plotting WA regresssion with CIs -----------------
  int N_pred = pred_lnWA.size();
  vector <Type> pred_lnSMSY_stream_CI(N_pred);
  vector <Type> pred_lnSMSY_ocean_CI(N_pred);
  vector <Type> pred_lnSREP_stream_CI(N_pred);
  vector <Type> pred_lnSREP_ocean_CI(N_pred);
  
  for (int i=0; i<N_pred; i++){
    pred_lnSMSY_stream_CI(i) = logDelta1 + exp(logDelta2) * pred_lnWA(i);
    pred_lnSMSY_ocean_CI(i) = logDelta1 + logDelta1_ocean + (exp(logDelta2) + Delta2_ocean) * pred_lnWA(i);
    pred_lnSREP_stream_CI(i) = logNu1 + exp(logNu2) * pred_lnWA(i);
    pred_lnSREP_ocean_CI(i) = logNu1 + logNu1_ocean + (exp(logNu2) + Nu2_ocean) * pred_lnWA(i);
  }
  
  //// Get predicted values for stream-type target stocks with CIs -------------
  int N_target_stream = target_lnWA_stream.size();
  vector <Type> target_lnSMSY_stream(N_target_stream);
  vector <Type> target_lnSREP_stream(N_target_stream);
  
  for (int i=0; i<N_target_stream; i++){
    target_lnSMSY_stream(i) = logDelta1 + exp(logDelta2) * target_lnWA_stream(i);
    target_lnSREP_stream(i) = logNu1 + exp(logNu2) * target_lnWA_stream(i);
  } 
  
  ///Get predicted values for ocean-type target stocks with CIs ----------------
  int N_target_ocean = target_lnWA_ocean.size();
  vector <Type> target_lnSMSY_ocean(N_target_ocean);
  vector <Type> target_lnSREP_ocean(N_target_ocean);
  
  for (int i=0; i<N_target_ocean; i++){
    target_lnSMSY_ocean(i) = logDelta1 + logDelta1_ocean + (exp(logDelta2) + Delta2_ocean) * target_lnWA_ocean(i);
    target_lnSREP_ocean(i) = logNu1 + logNu1_ocean + (exp(logNu2) + Nu2_ocean) * target_lnWA_ocean(i);
  }
  
  // ///Get predicted lines values ------------------------------------------------
  // int N_target_ocean = target_lnWA_ocean.size();
  // vector <Type> target_lnSMSY_ocean(N_target_ocean);
  // vector <Type> target_lnSREP_ocean(N_target_ocean);
  // 
  // for (int i=0; i<N_target_ocean; i++){
  //   target_lnSMSY_ocean(i) = logDelta1 + logDelta1_ocean + (exp(logDelta2) + Delta2_ocean) * target_lnWA_ocean(i);
  //   target_lnSREP_ocean(i) = logNu1 + logNu1_ocean + (exp(logNu2) + Nu2_ocean) * target_lnWA_ocean(i);
  // }
  
  
  
  // REPORTING -----------------------------------------------------------------
  vector <Type> lnSMSY = log(SMSY*scale); //This is taking a value on the real
  // scale * a constant scalar (ignore) and putting it on the log-scale
  // This shouldn't matter because it was calculated on the real-scale
  // It's only if you are exponentiating it back that it would be bad
  vector <Type> lnSREP = log(SREP*scale);
  
  ADREPORT(SMSY); // Removed _std
  ADREPORT(SREP); // Removed _std
  ADREPORT(logRS_pred);
  ADREPORT(pred_lnSMSY); // lnSMSY for synoptic set
  ADREPORT(pred_lnSREP);
  ADREPORT(lnSMSY);
  ADREPORT(lnSREP); 
  
  ADREPORT(pred_lnSMSY_stream_CI); // synoptic "predicted" lnSMSY along a line
  ADREPORT(pred_lnSMSY_ocean_CI);
  ADREPORT(pred_lnSREP_stream_CI);
  ADREPORT(pred_lnSREP_ocean_CI);
  
  ADREPORT(target_lnSMSY_ocean); // actual predicted lnSMSY and lnSREP from watershed model for WCVI examples
  ADREPORT(target_lnSREP_ocean);
  ADREPORT(target_lnSMSY_stream);
  ADREPORT(target_lnSREP_stream);
  
  REPORT(nLL); // Removed _std
  return ans;
}




