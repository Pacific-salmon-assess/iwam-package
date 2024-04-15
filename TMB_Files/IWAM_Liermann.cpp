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
  
// REAL CODE STARTS  
template<class Type>
Type objective_function<Type>:: operator() ()
{
  // REMOVE ALL _std objects
  DATA_VECTOR(S); // Org _std 
  DATA_VECTOR(logRS); // Org _std
  DATA_IVECTOR(stk); // stock number // Org _std
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
  DATA_INTEGER(SigRicPriorNorm);
  DATA_INTEGER(SigRicPriorGamma);
  DATA_INTEGER(SigRicPriorCauchy);
  DATA_INTEGER(biasCor);
  DATA_INTEGER(SigDeltaPriorNorm);
  DATA_INTEGER(SigDeltaPriorGamma);
  DATA_INTEGER(SigDeltaPriorCauchy);
  DATA_SCALAR(Tau_dist);
  DATA_SCALAR(Tau_D_dist);
  //DATA_SCALAR(logDeltaSigma);
  //DATA_SCALAR(logNuSigma);
  
  DATA_SCALAR(SigDelta_mean);
  DATA_SCALAR(SigDelta_sig);
  DATA_SCALAR(SigNu_mean);
  DATA_SCALAR(SigNu_sig);
  
  DATA_VECTOR(pred_lnWA);
  DATA_VECTOR(target_lnWA_ocean);
  DATA_VECTOR(target_lnWA_stream);
  
  // Remove all _std
  PARAMETER_VECTOR(logA); // Org _std
  PARAMETER_VECTOR(logB); // Org _std
  PARAMETER_VECTOR(logSigma); // Org _std
  
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
  int N_stks = scale.size(); // Removed _std
  
  

  vector <Type> logRS_pred(N_Obs); // Removed _std
  vector <Type> sigma = exp(logSigma); // Removed _std
  Type sigmaA = exp(logSigmaA);
  vector <Type> nLL(N_Obs); // negative log-likelihood to calculate AIC - not need for est.

  // Standard Ricker model: 
  for (int i = 0; i<N_Obs; i++){
    if(biasCor == 0) {
      logRS_pred(i) = logA(stk(i)) - exp(logB(stk(i))) * S(i); // S is: Sp/scale
    }
    if(biasCor == 1) { // correcting for the back-calculation bias - from log transform to raw
      logRS_pred(i) = logA(stk(i)) - exp(logB(stk(i))) * S(i) - pow(sigma(stk(i)),2) / Type(2);
    } // power function - squared sigma / 2
    // look up TMB pow
    ans += -dnorm(logRS_pred(i), logRS(i),  sigma(stk(i)), true);    
    nLL(i) = -dnorm(logRS_pred(i), logRS(i),  sigma(stk(i)), true);
  }
  
  // Add hierarchical structure to A:
  for(int i=0; i<N_stks; i++){
    // add prior on logA, 
    ans += -dnorm(logA(i), logMuA_stream + logMuA_ocean * lifehist(i), sigmaA, true );
     // add prior on sigma, 
    if (SigRicPriorGamma == 1) {
       ans += -dgamma(pow(sigma(i),-2), Tau_dist, 1/Tau_dist, true);
    }
    if (SigRicPriorNorm == 1) {
      //ans += -abs( dnorm( sigma(i), HalfNormMean, HalfNormSig, true) );
      //3 June 2021. abs() function no longer works with TMB, so have removed
      ans += -( dnorm( sigma(i), HalfNormMean, HalfNormSig, true) );
    }
    if (SigRicPriorCauchy == 1) {
      //ans += - abs( dt( sigma(i), Type(1), true ));
      //3 June 2021. abs() function no longer works with TMB, so have removed
      ans += - ( dt( sigma(i), Type(1), true ));
    }
  }
  
  // Add priors for hyperpars:
  // MuA prior for stream type
  ans += -dnorm(logMuA_stream, logMuA_stream_mean, logMuA_stream_sig, true);
  // MuA prior for ocean type
  ans += -dnorm(logMuA_ocean, logMuA_ocean_mean, logMuA_ocean_sig, true);
  // sigmaA prior
  if (SigRicPriorGamma == 1) {
    ans += -dgamma(pow(sigmaA,-2), Tau_dist, 1/Tau_dist, true);
  }
  if (SigRicPriorNorm == 1) {
    //3June 2021. abs() functin no longer works in TMB
    //ans += -abs( dnorm( sigmaA, HalfNormMeanA, HalfNormSigA, true) );
    ans += -( dnorm( sigmaA, HalfNormMeanA, HalfNormSigA, true) );
  }
  if (SigRicPriorCauchy == 1) {
    //ans += - abs(dt( sigmaA, Type(1), true));
    ans += - (dt( sigmaA, Type(1), true));
  }
  
  
  
  //Calculate SMSY and SREP
  vector <Type> SMSY(N_stks); // Removed _std 
  vector <Type> SREP(N_stks); // Removed _std
  
  //For SMSY calculation, 
  for(int i=0; i<N_stks; i++){
    SMSY(i) =  (1 - LambertW( exp (1- logA(i)) ) ) / exp(logB(i)) ; //
  }
  SREP = logA / exp(logB);
 // SMSY(i) =  (1 - LambertW( exp (1- logA(i) - sigma^2/2) ) ) / exp(logB(i)) ; // non-bias correction version
  
  // Liermann's model with both stream and ocean type =================
  vector <Type> pred_lnSMSY(N_stks);
  vector <Type> pred_lnSREP(N_stks);
  Type sigma_delta = exp(logDeltaSigma);
  Type sigma_nu = exp(logNuSigma);
  
  for (int i=0; i<N_stks; i++){ // THE ACTUAL WATERSHED MODEL
    pred_lnSMSY(i) = logDelta1 + logDelta1_ocean * lifehist(i) + ( exp(logDelta2) + Delta2_ocean * lifehist(i) ) * log(WAbase(i)) ;
      // Confusion about log-space vs non log-space
      // From Parken model (allometric equation)
    ans += -dnorm( pred_lnSMSY(i), log(SMSY(i) * scale(i) ),  sigma_delta, true);
    pred_lnSREP(i) = logNu1 + logNu1_ocean * lifehist(i) + ( exp(logNu2) + Nu2_ocean * lifehist(i) ) * log(WAbase(i)) ;
    ans += -dnorm( pred_lnSREP(i), log(SREP(i) * scale(i) ),  sigma_nu, true);
  }
  // Stream-type is the base and deviation for the ocean
  
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
  
  // Half cauchy prior on sigma_delta and sigma_nu
  if (SigDeltaPriorCauchy == 1) {
    //3 June 2021. abs() no longer works in TMB
    //ans += -abs( dt( sigma_delta, Type(1), true));
    //ans += - abs( dt( sigma_nu, Type(1), true ));
    ans += -( dt( sigma_delta, Type(1), true));
    ans += -( dt( sigma_nu, Type(1), true ));
  }
  
    // Get predicted values for plotting WA regresssion with CIs
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
  
  
  
  //// Get predicted values for stream-type target stocks with CIs
  int N_target_stream = target_lnWA_stream.size();
  vector <Type> target_lnSMSY_stream(N_target_stream);
  vector <Type> target_lnSREP_stream(N_target_stream);

  for (int i=0; i<N_target_stream; i++){
    target_lnSMSY_stream(i) = logDelta1 + exp(logDelta2) * target_lnWA_stream(i);
    target_lnSREP_stream(i) = logNu1 + exp(logNu2) * target_lnWA_stream(i);
  } 
  
  ///Get predicted values for ocean-type target stocks with CIs
  int N_target_ocean = target_lnWA_ocean.size();
  vector <Type> target_lnSMSY_ocean(N_target_ocean);    
  vector <Type> target_lnSREP_ocean(N_target_ocean);

  for (int i=0; i<N_target_ocean; i++){
    target_lnSMSY_ocean(i) = logDelta1 + logDelta1_ocean + (exp(logDelta2) + Delta2_ocean) * target_lnWA_ocean(i);
    target_lnSREP_ocean(i) = logNu1 + logNu1_ocean + (exp(logNu2) + Nu2_ocean) * target_lnWA_ocean(i);
  }

  
  vector <Type> lnSMSY = log(SMSY*scale); //This is taking a value on the real
  // scale * a constant scalar (ignore) and putting it on the log-scale
  // This shouldn't matter because it was calculated on the real-scale
  // It's only if you are exponentiating it back that it would be bad
  vector <Type> lnSREP = log(SREP*scale);
  
  
  ADREPORT(SMSY); // Removed _std
  ADREPORT(SREP); // Removed _std
  ADREPORT(logRS_pred);
  ADREPORT(pred_lnSMSY);
  ADREPORT(pred_lnSREP);
  ADREPORT(lnSMSY);
  ADREPORT(lnSREP);
  ADREPORT(pred_lnSMSY_stream_CI);
  ADREPORT(pred_lnSMSY_ocean_CI);
  ADREPORT(pred_lnSREP_stream_CI);
  ADREPORT(pred_lnSREP_ocean_CI);
  ADREPORT(target_lnSMSY_ocean);
  ADREPORT(target_lnSREP_ocean);
  ADREPORT(target_lnSMSY_stream);
  ADREPORT(target_lnSREP_stream);
  REPORT(nLL); // Removed _std
  return ans;
}

