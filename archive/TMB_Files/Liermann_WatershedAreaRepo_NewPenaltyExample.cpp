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
  
  
template<class Type>
Type objective_function<Type>:: operator() ()
{
  DATA_VECTOR(S_std);
  DATA_VECTOR(logRS_std);
  DATA_IVECTOR(stk_std);
  DATA_IVECTOR(yr_std);
  DATA_SCALAR(logMuAs_mean);
  DATA_SCALAR(logMuAs_sig);
  DATA_SCALAR(logMuAo_mean);
  DATA_SCALAR(logMuAo_sig);
  DATA_SCALAR(HalfNormMean);
  DATA_SCALAR(HalfNormSig);
  DATA_SCALAR(HalfNormMeanA);
  DATA_SCALAR(HalfNormSigA);
  DATA_VECTOR(WA); 
  DATA_VECTOR(Scale);
  DATA_IVECTOR(Stream);
  DATA_INTEGER(SigRicPriorNorm);
  DATA_INTEGER(SigRicPriorGamma);
  DATA_INTEGER(SigRicPriorCauchy);
  DATA_INTEGER(biasCor);
  DATA_INTEGER(SigDeltaPriorNorm);
  DATA_INTEGER(SigDeltaPriorGamma);
  DATA_INTEGER(SigDeltaPriorCauchy);
  //DATA_SCALAR(Tau_dist);
  //DATA_SCALAR(Tau_D_dist);
  //DATA_SCALAR(logDeltaSigma);
  //DATA_SCALAR(logNuSigma);
  
  DATA_SCALAR(sigDelta_mean);
  DATA_SCALAR(sigDelta_sig);
  DATA_SCALAR(sigNu_mean);
  DATA_SCALAR(sigNu_sig);
  
  DATA_VECTOR(PredlnWA);
  DATA_VECTOR(TestlnWAo);
  //DATA_VECTOR(TestlnWAs);
  

  PARAMETER_VECTOR(logA_std);
  PARAMETER_VECTOR(logB_std);
  PARAMETER_VECTOR(logSigma_std);
  PARAMETER(logMuAs);
  PARAMETER(logSigmaA);
  PARAMETER(logMuAo);
  PARAMETER(logDelta1);
  PARAMETER(logDelta1ocean);
  PARAMETER(logDelta2);
  PARAMETER(Delta2ocean);
  PARAMETER(logDeltaSigma);
  
  PARAMETER(logNu1);
  PARAMETER(logNu1ocean);
  PARAMETER(logNu2);
  PARAMETER(Nu2ocean);
  PARAMETER(logNuSigma);
  

  
  Type ans=0.0;
  int N_Obs = S_std.size(); 
  int N_stks_std = Scale.size();
  
  

  vector <Type> LogRS_Pred_std(N_Obs);
  vector <Type> sigma_std = exp(logSigma_std);
  Type sigmaA = exp(logSigmaA);
  vector <Type> nLL_std(N_Obs);

  // Standard Ricker model: 
  for (int i = 0; i<N_Obs; i++){
    if(biasCor == 0) {
      LogRS_Pred_std(i) = logA_std(stk_std(i)) - exp(logB_std(stk_std(i))) * S_std(i);
      }
    if(biasCor == 1) {
      LogRS_Pred_std(i) = logA_std(stk_std(i)) - exp(logB_std(stk_std(i))) * S_std(i) - pow(sigma_std(stk_std(i)),2) / Type(2);
    }
    ans += -dnorm(LogRS_Pred_std(i), logRS_std(i),  sigma_std(stk_std(i)), true);    
    nLL_std(i) = -dnorm(LogRS_Pred_std(i), logRS_std(i),  sigma_std(stk_std(i)), true);
  }
  
  // Add hierarchical structure to A ==============
  for(int i=0; i<N_stks_std; i++){
    // add prior on logA_std, 
    ans += -dnorm(logA_std(i), logMuAs + logMuAo * Stream(i), sigmaA, true );
     // add prior on sigma 
    if (SigRicPriorGamma == 1) {
       //ans += -dgamma(pow(sigma_std(i),-2), Tau_dist, 1/Tau_dist, true);
      ans += -dgamma(sigma_std(i), Type(7.5), Type(0.1), true);
    }
    if (SigRicPriorNorm == 1) {
      //ans += -abs( dnorm( sigma_std(i), HalfNormMean, HalfNormSig, true) );
      //3 June 2021. abs() function no longer works with TMB, so have removed
      ans += -( dnorm( sigma_std(i), HalfNormMean, HalfNormSig, true) );
    }
    if (SigRicPriorCauchy == 1) {
      //ans += - abs( dt( sigma_std(i), Type(1), true ));
      //3 June 2021. abs() function no longer works with TMB, so have removed
      ans += - ( dt( sigma_std(i), Type(1), true ));
    }
  }
  
  // Add priors for hyperpars ====================
  // MuA prior for stream type
  ans += -dnorm(logMuAs, logMuAs_mean, logMuAs_sig, true);
  // MuA prior for ocean type
  ans += -dnorm(logMuAo, logMuAo_mean, logMuAo_sig, true);
  // sigmaA prior
  if (SigRicPriorGamma == 1) {
    //ans += -dgamma(pow(sigmaA,-2), Tau_dist, 1/Tau_dist, true);
    //ans += -dgamma(pow(sigmaA,-2), Type(3), Type(0.2), true);
    ans += -dgamma(sigmaA, Type(7.5), Type(0.1), true);//shape=2, rate=3, scale=1/3
    
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
  vector <Type> SMSY_std(N_stks_std);  
  vector <Type> SREP_std(N_stks_std);  
  
  //For SMSY calculation, 
  for(int i=0; i<N_stks_std; i++){
    SMSY_std(i) =  (1 - LambertW( exp (1- logA_std(i)) ) ) / exp(logB_std(i)) ;
  }
  SREP_std = logA_std / exp(logB_std);
  
  
  //Liermann's model with both stream and ocean type=================
  vector <Type> PredlnSMSY(N_stks_std);
  vector <Type> PredlnSREP(N_stks_std);
  Type sigma_delta = exp(logDeltaSigma);
  Type sigma_nu = exp(logNuSigma);
  
  for (int i=0; i<N_stks_std; i++){
    if(biasCor == 0) {
      PredlnSMSY(i) = logDelta1 + logDelta1ocean * Stream(i) + ( exp(logDelta2) + Delta2ocean * Stream(i) ) * log(WA(i)) ;
    }
    if(biasCor == 1) {
      PredlnSMSY(i) = logDelta1 + logDelta1ocean * Stream(i) + ( exp(logDelta2) + Delta2ocean * Stream(i) ) * log(WA(i)) - pow(sigma_delta,2) / Type(2);
    }
    ans += -dnorm( PredlnSMSY(i), log(SMSY_std(i) * Scale(i) ),  sigma_delta, true);
    
    if(biasCor == 0) {
      PredlnSREP(i) = logNu1 + logNu1ocean * Stream(i) + ( exp(logNu2) + Nu2ocean * Stream(i) ) * log(WA(i)) ;
    }
    if(biasCor == 1) {
      PredlnSREP(i) = logNu1 + logNu1ocean * Stream(i) + ( exp(logNu2) + Nu2ocean * Stream(i) ) * log(WA(i))  - pow(sigma_nu,2) / Type(2);
    }
    ans += -dnorm( PredlnSREP(i), log(SREP_std(i) * Scale(i) ),  sigma_nu, true);
  }
  
  // Normal prior on sigma_delta and sigma_nu
  if (SigDeltaPriorNorm == 1) {
    ans += -dnorm(sigma_delta, sigDelta_mean, sigDelta_sig, true);
    ans += -dnorm(sigma_nu, sigNu_mean, sigNu_sig, true);
  }
  
  // Inverse gamma prior on sigma_delta and sigma_nu
  if (SigDeltaPriorGamma == 1) {
    //ans += -dgamma(pow(sigma_delta,-2), Tau_D_dist, 1/Tau_D_dist, true);
    //ans += -dgamma(pow(sigma_delta,-2), Type(3), Type(1), true);
    ans += -dgamma(pow(sigma_delta,-2), Type(0.75), Type(1), true);
    ans += Type(2)*log(pow(sigma_delta,2));//Jacobian adjustment when invgamma penalty on variance
    //ans += -dgamma(pow(sigma_nu,-2), Tau_D_dist, 1/Tau_D_dist, true);
    //ans += -dgamma(pow(sigma_nu,-2), Type(3), Type(1), true);
    ans += -dgamma(pow(sigma_nu,-2), Type(0.75), Type(1), true);
    ans += Type(2)*log(pow(sigma_nu,2));///Jacobian adjustment when invgamma penalty on variance
  }
  
  // Half cauchy prior on sigma_delta and sigma_nu
  if (SigDeltaPriorCauchy == 1) {
    //3 June 2021. abs() no longer works in TMB
    //ans += -abs( dt( sigma_delta, Type(1), true));
    //ans += - abs( dt( sigma_nu, Type(1), true ));
    ans += -( dt( sigma_delta, Type(1), true));
    ans += -( dt( sigma_nu, Type(1), true ));
  }
  
    // Get predicted values for plotting  WA regresssion with CIs
  int N_pred = PredlnWA.size();
  vector <Type> PredlnSMSYs_CI(N_pred);
  vector <Type> PredlnSMSYo_CI(N_pred);
  vector <Type> PredlnSREPs_CI(N_pred);
  vector <Type> PredlnSREPo_CI(N_pred);
  
  for (int i=0; i<N_pred; i++){
    PredlnSMSYs_CI(i) = logDelta1 + exp(logDelta2) * PredlnWA(i);
    PredlnSMSYo_CI(i) = logDelta1 + logDelta1ocean + (exp(logDelta2) + Delta2ocean) * PredlnWA(i);
    PredlnSREPs_CI(i) = logNu1 + exp(logNu2) * PredlnWA(i);
    PredlnSREPo_CI(i) = logNu1 + logNu1ocean + (exp(logNu2) + Nu2ocean) * PredlnWA(i);
  }
  
  //// Get predicted values for stream-type Test stocks with CIs
  //int N_tests = TestlnWAs.size();
  //vector <Type> TestlnSMSYs(N_tests);
  //vector <Type> TestlnSREPs(N_tests);
  
  //for (int i=0; i<N_tests; i++){
  //  TestlnSMSYs(i) = logDelta1 + exp(logDelta2) * TestlnWAs(i);
  //  TestlnSREPs(i) = logNu1 + exp(logNu2) * TestlnWAs(i);
  //}
  
  
  ///Get predicted values for ocean-type Test stocks with CIs
  int N_testo = TestlnWAo.size();
  vector <Type> TestlnSMSYo(N_testo);
  vector <Type> TestlnSREPo(N_testo);
  
  for (int i=0; i<N_testo; i++){
    TestlnSMSYo(i) = logDelta1 + logDelta1ocean + (exp(logDelta2) + Delta2ocean) * TestlnWAo(i);
    TestlnSREPo(i) = logNu1 + logNu1ocean + (exp(logNu2) + Nu2ocean) * TestlnWAo(i);
  }
  
  vector <Type> lnSMSY = log(SMSY_std*Scale);
  vector <Type> lnSREP = log(SREP_std*Scale);
  
  
  ADREPORT(SMSY_std);
  ADREPORT(SREP_std);
  ADREPORT(LogRS_Pred_std);
  ADREPORT(PredlnSMSY);
  ADREPORT(PredlnSREP);
  ADREPORT(lnSMSY);
  ADREPORT(lnSREP);
  ADREPORT(PredlnSMSYs_CI);
  ADREPORT(PredlnSMSYo_CI);
  ADREPORT(PredlnSREPs_CI);
  ADREPORT(PredlnSREPo_CI);
  ADREPORT(TestlnSMSYo);
  ADREPORT(TestlnSREPo);
  //ADREPORT(TestlnSMSYs);
  //ADREPORT(TestlnSREPs);
  REPORT(nLL_std);
  return ans;
}

