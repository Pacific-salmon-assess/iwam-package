#include <TMB.hpp >

// Based on IWAM_FixedSep_Ricstd.cpp from the watershed-area-model Repository

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
Type LambertW(Type x)
  {
    CppAD::vector<Type> tx(1);
    tx[0] = x;
    return LambertW(tx)[0];
  }
  
  
template <class Type> 
vector <Type> minus_one_to_one(vector <Type> x) 
  { 
    return Type(2) * invlogit(x) - Type(1); 
  } 
  
template<class Type>
Type objective_function<Type>:: operator() ()
{
  DATA_VECTOR(S);
  DATA_VECTOR(logRS);
  
  DATA_IVECTOR(stk);
  DATA_IVECTOR(yr);

  DATA_VECTOR(WAbase);
  DATA_VECTOR(scale);
  DATA_IVECTOR(lifehist);
  DATA_INTEGER(biasCor);

  // DATA_VECTOR(pred_lnWA);
  // DATA_VECTOR(target_lnWA_ocean);
  // DATA_VECTOR(target_lnWA_stream);
  
  PARAMETER_VECTOR(logA);
  PARAMETER_VECTOR(logB);
  PARAMETER_VECTOR(logSigma);
  // PARAMETER(logDelta1);
  // PARAMETER(logDelta1_ocean);
  // PARAMETER(logDelta2);
  // PARAMETER(Delta2_ocean);
  // PARAMETER(logDeltaSigma);
  
  // PARAMETER(logNu1);
  // PARAMETER(logNu1_ocean);
  // PARAMETER(logNu2);
  // PARAMETER(Nu2_ocean);
  // PARAMETER(logNuSigma);

  Type ans=0.0;
  int N_Obs = S.size(); 
  int N_stks = scale.size(); 

  vector <Type> logRS_pred(N_Obs);
  vector <Type> sigma = exp(logSigma);
  // Type sigmaA = exp(logSigmaA);
  vector <Type> nLL(N_Obs);
  
  // Standard Ricker model
  for (int i = 0; i<N_Obs; i++){
    if(biasCor == 0 ){
      logRS_pred(i) = logA(stk(i)) - exp(logB(stk(i))) * S(i);
      }
    if(biasCor == 1 ){
      logRS_pred(i) = logA(stk(i)) - exp(logB(stk(i))) * S(i) - pow(sigma(stk(i)), 2) / Type(2);
    }
    ans += -dnorm(logRS_pred(i), logRS(i),  sigma(stk(i)), true);
    nLL(i) = -dnorm(logRS_pred(i), logRS(i),  sigma(stk(i)), true);
  }

  // Code for SMSY SREP 
  vector <Type> SMSY(N_stks);  
  vector <Type> SREP(N_stks);

  for(int i=0; i<N_stks; i++){
    SMSY(i) =  (1 - LambertW( exp (1- logA(i)) ) ) / exp(logB(i)) ; //
  }
  SREP = logA / exp(logB);
  
  
  
  //Liermann's model with both stream and ocean type=================
  // vector <Type> pred_lnSMSY(N_stks);
  // vector <Type> pred_lnSREP(N_stks);
  // 
  // Type sigma_delta = exp(logDeltaSigma);
  // Type sigma_nu = exp(logNuSigma);
  // 
  // for (int i=0; i<N_stks; i++){
  //   pred_lnSMSY(i) = logDelta1 + logDelta1_ocean * lifehist(i) + ( exp(logDelta2) + Delta2_ocean * lifehist(i) ) * log(WAbase(i)) ;
  //   ans += -dnorm(pred_lnSMSY(i), log(SMSY(i) * scale(i) ),  sigma_delta, true);
  // }
  // for (int i=0; i<N_stks; i++){
  //   if(biasCor == 0) {
  //     pred_lnSMSY(i) = logDelta1 + logDelta1_ocean * lifehist(i) + ( exp(logDelta2) + Delta2_ocean * lifehist(i) ) * log(WAbase(i)) ;
  //   }
  //   if(biasCor == 1) {
  //     pred_lnSMSY(i) = logDelta1 + logDelta1_ocean * lifehist(i) + ( exp(logDelta2) + Delta2_ocean * lifehist(i) ) * log(WAbase(i)) - pow(sigma_delta,2) / Type(2);
  //   }
  //   ans += -dnorm( pred_lnSMSY(i), log(SMSY(i) * scale(i) ),  sigma_delta, true);
  //   
  //   if(biasCor == 0) {
  //     pred_lnSREP(i) = logNu1 + logNu1_ocean * lifehist(i) + ( exp(logNu2) + Nu2_ocean * lifehist(i) ) * log(WAbase(i)) ;
  //   }
  //   if(biasCor == 1) {
  //     pred_lnSREP(i) = logNu1 + logNu1_ocean * lifehist(i) + ( exp(logNu2) + Nu2_ocean * lifehist(i) ) * log(WAbase(i))  - pow(sigma_nu,2) / Type(2);
  //   }
  //   
  //   ans += -dnorm( pred_lnSREP(i), log(SREP(i) * scale(i) ),  sigma_nu, true);
  // }

    
// Get predicted values for plotting  WA regression with CIs
// int N_pred = pred_lnWA.size();
// vector <Type> pred_lnSMSY_stream_CI(N_pred);
// vector <Type> pred_lnSMSY_ocean_CI(N_pred);
// vector <Type> pred_lnSREP_stream_CI(N_pred);
// vector <Type> pred_lnSREP_ocean_CI(N_pred);
// 
// for (int i=0; i<N_pred; i++){
//   pred_lnSMSY_stream_CI(i) = logDelta1 + exp(logDelta2) * pred_lnWA(i);
//   pred_lnSMSY_ocean_CI(i) = logDelta1 + logDelta1_ocean + (exp(logDelta2) + Delta2_ocean) * pred_lnWA(i);
// //   pred_lnSREP_stream_CI(i) = logNu1 + exp(logNu2) * pred_lnWA(i);
// //   pred_lnSREP_ocean_CI(i) = logNu1 + logNu1_ocean + (exp(logNu2) + Nu2_ocean) * pred_lnWA(i);
// }
// 
// //// Get predicted values for stream-type target stocks with CIs
// int N_target_stream = target_lnWA_stream.size();
// vector <Type> target_lnSMSY_stream(N_target_stream);
// vector <Type> target_lnSREP_stream(N_target_stream);
// 
// for (int i=0; i<N_target_stream; i++){
//   target_lnSMSY_stream(i) = logDelta1 + exp(logDelta2) * target_lnWA_stream(i);
//   target_lnSREP_stream(i) = logNu1 + exp(logNu2) * target_lnWA_stream(i);
// }
// 
// ///Get predicted values for ocean-type target stocks with CIs
// int N_target_ocean = target_lnWA_ocean.size();
// vector <Type> target_lnSMSY_ocean(N_target_ocean);
// vector <Type> target_lnSREP_ocean(N_target_ocean);
// 
// for (int i=0; i<N_target_ocean; i++){
//   target_lnSMSY_ocean(i) = logDelta1 + logDelta1_ocean + (exp(logDelta2) + Delta2_ocean) * target_lnWA_ocean(i);
//   target_lnSREP_ocean(i) = logNu1 + logNu1_ocean + (exp(logNu2) + Nu2_ocean) * target_lnWA_ocean(i);
// }

  vector <Type> lnSMSY = log(SMSY*scale);
  vector <Type> lnSREP = log(SREP*scale);
  
  ADREPORT(SMSY);
  ADREPORT(SREP)
  
  ADREPORT(logRS_pred);
  // ADREPORT(pred_lnSMSY); // lnSMSY for synoptic set
  // ADREPORT(pred_lnSREP);
  ADREPORT(lnSMSY);
  ADREPORT(lnSREP); 

  // ADREPORT(pred_lnSMSY_stream_CI); // synoptic "predicted" lnSMSY along a line
  // ADREPORT(pred_lnSMSY_ocean_CI);
  // ADREPORT(pred_lnSREP_stream_CI);
  // ADREPORT(pred_lnSREP_ocean_CI);
  
  // ADREPORT(target_lnSMSY_ocean); // actual predicted lnSMSY and lnSREP from watershed model for WCVI examples
  // ADREPORT(target_lnSREP_ocean);
  // ADREPORT(target_lnSMSY_stream);
  // ADREPORT(target_lnSREP_stream);
  
  REPORT(nLL);
  return ans;
}
