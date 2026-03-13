#include <TMB.hpp>

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
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(N_Stks);
  DATA_VECTOR(LM_Agg_Abund);
  DATA_VECTOR(N_Above_BM);
  DATA_VECTOR(Pred_Abund);
  DATA_SCALAR(p);
  DATA_SCALAR(B_penalty_mu);
  DATA_SCALAR(B_penalty_sig);
  DATA_INTEGER(Penalty);
  DATA_INTEGER(Bern_logistic);
  
  
  PARAMETER(B_0);
  PARAMETER(B_1);
  
  Type ans=0.0;
  
  // get number of years
  int Logistic_Mod_Yrs = LM_Agg_Abund.size();

 
  // vectors of number of "trials" for logistic mods
  vector<Type> N_bin(Logistic_Mod_Yrs);
  vector<Type> N_bern(Logistic_Mod_Yrs);
  // vector for logistic likelihood
  vector<Type> LogitP(Logistic_Mod_Yrs);
  // Also version with 0/1's if doing bernoulli
  vector <Type> All_Above_BM(Logistic_Mod_Yrs);
  // set to 0
  All_Above_BM.setZero();
  
  // create logistic likelihood and
  // Fill Ns with values
  for(int i=0; i<Logistic_Mod_Yrs; ++i){
    LogitP(i) = B_0 + B_1*LM_Agg_Abund(i);
    N_bin(i) = N_Stks;
    N_bern(i) = 1;
    // Also fill in bernoulli version with 1s when all above LRP
    if(N_Above_BM(i) == N_Stks){
      All_Above_BM(i) = 1;
    }
    
  }
  
  // Fit logistic model
  
  if(Bern_logistic == 1){
    ans += -sum(dbinom_robust(All_Above_BM, N_bern, LogitP, true));
  } else if(Bern_logistic == 0) {
    ans += -sum(dbinom_robust(N_Above_BM, N_bin, LogitP, true));
  }
  
  
  
  // Add log-normal penalty on the aggregate abundance at p=0.01. AggAbun= ( log(p/(1-p) - B_0)/B_1
  Type p_min = 0.01;
  if(Penalty == 1) {
    ans += -dnorm(  (log(p_min/(1-p_min)) - B_0) /B_1, B_penalty_mu , B_penalty_sig, true );
  }
  
  // Alternative formulations of prior
  // B_penalty_mu and sig are the mean and sig in log-space
  //ans += -dlnorm( (log(p_min/(1-p_min)) - B_0) /B_1, B_penalty_mu, B_penalty_sig, true );
  // dlnorm not incluced in TMB.hpp!
  
  // OR Add normal penalty on the log(aggregate abundance) at p=0.01. AggAbun= ( log(p/(1-p) - B_0)/B_1
  //ans += -dnorm( log( (log(p_min/(1-p_min)) - B_0) /B_1), B_penalty_mu, B_penalty_sig, true );
  
  // OR Add normal penalty on the aggregate abundance at p=0.01. AggAbun= ( log(p/(1-p) - B_0)/B_1
  //Type raw_B_mu = 0;//exp( B_penalty_mu + B_penalty_sig^2/2 );
  //Type raw_B_sig = 10;//sqrt( (exp ( B_penalty_sig^2) - 1) * exp( 2 * B_penalty_mu  + B_penalty_sig^2) );
  
  

    //Get final BM
  Type Agg_LRP = (log(p/(1-p)) - B_0)/(B_1);
   
   // Get estimates for plotting CIs
  int N_Preds = Pred_Abund.size();
  vector<Type> Logit_Preds(N_Preds);
  
  Logit_Preds = B_0 + B_1*Pred_Abund;
  
  REPORT(N_Above_BM);
  REPORT(All_Above_BM);
  REPORT(ans);
  ADREPORT(Agg_LRP);
  ADREPORT(Logit_Preds);
  
  return ans;
  
}
