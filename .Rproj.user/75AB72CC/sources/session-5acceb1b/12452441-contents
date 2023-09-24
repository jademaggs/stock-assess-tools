// Weight at length - maximum likelihood estimation

#include <TMB.hpp>

template<class Type>

Type objective_function<Type>::operator() ()
{
  // Data
  DATA_VECTOR(len);
  DATA_VECTOR(wgt);
  int n = wgt.size(); // get number of data points to loop over
  
  // Parameters
  PARAMETER(a);
  PARAMETER(b);
  PARAMETER(logSigma);
  
  Type sigma = exp(logSigma);
  vector<Type> wgt_pred(n);
  
  // Report transformed variables
  ADREPORT(sigma);
  
  // Fit the model
  for(int i = 0; i < n; i++){
    wgt_pred(i) = a * pow(len(i), b);
  }
  
  // Initialise negative log-likelihood
  Type nll = 0.0;
  
  // Get negative log-likelihood
  nll = -sum(dnorm(wgt, wgt_pred, sigma, true));
  
  return nll;
} // remember to end the .cpp file with an additional line to avoid error message 
