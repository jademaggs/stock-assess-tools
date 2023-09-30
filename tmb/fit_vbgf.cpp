// von Bertalanffy Growth Function
#include <TMB.hpp>

template<class Type>

Type objective_function<Type>::operator() ()
  {

  // data:
  DATA_VECTOR(age_obs);
  DATA_VECTOR(len_obs);
  int n = age_obs.size(); // get number of data points to loop over


  // parameters:
  PARAMETER(linf); // asymptoptic length
  PARAMETER(k); // growth rate
  PARAMETER(t0); // Age at length 0
  PARAMETER(logsigma); // log(residual SD)
  // Fit sigma on a log scale to keep it > 0


  // procedures: (transformed parameters)
  Type sigma = exp(logsigma);
  vector<Type> len_pred(n);


  // Report transformed variables
  ADREPORT(sigma);


  // Initialize negative log-likelihood
  Type nll = 0.0;


  // Fit the model
  for(int i = 0; i < n; i++){ // C++ starts loops at 0!
  // get negative log-likelihood (last argument is log = TRUE)
  len_pred(i) = linf * (1.0-exp(-k*(age_obs(i) - t0))) ;
  }


  // Calculate the negative log-likelihood
  nll = -sum(dnorm(len_obs, len_pred, sigma, true));
  // NOTE: the negative log-likelihood can also be included in the for loop using the subtraction assignment '-=' where:
  // nll -= dnorm(Length, LengthPred, sigma, true)
  // This will subtract each log-likelihood from the initial 0 value


  return nll;
}
