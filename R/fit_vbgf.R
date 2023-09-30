#' Fit von Bertalanffy Growth Function
#'
#' Fit von Bertalanffy Growth Function to length at age data using maximum likelihood estimation.
#' @param len_obs The observed length of fish.
#' @param age_obs The observed ages of fish.
#' @param linf_strt The starting value for L-infinity (average length at maximum age).
#' @param k_strt The starting value for Kappa - growth coefficent.
#' @param t0_start The starting value for t0 (age at zero length).
#' @param logSigma The starting value for the logsigma parameter.
#'     If no value specified, then 0 will be used.
#' @return A summary of the estimated parameters and plot
#' @examples
#' #library(TMB)
#' fit_vbgf(len_obs, age_obs, linf_strt, k_strt, t0_start, logsigma = 0)
#' @export
fit_vbgf <- function (len_obs, age_obs, linf_strt, k_strt, t0_start, logsigma = 0){

  # Get min and max age for plots
  age_obs_min = min(age_obs)
  age_obs_max = max(age_obs)

  # Plot length-at-age
  plot(age_obs, len_obs, xlab = "Age (yr)", ylab = "Length")


  # Load model
  dyn.load(dynlib("tmb/fit_vbgf"))


  # Run the model
  obj = MakeADFun(
    data = list(len = len_obs, age = age_obs),
    parameters = list(linf = linf_strt, k = k_strt, t0 = t0_start,
                      logsigma = logsigma),
    DLL = "fit_vbgf") # Prepare model object, data, and initial parameter estimates for fitting


  # Fit model using nlminb optimiser
  opt1 = nlminb(obj$par,
                obj$fn,
                obj$gr)

  summ = summary(sdreport(obj)) # Get summary of parameter estimates
  summ

  # Summary of estimated parameters
  summ = summary(sdreport(obj))
  return(summ)

  # Plot the fit
  curve(summ[1,1] * (1 - exp(-summ[2,1] * (x - summ[3,1]))),
        from = 0,
        to = age_obs_max,
        col = 2,
        add = T)

}
