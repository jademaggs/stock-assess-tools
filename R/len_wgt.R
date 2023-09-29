#' Fit length-weight model
#'
#' Fit length-weight model using maximum likelihood estimation.
#' @param len_obs The observed length of fish.
#' @param wgt_obs The observed weight of fish.
#' @param a_strt The starting value for the a parameter.
#' @param b_strt The starting value for the b parameter.
#' @param logSigma The starting value for the logSigma parameter.
#'     If no value specified, then 0 will be used
#' @return A summary of the estimated parameters and plot
#' @examples
#' #library(TMB)
#' fit_lw(len_obs, wgt_obs, a_strt, b_strt, logSigma = 0)
#' @export
fit_lw <- function (len_obs, wgt_obs, a_strt, b_strt, logSigma = 0){

  # Get min and max lengths for plots
  len_obs_min = min(len_obs)
  len_obs_max = max(len_obs)

  # Plot observed data
  plot(len_obs, wgt_obs , xlab = "Length", ylab = "Weight")

  # Load model
  dyn.load(dynlib("tmb/length_weight"))

  # Run the model
  obj = MakeADFun(
    data = list(wgt = wgt_obs, len = len_obs),
    parameters = list(a = a_strt, b = b_strt, logSigma = logSigma),
    DLL = "length_weight"
  )

  # Optimisation
  opt = nlminb(obj$par, obj$fn, obj$gr)

  # Summary of estimated parameters
  summ = summary(sdreport(obj))
  return(summ)

  # Plot the predictions
  curve(summ[1,1] * x ^ summ[2,1],
        from = len_obs_min,
        to = len_obs_max,
        col = 2,
        add = TRUE,
        xlab = 'Length',
        ylab = 'Weight')
}
