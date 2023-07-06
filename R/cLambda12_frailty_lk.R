#' cumulative hazard function
#'
#' @description evaluate cumulative hazard for the 1-2 transition
#'
#' @param x value at which to evaluate cumulative hazard
#' @param xdata1 covariate data
#' @param omega1 frailty term
#' @param shape12 shape of Weibull distribution
#' @param scale12 scale of Weibull distribution
#' @param c12 frailty coefficient
#' @param beta12_1 covariate coefficient
#' 
#' @return cumulative hazard
#'
#' @examples
#' example(cLambda12_frailty_lk(x = 1, xdata1 = 0, omega1 = 1, 
#' shape12 = 1, scale12 = 1, c12 = 1, beta12_1 = 0))
cLambda12_frailty_lk = function(x, xdata1, omega1, shape12, scale12, c12, beta12_1) {
  return(-pweibull(x, scale = 1 / scale12, shape = shape12, lower = F, log = T) * exp(beta12_1 * xdata1 + c12 * omega1))
}
