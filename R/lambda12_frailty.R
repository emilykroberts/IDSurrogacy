#' hazard function
#'
#' @description evaluate hazard for the 1-2 transition
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
#' example(lambda12_frailty(x = 1, xdata1 = 0, omega1 = 1, shape12 = 1, 
#' scale12 = 1, c12 = 1, beta12_1 = 0))
lambda12_frailty = function(x, xdata1, omega1, scale12, shape12, c12, beta12_1) {
  return(shape12 * (scale12) ^ (shape12) * (x) ^ (shape12 - 1) * exp(beta12_1 * xdata1 + c12 * omega1))
}
