#' cumulative hazard function
#'
#' @description evaluate cumulative hazard for the 1-3 transition
#'
#' @param
#' x: value at which to evaluate cumulative hazard
#' xdata1: covariate data
#' omega2: frailty term
#' shape13: shape of Weibull distribution
#' scale13: scale of Weibull distribution
#' c13: frailty coefficient
#' beta13_1: covariate coefficient
#' 
#' @return cumulative hazard
#'
#' @examples
#' example(cLambda13_frailty_lk(x = 1, xdata1 = 0, omega2 = 1, shape13 = 1, scale13 = 1, c13 = 1, beta13_1 = 0))
cLambda13_frailty_lk<-function(x, xdata1, omega2, scale13, shape13, c13, beta13_1) {
 return(-pweibull(x, scale=1/scale13, shape = shape13, lower = F, log = T)*exp(beta13_1*xdata1 + c13*omega2))
 }
