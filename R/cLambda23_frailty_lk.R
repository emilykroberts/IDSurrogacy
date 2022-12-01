#' cumulative hazard function
#'
#' @description evaluate cumulative hazard for the 2-3 transition
#'
#' @param
#' x: value at which to evaluate cumulative hazard
#' xdata1: covariate data
#' omega2: frailty term
#' shape23: shape of Weibull distribution
#' scale23: scale of Weibull distribution
#' c23: frailty coefficient
#' beta23_1: covariate coefficient
#' 
#' @return cumulative hazard
#'
#' @examples
#' example(cLambda123_frailty_lk(x = 1, xdata1 = 0, omega2 = 1, shape23 = 1, scale23 = 1, c23 = 1, beta23_1 = 0))
cLambda23_frailty_lk<-function(x, xdata1,v_predict, omega2, shape23, scale23, theta23, c23, beta23_1) {
 return(-pweibull(x, scale=1/scale23, shape = shape23, lower = F, log = T)*exp(-(beta23_1*xdata1+c23*omega2+theta23*v_predict)))}
