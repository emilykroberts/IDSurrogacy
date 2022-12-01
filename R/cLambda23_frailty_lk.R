#' cumulative hazard function
#'
#' @description evaluate cumulative hazard for the 2-3 transition
#'
#' @param x value at which to evaluate cumulative hazard
#' @param xdata1 covariate data
#' @param omega2 frailty term
#' @param shape23 shape of Weibull distribution
#' @param scale23 scale of Weibull distribution
#' @param c23 frailty coefficient
#' @param beta23_1 covariate coefficient
#' @param theta23 theta23 parameter in 2-3 transition
#' @param v_predict potential time that S occured to be integrated over
#' 
#' @return cumulative hazard
#'
#' @examples
#' example(cLambda123_frailty_lk(x = 1, xdata1 = 0, omega2 = 1, 
#' shape23 = 1, scale23 = 1, c23 = 1, beta23_1 = 0, theta23 = 0, v_predict = 1))
cLambda23_frailty_lk<-function(x, xdata1,v_predict, omega2, shape23, scale23, theta23, c23, beta23_1) {
 return(-pweibull(x, scale=1/scale23, shape = shape23, lower = F, log = T)*exp(-(beta23_1*xdata1+c23*omega2+theta23*v_predict)))}
