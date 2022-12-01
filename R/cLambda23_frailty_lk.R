#' likelihood function
#'
#' @description evaluate likelihood ...
#'
#' @param weibull model parameters
#'
#' @return likelihood
#'
#' @examples
#' example(cLambda123_frailty_lk(x = 1, xdata1 = 0, omega2 = 1, shape23 = 1, scale23 = 1, c23 = 1, beta23_1 = 0))
cLambda23_frailty_lk<-function(x,xdata1,v_predict, omega2, shape23, scale23, theta23, c23, beta23_1) {
 return(-pweibull(x, scale=1/scale23, shape = shape23, lower = F, log = T)*exp(-(beta23_1*xdata1+c23*omega2+theta23*v_predict)))}
