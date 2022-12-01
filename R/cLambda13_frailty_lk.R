#' likelihood function
#'
#' @description evaluate likelihood ...
#'
#' @param weibull model parameters
#'
#' @return likelihood
#'
#' @examples
#' example(cLambda13_frailty_lk(x = 1, xdata1 = 0, omega2 = 1, shape13 = 1, scale13 = 1, c13 = 1, beta13_1 = 0))
cLambda13_frailty_lk<-function(x,xdata1, omega2, scale13, shape13, c13, beta13_1) {
 return(-pweibull(x, scale=1/scale13, shape = shape13, lower = F, log = T)*exp(beta13_1*xdata1 + c13*omega2))}
