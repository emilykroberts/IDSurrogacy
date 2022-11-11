#' likelihood function
#'
#' @description evaluate likelihood ...
#'
#' @param weibull model parameters
#'
#' @return likelihood
#'
#' @examples likelihood for given parameter set
#' 
cLambda12_frailty_lk<-function(x,xdata1, omega1, shape12, scale12, c12, beta12_1) {
 return( -pweibull(x, scale=1/scale12, shape= shape12, lower = F, log = T)*exp(beta12_1*xdata1 + c12*omega1))
  }
