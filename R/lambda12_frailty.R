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
lambda12_frailty<-function(x,xdata1, omega1, scale12, shape12, c12, beta12_1) {
 return(shape12*(scale12)^(shape12)*(x)^(shape12-1)*exp(beta12_1*xdata1 + c12*omega1))}
