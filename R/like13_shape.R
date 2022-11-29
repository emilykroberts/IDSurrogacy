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
like13_shape = function(shape13_0, dat0, c13_0, scale13_0, omega13_z0){

b = (((1/shape13_0) - 1) * (sum(log(dat0$y13 ^ dat0$s13), na.rm = T)) - 
 scale13_0 * shape13_0 * sum(dat0$y13 ^ (1/shape13_0) * exp(c13_0 * omega13_z0
 ), na.rm = T))
 
return(b)
}
