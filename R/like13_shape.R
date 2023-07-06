#' likelihood function
#'
#' @description evaluate likelihood ...
#'
#' @param omega13_z0 frailty term
#' @param shape13_0 shape of Weibull distribution
#' @param scale13_0 scale of Weibull distribution
#' @param c13_0 frailty coefficient
#' @param dat0 data
#'
#' @return likelihood
#'
#' @examples 
#' example(like13_shape(shape13_0 = 1, dat0 = 1, c13_0 = 1, 
#' scale13_0 = 1, omega13_z0 = 1))
like13_shape = function(shape13_0, dat0, c13_0, scale13_0, omega13_z0){
  
  b = (((1 / shape13_0) - 1) * (sum(log(dat0$y13 ^ dat0$s13), na.rm = T)) - 
         scale13_0 * shape13_0 * sum(dat0$y13 ^ (1 / shape13_0) * exp(c13_0 * omega13_z0
         ), na.rm = T))
  
  return(b)
}
